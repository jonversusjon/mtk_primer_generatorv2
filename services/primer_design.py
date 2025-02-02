import os
import time
from typing import Dict, List, Optional, Any, Tuple
import numpy as np
from tqdm import tqdm
import primer3

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

from .mutation import find_codon_replacements_in_range
from primer_select import is_overhang_compatible


def find_off_targets_for_primer(binding_seq: str,
                                template_seq: Optional[str] = None,
                                tm_threshold: float = 45.0,  
                                min_3p_match: int = 10,  
                                max_mismatches: int = 1,  
                                mv_conc: float = 50.0,
                                dv_conc: float = 1.5,
                                dntp_conc: float = 0.2,
                                dna_conc: float = 250.0,
                                verbose: bool = True) -> List[Dict[str, float]]:
    """
    Finds potential off-target binding sites using primer3, enforcing SnapGene-like rules.
    
    Returns:
        List of off-target matches exceeding the filtering criteria.
    """
    if not isinstance(binding_seq, str):
        raise TypeError(f"Expected a string for 'binding_seq', but got {type(binding_seq)} instead.")

    if not template_seq:
        if verbose:
            print(f"⚠️ Skipping off-target analysis for primer: {binding_seq} (No template provided)")
        return []

    template_seq = template_seq.upper()
    binding_seq = binding_seq.upper()

    off_targets = []
    primer_len = len(binding_seq)
    template_len = len(template_seq)

    num_off_target_checks = template_len - primer_len + 1  

    for i in range(num_off_target_checks):
        window = template_seq[i : i + primer_len]

        match_3p = sum(1 for a, b in zip(binding_seq[-min_3p_match:], window[-min_3p_match:]) if a == b)
        if match_3p < min_3p_match:
            continue  

        result = primer3.bindings.calc_heterodimer(
            seq1=binding_seq,
            seq2=window,
            mv_conc=mv_conc,
            dv_conc=dv_conc,
            dntp_conc=dntp_conc,
            dna_conc=dna_conc
        )

        if result.tm >= tm_threshold:
            mismatches = sum(1 for a, b in zip(binding_seq, window) if a != b)

            if mismatches <= max_mismatches:  
                off_targets.append({
                    "index": i,
                    "tm": round(result.tm, 2),
                    "mismatches": mismatches
                })

    if off_targets:
        print("\n🔍 **Off-Target Primer Found** 🔍")
        print(f"Primer Binding Sequence: {binding_seq}")
        print(f"Total Off-Targets: {len(off_targets)}\n")
        for off_target in off_targets:
            print(f"   📌 Off-Target at Index: {off_target['index']}")
            print(f"   🔥 Tm: {off_target['tm']}°C")
            print(f"   ❌ Mismatches: {off_target['mismatches']}\n")

    return off_targets


def find_binding_sequence(seq: Seq, start: int, end: int, min_tm: float = 57) -> str:
    """Finds the shortest binding sequence with Tm >= min_tm."""
    for i in range(15, 36):  # Typical primer length range
        binding_seq = seq[start:end + i]
        if mt.Tm_NN(binding_seq) >= min_tm:
            return binding_seq
    return ""


def construct_primer(seq: Seq, spacer: str, bsmbi_site: str, six_nuc_seq: str, binding_seq: str, reverse=False) -> Dict[str, str]:
    """
    Constructs a forward or reverse primer, returning its components as strings.
    """
    # Ensure the sequences are strings.
    bsmbi_site = str(bsmbi_site)
    spacer = str(spacer)
    six_nuc_seq = str(six_nuc_seq)
    binding_seq = str(binding_seq)

    if reverse:
        overhang_sequence = bsmbi_site + str(Seq(six_nuc_seq).reverse_complement())
        primer_sequence = spacer + overhang_sequence + str(binding_seq)
    else:
        overhang_sequence = bsmbi_site + six_nuc_seq
        primer_sequence = spacer + overhang_sequence + binding_seq

    # Force everything to be Python strings, not Seq
    return {
        "binding_sequence": str(binding_seq),
        "overhang_sequence": str(overhang_sequence),
        "primer_sequence": str(primer_sequence),
    }


def design_primers_for_mutation(
    seq: Seq,
    nucleotide_index: int,
    new_nucleotide: str,
    spacer: str,
    bsmbi_site: str,
    min_tm: float = 57,
    template_seq: Optional[str] = None,
    verbose: bool = False
) -> Dict[str, List[Dict]]:
    """
    Designs primers for introducing a specific mutation in a DNA sequence.
    If `template_seq` is provided, off-target analysis is performed.
    """

    if not isinstance(seq, (str, Seq)):
        raise TypeError(f"seq should be a string or Seq object, got {type(seq)}")

    if not (0 <= nucleotide_index < len(seq)):
        raise IndexError(f"Index {nucleotide_index} is out of bounds for sequence of length {len(seq)}")

    if not isinstance(new_nucleotide, str) or len(new_nucleotide) != 1:
        raise TypeError(f"new_nucleotide should be a single character string, got {new_nucleotide}")

    primers = {"forward": [], "reverse": []}
    total_primers = 0

    # Step 1: Generate Primers
    for shift in range(6):
        left = nucleotide_index - (6 - 1) + shift
        right = nucleotide_index + 1 + shift

        if not (0 <= left < len(seq)) or not (0 <= right <= len(seq)):
            if verbose:
                print(f"Skipping out-of-bounds shift {shift}. Left: {left}, Right: {right}")
            continue

        six_nuc_seq = seq[:nucleotide_index] + new_nucleotide + seq[nucleotide_index + 1:]
        six_nuc_seq = six_nuc_seq[left:right]

        forward_binding_seq = find_binding_sequence(seq, right, right, min_tm)
        reverse_binding_seq = find_binding_sequence(seq.reverse_complement(), left - 1, left - 1, min_tm)

        if not forward_binding_seq or not reverse_binding_seq:
            if verbose:
                print(f"Warning: No binding sequence found for position {nucleotide_index} with shift {shift}.")
            continue

        # Construct primers
        forward_primer = construct_primer(seq, spacer, bsmbi_site, six_nuc_seq, forward_binding_seq)
        reverse_primer = construct_primer(seq, spacer, bsmbi_site, six_nuc_seq, reverse_binding_seq, reverse=True)

        forward_site_count = forward_primer["overhang_sequence"].count(bsmbi_site)
        reverse_site_count = reverse_primer["overhang_sequence"].count(bsmbi_site)

        if verbose:
            print(f"Shift {shift}: Forward site count = {forward_site_count}, Reverse site count = {reverse_site_count}")

        if forward_site_count == 1 and reverse_site_count == 1:
            fwd_dict = {
                "primer_sequence": forward_primer["primer_sequence"],
                "binding_sequence": forward_primer["binding_sequence"],
                "overhang_sequence": forward_primer["overhang_sequence"],
                "tm": round(mt.Tm_NN(forward_primer["primer_sequence"]), 1)
            }

            rev_dict = {
                "primer_sequence": reverse_primer["primer_sequence"],
                "binding_sequence": reverse_primer["binding_sequence"],
                "overhang_sequence": reverse_primer["overhang_sequence"],
                "tm": round(mt.Tm_NN(reverse_primer["primer_sequence"]), 1)
            }

            primers["forward"].append(fwd_dict)
            primers["reverse"].append(rev_dict)
            total_primers += 2  # Counting both forward and reverse primers

    if verbose:
        print(f"Total primers designed: {total_primers}")

    # Step 2: Off-target analysis (Only if `template_seq` is provided)
    if template_seq:
        for primer_list in ["forward", "reverse"]:
            for primer_dict in primers[primer_list]:
                try:
                    primer_dict["off_target_analysis"] = find_off_targets_for_primer(
                        primer_dict["binding_sequence"],  # ✅ Only pass binding sequence, NOT full dictionary
                        template_seq=template_seq
                    )
                except Exception as e:
                    if verbose:
                        print(f"Error in off-target analysis: {e}")
                    primer_dict["off_target_analysis"] = []
    else:
        if verbose:
            print("⚠️ Skipping off-target analysis (no template sequence provided).")

    return primers


def backtrack(
    seq: Seq,
    site_index: int,
    selected_primers: List[Tuple[Dict, Dict]],
    sites_to_mutate_list: List[Dict],
    codon_usage_dict: Dict,
    spacer: str,
    bsmbi_site: str,
    min_tm: float,
    template_seq: Optional[str],
    verbose: bool,
) -> List[Tuple[Dict, Dict]]:

    num_sites_to_mutate = len(sites_to_mutate_list)

    if site_index == num_sites_to_mutate:
        return selected_primers

    site_data = sites_to_mutate_list[site_index]  # Should be a dictionary
    # print(f"DEBUG: site_data = {site_data}")  # Check what is in site_data

    site_start = int(site_data["mutation"]["index"])  # Ensure this is an integer
    site_sequence = str(site_data["sticky_end"])  # Ensure it's a string

    print(f"🔬 Designing primers for site {site_index + 1} of {num_sites_to_mutate} (Position: {site_start})")

    mutations = get_mutations_for_site(
        seq, site_start, site_start + len(site_sequence), site_sequence, codon_usage_dict, verbose=verbose
    )

    if verbose:
        print(f"Processing restriction site at index {site_start}, {len(mutations)} possible mutations")

    for mutation in mutations:
        if not isinstance(mutation, dict):  # Ensure mutation is a dictionary
            print(f"⚠️ ERROR: Unexpected mutation format: {mutation}")
            continue  # Skip invalid mutation

        nucleotide_index = int(mutation["index"])
        new_nucleotide = mutation["new_nucleotide"]

        primers = design_primers_for_mutation(
            seq,
            nucleotide_index,
            new_nucleotide,
            spacer,
            bsmbi_site,
            min_tm,
            template_seq=template_seq,
            verbose=verbose,
        )

        if not primers["forward"] or not primers["reverse"]:
            continue  

        for forward_primer in primers["forward"]:
            for reverse_primer in primers["reverse"]:
                if selected_primers:
                    last_reverse_primer = selected_primers[-1][1]
                    if not is_overhang_compatible([
                        last_reverse_primer["overhang_sequence"],
                        forward_primer["overhang_sequence"]
                    ]):
                        continue  

                selected_primers.append((forward_primer, reverse_primer))
                result = backtrack(
                    seq,
                    site_index + 1,
                    selected_primers,
                    sites_to_mutate_list,
                    codon_usage_dict,
                    spacer,
                    bsmbi_site,
                    min_tm,
                    template_seq,
                    verbose
                )
                if result:
                    print(f"✅ Selected primers at site index {site_index}: {selected_primers}")
                    return result  

                selected_primers.pop()  

    print(f"⚠️ No valid primers found for mutation at site {site_index}, backtracking...")
    return None  



def generate_GG_edge_primers(seq, part_end_dict, part_num, primer_orientation, kozak, primer_name=None, verbose=False):
    """
    Generates a Golden Gate edge primer considering the specified part number and orientation.

    Args:
        seq (Seq or str): The DNA sequence for which primers are being designed.
        part_end_dict (dict): Dictionary containing part-specific sequences.
        part_num (str): Part number for which the primer is designed.
        primer_orientation (str): Orientation of the primer to generate ('forward' or 'reverse').
        kozak (str): Whether to yield MTK kozak or canonical kozak.
        primer_name (str, optional): Custom name for the primer. If None, defaults to MTK naming.
        verbose (bool): If True, prints debug information.

    Returns:
        tuple: A tuple containing primer name and sequence (primer_name, primer_sequence).
    """
    try:
        seq = Seq(seq)

        # Determine part-specific key and primer name suffix
        if part_num == '2' and primer_orientation == 'reverse' and kozak == 'canonical':
            part_specific_key = '2starreverse'
            primer_suffix = '2*'
        elif part_num in ['3', '3a'] and primer_orientation == 'forward' and kozak == 'canonical':
            part_specific_key = f'{part_num}starforward'
            primer_suffix = f'{part_num}*'
        else:
            part_specific_key = str(part_num) + primer_orientation
            primer_suffix = part_num

        # Retrieve the part-specific sequence
        if part_specific_key not in part_end_dict:
            raise ValueError(f"No sequence found for {part_specific_key} in part_end_dict.")
        part_specific_seq = part_end_dict[part_specific_key]

        # Generate primer sequence
        if primer_orientation == 'forward':
            n_R = 20  # Placeholder for primer length
            edge_primer = part_specific_seq + str(seq[:n_R])
            primer_name = f"MTK{primer_suffix}_{primer_name}_FW" if primer_name else f"MTK{primer_suffix}_FW"
        elif primer_orientation == 'reverse':
            n_L = 20  # Placeholder for primer length
            edge_primer = part_specific_seq + str(seq[-n_L:].reverse_complement())
            primer_name = f"MTK{primer_suffix}_{primer_name}_RV" if primer_name else f"MTK{primer_suffix}_RV"
        else:
            raise ValueError("primer_orientation must be 'forward' or 'reverse'.")

        # Optionally print debug information
        if verbose: print(f"Generated primer: {primer_name}, Sequence: {edge_primer}")

        return primer_name, edge_primer

    except Exception as e:
        # Gracefully handle errors
        print(f"Error generating primer for part {part_num} ({primer_orientation}): {e}")
        return None, None


# def count_restriction_sites(seq, recog_seq_f, recog_seq_r):
#     """
#     Counts the occurrences of a forward and reverse restriction site sequence.

#     Args:
#         seq (str): DNA sequence.
#         recog_seq_f (str): Forward recognition sequence.
#         recog_seq_r (str): Reverse recognition sequence.

#     Returns:
#         Tuple[int, int]: Counts of forward and reverse recognition sequences.
#     """
#     seq = seq.upper()
#     potential_recognition_sites = np.array([seq[i: i + len(recog_seq_f)] for i in range(len(seq) - len(recog_seq_f) + 1)])
#     forward_count = np.sum(potential_recognition_sites == recog_seq_f)
#     reverse_count = np.sum(potential_recognition_sites == recog_seq_r)
#     return forward_count, reverse_count

    
# def add_cut_sites(dna_seq, enzyme):
#     """
#     Adds restriction cut sites to a DNA sequence if missing.

#     Args:
#         dna_seq (str): The DNA sequence to modify.
#         enzyme (str): The enzyme for which cut sites are required (e.g., 'bsmbi', 'bsai').

#     Returns:
#         str: Modified DNA sequence with cut sites included.
#     """
#     recognition_sequences = {
#         'bsmbi': ('CGTCTC', 7, -5),
#         'bsai': ('GGTCTC', 7, -5)
#     }

#     if enzyme not in recognition_sequences:
#         raise ValueError(f"Unsupported enzyme: {enzyme}")

#     site, forward_cut, reverse_cut = recognition_sequences[enzyme]
    
#     # Ensure the recognition site is present in the sequence
#     if site not in dna_seq:
#         dna_seq = site + dna_seq  # Prepend recognition site if missing
#     return dna_seq
