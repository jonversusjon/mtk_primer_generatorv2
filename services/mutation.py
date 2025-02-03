from Bio.Seq import Seq
from typing import Dict, List, Optional, Tuple
import logging
from utils.utils import gc_content
from .primer_select import is_overhang_compatible

def find_alternative_codons(codon, codon_usage_dict, max_mutations=1, verbose=False):
    """
    Finds alternative codons for the amino acid encoded by the given codon,
    allowing up to a user-defined maximum number of mutations.

    Args:
        codon (str): The DNA codon to exclude (e.g., 'GAG').
        codon_usage_dict (dict): Codon usage dictionary (amino acids as keys).
        max_mutations (int): Maximum number of base changes allowed per codon (default: 1 for point mutations).
        verbose (bool): Whether to print debug information.

    Returns:
        List[Dict[str, float]]: Sorted list of dictionaries with alternative 
        codons and their usage frequencies, filtered by mutation count.
    """
    # Ensure the codon is in DNA format
    codon_dna = codon.replace('U', 'T')  # Handle RNA input
    codon_rna = codon_dna.replace('T', 'U')  # Convert to RNA

    # Determine the amino acid encoded by the codon
    amino_acid = str(Seq(codon_rna).translate())
    if verbose:
        print(f"Original Codon: {codon_dna}, Amino Acid: {amino_acid}")

    # Validate the amino acid in the codon usage dictionary
    if amino_acid not in codon_usage_dict:
        if verbose:
            print(f"Amino acid {amino_acid} not found in codon usage dictionary.")
        return []

    # Find and filter alternative codons by the allowed number of mutations
    def mutation_count(codon1, codon2):
        """Counts the number of base differences between two codons."""
        return sum(1 for a, b in zip(codon1, codon2) if a != b)

    alternative_codons = [
        {"codon": alt_codon.replace('U', 'T'), "frequency": usage}
        for alt_codon, usage in codon_usage_dict[amino_acid].items()
        if alt_codon.replace('U', 'T') != codon_dna and mutation_count(codon_dna, alt_codon.replace('U', 'T')) <= max_mutations
    ]

    # Sort alternative codons by usage frequency
    sorted_codons = sorted(alternative_codons, key=lambda x: x["frequency"], reverse=True)

    if verbose:
        print(f"Allowed Codons for {amino_acid} with ≤{max_mutations} mutations: {sorted_codons}")

    return sorted_codons


def find_codon_replacements_in_range(seq, start_idx, end_idx, codon_usage_dict, max_mutations, verbose=False):
    """
    Finds alternative codons for codons in a specified sequence range.

    Args:
        seq (str): The DNA sequence to analyze.
        start_idx (int): Start index of the range.
        end_idx (int): End index of the range.
        codon_usage_dict (dict): Codon usage dictionary.
        verbose (bool): Whether to print debug information.

    Returns:
        List[Dict]: List of proposed codon replacements.
    """
    proposed_replacements = []

    for i in range(start_idx, min(end_idx, len(seq) - 2), 3):
        codon = seq[i:i + 3]
        alternative_codons = find_alternative_codons(codon, codon_usage_dict, max_mutations=max_mutations, verbose=verbose)

        if alternative_codons:
            proposed_replacements.append({
                "original_codon": codon,
                "position": i,
                "alternative_codons": alternative_codons,
            })

    if verbose:
        print(f"Proposed codon replacements: {proposed_replacements}")

    return proposed_replacements


def gather_mutation_options(
    seq: Seq,
    sites_to_mutate: dict,
    codon_usage_dict: dict,
    spacer: str,
    bsmbi_site: str,
    min_tm: float = 57,
    template_seq: Optional[str] = None,
    verbose: bool = False,
) -> List[Tuple[Dict, Dict]]:
    """
    Designs primers with codon optimization, constraint-based design, and backtracking.
    """

    sites_to_mutate_list = []
    for enzyme, sites in sites_to_mutate.items():
        for site in sites:
            if "sequence" not in site or "position" not in site:
                raise ValueError(f"Site is missing required keys: {site}")

            site_sequence = str(site["sequence"])
            site_start = int(site["position"])
            site_end = site_start + len(site_sequence)
            sites_to_mutate_list.append((site_start, site_end, site_sequence))

    sites_to_mutate_list.sort(key=lambda x: x[0])

    # Step 1: Gather all mutation options
    mutation_sets = []
    for site_start, site_end, site_sequence in sites_to_mutate_list:
        mutations = get_mutations_for_site(
            seq, site_start, site_end, site_sequence, codon_usage_dict, verbose=verbose
        )
    
        # Ensure codon usage frequency is a number before sorting
        for mutation in mutations:
            mutation["codon_usage_frequency"] = float(mutation["codon_usage_frequency"])
    
        mutations.sort(key=lambda x: x["codon_usage_frequency"], reverse=True)
        mutation_sets.append(mutations)

    return mutation_sets


def find_best_mutation_set(mutation_options, verbose=False):
    """
    Iteratively tries mutation sets, eliminating used ones to avoid retries.
    
    Args:
        mutation_options (List[List[Dict]]): Ordered list of mutation sets (mutates in place).
        verbose (bool): Whether to print debug information.
    
    Returns:
        List[Dict]: The best valid mutation set found, or an empty list if none work.
    """
    while mutation_options:
        # Get and remove the highest-priority mutation set (first in list)
        mutation_set = mutation_options.pop(0)  

        if verbose:
            print(f"\n🔍 Attempting mutation set (Remaining: {len(mutation_options)})...")

        valid_mutation_subset = find_compatible_mutation_subsets(mutation_set, verbose=verbose)

        if valid_mutation_subset:
            if verbose:
                print("✅ Found a valid mutation set!")
            return valid_mutation_subset  # Return as soon as a valid set is found

    return []  # No valid mutation set found


def get_sticky_end_options(seq, mutation_index, sticky_end_length=4):
    """
    Determines potential sticky end sequences based on a mutation location.
    """
    sticky_end_options = []

    for shift in range(-sticky_end_length + 1, 1):  
        cut_site = mutation_index + shift  

        if cut_site < 0 or cut_site + sticky_end_length > len(seq):
            continue  

        sticky_end = seq[cut_site : cut_site + sticky_end_length]

        if 0.25 <= gc_content(sticky_end) <= 0.75:  
            sticky_end_options.append({"sticky_end": sticky_end, "cut_site": cut_site})

    return sticky_end_options


def get_mutations_for_site(seq, site_start, site_end, site_sequence, codon_usage_dict, max_mutations=1, verbose=False):
    """
    Finds and prioritizes mutations within a restriction site based on codon usage.
    Also generates potential sticky end cut sites based on valid mutations.
    """
    mutations = []
    
    # Use find_codon_replacements_in_range to get codon-level replacements
    codon_replacements = find_codon_replacements_in_range(
        seq, site_start, site_end, codon_usage_dict, max_mutations
    )
    
    for codon_data in codon_replacements:
        codon_start = codon_data["position"]
        original_codon = codon_data["original_codon"]

        for alt_codon_data in codon_data["alternative_codons"]:
            new_codon = alt_codon_data["codon"]

            # Calculate individual nucleotide changes
            for i in range(3):
                if original_codon[i] != new_codon[i]:
                    nucleotide_index = codon_start + i
                    original_nucleotide = original_codon[i]
                    new_nucleotide = new_codon[i]

                    # Check if this mutation disrupts the restriction site
                    mutated_seq = seq[:nucleotide_index] + new_nucleotide + seq[nucleotide_index + 1:]
                    mutated_site_region = mutated_seq[site_start:site_end]
                    if site_sequence in mutated_site_region:
                        print(f"⚠️ Rejected mutation {new_codon} at {codon_start} because {site_sequence} is still present at {site_start}.")
                        continue

                    # Get potential sticky end overhangs based on this mutation
                    sticky_end_candidates = get_sticky_end_options(seq, nucleotide_index)

                    # Ensure sticky_end_candidates is a **flat list**
                    sticky_end_candidates = [se for se in sticky_end_candidates if isinstance(se, dict)]
                    
                    # Filter sticky ends based on GC content and uniqueness
                    valid_sticky_ends = [
                        se for se in sticky_end_candidates
                        if 0.25 <= (sum(1 for nt in se["sticky_end"] if nt in "GC") / len(se["sticky_end"])) <= 0.75
                    ]

                    if not valid_sticky_ends:
                        continue  # Skip if no valid sticky ends

                    # Append mutation as a **flat dictionary**
                    mutations.append({
                        "index": nucleotide_index,
                        "original_nucleotide": original_nucleotide,
                        "new_nucleotide": new_nucleotide,
                        "original_codon": original_codon,
                        "new_codon": new_codon,
                        "codon_usage_frequency": alt_codon_data["frequency"],
                        "sticky_end_options": valid_sticky_ends,  # Store potential sticky ends
                    })

    # 🛠 Flatten any nested lists (prevention)
    mutations = [m for m in mutations if isinstance(m, dict)]

    # Sort mutations by codon usage frequency (descending)
    mutations.sort(key=lambda x: x["codon_usage_frequency"], reverse=True)

    return mutations


def find_compatible_mutation_subsets(mutation_sets, verbose=False):
    """
    Finds a subset of mutations where at least one mutation per restriction site has a compatible set of sticky ends.
    """
    if verbose:
        print("\n🔬 Searching for a valid subset of point mutations...")

    all_mutation_options = []
    
    for mutation_set in mutation_sets:
        # Ensure mutation_set is always treated as a list
        if isinstance(mutation_set, dict):  # If it's a dict, wrap it in a list
            mutation_set = [mutation_set]

        mutation_options = []
        for mutation in mutation_set:  # Now mutation_set is correctly treated as a list
            if not isinstance(mutation, dict) or "sticky_end_options" not in mutation:
                print(f"⚠️ Unexpected data format: {mutation}")
                continue  # Skip malformed entries
            
            if not isinstance(mutation["sticky_end_options"], list):
                print(f"⚠️ sticky_end_options should be a list but got: {type(mutation['sticky_end_options'])}")
                continue  # Skip malformed entries
            
            for sticky_end in mutation["sticky_end_options"]:
                if not isinstance(sticky_end, dict):
                    print(f"⚠️ Expected dict but got {type(sticky_end)}: {sticky_end}")
                    continue  # Skip malformed entries
                
                mutation_options.append({
                    "mutation": mutation,
                    "sticky_end": str(sticky_end.get("sticky_end", "UNKNOWN")),  # Ensure string type
                    "cut_site": sticky_end.get("cut_site", "UNKNOWN")
                })

        all_mutation_options.append(mutation_options)

    # Try to find a valid subset
    selected_mutations = []
    for mutation_group in all_mutation_options:
        for candidate in mutation_group:
            # Check if sticky ends are compatible with previously selected ones
            if not selected_mutations or all(
                is_overhang_compatible([prev["sticky_end"], candidate["sticky_end"]])
                for prev in selected_mutations
            ):
                selected_mutations.append(candidate)
                break  # Select only one mutation per site

    if verbose:
        print(f"✅ Found valid subset with {len(selected_mutations)} compatible mutations.")

    return selected_mutations if len(selected_mutations) == len(mutation_sets) else None
