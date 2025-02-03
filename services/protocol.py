from Bio.Seq import Seq
from typing import List, Dict, Optional
from .sequence_prep import (
    adjust_sequence_for_frame_and_codons,
    find_bsmbi_bsai_sites,
    summarize_bsmbi_bsai_sites
)
from .primer_design import generate_GG_edge_primers, get_all_possible_internal_primers
from .primer_select import select_best_internal_primers, format_primers_for_output
from .mutation import *
# from utils.utils import plot_primers_for_mutation_site, rank_and_print_mutation_sets


def create_gg_protocol(
    seq: List[str],
    part_end_dict: Dict[str, str],
    species_codon_usage: Dict[str, Dict[str, float]],
    part_num_left: List[str],
    part_num_right: List[str],
    max_mutations: int,
    primer_name: Optional[List[str]] = None,
    template_seq: Optional[str] = None,
    kozak: str = "MTK",
    verbose: bool = False,
    output_tsv_path: str = "designed_primers.tsv",
) -> List[List[str]]:
    """Main function to orchestrate the Golden Gate protocol creation."""

    print("Starting Golden Gate protocol creation...")
    primer_data = []

    for i, single_seq in enumerate(seq):
        single_seq = preprocess_sequence(single_seq)
        sites_to_mutate = find_and_summarize_sites(single_seq, i, verbose)
        
        num_sites_to_mutate = sum(len(sites) for sites in sites_to_mutate.values())
        print(f"Number of sites to mutate: {num_sites_to_mutate}")
        
        # Step 1: Gather mutation options sorted by codon usage frequency
        mutation_options = gather_mutation_options(
            seq=single_seq,
            sites_to_mutate=sites_to_mutate,
            codon_usage_dict=species_codon_usage,
            spacer="GAA",
            bsmbi_site="CGTCTC",
            min_tm=57,
            template_seq=template_seq,
            verbose=verbose
        )
        # rank_and_print_mutation_sets(mutation_options, verbose)
        print(f'mutation_options: {mutation_options}')
        # Step 2: Iterate while mutation options exist
        selected_mutations = []
        while mutation_options and not selected_mutations:
            selected_mutations = find_best_mutation_set(mutation_options, verbose)
        
        # Step 3: Proceed with primer design
        if selected_mutations:
            primer_sets = get_all_possible_internal_primers(
                seq=single_seq,
                seq_index=i,
                proposed_mutations=selected_mutations,
                part_num_left=part_num_left[i],
                part_num_right=part_num_right[i],
                primer_name=primer_name[i],
                verbose=verbose
            )
            internal_primers = select_best_internal_primers(primer_sets, template_seq)
            primer_data += format_primers_for_output(internal_primers, primer_name[i])

        else:
            print(f"❌ No valid mutation subset found for sequence {i + 1}. Skipping primer design.")

       # Step 4: Design edge primers
        primer_data += get_edge_primers(
            single_seq, i, part_end_dict, part_num_left, part_num_right, kozak, primer_name
        )

        # Save and return results
        save_primers_to_tsv(primer_data, output_tsv_path)
        
        return primer_data

def preprocess_sequence(sequence: str) -> Seq:
    """Converts sequence to uppercase and adjusts for frame and codons."""
    sequence = Seq(sequence.upper())
    return adjust_sequence_for_frame_and_codons(sequence)

def find_and_summarize_sites(sequence: Seq, index: int, verbose: bool) -> Dict:
    """Finds and summarizes restriction sites needing mutation."""
    sites_to_mutate = find_bsmbi_bsai_sites(index, sequence, verbose=verbose)
    if sites_to_mutate:
        summarize_bsmbi_bsai_sites(sites_to_mutate)
    return sites_to_mutate


def get_edge_primers(
    sequence: Seq,
    index: int,
    part_end_dict: Dict[str, str],
    part_num_left: List[str],
    part_num_right: List[str],
    kozak: str,
    primer_name: Optional[List[str]]
) -> List[List[str]]:
    """Fetches edge primers for a given sequence and returns primer data."""
    primer_data = []
    left_name, left_primer = generate_GG_edge_primers(
        sequence, part_end_dict, part_num_left[index], "forward", kozak, primer_name
    )
    right_name, right_primer = generate_GG_edge_primers(
        sequence, part_end_dict, part_num_right[index], "reverse", kozak, primer_name
    )
    primer_data.append([left_name, left_primer, f"Amplicon_{index + 1}"])
    primer_data.append([right_name, right_primer, f"Amplicon_{index + 1}"])
    
    return primer_data


def save_primers_to_tsv(primer_data: List[List[str]], output_tsv_path: str):
    """Saves primer data to a TSV file."""
    if not primer_data:
        print("Warning: No primer data to save.")
        return  # Exit early if no data

    try:
        with open(output_tsv_path, "w") as tsv_file:
            tsv_file.write("Primer Name\tSequence\tAmplicon\n")
            for row in primer_data:
                tsv_file.write("\t".join(map(str, row)) + "\n")
    except IOError as e:
        print(f"Error writing to file {output_tsv_path}: {e}")


