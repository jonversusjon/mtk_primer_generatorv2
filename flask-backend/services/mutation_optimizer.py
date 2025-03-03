from typing import Dict, List
from Bio.Seq import Seq
import numpy as np
from itertools import product
from .utils import GoldenGateUtils
import random
from tqdm import tqdm
import logging
from config.logging_config import logger
from services.base import debug_context

class MutationOptimizer:
    """
    Optimizes mutations for Golden Gate assembly by balancing codon usage,
    sequence stability, and restriction site compatibility.
    """

    def __init__(self, verbose: bool = False):
        self.logger = logger.getChild("MutationOptimizer")
        self.utils = GoldenGateUtils()
        
        """
        Initialize the optimizer with a precomputed compatibility table.
        
        Args:
            verbose (bool): If True, provide more user-friendly logs in production.
        """
        self.verbose = verbose

        # Load compatibility table
        self.compatibility_table_path = 'static/data/compatibility_table.bin'
        self.compatibility_table = self.utils._load_compatibility_table(self.compatibility_table_path)

        self.logger.debug("MutationOptimizer initialized with compatibility table.")

        if self.verbose:
            self.logger.info("MutationOptimizer is running in verbose mode.")


    def optimize_mutations(
        self,
        sequence: str,
        mutation_options: Dict
    ) -> List[Dict]:
        """
        Find the optimal combination of mutations that are compatible according to the BsmBI sticky-end strategy.

        Args:
            sequence (str): The full sequence being analyzed.
            sites_to_mutate (Dict): Dictionary of restriction sites to mutate.
            mutation_options (Dict): Mutation site properties and possible mutations.

        Returns:
            List[Dict]: A list of optimized mutation sets.
        """
        with debug_context("mutation_optimization"):
            logger.info("Step 1: Predicting BsmBI reassembly overhangs for each mutation option...")
            mutation_options_with_overhangs = self.add_predicted_overhangs(sequence, mutation_options)
            
            logger.info("Step 2: Generating possible mutation sets...")
            mutation_sets = self.generate_mutation_sets(mutation_options_with_overhangs)
            logger.debug(f"Generated {len(mutation_sets)} mutation sets.")
            
            logger.info("Step 3: Computing compatibility matrices...")
            compatibility_matrices = self.create_compatibility_matrices(mutation_sets)

            logger.info("Step 4: Filtering mutation sets based on compatibility...")
            optimized_mutations = self.filter_compatible_mutations(mutation_sets, compatibility_matrices)
            
            logger.debug(f"\nFinal number of optimized mutation sets: {len(optimized_mutations)}")

            return optimized_mutations, compatibility_matrices


    def add_predicted_overhangs(self, sequence: str, mutation_options: Dict) -> Dict:
        """
        Captures four possible BsmBI reassembly overhangs for each proposed alternative codon.
        Each alternative codon contributes four different overhangs based on where the mutated base is positioned.

        This version extracts the exact 4-bp overhang while also tracking an extended 6-bp sequence
        that includes one base on either side of the overhang.
        """
        
        def calculate_overhangs(sequence, position, offset, utils):
            """
            Extracts predicted overhangs from a given sequence at a specified mutation position.
            Ensures correct 4-bp overhang determination and maintains top/bottom strand consistency.

            Args:
                sequence (str): The full sequence being analyzed.
                position (int): The base position of the codon in the sequence.
                offset (int): The relative index of the mutation within the codon.
                utils: Utility functions for sequence manipulation.

            Returns:
                Dict: Contains extracted overhangs and extended sequences for both strands.
            """
            mutation_position = position + offset

            # Define base slicing parameters
            OVERHANG_LENGTH = 4
            EXTENDED_LENGTH = 6  # Overhang + 1 nt on either side

            # Ensure valid slicing boundaries
            start_index = max(0, mutation_position - 2)
            end_index = min(mutation_position, len(sequence) - OVERHANG_LENGTH)

            if start_index < 1:
                start_index = 1

            if end_index < start_index:
                return {}

            # Storage for different overhang possibilities
            top_strand_overhangs = []
            bottom_strand_overhangs = []
            top_extended_sequences = []
            bottom_extended_sequences = []

            # Loop over possible positions to extract overhangs
            for i in range(start_index, end_index + 1):
                chunk_start = i - 1  # 1 nt before the overhang
                chunk_end = i + OVERHANG_LENGTH + 1  # 1 nt after the overhang
                full_top = sequence[chunk_start:chunk_end]

                # Skip invalid sequences (boundary issues)
                if len(full_top) != EXTENDED_LENGTH:
                    continue

                # Reverse complement for bottom strand
                full_bottom = utils.reverse_complement(full_top)

                # Slice out the overhang and extended sequences
                t_overhang = full_top[1:1 + OVERHANG_LENGTH]  # Middle 4 bases
                t_extended = full_top  # Full 6 bases

                b_overhang = full_bottom[1:1 + OVERHANG_LENGTH]
                b_extended = full_bottom

                # Append results
                top_strand_overhangs.append(t_overhang)
                bottom_strand_overhangs.append(b_overhang)
                top_extended_sequences.append(t_extended)
                bottom_extended_sequences.append(b_extended)

            return {
                "top_strand_overhangs": top_strand_overhangs,
                "bottom_strand_overhangs": bottom_strand_overhangs,
                "top_extended_sequences": top_extended_sequences,
                "bottom_extended_sequences": bottom_extended_sequences
            }

        # -- Main logic: iterate through mutation sites, codons, and alternative codons --
        for site_key, site_data in mutation_options.items():
            for codon in site_data["codons"]:
                position = codon["position"]
                for idx, alternative in enumerate(codon["alternative_codons"]):
                    try:
                        offset = alternative["mutations"].index(1)
                    except ValueError:
                        continue

                    overhangs = calculate_overhangs(sequence, position, offset, self.utils)
                    alternative["overhangs"] = overhangs

        return mutation_options


    def generate_mutation_sets(self, mutation_options: Dict) -> List[Dict]:
        """
        Generates all possible mutation sets by selecting exactly one alternative codon per restriction site.
        Each mutation includes a predicted 4-nt reassembly overhang.

        Args:
            seq (str): The full sequence being analyzed.
            mutation_options (Dict): Mutation site properties and possible mutations.
            overhangs_per_site (Dict): Precomputed overhangs for each mutation option.

        Returns:
            List[Dict]: List of mutation sets, where each set is a dictionary containing one mutation per site.
        """

        mutation_choices_per_site = []

        for site, site_data in mutation_options.items():
            site_mutation_choices = []

            for codon in site_data["codons"]:
                for alternative in codon["alternative_codons"]:
                    mutation_entry = {
                        "site": site,
                        "position": site_data["position"],
                        "original_sequence": codon["original_sequence"],
                        "alternative_sequence": alternative["seq"],
                        "mutated_base_index": (codon["position"] - site_data["position"] + site_data["frame"]) % 6,
                        "usage": alternative["usage"],
                        "overhangs": alternative["overhangs"],
                    }
                    site_mutation_choices.append(mutation_entry)

            mutation_choices_per_site.append(site_mutation_choices)

        # Generate mutation sets: One mutation per site
        all_mutation_sets = [
            {mutation["site"]: mutation for mutation in mutation_combination}
            for mutation_combination in product(*mutation_choices_per_site)
        ]

        print(f"Created {len(all_mutation_sets)} mutation sets")

        return all_mutation_sets


    def create_compatibility_matrices(self, mutation_sets: List[Dict]) -> List[np.ndarray]:
        """
        Generates compatibility matrices for a set of mutation sites and their respective overhangs.

        Parameters:
        -----------
        mutation_sets : List[Dict]
            A list of dictionaries, where each dictionary represents a set of mutation sites.
            Each mutation set contains position keys with corresponding overhang lists under:
            `mutation_set[pos]["overhangs"]["+"]`.

        Returns:
        --------
        List[np.ndarray]
            A list of compatibility matrices (numpy arrays), each indicating valid overhang combinations.
            Each matrix has a shape `(4, 4, ..., 4)` based on the number of mutation sites.

        Debugging:
        ----------
        - Randomly selects a few overhang combinations to report as compatible or not compatible.
        - If an overhang combination is incompatible, prints a reason why.
        """

        compatibility_matrices = []
        
        for mutation_set_idx, mutation_set in enumerate(tqdm(mutation_sets, desc="Processing Mutation Sets", unit="set")):
            
            position_keys = list(mutation_set.keys())  # Identify mutation site positions
            overhang_lists = [mutation_set[pos]["overhangs"]["top_strand_overhangs"] for pos in position_keys]

            # Define the shape of the compatibility matrix
            matrix_shape = (4,) * len(position_keys)
            compatibility_matrix = np.zeros(matrix_shape, dtype=int)
            
            all_ones = 0  # Counter for valid combinations
            tested_combos = []  # Store some tested combinations for debug output

            for combo in product(*overhang_lists):
                # Extract the position of each element within its respective group
                try:
                    positions = tuple(group.index(combo[i]) for i, group in enumerate(overhang_lists))
                except ValueError as e:
                    print(f"Error: {combo[i]} not found in its corresponding overhang group! Overhang lists: {overhang_lists}")
                    raise e

                
                # Check pairwise compatibility and exit early if any pair is incompatible
                compatible = True
                n_sites = len(combo)
                for i in range(n_sites):
                    for j in range(i + 1, n_sites):
                        idx1 = self.utils.seq_to_index(combo[i])
                        idx2 = self.utils.seq_to_index(combo[j])
                        if self.compatibility_table[idx1][idx2] != 1:
                            compatible = False
                            break
                    if not compatible:
                        break

                if compatible:
                    all_ones += 1
                    if all(0 <= pos < 4 for pos in positions):  # Ensure all indices are in range (0-3)
                        compatibility_matrix[positions] = 1
                    else:
                        print(f"Warning: Computed index {positions} is out of range for compatibility_matrix with shape {matrix_shape}.")

                else:
                    tested_combos.append((combo, self.analyze_incompatibility_reason(combo)))

            # Append the matrix for this mutation set to the list
            compatibility_matrices.append(compatibility_matrix)
            
        # Debugging: Count how many matrices are composed entirely of zeroes
        zero_matrices_count = sum(1 for matrix in compatibility_matrices if np.all(matrix == 0))
        logger.debug(f"\n ⚠️ {zero_matrices_count} out of {len(compatibility_matrices)} compatibility matrices are all zeroes. ⚠️")
        
        if zero_matrices_count == len(compatibility_matrices):
            logger.warning("❗ No valid combinatations found.")


        return compatibility_matrices


    def filter_compatible_mutations(self, mutation_sets: list, compatibility_matrices: list) -> list:
        """
        Filters mutation sets by removing those whose compatibility matrices contain only zeros.
        Processes in reverse order to avoid index shifting issues.

        Args:
            mutation_sets (list): List of mutation sets, each containing mutation sites and precomputed overhangs.
            compatibility_matrices (list): List of computed compatibility matrices corresponding to mutation_sets.

        Returns:
            list: Filtered list of mutation sets with at least one valid compatibility.
        """
        print(f"self.verbose: {self.verbose}")
        # Identify indices of mutation sets where the compatibility matrix is entirely zero
        indices_to_remove = [i for i, matrix in enumerate(compatibility_matrices) if np.all(matrix == 0)]

        # Process in reverse order to avoid index shifting while removing elements
        for idx in reversed(indices_to_remove):
            del mutation_sets[idx]
            del compatibility_matrices[idx]  # Keep compatibility_matrices in sync if needed later
        
        if self.verbose:
            if len(indices_to_remove) == 0:
                logger.info(f"All mutation sets have at least 1 valid overhang combinations.")
            else:    
                logger.info(f"Removed {len(indices_to_remove)} mutation sets that had no valid overhang combinations.")
                
        return mutation_sets



    def analyze_incompatibility_reason(self, combo):
        """
        Analyzes why a given combination of overhangs is incompatible.

        Parameters:
        -----------
        combo : tuple of str
            A combination of overhang sequences.

        Returns:
        --------
        str
            A reason explaining why the overhang combination is incompatible.
        """
        # Rule 1: Check if any overhang has more than 3 consecutive identical bases
        for seq in combo:
            for i in range(len(seq) - 3):
                if seq[i] == seq[i + 1] == seq[i + 2] == seq[i + 3]:
                    return f"Overhang {seq} has more than 3 consecutive identical bases."

        # Rule 2: Check if any two overhangs share the same 3 consecutive bases
        seen_triplets = set()
        for seq in combo:
            for i in range(len(seq) - 2):
                triplet = seq[i:i + 3]
                if triplet in seen_triplets:
                    return f"Multiple overhangs share the triplet '{triplet}'."
                seen_triplets.add(triplet)

        # Rule 3: Check if any overhang has 0% or 100% GC content
        for seq in combo:
            gc_content = sum(1 for base in seq if base in "GC") / len(seq) * 100
            if gc_content == 0:
                return f"Overhang {seq} has 0% GC content (all A/T)."
            elif gc_content == 100:
                return f"Overhang {seq} has 100% GC content (all G/C)."

        return "Unknown reason (should not happen)."

