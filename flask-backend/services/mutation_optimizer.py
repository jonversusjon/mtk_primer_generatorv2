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
from services.debug.debug_utils import MutationDebugger, debug_function, visualize_matrix, visualize_overhang_compatibility
from functools import wraps


class MutationOptimizer:
    """
    Optimizes mutations for Golden Gate assembly by balancing codon usage,
    sequence stability, and restriction site compatibility.
    """

    def __init__(self, verbose: bool = False, debug: bool = False):
        """
        Initialize the optimizer with a precomputed compatibility table.

        Args:
            verbose (bool): If True, provide more user-friendly logs in production.
            debug (bool): If True, enables comprehensive debugging output.
        """
        self.logger = logger.getChild("MutationOptimizer")
        self.utils = GoldenGateUtils()

        self.verbose = verbose
        self.debug = debug

        # Initialize debugger if debug mode is enabled
        if self.debug:
            self.debugger = MutationDebugger()
            self.logger.info(
                "🔍 Debug mode enabled with detailed logging and validation")

            # Replace our debug_function decorator with an actual implementation
            def apply_debug_wrapper(method):
                @wraps(method)
                def wrapped(self, *args, **kwargs):
                    # Log function start
                    func_name = method.__name__
                    self.debugger.log_function_start(func_name, args, kwargs)

                    # Run the function
                    try:
                        result = method(self, *args, **kwargs)

                        # Basic validation
                        if result is not None:
                            if isinstance(result, tuple):
                                for i, item in enumerate(result):
                                    self.debugger.validate(
                                        item is not None,
                                        f"Return value {i+1} is not None"
                                    )
                            else:
                                self.debugger.validate(
                                    result is not None, "Return value is not None")

                        # Log function end
                        self.debugger.log_function_end(func_name, result)
                        return result

                    except Exception as e:
                        self.debugger.log_error(
                            f"Exception in {func_name}: {str(e)}")
                        raise

                return wrapped

            # Apply debug wrapper to all methods decorated with @debug_function
            for attr_name in dir(self.__class__):
                if attr_name.startswith('__'):
                    continue

                attr = getattr(self.__class__, attr_name)
                if callable(attr) and hasattr(attr, '__wrapped__'):
                    setattr(self, attr_name, apply_debug_wrapper(
                        attr).__get__(self, self.__class__))

        # Load compatibility table
        self.compatibility_table_path = 'static/data/compatibility_table.bin'
        self.compatibility_table = self.utils._load_compatibility_table(
            self.compatibility_table_path)

        self.logger.debug(
            "MutationOptimizer initialized with compatibility table.")

        if self.verbose:
            self.logger.info("MutationOptimizer is running in verbose mode.")

        # Validate compatibility table if debugging
        if self.debug:
            self.debugger.validate(
                self.compatibility_table is not None,
                "Compatibility table loaded successfully"
            )
            self.debugger.validate(
                isinstance(self.compatibility_table, np.ndarray),
                "Compatibility table is a numpy array"
            )

    @debug_function
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
        if not self.debug:
            with debug_context("mutation_optimization"):
                logger.info(
                    "Step 1: Predicting BsmBI reassembly overhangs for each mutation option...")
                mutation_options_with_overhangs = self.add_predicted_overhangs(
                    sequence, mutation_options)

                logger.info("Step 2: Generating possible mutation sets...")
                mutation_sets = self.generate_mutation_sets(
                    mutation_options_with_overhangs)
                logger.debug(f"Generated {len(mutation_sets)} mutation sets.")

                logger.info("Step 3: Computing compatibility matrices...")
                compatibility_matrices = self.create_compatibility_matrices(
                    mutation_sets)

                logger.info(
                    "Step 4: Filtering mutation sets based on compatibility...")
                optimized_mutations = self.filter_compatible_mutations(
                    mutation_sets, compatibility_matrices)

                logger.debug(
                    f"\nFinal number of optimized mutation sets: {len(optimized_mutations)}")

                return optimized_mutations, compatibility_matrices
        else:
            # Validate inputs
            self.debugger.validate(
                isinstance(sequence, str) and len(sequence) > 0,
                "Input sequence is valid",
                sequence[:50] + "..." if len(sequence) > 50 else sequence
            )

            self.debugger.validate(
                isinstance(mutation_options, dict) and len(
                    mutation_options) > 0,
                "Mutation options are valid",
                {k: v for k, v in list(mutation_options.items())[:2]} if len(
                    mutation_options) > 2 else mutation_options
            )

            # STEP 1: Predict BsmBI reassembly overhangs
            self.debugger.log_step(
                "Predict Overhangs",
                "Predicting BsmBI reassembly overhangs for each mutation option"
            )
            mutation_options_with_overhangs = self.add_predicted_overhangs(
                sequence, mutation_options)

            # Validate step 1 output
            overhang_count = sum(
                1 for site in mutation_options_with_overhangs.values()
                for codon in site["codons"]
                for alt in codon["alternative_codons"]
                if "overhangs" in alt
            )

            total_alts = sum(
                len(codon["alternative_codons"])
                for site in mutation_options_with_overhangs.values()
                for codon in site["codons"]
            )

            self.debugger.validate(
                overhang_count > 0,
                f"Generated overhangs for {overhang_count}/{total_alts} alternative codons",
                {"success_rate": f"{overhang_count/total_alts:.2%}" if total_alts > 0 else "N/A"}
            )

            # STEP 2: Generate possible mutation sets
            self.debugger.log_step(
                "Generate Mutation Sets",
                "Creating all possible combinations of mutations"
            )
            mutation_sets = self.generate_mutation_sets(
                mutation_options_with_overhangs)

            # Validate step 2 output
            self.debugger.validate(
                isinstance(mutation_sets, list),
                f"Generated {len(mutation_sets)} mutation sets",
                {"sample": mutation_sets[0] if mutation_sets else None}
            )

            # STEP 3: Compute compatibility matrices
            self.debugger.log_step(
                "Compute Compatibility",
                "Creating compatibility matrices for mutation sets"
            )
            compatibility_matrices = self.create_compatibility_matrices(
                mutation_sets)

            # Validate step 3 output
            nonzero_matrices = sum(
                1 for m in compatibility_matrices if np.any(m != 0))
            self.debugger.validate(
                len(compatibility_matrices) == len(mutation_sets),
                f"Created {len(compatibility_matrices)} compatibility matrices, {nonzero_matrices} have valid combinations",
                {"nonzero_rate":
                    f"{nonzero_matrices/len(compatibility_matrices):.2%}" if compatibility_matrices else "N/A"}
            )

            # Sample a matrix for visualization if any exist
            if compatibility_matrices and nonzero_matrices > 0:
                # Find first non-zero matrix
                for i, matrix in enumerate(compatibility_matrices):
                    if np.any(matrix != 0):
                        self.debugger.log_step(
                            "Matrix Visualization",
                            f"Sample compatibility matrix (set #{i+1})",
                            visualize_matrix(matrix)
                        )
                        break

            # STEP 4: Filter compatible mutations
            self.debugger.log_step(
                "Filter Compatible Mutations",
                "Removing mutation sets with no valid compatibility"
            )
            optimized_mutations = self.filter_compatible_mutations(
                mutation_sets, compatibility_matrices)

            # Validate step 4 output
            self.debugger.validate(
                isinstance(optimized_mutations, list),
                f"Final result: {len(optimized_mutations)} optimized mutation sets",
                {"reduction": f"{(1 - len(optimized_mutations)/len(mutation_sets)):.2%}" if mutation_sets else "N/A"}
            )

            # Summarize all validations at the end
            if self.debug:
                self.debugger.summarize_validations()

            return optimized_mutations, compatibility_matrices

    @debug_function
    def add_predicted_overhangs(self, sequence: str, mutation_options: Dict) -> Dict:
        """
        Captures BsmBI reassembly overhangs for each proposed alternative codon,
        ensuring the mutation is positioned within the overhang.
        """
        if self.debug:
            # Input validation
            self.debugger.validate(
                isinstance(sequence, str) and len(sequence) > 0,
                "Input sequence is valid",
                sequence[:50] + "..." if len(sequence) > 50 else sequence
            )

            self.debugger.validate(
                isinstance(mutation_options, dict) and len(
                    mutation_options) > 0,
                "Mutation options are valid",
                {k: f"{len(v['codons'])} codons" for k,
                 v in mutation_options.items()}
            )

            # Track all sites processed
            site_count = 0
            codon_count = 0
            alt_codon_count = 0
            successful_overhangs = 0

        def calculate_overhangs(
            original_seq: str,
            mutated_seq: str,
            mutation_start: int,
            mutation_length: int,
            utils
        ) -> Dict:
            """
            Extract predicted overhangs from the mutated sequence, ensuring the mutation
            is placed within the 4-nt BsmBI overhang.

            Args:
                original_seq (str): The original sequence (unused here, but kept for consistency).
                mutated_seq (str): The sequence after applying the mutation.
                mutation_start (int): 0-based index of where the mutation begins in mutated_seq.
                mutation_length (int): Number of nucleotides changed (could be 1 for a point mutation).
                utils: Utility class with a reverse_complement function, etc.

            Returns:
                Dict: Contains lists of possible top/bottom overhangs and extended sequences.
                    Keys:
                        - top_strand_overhangs
                        - bottom_strand_overhangs
                        - top_extended_sequences
                        - bottom_extended_sequences
            """
            OVERHANG_LENGTH = 4
            EXTENDED_LENGTH = 6  # 4-nt overhang + 1 extra nt on each side

            top_strand_overhangs = []
            bottom_strand_overhangs = []
            top_extended_sequences = []
            bottom_extended_sequences = []
            positions_checked = []

            # First, identify which specific bases are actually mutated
            mutated_positions = []
            for i in range(mutation_length):
                pos = mutation_start + i
                if pos >= len(original_seq) or pos >= len(mutated_seq):
                    continue
                if original_seq[pos] != mutated_seq[pos]:
                    mutated_positions.append(pos)

            # Debug log
            if hasattr(self, 'debug') and self.debug:
                self.debugger.log_step(
                    "Mutation Analysis",
                    f"Identified exact mutated bases",
                    {
                        "mutation_region": f"[{mutation_start}, {mutation_start+mutation_length})",
                        "mutated_positions": mutated_positions,
                        "original_region": original_seq[mutation_start:mutation_start+mutation_length],
                        "mutated_region": mutated_seq[mutation_start:mutation_start+mutation_length]
                    }
                )

            # If no actual mutations were found, return empty result
            if not mutated_positions:
                if hasattr(self, 'debug') and self.debug:
                    self.debugger.log_warning(
                        "No mutated bases identified in the specified region")
                return {
                    "top_strand_overhangs": [],
                    "bottom_strand_overhangs": [],
                    "top_extended_sequences": [],
                    "bottom_extended_sequences": []
                }

            unique_overhangs = {}
            for mut_pos in mutated_positions:
                # The 4 possible positions are:
                # 1. Mutated base is at position 1 of overhang: i = mut_pos
                # 2. Mutated base is at position 2 of overhang: i = mut_pos - 1
                # 3. Mutated base is at position 3 of overhang: i = mut_pos - 2
                # 4. Mutated base is at position 4 of overhang: i = mut_pos - 3

                possible_starts = [mut_pos - 3,
                                   mut_pos - 2, mut_pos - 1, mut_pos]

                for i in possible_starts:
                    # Skip invalid positions
                    if i < 0 or i + OVERHANG_LENGTH > len(mutated_seq):
                        continue

                    positions_checked.append(i)

                    # Now define the 6-nt "extended" region: 1 base before + 4-nt overhang + 1 base after
                    chunk_start = i - 1
                    chunk_end = i + OVERHANG_LENGTH + 1  # exclusive

                    # Check boundaries
                    if chunk_start < 0 or chunk_end > len(mutated_seq):
                        continue

                    full_top = mutated_seq[chunk_start:chunk_end]
                    if len(full_top) != EXTENDED_LENGTH:
                        continue

                    # Reverse complement for the bottom strand
                    full_bottom = utils.reverse_complement(full_top)

                    # Overhang is the middle 4 bases, e.g. full_top[1:5]
                    t_overhang = full_top[1:1 + OVERHANG_LENGTH]
                    b_overhang = full_bottom[1:1 + OVERHANG_LENGTH]

                    # Make sure we don't add duplicates
                    if t_overhang in unique_overhangs:
                        if self.debug:
                            self.debugger.log_step(
                                "Duplicate Overhang",
                                f"Duplicate overhang at position {i}",
                                {
                                    "overhang": t_overhang,
                                    "previous_position": unique_overhangs[t_overhang]
                                },
                                level=logging.DEBUG
                            )
                        continue

                    # Store the position with this overhang
                    unique_overhangs[t_overhang] = i

                    # Append to lists
                    top_strand_overhangs.append(t_overhang)
                    bottom_strand_overhangs.append(b_overhang)
                    top_extended_sequences.append(full_top)
                    bottom_extended_sequences.append(full_bottom)

                    if self.debug:
                        self.debugger.log_step(
                            "Found Overhang",
                            f"Valid overhang at position {i}",
                            {
                                "mutated_base_position": mut_pos,
                                "position_in_overhang": mut_pos - i + 1,  # 1-based position in overhang
                                "top_overhang": t_overhang,
                                "bottom_overhang": b_overhang,
                                "top_extended": full_top,
                                "bottom_extended": full_bottom
                            },
                            level=logging.DEBUG
                        )

            # In debug mode, log the results
            if self.debug:
                self.debugger.log_step(
                    "Position Analysis",
                    f"Checked {len(positions_checked)} positions",
                    {
                        "positions_checked": sorted(positions_checked),
                        "unique_positions": sorted(list(set(positions_checked))),
                        "mutated_positions": mutated_positions
                    }
                )

                self.debugger.validate(
                    len(top_strand_overhangs) <= 4 * len(mutated_positions),
                    f"Found {len(top_strand_overhangs)} valid overhangs (expected ≤ {4 * len(mutated_positions)})",
                    {"overhangs": top_strand_overhangs}
                )

            return {
                "top_strand_overhangs": top_strand_overhangs,
                "bottom_strand_overhangs": bottom_strand_overhangs,
                "top_extended_sequences": top_extended_sequences,
                "bottom_extended_sequences": bottom_extended_sequences
            }

        # -- Main logic: iterate through mutation sites, codons, and alternative codons --
        for site_key, site_data in mutation_options.items():
            if self.debug:
                site_count += 1
                self.debugger.log_step(
                    "Process Site",
                    f"Processing site {site_key} at position {site_data['position']}",
                    {"frame": site_data["frame"]
                        if "frame" in site_data else "N/A"}
                )

            for codon in site_data["codons"]:
                position = codon["position"]
                original_codon = codon["original_sequence"]

                if self.debug:
                    codon_count += 1
                    self.debugger.log_step(
                        "Process Codon",
                        f"Codon at position {position}: {original_codon}",
                        {"alternatives": len(codon["alternative_codons"])},
                        level=logging.DEBUG
                    )

                for idx, alternative in enumerate(codon["alternative_codons"]):
                    if self.debug:
                        alt_codon_count += 1

                    # Create the mutated sequence
                    alt_codon = alternative["seq"]
                    mutation_length = len(alt_codon)

                    # Apply the mutation to the sequence
                    mutated_seq = sequence[:position] + alt_codon + \
                        sequence[position + len(original_codon):]

                    # Try to find where in the codon the mutation occurs
                    try:
                        # Get indices of nucleotides that differ between original and alternative
                        mutations = []
                        for i in range(min(len(original_codon), len(alt_codon))):
                            if i < len(original_codon) and i < len(alt_codon):
                                if original_codon[i] != alt_codon[i]:
                                    mutations.append(i)

                        # If no mutations found, continue to next alternative
                        if not mutations:
                            continue

                        # Generate overhangs ensuring mutation is in the overhang
                        overhangs = calculate_overhangs(
                            sequence,
                            mutated_seq,
                            position,
                            mutation_length,
                            self.utils
                        )

                        # Only store if we found valid overhangs
                        if overhangs and all(key in overhangs for key in ["top_strand_overhangs", "bottom_strand_overhangs"]):
                            if len(overhangs["top_strand_overhangs"]) > 0:
                                alternative["overhangs"] = overhangs
                                if self.debug:
                                    successful_overhangs += 1

                    except Exception as e:
                        error_msg = f"Error calculating overhangs for {site_key}, {alt_codon}: {str(e)}"
                        self.logger.error(error_msg)
                        if self.debug:
                            self.debugger.log_error(error_msg)
                        continue

        # Final validation in debug mode
        if self.debug:
            self.debugger.validate(
                successful_overhangs > 0,
                f"Generated overhangs for {successful_overhangs}/{alt_codon_count} alternative codons",
                {
                    "sites_processed": site_count,
                    "codons_processed": codon_count,
                    "alternatives_processed": alt_codon_count,
                    "success_rate": f"{successful_overhangs/alt_codon_count:.2%}" if alt_codon_count > 0 else "N/A"
                }
            )

        return mutation_options

    @debug_function
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
        if self.debug:
            # Input validation
            self.debugger.validate(
                isinstance(mutation_options, dict) and len(
                    mutation_options) > 0,
                "Mutation options are valid",
                {k: f"{len(v['codons'])} codons" for k,
                 v in mutation_options.items()}
            )

            self.debugger.log_step(
                "Prepare Mutation Choices",
                "Organizing mutation choices per site"
            )

        mutation_choices_per_site = []

        for site, site_data in mutation_options.items():
            site_mutation_choices = []

            if self.debug:
                self.debugger.log_step(
                    "Process Site",
                    f"Processing site {site} at position {site_data['position']}"
                )

            for codon in site_data["codons"]:
                for alternative in codon["alternative_codons"]:
                    # Skip alternatives without valid overhangs
                    if "overhangs" not in alternative:
                        continue

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

            if self.debug:
                self.debugger.validate(
                    len(site_mutation_choices) > 0,
                    f"Found {len(site_mutation_choices)} valid mutation options for site {site}"
                )

            mutation_choices_per_site.append(site_mutation_choices)

        # Generate mutation sets: One mutation per site
        all_mutation_sets = [
            {mutation["site"]: mutation for mutation in mutation_combination}
            for mutation_combination in product(*mutation_choices_per_site)
        ]

        if self.debug:
            expected_count = np.prod([len(choices)
                                     for choices in mutation_choices_per_site])
            self.debugger.validate(
                len(all_mutation_sets) == expected_count,
                f"Created {len(all_mutation_sets)} mutation sets (expected: {expected_count})",
                {"sites": len(mutation_choices_per_site)}
            )

            # Display a sample mutation set
            if all_mutation_sets:
                sample_set = all_mutation_sets[0]
                sample_info = {
                    "sites": list(sample_set.keys()),
                    "sample_mutation": {
                        "original": next(iter(sample_set.values()))["original_sequence"],
                        "alternative": next(iter(sample_set.values()))["alternative_sequence"],
                        "overhangs_count": len(next(iter(sample_set.values()))["overhangs"]["top_strand_overhangs"])
                    }
                }
                self.debugger.log_step(
                    "Sample Mutation Set",
                    f"Example of generated mutation set",
                    sample_info
                )

        if self.verbose or self.debug:
            print(f"Created {len(all_mutation_sets)} mutation sets")

        return all_mutation_sets

    @debug_function
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
        if self.debug:
            # Input validation
            self.debugger.validate(
                isinstance(mutation_sets, list) and len(mutation_sets) > 0,
                f"Processing {len(mutation_sets)} mutation sets",
                {"first_set_sites": list(
                    mutation_sets[0].keys()) if mutation_sets else None}
            )

            self.debugger.log_step(
                "Matrix Creation",
                "Creating compatibility matrices for each mutation set"
            )

            # Track statistics
            total_combinations = 0
            compatible_combinations = 0
            zero_matrices = 0

        compatibility_matrices = []

        for mutation_set_idx, mutation_set in enumerate(tqdm(mutation_sets, desc="Processing Mutation Sets", unit="set")):
            if self.debug and mutation_set_idx < 3:  # Only log detailed info for first 3 sets
                self.debugger.log_step(
                    "Process Set",
                    f"Creating matrix for mutation set {mutation_set_idx+1}",
                    {"sites": list(mutation_set.keys())}
                )

            # Identify mutation site positions
            position_keys = list(mutation_set.keys())
            overhang_lists = [mutation_set[pos]["overhangs"]
                              ["top_strand_overhangs"] for pos in position_keys]

            # Define the shape of the compatibility matrix
            matrix_shape = tuple(len(overhangs)
                                 for overhangs in overhang_lists)

            if self.debug and mutation_set_idx < 3:  # Only log detailed info for first 3 sets
                self.debugger.log_step(
                    "Matrix Shape",
                    f"Creating matrix with shape {matrix_shape}",
                    {"overhangs_per_site": [len(o) for o in overhang_lists]}
                )

            compatibility_matrix = np.zeros(matrix_shape, dtype=int)

            all_ones = 0  # Counter for valid combinations
            tested_combos = []  # Store some tested combinations for debug output

            # Generate all possible combinations of overhangs
            for combo_indices in product(*[range(len(overhang_list)) for overhang_list in overhang_lists]):
                combo = tuple(overhang_lists[i][idx]
                              for i, idx in enumerate(combo_indices))

                if self.debug:
                    total_combinations += 1

                # Check pairwise compatibility
                compatible = True
                n_sites = len(combo)

                for i in range(n_sites):
                    for j in range(i + 1, n_sites):
                        idx1 = self.utils.seq_to_index(combo[i])
                        idx2 = self.utils.seq_to_index(combo[j])
                        if self.compatibility_table[idx1][idx2] != 1:
                            compatible = False
                            if self.debug and mutation_set_idx < 3 and len(tested_combos) < 5:
                                reason = self.analyze_incompatibility_reason(
                                    combo)
                                tested_combos.append((combo, reason))
                            break
                    if not compatible:
                        break

                if compatible:
                    all_ones += 1
                    if self.debug:
                        compatible_combinations += 1
                    compatibility_matrix[combo_indices] = 1

            # Check if the matrix has any valid combinations
            if np.all(compatibility_matrix == 0):
                if self.debug:
                    zero_matrices += 1

            # Append the matrix for this mutation set to the list
            compatibility_matrices.append(compatibility_matrix)

            # Log some sample test combinations for the first few matrices
            if self.debug and mutation_set_idx < 3 and tested_combos:
                for combo, reason in tested_combos:
                    self.debugger.log_step(
                        "Incompatible Combination",
                        f"Sample incompatible overhang set",
                        {"overhangs": combo, "reason": reason},
                        level=logging.DEBUG
                    )

                # Log matrix stats
                self.debugger.log_step(
                    "Matrix Stats",
                    f"Matrix {mutation_set_idx+1} has {all_ones} valid combinations out of {np.prod(matrix_shape)}",
                    {"valid_percentage": f"{all_ones/np.prod(matrix_shape):.2%}" if np.prod(
                        matrix_shape) > 0 else "N/A"}
                )

        # Debugging: Count how many matrices are composed entirely of zeroes
        zero_matrices_count = sum(
            1 for matrix in compatibility_matrices if np.all(matrix == 0))
        logger.debug(
            f"\n ⚠️ {zero_matrices_count} out of {len(compatibility_matrices)} compatibility matrices are all zeroes. ⚠️")

        if self.debug:
            self.debugger.validate(
                compatible_combinations > 0,
                f"Found {compatible_combinations} compatible combinations out of {total_combinations} total",
                {
                    "matrices_created": len(compatibility_matrices),
                    "zero_matrices": zero_matrices,
                    "success_rate": f"{compatible_combinations/total_combinations:.2%}" if total_combinations > 0 else "N/A"
                }
            )

        if zero_matrices_count == len(compatibility_matrices):
            logger.warning("❗ No valid combinatations found.")
            if self.debug:
                self.debugger.log_warning(
                    "No valid combinations found in any matrix!")

        return compatibility_matrices

    @debug_function
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
        if self.debug:
            # Input validation
            self.debugger.validate(
                isinstance(mutation_sets, list) and len(mutation_sets) > 0,
                f"Input: {len(mutation_sets)} mutation sets"
            )
            self.debugger.validate(
                isinstance(compatibility_matrices, list) and len(
                    compatibility_matrices) == len(mutation_sets),
                f"Input: {len(compatibility_matrices)} compatibility matrices"
            )

            self.debugger.log_step(
                "Filter Sets",
                "Identifying mutation sets with no valid combinations"
            )

        # Identify indices of mutation sets where the compatibility matrix is entirely zero
        indices_to_remove = [i for i, matrix in enumerate(
            compatibility_matrices) if np.all(matrix == 0)]

        if self.debug:
            self.debugger.log_step(
                "Remove Sets",
                f"Will remove {len(indices_to_remove)} sets with all-zero matrices",
                {"percentage": f"{len(indices_to_remove)/len(mutation_sets):.2%}" if mutation_sets else "N/A"}
            )

        # Process in reverse order to avoid index shifting while removing elements
        for idx in reversed(indices_to_remove):
            del mutation_sets[idx]
            # Keep compatibility_matrices in sync if needed later
            del compatibility_matrices[idx]

        if self.verbose:
            if len(indices_to_remove) == 0:
                logger.info(
                    f"All mutation sets have at least 1 valid overhang combinations.")
            else:
                logger.info(
                    f"Removed {len(indices_to_remove)} mutation sets that had no valid overhang combinations.")

        if self.debug:
            self.debugger.validate(
                len(mutation_sets) + len(indices_to_remove) == len(
                    compatibility_matrices) + len(indices_to_remove),
                f"Successfully filtered: {len(mutation_sets)} sets remaining"
            )

            # Sample a remaining set if available
            if mutation_sets:
                sample_set = mutation_sets[0]
                sample_info = {
                    "sites": list(sample_set.keys()),
                    "mutations": {key: value["alternative_sequence"] for key, value in sample_set.items()}
                }
                self.debugger.log_step(
                    "Sample Result",
                    f"Example of remaining valid mutation set",
                    sample_info
                )

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
            gc_content = sum(
                1 for base in seq if base in "GC") / len(seq) * 100
            if gc_content == 0:
                return f"Overhang {seq} has 0% GC content (all A/T)."
            elif gc_content == 100:
                return f"Overhang {seq} has 100% GC content (all G/C)."

        return "Unknown reason (should not happen)."
