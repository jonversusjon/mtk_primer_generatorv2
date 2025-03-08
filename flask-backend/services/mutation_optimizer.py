from typing import Dict, List
from Bio.Seq import Seq
import numpy as np
from itertools import product
from tqdm import tqdm
import logging
from functools import wraps
import random

from .utils import GoldenGateUtils
from config.logging_config import logger
from services.base import debug_context
from services.debug.debug_mixin import DebugMixin
from services.debug.debug_utils import MutationDebugger, visualize_matrix, visualize_overhang_compatibility


class MutationOptimizer(DebugMixin):
    """
    Optimizes mutations for Golden Gate assembly by balancing codon usage,
    sequence stability, and restriction site compatibility.
    """

    def __init__(self, verbose: bool = False, debug: bool = False):
        """
        Initialize the optimizer with a precomputed compatibility table.
        """
        self.logger = logger.getChild("MutationOptimizer")
        self.utils = GoldenGateUtils()

        self.verbose = verbose
        self.debug = debug
        self.debugger = MutationDebugger() if self.debug else None

        if self.debug:
            self.logger.info(
                "🔍 Debug mode enabled with detailed logging and validation")

        # Load compatibility table
        self.compatibility_table_path = 'static/data/compatibility_table.bin'
        self.compatibility_table = self.utils._load_compatibility_table(
            self.compatibility_table_path)
        self.logger.debug(
            "MutationOptimizer initialized with compatibility table.")

        if self.verbose:
            self.logger.info("MutationOptimizer is running in verbose mode.")

        if self.debug:
            self.validate(self.compatibility_table is not None,
                          "Compatibility table loaded successfully")
            self.validate(isinstance(self.compatibility_table,
                          np.ndarray), "Compatibility table is a numpy array")

    @DebugMixin.debug_wrapper
    def optimize_mutations(self, mutation_options: Dict) -> List[Dict]:
        """
        Find the optimal combination of mutations that are compatible according to the BsmBI sticky-end strategy.
        """
        if not self.debug:
            with debug_context("mutation_optimization"):
                logger.info("Step 1: Generating possible mutation sets...")
                mutation_sets = self.generate_mutation_sets(mutation_options)
                logger.debug(f"Generated {len(mutation_sets)} mutation sets.")

                logger.info("Step 2: Computing compatibility matrices...")
                compatibility_matrices = self.create_compatibility_matrices(
                    mutation_sets)

                logger.info(
                    "Step 3: Filtering mutation sets based on compatibility...")
                optimized_mutations = self.filter_compatible_mutations(
                    mutation_sets, compatibility_matrices)

                logger.debug(
                    f"Final number of optimized mutation sets: {len(optimized_mutations)}")
                return optimized_mutations, compatibility_matrices
        else:
            self.validate(isinstance(mutation_options, dict) and len(mutation_options) > 0,
                          "Mutation options are valid",
                          {k: v for k, v in list(mutation_options.items())[:2]})
            self.log_step("Generate Mutation Sets",
                          "Creating all possible combinations of mutations")
            mutation_sets = self.generate_mutation_sets(mutation_options)
            self.validate(isinstance(mutation_sets, list),
                          f"Generated {len(mutation_sets)} mutation sets",
                          {"sample": mutation_sets[0] if mutation_sets else None})
            self.log_step("Compute Compatibility",
                          "Creating compatibility matrices for mutation sets")
            compatibility_matrices = self.create_compatibility_matrices(
                mutation_sets)
            nonzero_matrices = sum(
                1 for m in compatibility_matrices if np.any(m != 0))
            self.validate(len(compatibility_matrices) == len(mutation_sets),
                          f"Created {len(compatibility_matrices)} compatibility matrices, {nonzero_matrices} have valid combinations",
                          {"nonzero_rate": f"{nonzero_matrices/len(compatibility_matrices):.2%}" if compatibility_matrices else "N/A"})
            if compatibility_matrices and nonzero_matrices > 0:
                for i, matrix in enumerate(compatibility_matrices):
                    if np.any(matrix != 0):
                        self.log_step(
                            "Matrix Visualization", f"Sample compatibility matrix (set #{i+1})", visualize_matrix(matrix))
                        break
            self.log_step("Filter Compatible Mutations",
                          "Removing mutation sets with no valid compatibility")
            optimized_mutations = self.filter_compatible_mutations(
                mutation_sets, compatibility_matrices)
            self.validate(isinstance(optimized_mutations, list),
                          f"Final result: {len(optimized_mutations)} optimized mutation sets",
                          {"reduction": f"{(1 - len(optimized_mutations)/len(mutation_sets)):.2%}" if mutation_sets else "N/A"})
            if self.debug:
                self.debugger.summarize_validations()
            return optimized_mutations, compatibility_matrices

    @DebugMixin.debug_wrapper
    def generate_mutation_sets(self, mutation_options: Dict) -> List[Dict]:
        """
        Generates all possible mutation sets by selecting exactly one alternative codon per restriction site.
        """
        if self.debug:
            self.validate(isinstance(mutation_options, dict) and len(mutation_options) > 0,
                          "Mutation options are valid",
                          {k: f"{len(v['codons'])} codons" for k, v in mutation_options.items()})
            self.log_step("Prepare Mutation Choices",
                          "Organizing mutation choices per site")

        mutation_choices_per_site = []

        for site, site_data in mutation_options.items():
            site_mutation_choices = []
            if self.debug:
                self.log_step(
                    "Process Site", f"Processing site {site} at position {site_data['position']}")
            for codon in site_data["codons"]:
                for alternative in codon["alternative_codons"]:
                    if "sticky_ends" not in alternative:
                        continue

                    # This unified list will store all overhang options for this alternative
                    overhang_options = []

                    context = site_data.get("context", "")
                    context_mutated_indices = site_data.get(
                        "context_mutated_indices", [])
                    if context and len(context) > 0:
                        codon_pos = codon["position"]
                        # Assume the first index in context_mutated_indices is the site start in context.
                        site_start_in_context = context_mutated_indices[0]
                        codon_pos_in_context = site_start_in_context + \
                            (codon_pos - site_data["position"])
                        # Update context with the alternative sequence (the mutation)
                        mutated_context = list(context)
                        for i in range(len(alternative["seq"])):
                            if codon_pos_in_context + i < len(mutated_context):
                                mutated_context[codon_pos_in_context +
                                                i] = alternative["seq"][i]
                        mutated_context = ''.join(mutated_context)

                        # Loop over each sticky end specification (each pos_key)
                        for pos_key in alternative["sticky_ends"]:
                            try:
                                pos = int(pos_key.split('_')[1])
                            except (IndexError, ValueError):
                                continue
                            pos_in_context = codon_pos_in_context + pos

                            # If top strand information is provided, iterate over its options
                            if "top_strand" in alternative["sticky_ends"][pos_key]:
                                for i, top_overhang in enumerate(alternative["sticky_ends"][pos_key]["top_strand"]):
                                    # Determine extended region boundaries based on the option index
                                    if i == 0:
                                        start, end = pos_in_context - 4, pos_in_context + 2
                                    elif i == 1:
                                        start, end = pos_in_context - 3, pos_in_context + 3
                                    elif i == 2:
                                        start, end = pos_in_context - 2, pos_in_context + 4
                                    elif i == 3:
                                        start, end = pos_in_context - 1, pos_in_context + 5
                                    else:
                                        continue

                                    start = max(0, start)
                                    end = min(len(mutated_context), end)
                                    if end - start >= 6:
                                        top_extended = mutated_context[start:end]
                                    else:
                                        top_extended = mutated_context[start:end]
                                        if len(top_extended) < 6:
                                            padding_needed = 6 - \
                                                len(top_extended)
                                            if start > 0:
                                                prefix = mutated_context[max(
                                                    0, start - padding_needed):start]
                                                top_extended = prefix + top_extended
                                            if len(top_extended) < 6 and end < len(mutated_context):
                                                suffix = mutated_context[end:min(
                                                    len(mutated_context), end + padding_needed)]
                                                top_extended = top_extended + suffix
                                        top_extended = top_extended[:6]

                                    # Get the corresponding bottom strand overhang (if available)
                                    bottom_overhang = None
                                    if "bottom_strand" in alternative["sticky_ends"][pos_key]:
                                        try:
                                            bottom_overhang = alternative["sticky_ends"][pos_key]["bottom_strand"][i]
                                        except IndexError:
                                            bottom_overhang = None

                                    # Compute bottom extended as the reverse complement of the top extended sequence
                                    bottom_extended = self.utils.reverse_complement(
                                        top_extended)

                                    # Append the unified overhang option
                                    overhang_options.append({
                                        "top_overhang": top_overhang,
                                        "bottom_overhang": bottom_overhang,
                                        "top_extended": top_extended,
                                        "bottom_extended": bottom_extended
                                    })
                    else:
                        # If no context is provided, use a default extended sequence method.
                        for pos_key in alternative["sticky_ends"]:
                            if "top_strand" in alternative["sticky_ends"][pos_key]:
                                for i, top_overhang in enumerate(alternative["sticky_ends"][pos_key]["top_strand"]):
                                    top_extended = "A" + top_overhang + "T"
                                    bottom_overhang = None
                                    if "bottom_strand" in alternative["sticky_ends"][pos_key]:
                                        try:
                                            bottom_overhang = alternative["sticky_ends"][pos_key]["bottom_strand"][i]
                                        except IndexError:
                                            bottom_overhang = None
                                    bottom_extended = "A" + bottom_overhang + "T" if bottom_overhang else None
                                    overhang_options.append({
                                        "top_overhang": top_overhang,
                                        "bottom_overhang": bottom_overhang,
                                        "top_extended": top_extended,
                                        "bottom_extended": bottom_extended
                                    })

                    # Now, instead of storing four separate lists, we store a single unified structure:
                    overhangs = {"overhang_options": overhang_options}

                    # Use this 'overhangs' dictionary in your mutation_entry:
                    original_seq = codon.get("original_sequence", codon.get(
                        "codon_seq", alternative.get("seq", "")))
                    mutation_entry = {
                        "site": site,
                        "position": site_data["position"],
                        "original_sequence": original_seq,
                        "alternative_sequence": alternative["seq"],
                        "mutated_base_index": (codon["position"] - site_data["position"] + site_data["frame"]) % 6,
                        "usage": alternative["usage"],
                        "overhangs": overhangs
                    }
                    site_mutation_choices.append(mutation_entry)

            if self.debug:
                self.validate(len(site_mutation_choices) > 0,
                              f"Found {len(site_mutation_choices)} valid mutation options for site {site}")
            mutation_choices_per_site.append(site_mutation_choices)

        all_mutation_sets = [
            {mutation["site"]: mutation for mutation in mutation_combination}
            for mutation_combination in product(*mutation_choices_per_site)
        ]
        if self.debug:
            expected_count = np.prod([len(choices)
                                      for choices in mutation_choices_per_site])
            self.validate(len(all_mutation_sets) == expected_count,
                          f"Created {len(all_mutation_sets)} mutation sets (expected: {expected_count})",
                          {"sites": len(mutation_choices_per_site)})
            if all_mutation_sets:
                sample_set = all_mutation_sets[0]
                sample_info = {
                    "sites": list(sample_set.keys()),
                    "sample_mutation": {
                        "original": next(iter(sample_set.values()))["original_sequence"],
                        "alternative": next(iter(sample_set.values()))["alternative_sequence"],
                        "overhangs_count": len(next(iter(sample_set.values()))["overhangs"]["overhang_options"])
                    }
                }
                self.debugger.log_step(
                    "Sample Mutation Set", "Example of generated mutation set", sample_info)
        if self.verbose or self.debug:
            print(f"Created {len(all_mutation_sets)} mutation sets")
        return all_mutation_sets

    @DebugMixin.debug_wrapper
    def create_compatibility_matrices(self, mutation_sets: List[Dict]) -> List[np.ndarray]:
        if self.debug:
            self.validate(isinstance(mutation_sets, list) and len(mutation_sets) > 0,
                          f"Processing {len(mutation_sets)} mutation sets",
                          {"first_set_sites": list(mutation_sets[0].keys()) if mutation_sets else None})
            self.log_step(
                "Matrix Creation", "Creating compatibility matrices for each mutation set")
            total_combinations = 0
            compatible_combinations = 0
            zero_matrices = 0
        compatibility_matrices = []
        for mutation_set_idx, mutation_set in enumerate(tqdm(mutation_sets, desc="Processing Mutation Sets", unit="set")):
            if self.debug and mutation_set_idx < 3:
                self.log_step("Process Set", f"Creating matrix for mutation set {mutation_set_idx+1}",
                              {"sites": list(mutation_set.keys())})
            position_keys = list(mutation_set.keys())
            overhang_lists = [
                [option["top_overhang"] for option in mutation_set[pos]["overhangs"]["overhang_options"]]
                for pos in position_keys
            ]

            matrix_shape = tuple(len(overhangs)
                                 for overhangs in overhang_lists)
            if self.debug and mutation_set_idx < 3:
                self.log_step("Matrix Shape", f"Creating matrix with shape {matrix_shape}",
                              {"overhangs_per_site": [len(o) for o in overhang_lists]})
            compatibility_matrix = np.zeros(matrix_shape, dtype=int)
            all_ones = 0
            tested_combos = []
            for combo_indices in product(*[range(len(overhang_list)) for overhang_list in overhang_lists]):
                combo = tuple(overhang_lists[i][idx]
                              for i, idx in enumerate(combo_indices))
                if self.debug:
                    total_combinations += 1
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
            if np.all(compatibility_matrix == 0):
                if self.debug:
                    zero_matrices += 1
            compatibility_matrices.append(compatibility_matrix)
            if self.debug and mutation_set_idx < 3 and tested_combos:
                for combo, reason in tested_combos:
                    self.debugger.log_step("Incompatible Combination", "Sample incompatible overhang set",
                                           {"overhangs": combo, "reason": reason}, level=logging.DEBUG)
                self.debugger.log_step("Matrix Stats", f"Matrix {mutation_set_idx+1} has {all_ones} valid combinations out of {np.prod(matrix_shape)}",
                                       {"valid_percentage": f"{all_ones/np.prod(matrix_shape):.2%}" if np.prod(matrix_shape) > 0 else "N/A"})
        zero_matrices_count = sum(
            1 for matrix in compatibility_matrices if np.all(matrix == 0))
        logger.debug(
            f"\n ⚠️ {zero_matrices_count} out of {len(compatibility_matrices)} compatibility matrices are all zeroes. ⚠️")
        if self.debug:
            self.validate(compatible_combinations > 0,
                          f"Found {compatible_combinations} compatible combinations out of {total_combinations} total",
                          {"matrices_created": len(compatibility_matrices),
                           "zero_matrices": zero_matrices,
                           "success_rate": f"{compatible_combinations/total_combinations:.2%}" if total_combinations > 0 else "N/A"})
        if zero_matrices_count == len(compatibility_matrices):
            logger.warning("❗ No valid combinatations found.")
            if self.debug:
                self.debugger.log_warning(
                    "No valid combinations found in any matrix!")
        return compatibility_matrices

    @DebugMixin.debug_wrapper
    def filter_compatible_mutations(self, mutation_sets: list, compatibility_matrices: list) -> list:
        if self.debug:
            self.validate(isinstance(mutation_sets, list) and len(mutation_sets) > 0,
                          f"Input: {len(mutation_sets)} mutation sets")
            self.validate(isinstance(compatibility_matrices, list) and len(compatibility_matrices) == len(mutation_sets),
                          f"Input: {len(compatibility_matrices)} compatibility matrices")
            self.log_step(
                "Filter Sets", "Identifying mutation sets with no valid combinations")
        indices_to_remove = [i for i, matrix in enumerate(
            compatibility_matrices) if np.all(matrix == 0)]
        if self.debug:
            self.log_step("Remove Sets", f"Will remove {len(indices_to_remove)} sets with all-zero matrices",
                          {"percentage": f"{len(indices_to_remove)/len(mutation_sets):.2%}" if mutation_sets else "N/A"})
        for idx in reversed(indices_to_remove):
            del mutation_sets[idx]
            del compatibility_matrices[idx]
        if self.verbose:
            if len(indices_to_remove) == 0:
                logger.info(
                    "All mutation sets have at least 1 valid overhang combinations.")
            else:
                logger.info(
                    f"Removed {len(indices_to_remove)} mutation sets that had no valid overhang combinations.")
        if self.debug:
            self.validate(len(mutation_sets) + len(indices_to_remove) == len(compatibility_matrices) + len(indices_to_remove),
                          f"Successfully filtered: {len(mutation_sets)} sets remaining")
            if mutation_sets:
                sample_set = mutation_sets[0]
                sample_info = {
                    "sites": list(sample_set.keys()),
                    "mutations": {key: value["alternative_sequence"] for key, value in sample_set.items()}
                }
                self.debugger.log_step(
                    "Sample Result", "Example of remaining valid mutation set", sample_info)
        return mutation_sets

    def analyze_incompatibility_reason(self, combo):
        """
        Analyzes why a given combination of overhangs is incompatible.
        """
        for seq in combo:
            for i in range(len(seq) - 3):
                if seq[i] == seq[i + 1] == seq[i + 2] == seq[i + 3]:
                    return f"Overhang {seq} has more than 3 consecutive identical bases."
        seen_triplets = set()
        for seq in combo:
            for i in range(len(seq) - 2):
                triplet = seq[i:i + 3]
                if triplet in seen_triplets:
                    return f"Multiple overhangs share the triplet '{triplet}'."
                seen_triplets.add(triplet)
        for seq in combo:
            gc_content = sum(
                1 for base in seq if base in "GC") / len(seq) * 100
            if gc_content == 0:
                return f"Overhang {seq} has 0% GC content (all A/T)."
            elif gc_content == 100:
                return f"Overhang {seq} has 100% GC content (all G/C)."
        return "Unknown reason (should not happen)."
