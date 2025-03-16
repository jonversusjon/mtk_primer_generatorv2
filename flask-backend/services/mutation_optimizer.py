from typing import Dict, List, Tuple, Any
import numpy as np
from itertools import product
from tqdm import tqdm

from services.debug import DebugMixin, MutationDebugger, debug_context
from config.logging_config import logger
from .utils import GoldenGateUtils


class MutationOptimizer(DebugMixin):
    """
    MutationOptimizer Module

    This module optimizes mutation combinations generated by the MutationAnalyzer module
    for Golden Gate assembly compatibility. It evaluates alternative codon mutations based on
    sticky-end compatibility, selecting sets of mutations that satisfy the Golden Gate cloning criteria.
    """

    def __init__(self, verbose: bool = False, debug: bool = False):
        self.logger = logger.getChild("MutationOptimizer")
        self.utils = GoldenGateUtils()
        self.verbose = verbose
        self.debug = debug
        self.debugger = MutationDebugger() if self.debug else None

        if self.debug:
            self.logger.info(
                "🔍 Debug mode enabled with detailed logging and validation")

        self.compatibility_table = self.utils.load_compatibility_table(
            'static/data/compatibility_table.bin')
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
    def optimize_mutations(self, mutation_options: Dict) -> Tuple[List[Dict], Any]:
        with debug_context("mutation_optimization"):
            self.log_step("Start Optimization",
                          "Beginning mutation optimization process.")
            # 1. Generate all possible mutation sets
            self.log_step("Generate Mutation Sets",
                          "Creating all possible mutation combinations.")
            mutation_sets = self.generate_mutation_sets(mutation_options)
            self.log_step("Generate Mutation Sets",
                          f"Total mutation sets generated: {len(mutation_sets)}")

            # 2. Create compatibility matrices for each mutation set
            self.log_step("Compute Compatibility",
                          "Creating compatibility matrices for each mutation set.")
            matrices = self.create_compatibility_matrices(mutation_sets)
            self.log_step("Compute Compatibility",
                          f"Generated {len(matrices)} compatibility matrices.")

            # 3. Filter out mutation sets without any compatible combination
            self.log_step("Filter Compatible Mutations",
                          "Filtering out mutation sets with no compatible combinations.")
            optimized = self.filter_compatible_mutations(
                mutation_sets, matrices)
            self.log_step("Filter Compatible Mutations",
                          f"Remaining mutation sets after filtering: {len(optimized)}")

            # Pydantic validation of optimized mutation sets
            from models.mutations import MutationSet
            validated_sets = []
            self.log_step("Validate Mutation Sets",
                          "Starting Pydantic validation of each mutation set.")
            for idx, m_set in enumerate(optimized):
                try:
                    validated_set = MutationSet.model_validate(
                        m_set)  # Pydantic v2 method
                    validated_sets.append(validated_set.model_dump())
                    self.log_step("Validate Mutation Set",
                                  f"Mutation set {idx} validated successfully.")
                except Exception as e:
                    self.logger.error(
                        f"Validation error in mutation set at index {idx}: {e}")
                    raise e
            self.log_step("Optimization Complete",
                          "Mutation optimization process completed successfully.")
            return validated_sets, matrices

    @DebugMixin.debug_wrapper
    def generate_mutation_sets(self, mutation_options: dict) -> list:
        self.validate(mutation_options and isinstance(mutation_options, dict),
                      f"Received {len(mutation_options)} mutation option sites")
        mutation_choices = []
        self.log_step("Generate Mutation Sets",
                      "Starting generation of mutation entries per site.")
        for site_key, site_data in mutation_options.items():
            self.log_step(
                "Process Site", f"Processing site {site_key} at position {site_data['position']}")
            mutation_entries = []
            # Iterate through each codon in the site.
            for codon in site_data.get("codons", []):
                for alt in codon.get("alternative_codons", []):
                    overhang_options = []
                    for sticky_key, sticky_val in alt.get("sticky_ends", {}).items():
                        top_list = sticky_val.get("top_strand", [])
                        bottom_list = sticky_val.get("bottom_strand", [])
                        for idx in range(len(top_list)):
                            top_option = top_list[idx]
                            bottom_option = bottom_list[idx] if idx < len(
                                bottom_list) else {}
                            overhang_options.append({
                                "top_overhang": top_option,
                                "bottom_overhang": bottom_option,
                                "overhang_start_index": top_option.get("overhang_start_index")
                            })
                    mutation_entry = {
                        "site": site_key,
                        "position": site_data.get("position"),
                        # Ensure upstream renaming
                        "codon_sequence": codon.get("codon_sequence"),
                        "alternative_codon_sequence": alt.get("seq"),
                        "mutated_base_index": alt.get("changes_in_site", [None])[0],
                        "overhangs": {"overhang_options": overhang_options},
                        "mutated_context": alt.get("mutated_context"),
                        "mutation_positions_in_context": alt.get("mutation_positions_in_context")
                    }
                    mutation_entries.append(mutation_entry)
            self.validate(
                mutation_entries, f"Found {len(mutation_entries)} mutation entries for {site_key}")
            mutation_choices.append(mutation_entries)
            self.log_step(
                "Process Site", f"Site {site_key} produced {len(mutation_entries)} mutation entries.")
        all_sets = [list(combo) for combo in product(*mutation_choices)]
        self.log_step("Cartesian Product",
                      f"Total Cartesian product sets: {len(all_sets)}")
        if self.debug:
            expected = np.prod([len(choices) for choices in mutation_choices])
            self.validate(len(all_sets) == expected,
                          f"Created {len(all_sets)} mutation sets (expected: {expected})",
                          {"sites": len(mutation_choices)})
        return all_sets

    @DebugMixin.debug_wrapper
    def create_compatibility_matrices(self, mutation_sets: List[Dict]) -> List[np.ndarray]:
        compatibility_matrices = []
        self.log_step("Compatibility Matrices",
                      f"Starting creation for {len(mutation_sets)} mutation sets.")
        for idx, mutation_set in enumerate(tqdm(mutation_sets, desc="Processing Mutation Sets", unit="set")):
            overhang_lists = [
                [opt["overhangs"]["overhang_options"][i]
                    for i in range(len(opt["overhangs"]["overhang_options"]))]
                for opt in mutation_set
            ]
            shape = tuple(len(lst) for lst in overhang_lists)
            matrix = np.zeros(shape, dtype=int)
            for combo_indices in product(*[range(len(lst)) for lst in overhang_lists]):
                combo = tuple(overhang_lists[i][idx]
                              for i, idx in enumerate(combo_indices))
                n = len(combo)
                if all(
                    self.compatibility_table[self.utils.seq_to_index(
                        combo[i]["top_overhang"]["seq"])]
                    [self.utils.seq_to_index(
                        combo[j]["top_overhang"]["seq"])] == 1
                    for i in range(n) for j in range(i + 1, n)
                ):
                    matrix[combo_indices] = 1
            compatibility_matrices.append(matrix)
            self.log_step("Compatibility Matrix",
                          f"Matrix for mutation set {idx} created with shape {matrix.shape}.")
        self.log_step("Compatibility Matrices",
                      "Completed creation of compatibility matrices.")
        return compatibility_matrices

    @DebugMixin.debug_wrapper
    def filter_compatible_mutations(self, mutation_sets: list, compatibility_matrices: list) -> list:
        self.log_step("Filter Mutations",
                      "Starting filtering of mutation sets based on compatibility matrices.")
        indices = [i for i, mat in enumerate(
            compatibility_matrices) if np.all(mat == 0)]
        if indices:
            self.log_step(
                "Filter Mutations", f"Found {len(indices)} mutation sets with zero compatibility. Removing them.")
        for i in sorted(indices, reverse=True):
            del mutation_sets[i]
            del compatibility_matrices[i]
        self.log_step("Filter Mutations",
                      f"Remaining mutation sets after filtering: {len(mutation_sets)}")
        return mutation_sets
