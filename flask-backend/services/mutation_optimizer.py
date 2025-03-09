from typing import Dict, List
from Bio.Seq import Seq
import numpy as np
from itertools import product
from tqdm import tqdm
import logging
from services.base import debug_context
from services.debug.debug_mixin import DebugMixin
from services.debug.debug_utils import MutationDebugger, visualize_matrix, visualize_overhang_compatibility
from config.logging_config import logger
from .utils import GoldenGateUtils


class MutationOptimizer(DebugMixin):
    """
    Optimizes mutations for Golden Gate assembly by balancing codon usage,
    sequence stability, and restriction site compatibility.
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

    # --- Sequence & Context Utilities ---

    def apply_mutation_to_context(self, context: str, codon_pos_in_context: int, alt_seq: str) -> str:
        """
        Applies the alternative sequence to the context at the given position,
        returning a new mutated context.
        """
        context_list = list(context)
        for i, nucleotide in enumerate(alt_seq):
            if codon_pos_in_context + i < len(context_list):
                context_list[codon_pos_in_context + i] = nucleotide
        return "".join(context_list)

    def extract_primer_context(self, mutated_context: str, mutation_start: int, mutation_end: int, flank: int = 30) -> str:
        """
        Extracts a primer design region from the mutated context.
        The region includes 'flank' number of bases upstream and downstream of the mutation.
        """
        start = max(0, mutation_start - flank)
        end = min(len(mutated_context), mutation_end + flank)
        return mutated_context[start:end]

    def process_overhangs(self, alternative: dict, codon_pos_in_context: int, context: str = None) -> list:
        """
        Processes sticky end information for an alternative mutation.
        If the alternative already contains precomputed overhangs (under "overhangs"/"overhang_options"),
        we update missing overhang_start_index values using the current codon position.
        Otherwise, if raw sticky-end data (under "sticky_ends") is available, we compute them.
        """
        # If the alternative already includes overhang options, refresh missing start indices.
        if "overhangs" in alternative and "overhang_options" in alternative["overhangs"]:
            options = alternative["overhangs"]["overhang_options"]
            for option in options:
                if option.get("overhang_start_index") is None:
                    # For now, set the start index to codon_pos_in_context.
                    option["overhang_start_index"] = codon_pos_in_context
                    logger.info(f"Updated overhang_start_index for site {alternative.get('site', 'unknown')} "
                                f"to {codon_pos_in_context}")
            return options

        # Otherwise, process raw sticky end data (if present)
        overhang_options = []
        for pos_key, overhang_data in alternative.get("sticky_ends", {}).items():
            try:
                pos = int(pos_key.split('_')[1])
            except (IndexError, ValueError) as e:
                logger.warning(
                    f"Skipping invalid sticky end key {pos_key}: {e}")
                continue

            if context is not None:
                overhang_start_index = codon_pos_in_context + pos
                if overhang_start_index < 0 or overhang_start_index >= len(context):
                    logger.warning(f"Overhang start index out of bounds: {overhang_start_index} "
                                   f"for mutation at {codon_pos_in_context}, skipping.")
                    continue
            else:
                overhang_start_index = None

            top_strands = overhang_data.get("top_strand", [])
            bottom_strands = overhang_data.get("bottom_strand", [])
            for i, top_overhang in enumerate(top_strands):
                bottom_overhang = bottom_strands[i] if i < len(
                    bottom_strands) else None
                if not top_overhang or not bottom_overhang:
                    logger.warning(
                        f"Skipping overhang at {overhang_start_index}: missing strands.")
                    continue
                entry = {
                    "top_overhang": top_overhang,
                    "bottom_overhang": bottom_overhang,
                    "overhang_start_index": overhang_start_index
                }
                overhang_options.append(entry)
                logger.info(f"Processed overhang: {entry}")
        return overhang_options

    # --- Mutation Set Generation ---

    @DebugMixin.debug_wrapper
    def generate_mutation_sets(self, mutation_options: dict) -> list:
        """
        Generates all possible mutation sets by selecting exactly one alternative codon per site.
        Each mutation entry includes a primer design context and updated overhang start indices.
        """
        if self.debug:
            self.validate(isinstance(mutation_options, dict) and len(mutation_options) > 0,
                          "Mutation options are valid",
                          {k: f"{len(v['codons'])} codons" for k, v in mutation_options.items()})
            self.log_step("Prepare Mutation Choices",
                          "Organizing mutation choices per site")

        mutation_choices_per_site = []
        for site, site_data in mutation_options.items():
            if self.debug:
                self.log_step(
                    "Process Site", f"Processing site {site} at position {site_data['position']}")
            site_mutation_choices = self.process_site_mutation_options(
                site, site_data)
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

    def process_site_mutation_options(self, site: str, site_data: dict) -> list:
        """
        Processes mutation options for a single site.
        If context is available, applies the mutation to generate a mutated context,
        extracts a primer region, and updates overhang start indices.
        """
        mutation_entries = []
        context = site_data.get("context", "")
        context_indices = site_data.get("context_mutated_indices", [])
        print(f"site_data: {site_data}")
        for codon in site_data.get("codons", []):
            for alternative in codon.get("alternative_codons", []):
                # If precomputed overhangs exist, use them (and update missing start indices).
                if "overhangs" in alternative and "overhang_options" in alternative["overhangs"]:
                    overhang_options = alternative["overhangs"]["overhang_options"]
                    if context:
                        site_start_in_context = context_indices[0] if context_indices else 0
                        codon_pos_in_context = site_start_in_context + \
                            (codon["position"] - site_data["position"])
                        # Update any missing overhang_start_index using codon_pos_in_context.
                        for option in overhang_options:
                            if option.get("overhang_start_index") is None:
                                option["overhang_start_index"] = codon_pos_in_context
                                logger.info(
                                    f"Setting overhang_start_index for site {site} to {codon_pos_in_context}")
                    final_overhang_options = overhang_options
                else:
                    # Fallback: compute overhangs from raw sticky-end data.
                    if context:
                        site_start_in_context = context_indices[0] if context_indices else 0
                        codon_pos_in_context = site_start_in_context + \
                            (codon["position"] - site_data["position"])
                        mutated_context = self.apply_mutation_to_context(
                            context, codon_pos_in_context, alternative["seq"])
                        primer_context = self.extract_primer_context(
                            mutated_context, codon_pos_in_context, codon_pos_in_context + len(alternative["seq"]), flank=30)
                        final_overhang_options = self.process_overhangs(
                            alternative, codon_pos_in_context, context=mutated_context)
                    else:
                        primer_context = ""
                        final_overhang_options = self.process_overhangs(
                            alternative, 0, context=None)

                mutated_base_index = (
                    codon["position"] - site_data["position"] + site_data["frame"]) % 6
                original_seq = codon.get("original_sequence", codon.get(
                    "codon_seq", alternative.get("seq", "")))
                mutation_entry = {
                    "site": site,
                    "position": site_data["position"],
                    "original_sequence": original_seq,
                    "alternative_sequence": alternative["seq"],
                    "mutated_base_index": mutated_base_index,
                    "overhangs": {"overhang_options": final_overhang_options},
                    "primer_context": primer_context if context else ""
                }
                mutation_entries.append(mutation_entry)
        return mutation_entries

    # --- Compatibility Matrix & Filtering ---

    @DebugMixin.debug_wrapper
    def create_compatibility_matrices(self, mutation_sets: List[Dict]) -> List[np.ndarray]:
        if self.debug:
            self.validate(isinstance(mutation_sets, list) and len(mutation_sets) > 0,
                          f"Processing {len(mutation_sets)} mutation sets",
                          {"first_set_sites": list(mutation_sets[0].keys()) if mutation_sets else None})
            self.log_step(
                "Matrix Creation", "Creating compatibility matrices for each mutation set")
        compatibility_matrices = []
        for mutation_set_idx, mutation_set in enumerate(tqdm(mutation_sets, desc="Processing Mutation Sets", unit="set")):
            position_keys = list(mutation_set.keys())
            overhang_lists = [
                [option["top_overhang"]
                    for option in mutation_set[pos]["overhangs"]["overhang_options"]]
                for pos in position_keys
            ]
            matrix_shape = tuple(len(overhangs)
                                 for overhangs in overhang_lists)
            if self.debug and mutation_set_idx < 3:
                self.log_step("Matrix Shape", f"Creating matrix with shape {matrix_shape}",
                              {"overhangs_per_site": [len(o) for o in overhang_lists]})
            compatibility_matrix = np.zeros(matrix_shape, dtype=int)
            for combo_indices in product(*[range(len(overhang_list)) for overhang_list in overhang_lists]):
                combo = tuple(overhang_lists[i][idx]
                              for i, idx in enumerate(combo_indices))
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
                    compatibility_matrix[combo_indices] = 1
            compatibility_matrices.append(compatibility_matrix)
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
        for idx in reversed(indices_to_remove):
            del mutation_sets[idx]
            del compatibility_matrices[idx]
        return mutation_sets

    @DebugMixin.debug_wrapper
    def optimize_mutations(self, mutation_options: Dict) -> List[Dict]:
        """
        Finds the optimal combination of mutations that are compatible according to the BsmBI sticky-end strategy.
        """
        if not self.debug:
            with debug_context("mutation_optimization"):
                logger.info("Generating possible mutation sets...")
                mutation_sets = self.generate_mutation_sets(mutation_options)
                logger.info("Computing compatibility matrices...")
                compatibility_matrices = self.create_compatibility_matrices(
                    mutation_sets)
                logger.info(
                    "Filtering mutation sets based on compatibility...")
                optimized_mutations = self.filter_compatible_mutations(
                    mutation_sets, compatibility_matrices)
                logger.debug(
                    f"Optimized mutation sets: {len(optimized_mutations)}")
                return optimized_mutations, compatibility_matrices
        else:
            self.validate(isinstance(mutation_options, dict) and len(mutation_options) > 0,
                          "Mutation options are valid",
                          {k: v for k, v in list(mutation_options.items())[:2]})
            self.log_step("Generate Mutation Sets",
                          "Creating all possible combinations of mutations")
            mutation_sets = self.generate_mutation_sets(mutation_options)
            self.log_step("Compute Compatibility",
                          "Creating compatibility matrices for mutation sets")
            compatibility_matrices = self.create_compatibility_matrices(
                mutation_sets)
            self.log_step("Filter Compatible Mutations",
                          "Removing mutation sets with no valid compatibility")
            optimized_mutations = self.filter_compatible_mutations(
                mutation_sets, compatibility_matrices)
            return optimized_mutations, compatibility_matrices
