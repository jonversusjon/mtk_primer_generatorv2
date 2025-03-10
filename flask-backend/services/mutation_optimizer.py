from typing import Dict, List
from Bio.Seq import Seq
import numpy as np
from itertools import product
from tqdm import tqdm
import logging

from services.base import debug_context
from services.debug.debug_mixin import DebugMixin
from services.debug.debug_utils import MutationDebugger
from config.logging_config import logger
from .utils import GoldenGateUtils


class MutationOptimizer(DebugMixin):
    """
    MutationOptimizer Module

    This module optimizes mutation combinations generated by the MutationAnalyzer module
    for Golden Gate assembly compatibility. It evaluates alternative codon mutations based on
    sticky-end compatibility, selecting sets of mutations that satisfy the Golden Gate cloning criteria.

    Input:
        - mutation_options (Dict): Output from the MutationAnalyzer.get_all_mutations function.
        Structure:
            {
            "mutation_{recognition_site_position}": {
                "position": int,
                "sequence": str,
                "frame": int,
                "strand": str,
                "enzyme": str,
                "context_sequence": str,
                "codons": [
                    {
                        "original_codon_sequence": str,
                        "position": int,
                        "amino_acid": str,
                        "alternative_codons": [
                            {
                                "seq": str,
                                "usage": float,
                                "mutations": tuple(int, int, int),
                                "changes_in_site": List[int],
                                "sticky_ends":  {
                                    "top_strand": {
                                        "seq": str,
                                        "start_index": int
                                    },
                                    "bottom_strand": {
                                        "seq": str,
                                        "start_index": int
                                    }
                                },
                                "mutation_positions_in_context": List[int]
                            }

    Output:
        - optimized_mutations (List[Dict]): A list of mutation sets, where each mutation set
        represents a combination of mutations that is fully compatible according to the Golden Gate
        sticky-end assembly rules. Each mutation entry within a set contains:
            • "site": recognition site identifier (e.g., "mutation_477").
            • "position": recognition site position.
            • "original_codon_sequence": the original codon sequence.
            • "alternative_codon_sequence": the selected alternative codon.
            • "mutated_base_index": position of the mutation within the recognition site.
            • "overhangs": Dict containing "overhang_options", a list with:
                    - "top_overhang": dict {"seq": str, "start_index": int}
                    - "bottom_overhang": dict {"seq": str, "start_index": int}
                    - "overhang_start_index": int
            • "primer_context": DNA context suitable for primer design around the mutation.

        - compatibility_matrices (List[np.ndarray]): List of N-dimensional numpy arrays corresponding
        to each mutation set, indicating compatible (1) and incompatible (0) overhang combinations.
        The compatibility matrix list is aligned with the mutation set list, where each matrix
        corresponds to the respective mutation set.

    Entry and Exit Point:
        - optimize_mutations is the orchestrator function that serves as the entry and exit point of the module,
        coordinating mutation set generation, compatibility analysis, and filtering steps.
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

        self.compatibility_table = self.utils._load_compatibility_table(
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

    # --- Context Utilities ---

    def apply_mutation_to_context(self, context: str, pos: int, alt_seq: str) -> str:
        return context[:pos] + alt_seq + context[pos+len(alt_seq):]

    def extract_primer_context(self, mutated_ctx: str, start: int, end: int, flank: int = 30) -> str:
        return mutated_ctx[max(0, start - flank):min(len(mutated_ctx), end + flank)]

    def process_overhangs(self, alternative: dict, codon_pos: int, context: str = None) -> list:
        overhang_options = []
        for key, overhang_data in alternative.get("sticky_ends", {}).items():
            # TODO: whatever used the position_ notation is no longer used and needs to be removed
            if not key.startswith("position_"):
                logger.warning(f"Skipping invalid sticky end key {key}")
                continue
            pos = int(key.split("_")[1])
            overhang_start = codon_pos + pos if context else None
            if context and (overhang_start < 0 or overhang_start >= len(context)):
                logger.warning(
                    f"Overhang start index {overhang_start} out of bounds for mutation at {codon_pos}")
                continue
            top_strands = overhang_data.get("top_strand", [])
            bottom_strands = overhang_data.get("bottom_strand", [])
            for top, bottom in zip(top_strands, bottom_strands):
                if not top or not bottom:
                    logger.warning(
                        f"Missing overhang strands at {overhang_start}")
                    continue
                entry = {"top_overhang": top, "bottom_overhang": bottom,
                         "overhang_start_index": overhang_start}
                overhang_options.append(entry)
                logger.info(f"Processed overhang: {entry}")
        return overhang_options

    # --- Mutation Set Generation ---

    @DebugMixin.debug_wrapper
    def generate_mutation_sets(self, mutation_options: dict) -> list:
        self.validate(mutation_options and isinstance(mutation_options, dict),
                      f"Received {len(mutation_options)} mutation option sites")
        mutation_choices = []
        for site_key, site_data in mutation_options.items():
            self.log_step(
                "Process Site", f"Processing site {site_key} at position {site_data['position']}")
            site_choices = self.process_site_mutation_options(
                site_key, site_data)
            self.validate(
                site_choices, f"Found {len(site_choices)} valid mutation options for {site_key}")
            mutation_choices.append(site_choices)
        all_sets = [{mutation["site"]: mutation for mutation in combo}
                    for combo in product(*mutation_choices)]
        if self.debug:
            expected = np.prod([len(choices) for choices in mutation_choices])
            self.validate(len(all_sets) == expected,
                          f"Created {len(all_sets)} mutation sets (expected: {expected})",
                          {"sites": len(mutation_choices)})
        return all_sets

    def process_site_mutation_options(self, site: str, site_data: dict) -> list:
        self.log_step("Process Site Mutations",
                      f"Processing mutations for site {site} with site data {site_data}")
        entries = []
        context_sequence = site_data.get("context_sequence", "")
        indices = site_data.get("context_recognition_site_indices", [])
        for codon in site_data.get("codons", []):
            for alt in codon.get("alternative_codons", []):
                if context_sequence:
                    site_start = indices[0] if indices else 0
                    codon_pos = site_start + \
                        (codon["position"] - site_data["position"])
                    mutated_ctx = self.apply_mutation_to_context(
                        context_sequence, codon_pos, alt["seq"])
                    primer_ctx = self.extract_primer_context(
                        mutated_ctx, codon_pos, codon_pos + len(alt["seq"]))
                    overhang_opts = self.process_overhangs(
                        alt, codon_pos, context=mutated_ctx)
                else:
                    primer_ctx = ""
                    overhang_opts = self.process_overhangs(
                        alt, 0, context=None)
                m_index = (codon["position"] -
                           site_data["position"] + site_data["frame"]) % 6
                original_seq = codon.get("original_codon_sequence", codon.get(
                    "codon_seq", alt.get("seq", "")))
                entry = {
                    "site": site,
                    "position": site_data["position"],
                    "original_codon_sequence": original_seq,
                    "alternative_codon_sequence": alt["seq"],
                    "mutated_base_index": m_index,
                    "overhangs": {"overhang_options": overhang_opts},
                    "primer_context": primer_ctx
                }
                entries.append(entry)
        return entries

    # --- Compatibility Matrix & Filtering ---

    @DebugMixin.debug_wrapper
    def create_compatibility_matrices(self, mutation_sets: List[Dict]) -> List[np.ndarray]:
        compatibility_matrices = []
        for mutation_set in tqdm(mutation_sets, desc="Processing Mutation Sets", unit="set"):
            keys = list(mutation_set.keys())
            overhang_lists = [
                [opt["top_overhang"]
                    for opt in mutation_set[k]["overhangs"]["overhang_options"]]
                for k in keys
            ]
            shape = tuple(len(lst) for lst in overhang_lists)
            matrix = np.zeros(shape, dtype=int)
            for combo_indices in product(*[range(len(lst)) for lst in overhang_lists]):
                combo = tuple(overhang_lists[i][idx]
                              for i, idx in enumerate(combo_indices))
                n = len(combo)
                if all(self.compatibility_table[self.utils.seq_to_index(combo[i]["seq"])][self.utils.seq_to_index(combo[j]["seq"])] == 1

                       for i in range(n) for j in range(i+1, n)):
                    matrix[combo_indices] = 1
            compatibility_matrices.append(matrix)
        return compatibility_matrices

    @DebugMixin.debug_wrapper
    def filter_compatible_mutations(self, mutation_sets: list, compatibility_matrices: list) -> list:
        indices = [i for i, mat in enumerate(
            compatibility_matrices) if np.all(mat == 0)]
        for i in sorted(indices, reverse=True):
            del mutation_sets[i]
            del compatibility_matrices[i]
        return mutation_sets

    @DebugMixin.debug_wrapper
    def optimize_mutations(self, mutation_options: Dict) -> List[Dict]:
        with debug_context("mutation_optimization"):
            self.log_step("Generate Mutation Sets",
                          "Creating all possible mutation combinations")
            mutation_sets = self.generate_mutation_sets(mutation_options)
            self.log_step("Compute Compatibility",
                          "Creating compatibility matrices")
            matrices = self.create_compatibility_matrices(mutation_sets)
            self.log_step("Filter Compatible Mutations",
                          "Removing incompatible mutation sets")
            optimized = self.filter_compatible_mutations(
                mutation_sets, matrices)

            return optimized, matrices
