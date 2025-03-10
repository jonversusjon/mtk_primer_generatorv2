import logging
from typing import Dict, List, Optional, Any
from Bio.Seq import Seq
from .utils import GoldenGateUtils
from .mutation_optimizer import MutationOptimizer
from itertools import product
from config.logging_config import logger
from services.base import debug_context
from services.debug.debug_utils import MutationDebugger
from functools import wraps
from services.debug.debug_mixin import DebugMixin


class MutationAnalyzer(DebugMixin):
    """
    MutationAnalyzer Module

    This module provides functionality for analyzing DNA sequences to generate mutation options 
    for Golden Gate assembly. The main entry and exit point is the get_all_mutations function.

    Input Parameters:
        - sequence to domesticate (str): The full DNA sequence that will be domesticated.
        - recognition_sites_to_mutate (List[Dict]): A list of dictionaries, each representing a recognition site to mutate.
        Each recognition site dictionary must include the following keys:
            • "position" (int): The start position of the recognition site within the sequence to domesticate.
            • "sequence" (str): The recognition site sequence.
            • "enzyme" (str): The name of the restriction enzyme.
            • "codons" (List[Dict]): A list of codon dictionaries for the recognition site, where each codon dictionary includes:
                    - "codon_seq" (str): The original codon sequence (e.g., "GAG").
                    - "position" (int): The position of the codon within the sequence to domesticate.
                    - "amino_acid" (str): The amino acid encoded by the codon.
            • "frame" (int): The reading frame of the recognition site.
            • "strand" (str): The DNA strand ("+" or "-") on which the recognition site is located.
        During processing, each recognition site dictionary is enriched with:
            • "context_sequence": A subsequence spanning 30 bp upstream and 30 bp downstream of the recognition site.
            • "overlapping_indices": Computed indices (relative to codon positions 0, 1, 2) that fall within 
                the recognition region (the recognition sequence is assumed to start at index 30 within the context).

    Additional Initialization Parameter:
        - codon_usage_dict (Dict[str, Dict[str, float]]): A dictionary providing codon usage frequencies for each amino acid,
        organized by organism.

    Output Data Structure:
        The get_all_mutations function returns a dictionary keyed by "mutation_{recognition_site_position}", where each value is a dictionary with:
            - "position": The recognition site position.
            - "sequence": The recognition site sequence.
            - "frame": The reading frame.
            - "strand": The DNA strand.
            - "enzyme": The enzyme name.
            - "codons": A list of dictionaries for each codon, each containing:
                    • "original_codon_sequence": The original codon sequence.
                    • "position": The codon position.
                    • "amino_acid": The encoded amino acid.
                    • "alternative_codons": A list of dictionaries for each valid alternative codon, each containing:
                            - "seq": The alternative codon sequence.
                            - "usage": The codon usage frequency.
                            - "mutations": A 3-tuple of integers indicating changes (1 for a mutation, 0 otherwise).
                            - "changes_in_site": A list of codon index positions (0, 1, or 2) where the candidate alters the recognition region.
                            - "sticky_ends": A dictionary keyed by "position_{codon_index}" that contains:
                                    • "top_strand": A list of dictionaries, each with:
                                            - "seq": The overhang sequence.
                                            - "start_index": The starting index in the context.
                                    • "bottom_strand": A list of dictionaries with corresponding reverse-complement data.
                            - "mutation_positions_in_context": A list of indices in the context where the mutation occurs.
    """

    def __init__(self, codon_usage_dict: Dict[str, Dict[str, float]], max_mutations: int = 1, verbose: bool = False, debug: bool = False):
        self.logger = logger.getChild("MutationAnalyzer")
        self.utils = GoldenGateUtils()
        self.state = {'current_codon': '',
                      'current_position': 0, 'mutations_found': []}
        self.codon_usage_dict = codon_usage_dict
        self.max_mutations = max_mutations
        self.verbose = verbose
        self.debug = debug
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        if self.debug:
            self.debugger = MutationDebugger(
                parent_logger=logger, use_custom_format=True)
            if hasattr(self.debugger.logger, 'propagate'):
                self.debugger.logger.propagate = False
            self.logger.info("Debug mode enabled for MutationAnalyzer")
            self.validate(codon_usage_dict and isinstance(codon_usage_dict, dict),
                          "Codon usage dictionary is valid",
                          {"organisms": list(codon_usage_dict.keys())})
            self.validate(isinstance(max_mutations, int) and max_mutations > 0,
                          f"Max mutations set to {max_mutations}",
                          {"valid_range": "1-3"})

    @DebugMixin.debug_wrapper
    def get_all_mutations(self, sites_to_mutate: List[Dict]) -> Dict[str, List[Dict]]:
        self.validate(sites_to_mutate and isinstance(sites_to_mutate, list),
                      f"Received {len(sites_to_mutate)} sites to mutate",
                      {"first_site": sites_to_mutate[0] if sites_to_mutate else None})
        self.log_step(
            "Analysis Start", f"Starting mutation analysis for {len(sites_to_mutate)} sites")
        mutation_options = {}
        try:
            with debug_context("mutation_analysis"):
                for site_idx, site in enumerate(sites_to_mutate):
                    print(f" ****** site ****** : {site}")
                    """
                    {
                    "position": 477,
                    "sequence": "GAGACC",
                    "frame": 0,
                    "codons": [
                        {
                        "codon_seq": "GAG",
                        "amino_acid": "E",
                        "position": 477
                        },
                        {
                        "codon_seq": "ACC",
                        "amino_acid": "T",
                        "position": 480
                        }
                    ],
                    "strand": "-",
                    "context_sequence": "GCCGAGCGCGCCGGGGTGCCCGCCTTCCTGGAGACCTCCGCGCCCCGCAACCTCCCCTTCTACGAG",
                    "context_recognition_site_indices": [
                        30,
                        31,
                        32,
                        33,
                        34,
                        35
                    ],
                    "enzyme": "BsaI"
                    }"""
                    self.log_step("Process Site",
                                  f"Analyzing site {site_idx+1}/{len(sites_to_mutate)}",
                                  {"position": site["position"],
                                   "sequence": site["sequence"],
                                   "enzyme": site["enzyme"],
                                   "codons": len(site["codons"])})

                    frame = site["frame"]
                    # # Initialize a streamlined mutations dict for the site.
                    # site_mutations = {"position": site["position"],
                    #                   "sequence": site["sequence"],
                    #                   "frame": frame,
                    #                   "strand": site["strand"],
                    #                   "enzyme": site["enzyme"],
                    #                   "codons": []}
                    # TODO: fixing the way context and mutated context are stored

                    for codon_idx, codon in enumerate(site["codons"]):
                        self.log_step("Process Codon",
                                      f"Analyzing codon {codon_idx+1}/{len(site['codons'])}",
                                      {"codon_seq": codon["codon_seq"],
                                       "position": codon["position"],
                                       "amino_acid": codon["amino_acid"]})
                        codon_offset = codon["position"] - site["position"]

                        # Inject overlapping_indices for this codon.
                        site_with_indices = {
                            **site, "overlapping_indices": self.utils.get_recognition_site_bases(frame, codon_idx)}
                        alternatives = self._find_alternative_codons(
                            site_with_indices, codon["codon_seq"], codon["amino_acid"])
                        self.validate(alternatives,
                                      f"Found {len(alternatives)} alternative codons",
                                      {"alternatives": alternatives})
                        # Enrich each alternative with sticky ends and context mutation positions.
                        for alt in alternatives:
                            candidate = alt["seq"]
                            # Determine which positions in the codon are mutated.
                            muts = [i for i in range(
                                3) if candidate[i] != codon["codon_seq"][i]]

                            # Calculate the mutated context and record mutation positions.
                #             mutation_info = self._calculate_mutated_context(site["context_sequence"],
                #                                                             codon_start_in_context,
                #                                                             codon["codon_seq"],
                #                                                             candidate)

                #             alt.update(mutation_info)

                #             alt["sticky_ends"] = self._calculate_sticky_ends_with_context(site["context_sequence"],
                #                                                                           mutation_info["mutated_context"],
                #                                                                           codon_start_in_context,
                #                                                                           muts)
                #         site_mutations["codons"].append({
                #             "original_codon_sequence": codon["codon_seq"],
                #             "position": codon["position"],
                #             "amino_acid": codon["amino_acid"],
                #             "alternative_codons": alternatives
                #         })
                #     mutation_options[f"mutation_{site_mutations['position']}"] = site_mutations
                #     self.log_step("Site Result",
                #                   f"Completed analysis for site at position {site['position']}",
                #                   {"alternatives_found": sum(len(c["alternative_codons"]) for c in site_mutations["codons"])})
                # if self.verbose:
                #     logger.info(
                #         f"Found {len(mutation_options)} sites with valid mutations")
                # self.validate(mutation_options,
                #               f"Generated mutation options for {len(mutation_options)} sites",
                #               {"site_keys": list(mutation_options.keys())})
                return mutation_options

        except Exception as e:
            error_msg = f"Error in mutation analysis: {str(e)}"
            logger.error(error_msg)
            if getattr(self, 'debugger', None):
                self.debugger.log_error(error_msg)
            return mutation_options

    @DebugMixin.debug_wrapper
    def _find_alternative_codons(self, site: Dict, original_codon: str, amino_acid: str) -> List[Dict]:
        for key in ["overlapping_indices", "sequence", "context_sequence"]:
            if key not in site:
                self.log_step("Site Validation",
                              f"Missing key: {key}", level=logging.ERROR)
                raise ValueError(f"Site dictionary missing key: {key}")
        self.log_step("Find Alternatives",
                      f"Finding alternatives for codon {original_codon} (AA: {amino_acid})",
                      {"max_mutations": self.max_mutations,
                       "site_sequence": site["sequence"],
                       "context_sequence": site["context_sequence"],
                       "overlapping_indices": site["overlapping_indices"]})
        alt_codons = self.utils.get_codons_for_amino_acid(amino_acid)
        self.log_step("Possible Codons",
                      f"Found {len(alt_codons)} synonymous codons for {amino_acid}",
                      {"codons": alt_codons})
        valid = []
        for cand in alt_codons:
            if cand == original_codon:
                continue
            muts = [i for i in range(3) if cand[i] != original_codon[i]]
            if len(muts) > self.max_mutations:
                self.log_step("Skip Candidate",
                              f"Candidate {cand} has {len(muts)} mutations (max: {self.max_mutations})",
                              {"mutations": muts}, level=logging.DEBUG)
                continue
            changes = [i for i in muts if i in site["overlapping_indices"]]
            if not changes:
                self.log_step("Skip Candidate",
                              f"Candidate {cand} doesn't alter recognition region",
                              {"mutations": muts,
                                  "overlapping_indices": site["overlapping_indices"]},
                              level=logging.DEBUG)
                continue
            mut_tuple = tuple(
                1 if cand[i] != original_codon[i] else 0 for i in range(3))
            usage = self.utils.get_codon_usage(
                cand, amino_acid, self.codon_usage_dict)
            candidate_info = {"seq": cand, "usage": usage,
                              "mutations": mut_tuple, "changes_in_site": changes}
            valid.append(candidate_info)
            self.log_step("Valid Alternative",
                          f"Candidate {cand} accepted",
                          {"mutations": mut_tuple, "changes_in_site": changes},
                          level=logging.DEBUG)
        self.validate(valid,
                      f"Generated {len(valid)} valid alternatives",
                      {"original": original_codon, "amino_acid": amino_acid})
        return valid

    def _calculate_sticky_ends_with_context(self, original_ctx: str, mutated_ctx: str, codon_start: int, mutation_positions: List[int]) -> Dict:
        sticky = {}
        for pos in mutation_positions:
            m_pos = codon_start + pos
            pos_sticky = {"top_strand": [], "bottom_strand": []}
            ranges = [range(m_pos - 3, m_pos + 1),
                      range(m_pos - 2, m_pos + 2),
                      range(m_pos - 1, m_pos + 3),
                      range(m_pos, m_pos + 4)]
            self.log_step("Calculate Sticky Ends",
                          f"Calculating sticky ends for mutation at position {pos}",
                          {"mutation_pos_in_context": m_pos,
                              "sticky_positions": [list(r) for r in ranges]},
                          level=logging.DEBUG)
            for r in ranges:
                if 0 <= min(r) and max(r) < len(mutated_ctx):
                    top = "".join(mutated_ctx[i] for i in r)
                    bottom = self.utils.reverse_complement(top)
                    pos_sticky["top_strand"].append(
                        {"seq": top, "start_index": r.start})
                    pos_sticky["bottom_strand"].append(
                        {"seq": bottom, "start_index": r.start})
                    self.log_step("Sticky End",
                                  f"Generated sticky end for range {list(r)}",
                                  {"top_strand": {"seq": top, "start_index": r.start},
                                   "bottom_strand": {"seq": bottom, "start_index": r.start}},
                                  level=logging.DEBUG)
            sticky[f"position_{pos}"] = pos_sticky
        return sticky

    @DebugMixin.debug_wrapper
    def _calculate_mutated_context(self, original_context: str, codon_start_in_context: int,
                                   original_codon: str, candidate_codon: str) -> Dict[str, Any]:
        """
        Calculates the mutated context for a candidate codon by performing only point mutations.
        The mutated context will retain the full length of the original context.

        Parameters:
            original_context (str): The full context string (e.g. 66 bp).
            codon_start_in_context (int): The starting index of the codon within the original_context.
            original_codon (str): The original codon sequence (length 3).
            candidate_codon (str): The candidate alternative codon (length 3).

        Returns:
            dict: Contains:
                - "mutated_context": The context string after applying the point mutation(s).
                - "mutation_positions_in_context": A list of indices (absolute in original_context)
                where the mutation(s) occurred.

        Raises:
            ValueError: If the codon indices are out of bounds or the mutated context length doesn't match.
        """
        # Validate input codon lengths.
        if len(original_codon) != 3 or len(candidate_codon) != 3:
            error_msg = "Both original and candidate codons must be of length 3."
            self.log_step("Input Error", error_msg, level=logging.ERROR)
            raise ValueError(error_msg)

        expected_length = len(original_context)

        # Check boundaries: ensure the codon (3 bp) fits within the original_context.
        if codon_start_in_context < 0 or codon_start_in_context + 3 > expected_length:
            error_msg = (f"Codon start index {codon_start_in_context} with codon length 3 "
                         f"is out of bounds for context of length {expected_length}.")
            self.log_step("Boundary Error", error_msg, level=logging.ERROR)
            raise ValueError(error_msg)

        # Perform point mutation: replace only the differing base(s) without changing overall length.
        mutated_ctx_list = list(original_context)
        mutation_positions = []
        for i in range(3):
            if candidate_codon[i] != original_codon[i]:
                pos = codon_start_in_context + i
                mutated_ctx_list[pos] = candidate_codon[i]
                mutation_positions.append(pos)

        mutated_context = "".join(mutated_ctx_list)

        # Verify the length remains unchanged.
        if len(mutated_context) != expected_length:
            error_msg = (f"Mutated context length ({len(mutated_context)}) does not match "
                         f"original context length ({expected_length}).")
            self.log_step("Length Mismatch", error_msg, level=logging.ERROR)
            raise ValueError(error_msg)

        result = {
            "mutated_context": mutated_context,
            "mutation_positions_in_context": mutation_positions
        }
        self.log_step("Mutated Context",
                      "Calculated mutated context for candidate codon",
                      result,
                      level=logging.DEBUG)
        return result

    @DebugMixin.debug_wrapper
    def _create_mutation_entry(self, site_sequence: str, site_start: int, codon_start: int, codon_pos: int,
                               original_codon: str, new_codon: str, frequency: float, enzyme: str, strand: str) -> Optional[Dict]:
        self.log_step("Create Mutation Entry",
                      f"Creating mutation entry for {original_codon} → {new_codon}",
                      {"site_sequence": site_sequence, "site_start": site_start, "codon_start": codon_start,
                       "codon_pos": codon_pos, "frequency": frequency, "enzyme": enzyme, "strand": strand})
        mutation_position = codon_start + codon_pos
        original_base = original_codon[codon_pos]
        new_base = new_codon[codon_pos]
        mutated_site = site_sequence[:codon_pos] + \
            new_base + site_sequence[codon_pos+1:]
        self.log_step("Mutation Check", "Checking if mutation disrupts restriction site",
                      {"original_site": site_sequence, "mutated_site": mutated_site,
                       "changed_base_pos": codon_pos, "original_base": original_base, "new_base": new_base})
        if mutated_site == site_sequence:
            if getattr(self, 'debugger', None):
                self.debugger.log_warning(
                    f"Mutation does not disrupt restriction site: {original_base} → {new_base} at position {mutation_position}")
            return None
        entry = {"position": mutation_position,
                 "original_base": original_base,
                 "new_base": new_base,
                 "original_codon": original_codon,
                 "new_codon": new_codon,
                 "frequency": frequency,
                 "site_start": site_start,
                 "enzyme": enzyme,
                 "strand": strand}
        mutation_id = f"{mutation_position}_{original_base}_{new_base}"
        if not hasattr(self, 'mutation_dict'):
            self.mutation_dict = {}
        self.mutation_dict[mutation_id] = entry
        self.validate(mutation_id in self.mutation_dict,
                      f"Successfully created mutation entry: {mutation_id}",
                      {"entry": entry})
        return entry
