import logging
from typing import Dict, List, Optional, Any
from config.logging_config import logger
from services.base import debug_context
from services.debug.debug_mixin import DebugMixin
from services.utils import GoldenGateUtils
from services.debug.debug_utils import MutationDebugger


class MutationAnalyzer(DebugMixin):
    """
    MutationAnalyzer Module

    This module provides functionality for analyzing DNA sequences to generate mutation options 
    for Golden Gate assembly. The main entry and exit point is the get_all_mutations function.

    [Documentation truncated for brevity…]
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
    def get_all_mutations(self, sites_to_mutate: List[Dict]) -> Dict[str, Dict]:
        # Log entry with full input summary
        self.logger.debug(
            f"Entering get_all_mutations with {len(sites_to_mutate) if sites_to_mutate else 0} sites to mutate.")
        self.validate(sites_to_mutate and isinstance(sites_to_mutate, list),
                      f"Received {len(sites_to_mutate)} sites to mutate",
                      {"first_site": sites_to_mutate[0] if sites_to_mutate else None})
        self.log_step("Analysis Start", f"Starting mutation analysis for {len(sites_to_mutate)} sites",
                      {"sites_overview": [site.get('position') for site in sites_to_mutate]})

        mutation_options = {}
        try:
            with debug_context("mutation_analysis"):
                for site_idx, site in enumerate(sites_to_mutate):
                    site_key = f"mutation_{site.get('position')}"
                    self.log_step("Process Site",
                                  f"Analyzing site {site_idx+1}/{len(sites_to_mutate)} at position {site.get('position')}",
                                  {"site_details": site})

                    # Prepare site-level mutation data and start detailed logging
                    site_mutations = {
                        "position": site.get("position"),
                        "sequence": site.get("sequence"),
                        "frame": site.get("frame"),
                        "strand": site.get("strand"),
                        "enzyme": site.get("enzyme"),
                        "codons": []
                    }
                    site_has_alternatives = False

                    codons_list = site.get("codons", [])
                    self.logger.debug(
                        f"Site {site.get('position')}: Found {len(codons_list)} codons for analysis.")
                    for codon_idx, codon in enumerate(codons_list):
                        self.log_step("Process Codon",
                                      f"Analyzing codon {codon_idx+1}/{len(codons_list)}",
                                      {"raw_codon_data": codon})
                        if codon.get("codon_sequence") is None:
                            self.logger.warning(
                                f"Site {site.get('position')}: Codon {codon_idx+1} missing codon_sequence. Skipping codon. Full codon data: {codon}")
                            continue

                        self.log_step("Process Codon",
                                      f"Codon {codon_idx+1}: Processing sequence {codon.get('codon_sequence')}",
                                      {"codon_details": {
                                          "codon_sequence": codon.get("codon_sequence"),
                                          "context_position": codon.get("context_position"),
                                          "amino_acid": codon.get("amino_acid")
                                      }})
                        # Inject overlapping indices into a copy of the site.
                        site_with_indices = {**site,
                                             "overlapping_indices": self.utils.get_recognition_site_bases(site.get("frame"), codon_idx)}
                        alternatives = self._find_alternative_codons(
                            site_with_indices,
                            codon["codon_sequence"],
                            codon.get("amino_acid"),
                            codon.get("context_position"))

                        if alternatives:
                            site_has_alternatives = True
                            self.log_step("Found Alternatives",
                                          f"Site {site.get('position')}, codon {codon.get('codon_sequence')}: Found {len(alternatives)} alternatives",
                                          {"alternatives": alternatives})
                            # Process each alternative codon
                            for alt in alternatives:
                                mutated_context, first_mut, last_mut = self.get_mutated_context(
                                    context_sequence=site.get(
                                        "context_sequence"),
                                    codon_context_position=codon.get(
                                        "context_position"),
                                    new_codon_sequence=alt.get("seq"),
                                    mutated_context_mutated_bases=alt.get("mutations"))
                                alt["mutated_context"] = mutated_context
                                alt["sticky_ends"] = self._calculate_sticky_ends_with_context(
                                    mutated_context, first_mut, last_mut)
                                self.logger.debug(
                                    f"Alternative processed for codon {codon.get('codon_sequence')}: {alt}")
                            site_mutations["codons"].append({
                                "codon_sequence": codon["codon_sequence"],
                                "context_position": codon.get("context_position"),
                                "amino_acid": codon.get("amino_acid"),
                                "alternative_codons": alternatives
                            })
                        else:
                            self.logger.debug(
                                f"Site {site.get('position')}, codon {codon.get('codon_sequence')}: No alternatives found.")

                    if site_has_alternatives:
                        try:
                            from models.mutations import MutationSite
                            validated_site = MutationSite.model_validate(
                                site_mutations)
                            mutation_options[site_key] = validated_site.model_dump(
                            )
                            self.log_step("Site Result",
                                          f"Completed analysis for site at position {site.get('position')}",
                                          {"alternatives_found": sum(len(c.get("alternative_codons", [])) for c in site_mutations["codons"]),
                                           "validated_site": mutation_options[site_key]})
                        except Exception as e:
                            self.logger.error(
                                f"Validation error for site {site.get('position')}: {e}. Site mutations data: {site_mutations}")
                            continue  # Skip this site if validation fails
                    else:
                        self.log_step("No Alternatives Found",
                                      f"Site {site.get('position')}: No alternative codons found",
                                      {"site": site.get("position")},
                                      level=logging.WARNING)
                        self.logger.debug(
                            f"Site {site.get('position')}: Mutation analysis skipped due to absence of alternatives.")

                self.logger.debug(
                    f"Mutation options collected: {mutation_options}")
                if self.verbose:
                    self.logger.info(
                        f"Found {len(mutation_options)} sites with valid mutations")
                self.validate(mutation_options,
                              f"Generated mutation options for {len(mutation_options)} sites",
                              {"site_keys": list(mutation_options.keys())})
                return mutation_options

        except Exception as e:
            self.logger.error(
                f"Critical error in mutation analysis: {e}", exc_info=True)
            raise e

    @DebugMixin.debug_wrapper
    def _find_alternative_codons(self, site: Dict, original_codon: str, amino_acid: str, codon_context_position: int) -> List[Dict]:
        # Verify required keys exist
        for key in ["overlapping_indices", "sequence", "context_sequence"]:
            if key not in site:
                self.log_step("Site Validation",
                              f"Missing key: {key}",
                              {"site_keys": list(site.keys())},
                              level=logging.ERROR)
                raise ValueError(f"Site dictionary missing key: {key}")

        # Log the initial parameters for finding alternatives
        self.log_step("Find Alternatives Start",
                      f"Finding alternatives for codon {original_codon} (AA: {amino_acid})",
                      {"max_mutations": self.max_mutations,
                       "site_sequence": site.get("sequence"),
                       "context_sequence": site.get("context_sequence"),
                       "overlapping_indices": site.get("overlapping_indices"),
                       "codon_context_position": codon_context_position})
        print(
            f"Finding alternatives for codon {original_codon} (AA: {amino_acid})")

        # Retrieve all synonymous codons for the given amino acid
        alt_codons = self.utils.get_codons_for_amino_acid(amino_acid)
        self.log_step("Possible Codons",
                      f"Found {len(alt_codons)} synonymous codons for amino acid {amino_acid}",
                      {"codons": alt_codons})

        valid = []
        for cand in alt_codons:
            self.log_step("Candidate Evaluation Start",
                          f"Evaluating candidate codon: {cand}",
                          {"original_codon": original_codon})
            # Skip candidate if it is identical to the original codon
            if cand == original_codon:
                self.log_step("Candidate Skipped",
                              f"Candidate {cand} is identical to the original codon.",
                              level=logging.DEBUG)
                continue

            # Determine which positions differ between candidate and original
            muts = [i for i in range(3) if cand[i] != original_codon[i]]
            self.log_step("Mutation Differences",
                          f"Candidate {cand} differs from {original_codon} at positions {muts}",
                          {"candidate": cand, "original": original_codon,
                              "differences": muts},
                          level=logging.DEBUG)

            # Skip candidate if number of mutations exceeds allowed max_mutations
            if len(muts) > self.max_mutations:
                self.log_step("Candidate Skipped",
                              f"Candidate {cand} has {len(muts)} mutations (allowed: {self.max_mutations}).",
                              {"mutations": muts},
                              level=logging.DEBUG)
                continue

            # Compute the positions in the context sequence where mutations occur
            mutation_positions_in_context = [
                codon_context_position + i for i in muts]
            self.log_step("Mutation Positions",
                          f"Candidate {cand} mutation positions in context: {mutation_positions_in_context}",
                          {"mutations": muts,
                              "codon_context_position": codon_context_position},
                          level=logging.DEBUG)

            # Determine which of these mutation positions are within the recognition (overlapping) region
            changes = [i for i in muts if i in site.get(
                "overlapping_indices", [])]
            self.log_step("Overlap Check",
                          f"Candidate {cand} changes within recognition region: {changes}",
                          {"overlapping_indices": site.get(
                              "overlapping_indices"), "mutations": muts},
                          level=logging.DEBUG)

            # If no changes occur in the overlapping (recognition) region, skip candidate
            if not changes:
                self.log_step("Candidate Skipped",
                              f"Candidate {cand} does not alter the recognition region.",
                              {"mutations": muts, "overlapping_indices": site.get(
                                  "overlapping_indices")},
                              level=logging.DEBUG)
                continue

            # Build a mutation tuple representing which positions are mutated (1) or not (0)
            mut_tuple = tuple(
                1 if cand[i] != original_codon[i] else 0 for i in range(3))
            # Retrieve codon usage score for logging purposes
            usage = self.utils.get_codon_usage(
                cand, amino_acid, self.codon_usage_dict)
            candidate_info = {
                "seq": cand,
                "usage": usage,
                "mutations": mut_tuple,
                "changes_in_site": changes,
                "mutation_positions_in_context": mutation_positions_in_context
            }
            valid.append(candidate_info)
            self.log_step("Candidate Accepted",
                          f"Candidate {cand} accepted.",
                          {"mutations": mut_tuple,
                              "changes_in_site": changes, "usage": usage},
                          level=logging.DEBUG)

        # Final validation log of how many valid candidates were generated
        self.validate(valid,
                      f"Generated {len(valid)} valid alternatives for codon {original_codon}",
                      {"original": original_codon, "amino_acid": amino_acid})
        self.log_step("Find Alternatives Complete",
                      f"Total valid alternatives for codon {original_codon}: {len(valid)}",
                      {"valid_candidates": valid})
        return valid

    def _calculate_sticky_ends_with_context(self, mutated_ctx: str, first_mut_idx: int, last_mut_idx: int) -> Dict:
        sticky = {}
        self.log_step("Calculate Sticky Ends",
                      f"Calculating sticky ends for boundaries: first_mut_idx={first_mut_idx}, last_mut_idx={last_mut_idx}",
                      {"mutation_context_preview": mutated_ctx[:50] +
                          "...", "context_length": len(mutated_ctx)},
                      level=logging.DEBUG)
        for pos in sorted({first_mut_idx, last_mut_idx}):
            pos_sticky = {"top_strand": [], "bottom_strand": []}
            ranges = [
                range(pos - 3, pos + 1),
                range(pos - 2, pos + 2),
                range(pos - 1, pos + 3),
                range(pos, pos + 4)
            ]
            self.log_step("Calculate Sticky Ends",
                          f"Mutation boundary at position {pos}",
                          {"sticky_ranges": [list(r) for r in ranges]},
                          level=logging.DEBUG)
            for r in ranges:
                if 0 <= min(r) and max(r) < len(mutated_ctx):
                    top = "".join(mutated_ctx[i] for i in r)
                    bottom = self.utils.reverse_complement(top)
                    pos_sticky["top_strand"].append(
                        {"seq": top, "overhang_start_index": r.start})
                    pos_sticky["bottom_strand"].append(
                        {"seq": bottom, "overhang_start_index": r.start})
                    self.log_step("Sticky End",
                                  f"Generated sticky end for range {list(r)}",
                                  {"top_strand": {"seq": top, "overhang_start_index": r.start},
                                   "bottom_strand": {"seq": bottom, "overhang_start_index": r.start}},
                                  level=logging.DEBUG)
            sticky[f"position_{pos}"] = pos_sticky
        return sticky

    @DebugMixin.debug_wrapper
    def get_mutated_context(self, context_sequence: str, codon_context_position: int, new_codon_sequence: str, mutated_context_mutated_bases: tuple) -> tuple:
        self.log_step("get_mutated_context",
                      f"Starting codon swap at position {codon_context_position} with new codon '{new_codon_sequence}'")
        if len(new_codon_sequence) != 3:
            raise ValueError("New codon must be 3 nucleotides long.")
        if codon_context_position < 0 or codon_context_position + 3 > len(context_sequence):
            raise ValueError(
                "Invalid codon_context_position or context_sequence too short for the swap.")
        mutated_context = context_sequence[:codon_context_position] + \
            new_codon_sequence + context_sequence[codon_context_position + 3:]
        try:
            first_mutation_index = mutated_context_mutated_bases.index(1)
            last_mutation_index = len(
                mutated_context_mutated_bases) - 1 - mutated_context_mutated_bases[::-1].index(1)
        except ValueError:
            first_mutation_index = None
            last_mutation_index = None
        self.validate(isinstance(first_mutation_index, int),
                      "Cannot calculate mutated context, no mutated base found.")
        num_mutations = mutated_context_mutated_bases.count(1)
        self.validate((first_mutation_index == last_mutation_index and num_mutations == 1) or
                      (first_mutation_index !=
                       last_mutation_index and num_mutations > 1),
                      "Mutation indices do not match the expected number of mutated bases")
        self.log_step("get_mutated_context", "Codon swap complete",
                      {"mutated_context": mutated_context,
                       "first_mutation_index": first_mutation_index,
                       "last_mutation_index": last_mutation_index})
        mutated_context_first_mutation_index = first_mutation_index + codon_context_position
        mutated_context_last_mutation_index = last_mutation_index + codon_context_position
        return mutated_context, mutated_context_first_mutation_index, mutated_context_last_mutation_index

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
