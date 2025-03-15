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
                    • "codon_original_sequence": The original codon sequence.
                    • "context_position": The codon position in the context sequence.
                    • "amino_acid": The encoded amino acid.
                    • "alternative_codons": A list of dictionaries for each valid alternative codon, each containing:
                            - "seq": The alternative codon sequence.
                            - "usage": The codon usage frequency.
                            - "mutations": A 3-tuple of integers indicating changes (1 for a mutation, 0 otherwise).
                            - "changes_in_site": A list of codon index positions (0, 1, or 2) where the candidate alters the recognition region.
                            - "sticky_ends": A dictionary keyed by "position_{codon_index}" that contains:
                                    • "top_strand": A list of dictionaries, each with:
                                            - "seq": The overhang sequence.
                                            - "overhang_start_index": The starting index in the context.
                                    • "bottom_strand": A list of dictionaries with corresponding reverse-complement data.
                            - "mutated_context": The context sequence with the codon replaced by the alternative codon.
                            - "mutation_positions_in_context": A list of indices in the context where the mutation occurs.
    """

    @DebugMixin.debug_wrapper
    def get_all_mutations(self, sites_to_mutate: List[Dict]) -> Dict[str, Dict]:
        self.validate(sites_to_mutate and isinstance(sites_to_mutate, list),
                      f"Received {len(sites_to_mutate)} sites to mutate",
                      {"first_site": sites_to_mutate[0] if sites_to_mutate else None})
        self.log_step(
            "Analysis Start", f"Starting mutation analysis for {len(sites_to_mutate)} sites")

        mutation_options = {}
        try:
            with debug_context("mutation_analysis"):
                for site_idx, site in enumerate(sites_to_mutate):
                    site_key = f"mutation_{site['position']}"
                    self.log_step("Process Site",
                                  f"Analyzing site {site_idx+1}/{len(sites_to_mutate)}",
                                  {"position": site["position"],
                                   "sequence": site["sequence"],
                                   "enzyme": site["enzyme"],
                                   "codons": len(site["codons"])})

                    # Assemble the site-level data.
                    site_mutations = {
                        "position": site["position"],
                        "sequence": site["sequence"],
                        "frame": site["frame"],
                        "strand": site["strand"],
                        "enzyme": site["enzyme"],
                        "codons": []
                    }

                    site_has_alternatives = False

                    for codon_idx, codon in enumerate(site["codons"]):
                        print(f" *** codon: {codon}")
                        # Re-key "codon_seq" to "codon_sequence"
                        codon["codon_sequence"] = codon.pop("codon_seq", None)

                        if codon["codon_sequence"] is None:
                            self.logger.warning(
                                f"Codon sequence missing at site {site['position']}. Skipping.")
                            continue  # Skip this codon if data is missing

                        self.log_step("Process Codon",
                                      f"Analyzing codon {codon_idx+1}/{len(site['codons'])}",
                                      {"codon_sequence": codon["codon_sequence"],
                                       "context_position": codon["context_position"],
                                       "amino_acid": codon["amino_acid"]})

                        # Inject overlapping indices
                        site_with_indices = {
                            **site,
                            "overlapping_indices": self.utils.get_recognition_site_bases(site["frame"], codon_idx)
                        }

                        alternatives = self._find_alternative_codons(
                            site_with_indices,
                            codon["codon_sequence"],
                            codon["amino_acid"],
                            codon["context_position"],
                        )

                        if alternatives:
                            site_has_alternatives = True
                            self.log_step("Found Alternatives",
                                          f"Found {len(alternatives)} alternative codons for {codon['codon_sequence']}",
                                          {"alternatives": alternatives})

                            for alt in alternatives:
                                mutated_context, mutated_context_first_mutation_index, mutated_context_last_mutation_index = self.get_mutated_context(
                                    context_sequence=site["context_sequence"],
                                    codon_context_position=codon["context_position"],
                                    new_codon_sequence=alt["seq"],
                                    mutated_context_mutated_bases=alt["mutations"],
                                )
                                alt["mutated_context"] = mutated_context
                                alt["sticky_ends"] = self._calculate_sticky_ends_with_context(
                                    mutated_context,
                                    mutated_context_first_mutation_index,
                                    mutated_context_last_mutation_index
                                )

                            site_mutations["codons"].append({
                                "original_codon_sequence": codon["codon_sequence"],
                                "context_position": codon["context_position"],
                                "amino_acid": codon["amino_acid"],
                                "alternative_codons": alternatives
                            })

                    if site_has_alternatives:
                        try:
                            from models.mtk import MutationSite
                            validated_site = MutationSite.model_validate(
                                site_mutations)
                            mutation_options[site_key] = validated_site.model_dump(
                            )
                        except Exception as e:
                            self.logger.error(
                                f"Validation error for site {site['position']}: {e}")
                            continue  # Log error and skip adding this site

                        self.log_step("Site Result",
                                      f"Completed analysis for site at position {site['position']}",
                                      {"alternatives_found": sum(len(c["alternative_codons"]) for c in site_mutations["codons"])})

                    else:
                        self.log_step("No Alternatives Found",
                                      f"No alternative codons found for site at position {site['position']}",
                                      {"site": site["position"]},
                                      level=logging.WARNING)

                if self.verbose:
                    logger.info(
                        f"Found {len(mutation_options)} sites with valid mutations")

                self.validate(mutation_options,
                              f"Generated mutation options for {len(mutation_options)} sites",
                              {"site_keys": list(mutation_options.keys())})

                return mutation_options

        except Exception as e:
            self.logger.error(f"Critical error in mutation analysis: {e}")
            raise e  # Let the error propagate so debugging can happen

    @DebugMixin.debug_wrapper
    def _find_alternative_codons(self, site: Dict, original_codon: str, amino_acid: str, codon_context_position: int) -> List[Dict]:
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

            # Always use the provided codon_context_position to calculate the mutation positions.
            mutation_positions_in_context = [
                codon_context_position + i for i in muts]

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
            candidate_info = {"seq": cand,
                              "usage": usage,
                              "mutations": mut_tuple,
                              "changes_in_site": changes,
                              "mutation_positions_in_context": mutation_positions_in_context}
            valid.append(candidate_info)
            self.log_step("Valid Alternative",
                          f"Candidate {cand} accepted",
                          {"mutations": mut_tuple, "changes_in_site": changes},
                          level=logging.DEBUG)
        self.validate(valid,
                      f"Generated {len(valid)} valid alternatives",
                      {"original": original_codon, "amino_acid": amino_acid})
        return valid

    def _calculate_sticky_ends_with_context(self, mutated_ctx: str, first_mut_idx: int, last_mut_idx: int) -> Dict:
        sticky = {}
        print(f"first_mut_idx: {first_mut_idx}, last_mut_idx: {last_mut_idx}")
        # Use a set to ensure unique positions (if only one mutation exists, for example)
        for pos in sorted({first_mut_idx, last_mut_idx}):
            pos_sticky = {"top_strand": [], "bottom_strand": []}
            # For each boundary, generate 4 possible 4nt windows (offsets mimic bsai and bsmbi sticky end options)
            ranges = [
                range(pos - 3, pos + 1),
                range(pos - 2, pos + 2),
                range(pos - 1, pos + 3),
                range(pos, pos + 4)
            ]
            self.log_step("Calculate Sticky Ends",
                          f"Calculating sticky ends for mutation boundary at position {pos}",
                          {"mutation_pos_in_context": pos,
                           "sticky_positions": [list(r) for r in ranges]},
                          level=logging.DEBUG)
            for r in ranges:
                # Check that the range is within bounds of mutated_ctx
                if 0 <= min(r) and max(r) < len(mutated_ctx):
                    # Since pos is a boundary of the mutated block, every window here includes at least one mutated base.
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
    def get_mutated_context(self,
                            context_sequence: str,
                            codon_context_position: int,
                            new_codon_sequence: str,
                            mutated_context_mutated_bases: tuple) -> tuple:
        """
        Replace the codon in the context_sequence at the specified codon_context_position with new_codon_sequence.

        Parameters:
            context_sequence (str): The original context sequence.
            codon_context_position (int): The starting index in the context_sequence where the codon is located.
            new_codon_sequence (str): The new codon sequence (must be exactly 3 nucleotides).
            mutated_context_mutated_bases (tuple): A tuple indicating mutation positions (1 for mutation, 0 for no mutation).

        Returns:
            tuple: The mutated context sequence, index of the first mutation, index of the last mutation.

        Raises:
            ValueError: If new_codon_sequence is not exactly 3 nucleotides long or if the position is invalid.
        """
        self.log_step("get_mutated_context",
                      f"Starting codon swap at position {codon_context_position} with new codon '{new_codon_sequence}'")

        if len(new_codon_sequence) != 3:
            raise ValueError("New codon must be 3 nucleotides long.")

        if codon_context_position < 0 or codon_context_position + 3 > len(context_sequence):
            raise ValueError(
                "Invalid codon_context_position or context_sequence too short for the swap.")

        mutated_context = (context_sequence[:codon_context_position] +
                           new_codon_sequence +
                           context_sequence[codon_context_position + 3:])

        # Compute mutation indices
        try:
            first_mutation_index = mutated_context_mutated_bases.index(
                1)
            last_mutation_index = len(
                mutated_context_mutated_bases) - 1 - mutated_context_mutated_bases[::-1].index(1)
        except ValueError:
            first_mutation_index = None
            last_mutation_index = None

        # Validate mutation indices
        self.validate(isinstance(first_mutation_index, int),
                      "cannot calculate mutated context, no mutated base")

        num_mutations = mutated_context_mutated_bases.count(1)
        self.validate(
            (first_mutation_index == last_mutation_index and num_mutations == 1) or
            (first_mutation_index !=
             last_mutation_index and num_mutations > 1),
            "Mutation indices do not match the expected number of mutated bases"
        )

        self.log_step("get_mutated_context", "Codon swap complete",
                      data={"mutated_context": mutated_context,
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
