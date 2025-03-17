import logging
from typing import Dict, List

from models import RestrictionSite, Codon, Mutation, OverhangOption
from .utils import GoldenGateUtils

from log_utils import logger

class MutationAnalyzer():
    """
    MutationAnalyzer Module

    This module provides functionality for analyzing DNA sequences to generate
    mutation options for Golden Gate assembly. The main entry and exit point
    is the get_all_mutations function.
    """

    def __init__(self,
                 codon_usage_dict: Dict[str, Dict[str, float]],
                 max_mutations: int = 1,
                 verbose: bool = False,
                 debug: bool = False):
        self.utils = GoldenGateUtils()
        self.state = {'current_codon': '',
                      'current_position': 0, 'mutations_found': []}
        self.codon_usage_dict = codon_usage_dict
        self.max_mutations = max_mutations
        self.verbose = verbose
        self.debug = debug

        logger.log_step("Initialization", f"Initializing MutationAnalyzer with verbose={verbose} and debug={debug}")
        if self.verbose:
            logger.log_step("Debug Mode", "Debug mode enabled for MutationAnalyzer")
            logger.validate(codon_usage_dict and isinstance(codon_usage_dict, dict),
                          "Codon usage dictionary is valid",
                          {"organisms": list(codon_usage_dict.keys())})
            logger.validate(isinstance(max_mutations, int) and max_mutations > 0,
                          f"Max mutations set to {max_mutations}",
                          {"valid_range": "1+"})
 
    def get_all_mutations(
        self,
        sites_to_mutate: List[RestrictionSite]
    ) -> Dict[str, Dict]:
        logger.log_step("Mutation Analysis", f"Starting mutation analysis for {len(sites_to_mutate)} site(s)")
        mutation_options = {}
        try:
            for site_idx, site in enumerate(sites_to_mutate):
                site_key = f"mutation_{site.position}"
                logger.log_step(
                    "Process Site",
                    f"Analyzing site {site_idx+1}/{len(sites_to_mutate)} at position {site.position}",
                    {"site_details": site}
                )
                site_has_alternatives = False
                valid_mutations = []

                for codon_idx, codon in enumerate(site.codons):
                    logger.log_step("Process Codon", f"Analyzing codon {codon_idx+1}/{len(site.codons)}: {codon.codon_sequence} at context position {codon.context_position}")
                    # Retrieve all synonymous codons for the given amino acid
                    alt_codon_seqs = self.utils.get_codon_seqs_for_amino_acid(codon.amino_acid)
                    for alt_codon_seq in alt_codon_seqs:
                        # Skip candidate if it is identical to the original codon
                        if alt_codon_seq == codon.codon_sequence:
                            logger.debug(f"Skipping alternative {alt_codon_seq} as it matches the original codon.")
                            continue

                        # Determine which positions differ between candidate and original
                        muts = [i for i in range(3) if alt_codon_seq[i] != codon.codon_sequence[i]]

                        # Skip candidate if number of mutations exceeds allowed max_mutations
                        if len(muts) > self.max_mutations:
                            logger.debug(f"Skipping alternative {alt_codon_seq} due to {len(muts)} mutations (max allowed: {self.max_mutations}).")
                            continue

                        # Reject the candidate if it only has mutations outside the recognition region
                        mutations_in_rs = list(set(muts) & set(codon.rs_overlap))
                        if not mutations_in_rs:
                            logger.debug(f"Skipping alternative {alt_codon_seq} as mutations {muts} fall outside recognition site bases {codon.rs_overlap}.")
                            continue

                        # If we're still here, it's a valid alternative
                        usage = self.utils.get_codon_usage(
                            alt_codon_seq, codon.amino_acid, self.codon_usage_dict
                        )
    
                        valid_alternative = Codon(
                            amino_acid=codon.amino_acid,
                            context_position=codon.context_position,
                            codon_sequence=alt_codon_seq,
                            rs_overlap=codon.rs_overlap,
                            usage=usage,
                        )
                        site_has_alternatives = True

                        mutated_context, first_mut, last_mut = self._get_mutated_context(
                            context_sequence=site.context_seq,
                            codon_context_position=site.context_first_base,
                            new_codon_sequence=valid_alternative.codon_sequence,
                            new_codon_mutated_bases=muts
                        )
                        logger.log_step("Mutated Context", f"Calculated mutated context with first mutation index {first_mut} and last mutation index {last_mut}")

                        overhang_options: List[OverhangOption] = self._calculate_sticky_ends_with_context(
                            mutated_context, first_mut, last_mut
                        )
                        logger.log_step("Sticky Ends", f"Calculated {len(overhang_options)} overhang option(s)")

                        valid_mutation = Mutation(
                            alt_codon=valid_alternative,
                            mut_indices_rs=mutations_in_rs,
                            mut_indices_codon=muts,
                            mut_context=mutated_context,
                            first_mut_idx=first_mut,
                            last_mut_idx=last_mut,
                            overhang_options=overhang_options,
                        )
                        valid_mutations.append(valid_mutation)
                        logger.log_step("Valid Mutation", f"Added valid mutation for codon {codon.codon_sequence}", {"alternative": alt_codon_seq, "mutations": muts})

                if site_has_alternatives:
                    mutation_options[site_key] = valid_mutations
                    logger.log_step("Site Completed", f"Site {site.position}: {len(valid_mutations)} valid mutation(s) found")
                else:
                    logger.log_step("No Alternatives Found", f"Site {site.position}: No alternative codons found", {"site": site.position}, level=logging.WARNING)
                    logger.debug(f"Site {site.position}: Mutation analysis skipped due to absence of alternatives.")

            logger.debug(f"Mutation options collected: {mutation_options}")
            if self.verbose:
                logger.log_step("Mutation Summary", f"Found {len(mutation_options)} site(s) with valid mutations")
            logger.validate(
                mutation_options,
                f"Generated mutation options for {len(mutation_options)} site(s)",
                {"site_keys": list(mutation_options.keys())}
            )
            return mutation_options

        except Exception as e:
            logger.error(f"Critical error in mutation analysis: {e}", exc_info=True)
            raise e

    def _calculate_sticky_ends_with_context(self, mutated_ctx: str,
                                            first_mut_idx: int,
                                            last_mut_idx: int) -> Dict:
        logger.log_step("Calculate Sticky Ends", "Calculating sticky ends for the mutated context")
        sticky = {}
        for pos in sorted({first_mut_idx, last_mut_idx}):
            pos_sticky = {"top_strand": [], "bottom_strand": []}
            ranges = [
                range(pos - 3, pos + 1),
                range(pos - 2, pos + 2),
                range(pos - 1, pos + 3),
                range(pos, pos + 4)
            ]
            for r in ranges:
                if 0 <= min(r) and max(r) < len(mutated_ctx):
                    top = "".join(mutated_ctx[i] for i in r)
                    bottom = self.utils.reverse_complement(top)
                    pos_sticky["top_strand"].append({
                        "seq": top,
                        "overhang_start_index": r.start
                    })
                    pos_sticky["bottom_strand"].append({
                        "seq": bottom,
                        "overhang_start_index": r.start
                    })
            sticky[f"position_{pos}"] = pos_sticky
            logger.log_step("Sticky Ends Calculated", f"Calculated sticky ends for position {pos} with {len(pos_sticky['top_strand'])} option(s)")
        return sticky

    def _get_mutated_context(
        self,
        context_sequence: str,
        codon_context_position: int,
        new_codon_sequence: str,
        new_codon_mutated_bases: tuple) -> tuple:
        logger.log_step("Mutated Context Start", f"Generating mutated context for codon swap at position {codon_context_position}")
        if len(new_codon_sequence) != 3:
            raise ValueError("New codon must be 3 nucleotides long.")
        if (codon_context_position < 0 or codon_context_position + 3 > len(context_sequence)):
            raise ValueError("Invalid codon_context_position or context_sequence too short for the swap.")
        mutated_context = (context_sequence[:codon_context_position] +
                           new_codon_sequence +
                           context_sequence[codon_context_position + 3:])
        try:
            first_mutation_index = new_codon_mutated_bases.index(1)
            last_mutation_index = (len(new_codon_mutated_bases) - 1 -
                                   new_codon_mutated_bases[::-1].index(1))
        except ValueError:
            first_mutation_index = None
            last_mutation_index = None

        mutated_context_first_mutation_index = first_mutation_index + codon_context_position if first_mutation_index is not None else None
        mutated_context_last_mutation_index = last_mutation_index + codon_context_position if last_mutation_index is not None else None
        logger.log_step("Mutated Context Completed", f"Mutated context generated with first mutation index {mutated_context_first_mutation_index} and last mutation index {mutated_context_last_mutation_index}")
        return (mutated_context, mutated_context_first_mutation_index, mutated_context_last_mutation_index)