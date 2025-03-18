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
                            "Codon usage dictionary is valid")
            logger.validate(isinstance(max_mutations, int) and max_mutations > 0,
                            f"Max mutations set to {max_mutations}",
                            {"valid_range": "1+"})

    def get_all_mutations(self,
                          sites_to_mutate: List[RestrictionSite]) -> Dict[str, Dict]:
        logger.log_step("Mutation Analysis", f"Starting mutation analysis for {len(sites_to_mutate)} site(s)")
        mutation_options = {}
        try:
            for site_idx, site in enumerate(sites_to_mutate):
                site_key = f"mutation_{site.position}"
                logger.log_step("Process Site",
                                f"Analyzing site {site_idx+1}/{len(sites_to_mutate)} at position {site.position}",
                                {"site_details": site})

                site_has_alternatives = False
                valid_mutations = []

                for codon_idx, codon in enumerate(site.codons):
                    logger.log_step("Process Codon",
                                    f"Analyzing codon {codon_idx+1}/{len(site.codons)}: {codon.codon_sequence} at context position {codon.context_position}")
                    # Retrieve all synonymous codon sequences for the given amino acid.
                    alt_codon_seqs = self.utils.get_codon_seqs_for_amino_acid(codon.amino_acid)
                    for alt_codon_seq in alt_codon_seqs:
                        # Skip candidate if identical to original codon
                        if alt_codon_seq == codon.codon_sequence:
                            continue

                        # Identify the positions where the alternative codon differs from the original.
                        muts = [i for i in range(3) if alt_codon_seq[i] != codon.codon_sequence[i]]

                        # Skip if the number of mutations exceeds allowed maximum.
                        if len(muts) > self.max_mutations:
                            logger.debug(f"Skipping alternative {alt_codon_seq} due to {len(muts)} mutations (max allowed: {self.max_mutations}).")
                            continue

                        # Ensure there is at least one mutation within the recognition site.
                        mutations_in_rs = list(set(muts) & set(codon.rs_overlap))
                        if not mutations_in_rs:
                            logger.debug(f"Skipping alternative {alt_codon_seq} as mutations {muts} fall outside recognition site bases {codon.rs_overlap}.")
                            continue

                        # Retrieve codon usage information.
                        usage = self.utils.get_codon_usage(alt_codon_seq, codon.amino_acid, self.codon_usage_dict)
                        valid_alternative = Codon(
                            amino_acid=codon.amino_acid,
                            context_position=codon.context_position,
                            codon_sequence=alt_codon_seq,
                            rs_overlap=codon.rs_overlap,
                            usage=usage,
                        )
                        site_has_alternatives = True

                        logger.debug(f"context_sequence: {site.context_seq}")
                        logger.debug(f"codon_context_position: {codon.context_position}")
                        logger.debug(f"new_codon_sequence: {valid_alternative.codon_sequence}")
                        logger.debug(f"new_codon_mutated_bases: {muts}")

                        # Calculate the mutated context and mutation indices.
                        mutated_context, first_mut, last_mut = self._get_mutated_context(
                            context_sequence=site.context_seq,
                            codon_context_position=codon.context_position,
                            new_codon_sequence=valid_alternative.codon_sequence,
                            new_codon_mutated_bases=muts
                        )
                        logger.log_step("Mutated Context", f"Calculated mutated context with first mutation index {first_mut} and last mutation index {last_mut}")

                        # Calculate sticky ends and map the returned dictionary into a list of OverhangOption objects.
                        overhang_options_raw = self._calculate_sticky_ends_with_context(mutated_context, first_mut, last_mut)
                        overhang_options = []
                        if isinstance(overhang_options_raw, dict):
                            # Iterate over each position's sticky options.
                            for pos_key, pos_options in overhang_options_raw.items():
                                top_options = pos_options.get("top_strand", [])
                                bottom_options = pos_options.get("bottom_strand", [])
                                # Use zip to pair each top and bottom option.
                                for top_option, bottom_option in zip(top_options, bottom_options):
                                    mapped_option = {
                                        "top_overhang": top_option.get("seq", ""),
                                        "bottom_overhang": bottom_option.get("seq", ""),
                                        "overhang_start_index": top_option.get("overhang_start_index")
                                    }
                                    if mapped_option["overhang_start_index"] is None:
                                        raise ValueError("overhang_start_index is missing in the overhang option")
                                    overhang_options.append(OverhangOption(**mapped_option))
                        elif isinstance(overhang_options_raw, list):
                            overhang_options = overhang_options_raw
                        else:
                            raise ValueError("Unexpected type returned for overhang options.")

                        logger.log_step("Sticky Ends", f"Calculated {len(overhang_options)} overhang option(s)")

                        # Create the Mutation instance ensuring that alt_codons is a list.
                        valid_mutation = Mutation(
                            alt_codons=[valid_alternative],
                            mut_indices_rs=mutations_in_rs,
                            mut_indices_codon=muts,
                            mut_context=mutated_context,
                            first_mut_idx=first_mut,
                            last_mut_idx=last_mut,
                            overhang_options=overhang_options,
                        )
                        valid_mutations.append(valid_mutation)
                        logger.log_step("Valid Mutation",
                                       f"Added valid mutation for codon {codon.codon_sequence}",
                                       {"alternative": alt_codon_seq, "mutations": muts})

                if site_has_alternatives:
                    mutation_options[site_key] = valid_mutations
                    logger.log_step("Site Completed", f"Site {site.position}: {len(valid_mutations)} valid mutation(s) found")
                else:
                    logger.log_step("No Alternatives Found",
                                    f"Site {site.position}: No alternative codons found",
                                    {"site": site.position}, level=logging.WARNING)
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

    def _calculate_sticky_ends_with_context(self, mutated_ctx: str, first_mut_idx: int, last_mut_idx: int) -> Dict:
        logger.log_step("Calculate Sticky Ends", "Calculating sticky ends for the mutated context")
        logger.debug(f"first_mut_idx: {first_mut_idx}")
        logger.debug(f"last_mut_idx: {last_mut_idx}")
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

    def _get_mutated_context(self,
                             context_sequence: str,
                             codon_context_position: int,
                             new_codon_sequence: str,
                             new_codon_mutated_bases: List[int]) -> tuple:
        """
        Generate a mutated context sequence by swapping the codon at the specified position and
        calculate the first and last mutation indices within the overall context.

        Args:
            context_sequence (str): The full context sequence.
            codon_context_position (int): The starting index of the codon in the context sequence.
            new_codon_sequence (str): The new codon sequence (must be exactly 3 nucleotides).
            new_codon_mutated_bases (List[int]): List of indices (0-indexed) of the codon bases that have been mutated.

        Returns:
            tuple: A tuple containing:
                - mutated_context (str): The updated context sequence.
                - mutated_context_first_mutation_index (int): Global position of the first mutated base.
                - mutated_context_last_mutation_index (int): Global position of the last mutated base.

        Raises:
            ValueError: If the new codon is not 3 nucleotides long, if the codon_context_position is invalid,
                        or if no mutated bases are provided.
        """
        logger.log_step("Mutated Context Start", f"Generating mutated context for codon swap at position {codon_context_position}")
        if len(new_codon_sequence) != 3:
            raise ValueError("New codon must be 3 nucleotides long.")
        if codon_context_position < 0 or codon_context_position + 3 > len(context_sequence):
            raise ValueError("Invalid codon_context_position or context_sequence too short for the swap.")

        mutated_context = (
            context_sequence[:codon_context_position] +
            new_codon_sequence +
            context_sequence[codon_context_position + 3:]
        )
        
        # Ensure that we have at least one mutation index.
        if not new_codon_mutated_bases:
            raise ValueError("No mutated bases provided in new_codon_mutated_bases.")
        
        # Calculate the first and last mutation indices relative to the codon.
        first_mutation_index = min(new_codon_mutated_bases)
        last_mutation_index = max(new_codon_mutated_bases)
        
        logger.debug(f"First mutation index (relative to codon): {first_mutation_index}")
        logger.debug(f"Last mutation index (relative to codon): {last_mutation_index}")
        
        # Convert these indices to positions in the full context sequence.
        mutated_context_first_mutation_index = first_mutation_index + codon_context_position
        mutated_context_last_mutation_index = last_mutation_index + codon_context_position
        logger.log_step("Mutated Context Completed",
                        f"Mutated context generated with first mutation index {mutated_context_first_mutation_index} and last mutation index {mutated_context_last_mutation_index}")
        return (mutated_context, mutated_context_first_mutation_index, mutated_context_last_mutation_index)
