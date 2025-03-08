import logging
from typing import Dict, List, Optional
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
    Analyzes DNA sequences and generates possible mutations for Golden Gate assembly.
    Focuses on finding and evaluating potential mutation sites.
    """

    def __init__(
        self,
        codon_usage_dict: Dict[str, Dict[str, float]],
        max_mutations: int = 1,
        verbose: bool = False,
        debug: bool = False
    ):
        """
        Initialize the MutationAnalyzer with given parameters and debugging options.
        """
        self.logger = logger.getChild("MutationAnalyzer")
        self.utils = GoldenGateUtils()
        self.state = {
            'current_codon': '',
            'current_position': 0,
            'mutations_found': []
        }
        self.codon_usage_dict = codon_usage_dict
        self.max_mutations = max_mutations
        self.verbose = verbose
        self.debug = debug

        if self.verbose:
            self.logger.setLevel(logging.DEBUG)

        if self.debug:
            self.debugger = MutationDebugger(
                parent_logger=logger,
                use_custom_format=True
            )
            # Avoid duplicate logs
            if hasattr(self.debugger.logger, 'propagate'):
                self.debugger.logger.propagate = False

            self.logger.info("Debug mode enabled for MutationAnalyzer")

            # Validate initialization parameters using mixin's validate
            self.validate(
                isinstance(codon_usage_dict, dict) and len(
                    codon_usage_dict) > 0,
                "Codon usage dictionary is valid",
                {"organisms": list(codon_usage_dict.keys())}
            )
            self.validate(
                isinstance(max_mutations, int) and max_mutations > 0,
                f"Max mutations set to {max_mutations}",
                {"valid_range": "1-3"}
            )

    @DebugMixin.debug_wrapper
    def get_all_mutations(
        self,
        sequence: str,
        sites_to_mutate: List[Dict],
    ) -> Dict[str, List[Dict]]:
        """
        Analyze and gather all possible mutation options for given restriction sites.
        """
        self.validate(
            isinstance(sites_to_mutate, list) and len(sites_to_mutate) > 0,
            f"Received {len(sites_to_mutate)} sites to mutate",
            {"first_site": sites_to_mutate[0] if sites_to_mutate else None}
        )
        self.log_step(
            "Analysis Start",
            f"Starting mutation analysis for {len(sites_to_mutate)} sites"
        )

        mutation_options = {}
        try:
            with debug_context("mutation_analysis"):
                for site_idx, site in enumerate(sites_to_mutate):
                    self.log_step(
                        "Process Site",
                        f"Analyzing site {site_idx+1}/{len(sites_to_mutate)}",
                        {
                            "position": site["position"],
                            "sequence": site["sequence"],
                            "enzyme": site["enzyme"],
                            "codons": len(site["codons"])
                        }
                    )

                    frame = site["frame"]
                    site_mutations = {
                        "position": site["position"],
                        "sequence": site["sequence"],
                        "frame": frame,
                        "strand": site["strand"],
                        "enzyme": site["enzyme"],
                        "codons": site["codons"]
                    }

                    processed_codons = []
                    for codon_idx, codon in enumerate(site['codons']):
                        self.log_step(
                            "Process Codon",
                            f"Analyzing codon {codon_idx+1}/{len(site['codons'])}",
                            {
                                "codon_seq": codon['codon_seq'],
                                "position": codon['position'],
                                "amino_acid": codon['amino_acid']
                            }
                        )
                        alternative_codons = self._find_alternative_codons_with_sequence(
                            site,
                            codon['codon_seq'],
                            codon['amino_acid'],
                            codon_idx,
                            frame,
                            codon['position'] - site["position"],
                        )
                        self.validate(
                            len(alternative_codons) > 0,
                            f"Found {len(alternative_codons)} alternative codons",
                            {"alternatives": alternative_codons}
                        )

                        processed_codons.append({
                            "original_sequence": str(codon['codon_seq']),
                            "position": codon['position'],
                            "amino_acid": codon['amino_acid'],
                            "alternative_codons": alternative_codons
                        })

                    site_mutations["codons"] = processed_codons
                    site_key = f"mutation_{site_mutations['position']}"
                    mutation_options[site_key] = site_mutations

                    self.log_step(
                        "Site Result",
                        f"Completed analysis for site at position {site['position']}",
                        {
                            "alternatives_found": sum(
                                len(codon["alternative_codons"])
                                for codon in processed_codons
                            )
                        }
                    )

                if self.verbose:
                    logger.info(
                        f"Found {len(mutation_options)} sites with valid mutations"
                    )
                self.validate(
                    len(mutation_options) > 0,
                    f"Generated mutation options for {len(mutation_options)} sites",
                    {"site_keys": list(mutation_options.keys())}
                )

                return mutation_options

        except Exception as e:
            error_msg = f"Error in mutation analysis: {str(e)}"
            logger.error(error_msg)
            if getattr(self, 'debugger', None):
                self.debugger.log_error(error_msg)
            return mutation_options  # Return what was processed to avoid silent failure

    def _find_alternative_codons_with_sequence(
        self,
        site: Dict,
        original_codon: str,
        amino_acid: str,
        codon_idx: int,
        frame: int,
        codon_offset: int
    ) -> List[Dict]:
        """
        Find alternative codons that disrupt the restriction site while staying within max_mutation.
        """
        context_sequence = site["context"]
        site_sequence = site["sequence"]
        
        # The site["context_mutated_indices"] currently represents the position of the 
        # recognition site in the context sequence, not the specific bases that will change.
        # We'll keep using it for locating the codon but rename it for clarity.
        recognition_site_indices = site["context_mutated_indices"]
        
        self.log_step(
            "Find Alternatives",
            f"Finding alternatives for codon {original_codon} (AA: {amino_acid})",
            {
                "codon_idx": codon_idx,
                "frame": frame,
                "max_mutations": self.max_mutations,
                "site_sequence": site_sequence,
                "codon_offset": codon_offset,
                "context_sequence": context_sequence,
                "recognition_site_indices": recognition_site_indices
            }
        )

        alternative_codons = self.utils.get_codons_for_amino_acid(amino_acid)

        self.log_step(
            "Possible Codons",
            f"Found {len(alternative_codons)} synonymous codons for {amino_acid}",
            {"codons": alternative_codons}
        )

        valid_alternatives = []
        overlapping_indices = self.utils.get_recognition_site_bases(
            frame, codon_idx)

        self.log_step(
            "Overlapping Indices",
            "Codon positions that overlap with recognition site",
            {"indices": overlapping_indices}
        )

        # Calculate the codon start position in the context
        if codon_offset >= len(site_sequence):
            codon_start_in_context = recognition_site_indices[-1] + (
                codon_offset - len(site_sequence) + 1)
        elif codon_offset < 0:
            codon_start_in_context = recognition_site_indices[0] + codon_offset
        else:
            codon_start_in_context = recognition_site_indices[codon_offset]

        for candidate in alternative_codons:
            if candidate == original_codon:
                continue

            # Find which positions differ between original and candidate codon
            mutations_list = [i for i in range(3) if candidate[i] != original_codon[i]]
            differences = len(mutations_list)

            if differences > self.max_mutations:
                self.log_step(
                    "Skip Candidate",
                    f"Candidate {candidate} has too many mutations",
                    {"differences": differences, "max_allowed": self.max_mutations},
                    level=logging.DEBUG
                )
                continue

            changes_in_site = [
                i for i in overlapping_indices if candidate[i] != original_codon[i]]
            if not changes_in_site:
                self.log_step(
                    "Skip Candidate",
                    f"Candidate {candidate} doesn't change the recognition site",
                    {"overlapping_indices": overlapping_indices},
                    level=logging.DEBUG
                )
                continue

            mutation_tuple = tuple(
                1 if candidate[i] != original_codon[i] else 0 for i in range(3)
            )
            usage = self.utils.get_codon_usage(
                str(candidate), amino_acid, self.codon_usage_dict)

            # Create a copy of the context sequence and apply the mutations
            mutated_context = list(context_sequence)
            
            # THIS IS THE KEY CHANGE - Calculate the actual mutated indices in the context
            # These are the indices that will change when swapping the codon
            context_mutated_indices = []
            for i in mutations_list:
                pos_in_context = codon_start_in_context + i
                if 0 <= pos_in_context < len(mutated_context):
                    mutated_context[pos_in_context] = candidate[i]
                    context_mutated_indices.append(pos_in_context)

            mutated_context = ''.join(mutated_context)
            
            # Log the specific mutation positions
            self.log_step(
                "Calculate Mutation Positions",
                f"Identified mutated indices for candidate {candidate}",
                {"context_mutated_indices": context_mutated_indices},
                level=logging.DEBUG
            )
            
            sticky_ends = self._calculate_sticky_ends_with_context(
                context_sequence,
                mutated_context,
                codon_start_in_context,
                mutations_list
            )

            valid_alternative = {
                "seq": candidate,
                "usage": usage,
                "mutations": mutation_tuple,
                "sticky_ends": sticky_ends,
                "mutation_positions_in_context": context_mutated_indices  # Use a new key to avoid confusion
            }

            valid_alternatives.append(valid_alternative)

            self.log_step(
                "Valid Alternative",
                f"Found valid alternative: {candidate}",
                {
                    "usage": usage,
                    "mutations": mutation_tuple,
                    "changes_in_site": changes_in_site,
                    "sticky_ends": sticky_ends,
                    "context_mutated_indices": context_mutated_indices
                },
                level=logging.DEBUG
            )

        if self.verbose:
            logger.debug(
                f"For input codon '{original_codon}' ({amino_acid}), found {len(valid_alternatives)} alternative codons"
            )
        self.validate(
            len(valid_alternatives) > 0,
            f"Generated {len(valid_alternatives)} valid alternatives",
            {"original": original_codon, "amino_acid": amino_acid}
        )
        return valid_alternatives
    def _calculate_sticky_ends_with_context(
        self,
        original_context: str,
        mutated_context: str,
        codon_start_in_context: int,
        mutation_positions: List[int]
    ) -> Dict:
        """
        Calculate all possible sticky end sequences for a given mutation using the full context.
        """
        sticky_ends = {}

        for pos in mutation_positions:
            mutation_pos_in_context = codon_start_in_context + pos
            position_sticky_ends = {
                "top_strand": [],
                "bottom_strand": []
            }

            sticky_positions = [
                range(mutation_pos_in_context - 3,
                      mutation_pos_in_context + 1),
                range(mutation_pos_in_context - 2,
                      mutation_pos_in_context + 2),
                range(mutation_pos_in_context - 1,
                      mutation_pos_in_context + 3),
                range(mutation_pos_in_context, mutation_pos_in_context + 4)
            ]

            self.log_step(
                "Calculate Sticky Ends",
                f"Calculating sticky ends for mutation at position {pos}",
                {
                    "mutation_pos_in_context": mutation_pos_in_context,
                    "sticky_positions": [list(r) for r in sticky_positions]
                },
                level=logging.DEBUG
            )

            for sticky_range in sticky_positions:
                if min(sticky_range) >= 0 and max(sticky_range) < len(mutated_context):
                    top_strand = ''.join(
                        mutated_context[i] for i in sticky_range)
                    bottom_strand = self.utils.reverse_complement(top_strand)
                    position_sticky_ends["top_strand"].append(top_strand)
                    position_sticky_ends["bottom_strand"].append(bottom_strand)

                    self.log_step(
                        "Sticky End",
                        f"Generated sticky end for range {list(sticky_range)}",
                        {"top_strand": top_strand, "bottom_strand": bottom_strand},
                        level=logging.DEBUG
                    )

            sticky_ends[f"position_{pos}"] = position_sticky_ends

        return sticky_ends

    @DebugMixin.debug_wrapper
    def _create_mutation_entry(
        self,
        site_sequence: str,
        site_start: int,
        codon_start: int,
        codon_pos: int,
        original_codon: str,
        new_codon: str,
        frequency: float,
        enzyme: str,
        strand: str
    ) -> Optional[Dict]:
        """
        Create a mutation entry and add it to the mutation dictionary.
        """
        self.log_step(
            "Create Mutation Entry",
            f"Creating mutation entry for {original_codon} → {new_codon}",
            {
                "site_sequence": site_sequence,
                "site_start": site_start,
                "codon_start": codon_start,
                "codon_pos": codon_pos,
                "frequency": frequency,
                "enzyme": enzyme,
                "strand": strand
            }
        )

        mutation_position = codon_start + codon_pos
        original_base = original_codon[codon_pos]
        new_base = new_codon[codon_pos]
        mutated_site_sequence = site_sequence[:codon_pos] + \
            new_base + site_sequence[codon_pos+1:]

        self.log_step(
            "Mutation Check",
            "Checking if mutation disrupts restriction site",
            {
                "original_site": site_sequence,
                "mutated_site": mutated_site_sequence,
                "changed_base_pos": codon_pos,
                "original_base": original_base,
                "new_base": new_base
            }
        )

        if mutated_site_sequence == site_sequence:
            if getattr(self, 'debugger', None):
                self.debugger.log_warning(
                    f"Mutation does not disrupt restriction site: {original_base} → {new_base} at position {mutation_position}"
                )
            return None

        mutation_entry = {
            "position": mutation_position,
            "original_base": original_base,
            "new_base": new_base,
            "original_codon": original_codon,
            "new_codon": new_codon,
            "frequency": frequency,
            "site_start": site_start,
            "enzyme": enzyme,
            "strand": strand
        }

        mutation_id = f"{mutation_position}_{original_base}_{new_base}"

        if not hasattr(self, 'mutation_dict'):
            self.mutation_dict = {}

        self.mutation_dict[mutation_id] = mutation_entry

        self.validate(
            mutation_id in self.mutation_dict,
            f"Successfully created mutation entry: {mutation_id}",
            {"entry": mutation_entry}
        )

        return mutation_entry
