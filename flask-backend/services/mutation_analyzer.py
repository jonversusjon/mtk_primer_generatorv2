import logging
from typing import Dict, List, Optional, Tuple, Union
from Bio.Seq import Seq
from .utils import GoldenGateUtils
from .mutation_optimizer import MutationOptimizer
from itertools import product
from config.logging_config import logger
from services.base import debug_context

class MutationAnalyzer:
    """
    Analyzes DNA sequences and generates possible mutations for Golden Gate assembly.
    Focuses on finding and evaluating potential mutation sites.
    """
    def __init__(
        self,
        sequence: Union[Seq, str],
        codon_usage_dict: Dict[str, Dict[str, float]],
        max_mutations: int = 1,
        verbose: bool = False
    ):
        self.logger = logger.getChild("MutationAnalyzer")
        self.utils = GoldenGateUtils()
        self.state = {
            'current_codon': '',
            'current_position': 0,
            'mutations_found': []
        }

        self.sequence = str(sequence) if isinstance(sequence, Seq) else sequence
        self.codon_usage_dict = codon_usage_dict
        self.max_mutations = max_mutations
        self.verbose = verbose

        if self.verbose:
            self.logger.setLevel(logging.DEBUG)

    def get_all_mutations(
        self,
        sites_to_mutate: List[Dict],
    ) -> Dict[str, List[Dict]]:  # Explicitly stating it returns a dict mapping to lists of dicts
        """
        Analyze and gather all possible mutation options for given restriction sites.

        Returns:
            Dictionary where keys are site identifiers and values are lists of mutation options.
        """
        try:
            mutation_options = {}
            with debug_context("mutation_analysis"):

                for site in sites_to_mutate:
                    frame = site["frame"]
                    site_mutations = {
                        "position": site["position"],
                        "sequence": site["sequence"],
                        "frame": frame,
                        "strand": site["strand"],
                        "enzyme": site["enzyme"],
                        "codons": site["codons"]
                    }

                    site_mutations["codons"] = [
                        {
                            "original_sequence": str(codon['codon_seq']),
                            "position": codon['position'],
                            "amino_acid": codon['amino_acid'],
                            "alternative_codons": self._find_alternative_codons(
                                codon['codon_seq'],
                                codon['amino_acid'],
                                codon_idx,
                                frame
                            )
                        }
                        for codon_idx, codon in enumerate(site['codons'])
                    ]

                    site_key = f"mutation_{site_mutations['position']}"
                    mutation_options[site_key] = site_mutations

                if self.verbose:
                    logger.info(f"Found {len(mutation_options)} sites with valid mutations")

                return mutation_options

        except Exception as e:
            logger.error(f"Error in mutation analysis: {str(e)}")
            return mutation_options  # Still return what was processed to avoid silent failure


    def _find_alternative_codons(
        self,
        original_codon: str,
        amino_acid: str,
        codon_idx: int,
        frame: int,
    ) -> List[str]:
        """
        Find alternative codons that disrupt the restriction site while staying within max_mutation.

        Rules:
        1. The candidate codon must not differ from the original codon in more than self.max_mutation bases.
        2. The candidate codon must change at least one base in the portion of the codon overlapping
            the restriction enzyme recognition site.

        Args:
            original_codon (str): The original codon (e.g., "GAG").
            amino_acid (str): The one-letter amino acid code (e.g., "E").
            codon_idx (int): The index of this codon within the recognition site (0 for first, etc.).
            frame (int): The reading frame of the site (0, 1, or 2).

        Returns:
            List[Dict[str, Any]]: A list of dictionaries containing alternative codon data.
        """

        alternative_codons = self.utils.get_codons_for_amino_acid(amino_acid)
        valid_alternatives = []

        # Get the indices (0, 1, or 2) within the codon that overlap with the recognition site.
        overlapping_indices = self.utils.get_recognition_site_bases(
            frame, codon_idx)

        for candidate in alternative_codons:
            # Skip candidate if it's identical to the original codon.
            if candidate == original_codon:
                continue

            # Rule 1: The candidate must not differ in more than self.max_mutation bases.
            differences = sum(1 for i in range(
                3) if candidate[i] != original_codon[i])
            if differences > self.max_mutations:
                continue

            # Rule 2: The candidate must change at least one base in the recognition site.
            if not any(candidate[i] != original_codon[i] for i in overlapping_indices):
                continue

            # Determine mutation tuple (0,1,0) format indicating changed bases.
            mutation_tuple = tuple(1 if candidate[i] != original_codon[i] else 0 for i in range(3))
            usage = self.utils.get_codon_usage(str(candidate), amino_acid, self.codon_usage_dict)
            
            valid_alternatives.append({
                "seq": candidate,
                "usage": usage,
                "mutations": mutation_tuple
            })

        if self.verbose:
            logger.debug(
                f"For input codon '{original_codon}' ({amino_acid}), found these alternative codons: {valid_alternatives}"
            )

        return valid_alternatives


    def _create_mutation_entry(self, site_sequence: str, site_start: int, codon_start: int, codon_pos: int,
                               original_codon: str, new_codon: str, frequency: float, enzyme: str, strand: str) -> Optional[Dict]:
        """
        Create a mutation entry and add it to the mutation dictionary.

        Returns:
            The created mutation entry if successful, None otherwise.
        """
        mutation_position = codon_start + codon_pos
        original_base = original_codon[codon_pos]
        new_base = new_codon[codon_pos]

        # Check if the mutation disrupts the restriction site
        mutated_site_sequence = site_sequence[:codon_pos] + \
            new_base + site_sequence[codon_pos+1:]
        if mutated_site_sequence == site_sequence:
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
        self.mutation_dict[mutation_id]=mutation_entry

        return mutation_entry
