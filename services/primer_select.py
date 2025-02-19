# services/primer_select.py
from Bio.Seq import Seq
import numpy as np
from typing import List, Dict, Optional, Union
import logging
from config.logging_config import logger
from services.base import debug_context


class PrimerSelector:
    """
    Selects the best primers based on scoring criteria.
    """

    def __init__(self, verbose: bool = False):
        self.logger = logger.getChild(
            "PrimerSelector")

        """
        Initialize the primer selection process.
        
        Args:
            verbose (bool): If True, provide more user-friendly logs in production.
        """
        self.verbose = verbose

        self.state = {
            'current_operation': '',
            'best_score': float('-inf'),
            'primers_evaluated': 0
        }

        self.logger.debug("PrimerSelector initialized.")

        if self.verbose:
            self.logger.info("PrimerSelector is running in verbose mode.")

    def format_primers_for_output(
        self,
        primer_sets: List[List[Seq]]
    ) -> List[List[str]]:
        """Formats primer sets for output in the final protocol."""
        with debug_context("format_primers"):
            primer_data = []
            try:
                for i, primers in enumerate(primer_sets):
                    for j, primer in enumerate(primers):
                        primer_entry = [
                            f"Primer_{i+1}_{j+1}",
                            str(primer),
                            f"Amplicon_{i+1}"
                        ]
                        primer_data.append(primer_entry)
                        logger.debug(f"Formatted primer: {primer_entry[0]}")

                return primer_data

            except Exception as e:
                logger.error(f"Error formatting primers: {str(e)}")
                raise

    def select_best_internal_primers(
        self,
        primer_sets: Dict,
        template_seq: Optional[str] = None
    ) -> Optional[Dict]:
        """Selects optimal internal primers based on multiple criteria."""
        with debug_context("select_best_primers"):
            if not primer_sets:
                logger.warning("No primer sets provided!")
                return None

            self.state['best_score'] = float('-inf')
            self.state['primers_evaluated'] = 0
            best_primers = None

            try:
                for site_index, codon_start_dict in primer_sets.items():
                    best_primers = self._evaluate_site_primers(
                        site_index,
                        codon_start_dict,
                        best_primers
                    )

                if best_primers:
                    logger.info(
                        f"Selected best primers with score {self.state['best_score']}"
                    )
                    return best_primers
                else:
                    logger.warning("No valid primers found after selection")
                    return None

            except Exception as e:
                logger.error(f"Error selecting primers: {str(e)}")
                raise

    def _evaluate_site_primers(
        self,
        site_index: int,
        codon_start_dict: Dict,
        current_best: Optional[Dict]
    ) -> Optional[Dict]:
        """Evaluates primers for a specific site."""
        for codon_start, primers in codon_start_dict.items():
            for mutation in primers:
                forward_primers = mutation.get(
                    "primers", {}).get("forward", [])
                reverse_primers = mutation.get(
                    "primers", {}).get("reverse", [])

                if not forward_primers or not reverse_primers:
                    logger.debug(
                        f"Skipping mutation at site {site_index}, "
                        f"codon {codon_start}: No primers found"
                    )
                    continue

                primer_pair = self._evaluate_primer_pairs(
                    forward_primers,
                    reverse_primers
                )

                if primer_pair:
                    current_best = primer_pair

        return current_best

    def _evaluate_primer_pairs(
        self,
        forward_primers: List[Dict],
        reverse_primers: List[Dict]
    ) -> Optional[Dict]:
        """Evaluates and scores primer pairs."""
        best_pair = None

        for fwd, rev in zip(forward_primers, reverse_primers):
            self.state['primers_evaluated'] += 1

            score = self._calculate_primer_score(fwd, rev)

            if score > self.state['best_score']:
                self.state['best_score'] = score
                best_pair = {"forward": fwd, "reverse": rev}

                logger.debug(
                    f"New best primer pair found with score {score}"
                )

        return best_pair

    def _calculate_primer_score(
        self,
        forward_primer: Dict,
        reverse_primer: Dict
    ) -> float:
        """Calculates a score for a primer pair based on various criteria."""
        tm_fwd = forward_primer.get("melting_temp", 0)
        tm_rev = reverse_primer.get("melting_temp", 0)
        gc_fwd = forward_primer.get("gc_content", 0)
        gc_rev = reverse_primer.get("gc_content", 0)

        # Scoring factors (can be adjusted based on importance)
        tm_score = (tm_fwd + tm_rev) / 2
        gc_score = (gc_fwd + gc_rev) / 2

        return tm_score + gc_score

    def is_overhang_compatible(
        self,
        overhangs: List[str]
    ) -> bool:
        """Checks compatibility between overhang sequences."""
        with debug_context("check_overhang_compatibility"):
            try:
                for i, overhang_i in enumerate(overhangs):
                    if not self._validate_single_overhang(overhang_i):
                        return False

                    # Check compatibility with other overhangs
                    for overhang_j in overhangs[i + 1:]:
                        if not self._check_overhang_pair_compatibility(
                            overhang_i,
                            overhang_j
                        ):
                            return False

                return True

            except Exception as e:
                logger.error(
                    f"Error checking overhang compatibility: {str(e)}")
                raise

    def _validate_single_overhang(self, overhang: str) -> bool:
        """Validates a single overhang sequence."""
        # Check GC content
        gc_content = sum(
            1 for base in overhang if base in "GC") / len(overhang)
        if gc_content in (0, 1):
            logger.debug(f"Invalid GC content ({gc_content}) for {overhang}")
            return False
        return True

    def _check_overhang_pair_compatibility(
        self,
        overhang1: str,
        overhang2: str
    ) -> bool:
        """Checks compatibility between two overhangs."""
        # Check for three consecutive matching bases
        for k in range(len(overhang1) - 2):
            if overhang1[k:k+3] in overhang2:
                logger.debug(
                    f"Found matching triplet {overhang1[k:k+3]} "
                    f"between {overhang1} and {overhang2}"
                )
                return False

        # Check for single base difference
        diff_count = sum(1 for a, b in zip(overhang1, overhang2) if a != b)
        if diff_count <= 1:
            logger.debug(
                f"Overhangs differ by only {diff_count} base: "
                f"{overhang1} vs {overhang2}"
            )
            return False

        return True
