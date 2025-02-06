# services/primer_select.py
from Bio.Seq import Seq
import numpy as np
from typing import List, Dict, Optional, Union
from .base import GoldenGateDesigner


class PrimerSelector(GoldenGateDesigner):
    def __init__(self, verbose: bool = False):
        super().__init__(verbose=verbose)
        self.state = {
            'current_operation': '',
            'best_score': float('-inf'),
            'primers_evaluated': 0
        }

    def format_primers_for_output(
        self,
        primer_sets: List[List[Seq]]
    ) -> List[List[str]]:
        """Formats primer sets for output in the final protocol."""
        with self.debug_context("format_primers"):
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
                        self.logger.debug(f"Formatted primer: {primer_entry[0]}")
                
                return primer_data
            
            except Exception as e:
                self.logger.error(f"Error formatting primers: {str(e)}")
                raise

    def select_best_internal_primers(
        self,
        primer_sets: Dict,
        template_seq: Optional[str] = None
    ) -> Optional[Dict]:
        """Selects optimal internal primers based on multiple criteria."""
        with self.debug_context("select_best_primers"):
            if not primer_sets:
                self.logger.warning("No primer sets provided!")
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
                    self.logger.info(
                        f"Selected best primers with score {self.state['best_score']}"
                    )
                    return best_primers
                else:
                    self.logger.warning("No valid primers found after selection")
                    return None

            except Exception as e:
                self.logger.error(f"Error selecting primers: {str(e)}")
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
                forward_primers = mutation.get("primers", {}).get("forward", [])
                reverse_primers = mutation.get("primers", {}).get("reverse", [])

                if not forward_primers or not reverse_primers:
                    self.logger.debug(
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
                
                self.logger.debug(
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
        with self.debug_context("check_overhang_compatibility"):
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
                self.logger.error(f"Error checking overhang compatibility: {str(e)}")
                raise

    def _validate_single_overhang(self, overhang: str) -> bool:
        """Validates a single overhang sequence."""
        # Check GC content
        gc_content = sum(1 for base in overhang if base in "GC") / len(overhang)
        if gc_content in (0, 1):
            self.logger.debug(f"Invalid GC content ({gc_content}) for {overhang}")
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
                self.logger.debug(
                    f"Found matching triplet {overhang1[k:k+3]} "
                    f"between {overhang1} and {overhang2}"
                )
                return False

        # Check for single base difference
        diff_count = sum(1 for a, b in zip(overhang1, overhang2) if a != b)
        if diff_count <= 1:
            self.logger.debug(
                f"Overhangs differ by only {diff_count} base: "
                f"{overhang1} vs {overhang2}"
            )
            return False

        return True

# from Bio.Seq import Seq
# import numpy as np
# from typing import List, Dict, Optional


# def format_primers_for_output(primer_sets: List[List[Seq]]) -> List[List[str]]:
#     """
#     Formats a list of primer sets for output in the final protocol.
#     """
#     primer_data = []
#     for i, primers in enumerate(primer_sets):
#         for j, primer in enumerate(primers):
#             primer_data.append([f"Primer_{i+1}_{j+1}", str(primer), f"Amplicon_{i+1}"])
#     return primer_data

# def select_best_internal_primers(primer_sets, template_seq):
#     """
#     Selects the best internal primers based on melting temperature, GC content, and off-target binding.
    
#     Arguments:
#     - primer_sets: A dictionary of designed primers (output from get_all_possible_internal_primers).
#     - template_seq: The reference sequence to check primer specificity.
    
#     Returns:
#     - The best primer set or None if no valid primers are found.
#     """
#     if not primer_sets:
#         print("❌ No primer sets provided!")
#         return None

#     best_primers = None
#     best_score = float('-inf')

#     for site_index, codon_start_dict in primer_sets.items():
#         for codon_start, primers in codon_start_dict.items():
#             for mutation in primers:
#                 # Extract forward and reverse primers
#                 forward_primers = mutation.get("primers", {}).get("forward", [])
#                 reverse_primers = mutation.get("primers", {}).get("reverse", [])

#                 if not forward_primers or not reverse_primers:
#                     print(f"❌ Skipping mutation at site {site_index}, codon {codon_start}: No primers found.")
#                     continue

#                 for fwd, rev in zip(forward_primers, reverse_primers):
#                     tm_fwd = fwd.get("melting_temp", 0)
#                     tm_rev = rev.get("melting_temp", 0)
#                     gc_fwd = fwd.get("gc_content", 0)
#                     gc_rev = rev.get("gc_content", 0)

#                     # Scoring: Prefer primers with balanced TM and GC content (adjust scoring criteria as needed)
#                     score = (tm_fwd + tm_rev) / 2 + (gc_fwd + gc_rev) / 2

#                     if score > best_score:
#                         best_score = score
#                         best_primers = {"forward": fwd, "reverse": rev}

#     if best_primers:
#         print(f"✅ Best primers selected: {best_primers}")
#         return best_primers
#     else:
#         print("❌ No valid primers found after selection.")
#         return None


# def is_overhang_compatible(overhangs):
#     """
#     Checks if a set of overhang sequences are compatible by:
#     1. Ensuring no three consecutive bases match across different overhangs.
#     2. Ensuring no two overhangs differ by only one base.
#     3. Ensuring each overhang has balanced GC content (not 0% or 100%).
#     """
#     num_overhangs = len(overhangs)

#     for i in range(num_overhangs):
#         overhang_i = overhangs[i]

#         for j in range(i + 1, num_overhangs):  # Only check pairs once
#             overhang_j = overhangs[j]

#             # Check for three consecutive bases being the same
#             for k in range(len(overhang_i) - 2):
#                 if overhang_i[k:k+3] in overhang_j:
#                     return False

#             # Check if they differ by only one base pair
#             diff_count = sum(1 for a, b in zip(overhang_i, overhang_j) if a != b)
#             if diff_count <= 1:
#                 return False

#         # Check GC content (must not be 0% or 100%)
#         gc_content = sum(1 for base in overhang_i if base in "GC") / len(overhang_i)
#         if gc_content == 0 or gc_content == 1:
#             return False

#     return True
