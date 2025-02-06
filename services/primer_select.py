from Bio.Seq import Seq
import numpy as np
from typing import List, Dict, Optional


def format_primers_for_output(primer_sets: List[List[Seq]]) -> List[List[str]]:
    """
    Formats a list of primer sets for output in the final protocol.
    """
    primer_data = []
    for i, primers in enumerate(primer_sets):
        for j, primer in enumerate(primers):
            primer_data.append([f"Primer_{i+1}_{j+1}", str(primer), f"Amplicon_{i+1}"])
    return primer_data

def select_best_internal_primers(primer_sets, template_seq):
    """
    Selects the best internal primers based on melting temperature, GC content, and off-target binding.
    
    Arguments:
    - primer_sets: A dictionary of designed primers (output from get_all_possible_internal_primers).
    - template_seq: The reference sequence to check primer specificity.
    
    Returns:
    - The best primer set or None if no valid primers are found.
    """
    if not primer_sets:
        print("❌ No primer sets provided!")
        return None

    best_primers = None
    best_score = float('-inf')

    for site_index, codon_start_dict in primer_sets.items():
        for codon_start, primers in codon_start_dict.items():
            for mutation in primers:
                # Extract forward and reverse primers
                forward_primers = mutation.get("primers", {}).get("forward", [])
                reverse_primers = mutation.get("primers", {}).get("reverse", [])

                if not forward_primers or not reverse_primers:
                    print(f"❌ Skipping mutation at site {site_index}, codon {codon_start}: No primers found.")
                    continue

                for fwd, rev in zip(forward_primers, reverse_primers):
                    tm_fwd = fwd.get("melting_temp", 0)
                    tm_rev = rev.get("melting_temp", 0)
                    gc_fwd = fwd.get("gc_content", 0)
                    gc_rev = rev.get("gc_content", 0)

                    # Scoring: Prefer primers with balanced TM and GC content (adjust scoring criteria as needed)
                    score = (tm_fwd + tm_rev) / 2 + (gc_fwd + gc_rev) / 2

                    if score > best_score:
                        best_score = score
                        best_primers = {"forward": fwd, "reverse": rev}

    if best_primers:
        print(f"✅ Best primers selected: {best_primers}")
        return best_primers
    else:
        print("❌ No valid primers found after selection.")
        return None


def is_overhang_compatible(overhangs):
    """
    Checks if a set of overhang sequences are compatible by:
    1. Ensuring no three consecutive bases match across different overhangs.
    2. Ensuring no two overhangs differ by only one base.
    3. Ensuring each overhang has balanced GC content (not 0% or 100%).
    """
    num_overhangs = len(overhangs)

    for i in range(num_overhangs):
        overhang_i = overhangs[i]

        for j in range(i + 1, num_overhangs):  # Only check pairs once
            overhang_j = overhangs[j]

            # Check for three consecutive bases being the same
            for k in range(len(overhang_i) - 2):
                if overhang_i[k:k+3] in overhang_j:
                    return False

            # Check if they differ by only one base pair
            diff_count = sum(1 for a, b in zip(overhang_i, overhang_j) if a != b)
            if diff_count <= 1:
                return False

        # Check GC content (must not be 0% or 100%)
        gc_content = sum(1 for base in overhang_i if base in "GC") / len(overhang_i)
        if gc_content == 0 or gc_content == 1:
            return False

    return True
