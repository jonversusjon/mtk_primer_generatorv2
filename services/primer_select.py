from Bio.Seq import Seq
import numpy as np
from typing import List, Dict, Optional


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
