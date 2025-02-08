from typing import Dict, List, Optional, Tuple, Set
import numpy as np
from dataclasses import dataclass
from services.base import GoldenGateDesigner
from itertools import product

@dataclass
class FragmentReassemblyRange:
    """
    Represents the range of valid fragment reassembly points for a restriction site.
    For any 6-nt restriction site, the region is defined as:
      - start: 24 bp upstream of the 5'-most base of the recognition site
      - end:   15 bp upstream of the 3'-most base of the recognition site
    """
    start: int  # 24 bp upstream of the 5'-most base
    end: int    # 15 bp upstream of the 3'-most base

    @property
    def range_size(self) -> int:
        """Number of possible reassembly points in the range"""
        return self.end - self.start + 1


def calculate_fragment_reassembly_range(recognition_site: Dict) -> FragmentReassemblyRange:
    """
    Calculate the valid fragment reassembly range based solely on recognition site coordinates.
    
    The range is defined as 24 bp upstream of the 5'-most base to 15 bp upstream of the 3'-most base.
    If the recognition site dictionary does not include an explicit end, it is assumed to be 6 nt wide.
    
    Args:
        recognition_site: Dict that should contain either:
            - 'start' and 'end' keys, or
            - 'start_index' (with a 6-nt site assumed).
            
    Returns:
        FragmentReassemblyRange with calculated start and end reassembly points.
    """
    # Determine the start of the recognition site.
    start = recognition_site.get('start', recognition_site.get('start_index'))
    # Determine the end of the recognition site. If not provided, assume a 6-nt site.
    end = recognition_site.get('end', start + 5)
    reassembly_start = start - 24  # 24 bp upstream of the first base
    reassembly_end = end - 15      # 15 bp upstream of the last base
    return FragmentReassemblyRange(start=reassembly_start, end=reassembly_end)

    
class MutationOptimizer(GoldenGateDesigner):
    """
    Optimizes mutation selection based on compatibility and codon usage frequency.
    Uses a precomputed binary compatibility table for fast lookups.
    Note: All references to "cut site" have been refactored to "fragment reassembly point".
    """
    
    def __init__(self, verbose: bool = False):
        """
        Initialize the optimizer with a precomputed compatibility table.
        
        Args:
            verbose: Whether to print detailed progress information
        """
        self.verbose = verbose
        self.compatibility_table_path = 'static/data/compatibility_table.bin'
        self.compatibility_matrix = self._load_compatibility_table(self.compatibility_table_path)

    def _load_compatibility_table(self, path: str) -> np.ndarray:
        """
        Loads the binary compatibility table into a numpy array.
        
        The binary file contains a 256×256 compatibility matrix where each element
        represents whether two 4-nucleotide sequences are compatible. The sequences
        are ordered alphabetically, so AA[AA] is at [0,0] and TT[TT] is at [255,255].
        
        Args:
            path: Path to the binary compatibility table file
                
        Returns:
            256x256 numpy array where element [i,j] indicates if sequence i is compatible with j
            
        Raises:
            FileNotFoundError: If compatibility table file cannot be found
            ValueError: If table dimensions or content are invalid
        """
        try:
            with open(path, 'rb') as f:
                binary_data = f.read()
            
            expected_size = 256 * 256 // 8  # 65,536 bits = 8,192 bytes
            if len(binary_data) != expected_size:
                raise ValueError(
                    f"Invalid compatibility table size. Expected {expected_size} bytes, "
                    f"got {len(binary_data)} bytes"
                )
                
            compatibility_bits = np.unpackbits(
                np.frombuffer(binary_data, dtype=np.uint8)
            )
            compatibility_matrix = compatibility_bits.reshape(256, 256)
            
            if self.verbose:
                print(f"Loaded compatibility matrix with shape: {compatibility_matrix.shape}")
                print(f"Number of compatible pairs: {np.sum(compatibility_matrix)}")
                print(f"Matrix density: {np.mean(compatibility_matrix):.2%}")
            
            return compatibility_matrix
            
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Compatibility table not found at {path}. "
                "Please ensure the binary file is in the correct location."
            )
        except Exception as e:
            raise ValueError(f"Error loading compatibility table: {str(e)}")
    


    def get_sticky_end(self, sequence: str, reassembly_point: int) -> str:
        """
        Get 4bp sticky end that would result from cutting at the given fragment reassembly point.
        
        Args:
            sequence: The full sequence to be cut.
            reassembly_point: The position at which the fragment reassembly point is defined.
            
        Returns:
            4-nucleotide sticky end sequence.
        """
        return sequence[reassembly_point:reassembly_point + 4]


    def optimize_mutations(
        self,
        mutation_result: Dict,
        compatibility_tensor: np.ndarray
    ) -> Tuple[List[Dict], List[int]]:
        """
        Find the optimal combination of mutations that are compatible according to the compatibility_tensor.
        
        For each restriction site, mutation options are extracted from mutation_result['restriction_sites'].
        The compatibility_tensor is assumed to have dimensions corresponding to the number of options per site.
        The optimal combination is defined as the one with the highest total usage score (sum of 'usage' across sites).
        
        Args:
            mutation_result: Dictionary containing mutation data (including a key 'restriction_sites')
            compatibility_tensor: A numpy array representing the compatibility of each combination.
        
        Returns:
            A tuple (best_mutations, valid_positions), where:
            - best_mutations is a list of chosen mutation dictionaries (one per site),
            - valid_positions is a list of indices (one per site) representing the position of the selected mutation.
        """
        # Extract mutation options per site.
        # Note: We now include the "site_index" key (using the site's start_index) for downstream compatibility.
        site_mutations = []
        for site in mutation_result.get("restriction_sites", []):
            site_muts = []
            for codon in site.get("codons", []):
                if "mutations" in codon:
                    for mutation in codon["mutations"]:
                        mut_copy = mutation.copy()
                        # Update the mutation info to include site_index (and other site-specific info)
                        mut_copy.update({
                            "site_index": site.get("start_index"),
                            "enzyme": site.get("enzyme"),
                            "strand": site.get("strand")
                        })
                        site_muts.append(mut_copy)
            if site_muts:
                site_mutations.append(site_muts)

        if not site_mutations:
            return ([], [])

        best_combo = None
        best_usage = -float("inf")
        best_indices = None

        # Iterate over every possible combination (one mutation per site)
        for indices in product(*(range(len(site)) for site in site_mutations)):
            combo = [site_mutations[i][j] for i, j in enumerate(indices)]
            # Check compatibility: compatibility_tensor is assumed to return 1 if compatible.
            if compatibility_tensor[indices] == 1:
                total_usage = sum(mut.get("usage", 0) for mut in combo)
                if total_usage > best_usage:
                    best_usage = total_usage
                    best_combo = combo
                    best_indices = indices

        if best_combo is None:
            return ([], [])
        else:
            return (best_combo, list(best_indices))


    def _gather_site_mutations(self, mutation_dict: Dict) -> List[List[Dict]]:
        """
        Gather possible mutations for each restriction site and attach the fragment reassembly range
        computed from the recognition site. (Note: mutable positions are no longer used.)
        
        Args:
            mutation_dict: Dictionary containing restriction sites and associated codon mutations.
            
        Returns:
            A list (per site) of mutation dictionaries.
        """
        site_mutations = []
        
        for site in mutation_dict['restriction_sites']:
            site_muts = []
            # Compute fragment reassembly range based on recognition site.
            # Use 'start' if available; otherwise, use 'start_index' (with a 6-nt assumption).
            frag_range = calculate_fragment_reassembly_range(site)
            
            for codon in site['codons']:
                for mutation in codon.get('mutations', []):
                    site_muts.append({
                        'site_index': site.get('start', site.get('start_index')),
                        'mutated_base': mutation['rs_index'],
                        'usage': mutation['usage'],
                        'mutation': mutation,
                        'fragment_reassembly_range': frag_range
                    })
            
            site_mutations.append(site_muts)
            print(f"Found {len(site_muts)} possible mutations for site {site.get('start', site.get('start_index'))}")
        
        return site_mutations

    def _filter_compatible_combinations(self, mutation_combinations, compatibility_matrix):
        compatible_sets = []
        
        for mut_set in mutation_combinations:
            cut_ranges = self._calculate_cut_ranges(mut_set)
            
            # Check compatibility matrix
            sliced_matrix = self._get_sliced_matrix(cut_ranges, compatibility_matrix)
            
            if np.any(sliced_matrix == 1):
                total_usage = sum(mut['usage'] for mut in mut_set)
                compatible_sets.append({
                    'mutations': mut_set,
                    'total_usage': total_usage,
                    'cut_ranges': cut_ranges
                })
        
        return compatible_sets

    def _calculate_cut_ranges(self, mut_set):
        cut_ranges = []
        
        for mut in mut_set:
            rs_index = mut['mutation']['rs_index']
            max_rs_index = max(mut['mutable_positions'])
            offset = max_rs_index - rs_index
            cut_ranges.append((offset, offset + 9))
        
        return cut_ranges

    def _get_sliced_matrix(self, cut_ranges, compatibility_matrix):
        sliced_matrix = compatibility_matrix[
            cut_ranges[0][0]:cut_ranges[0][1],
            cut_ranges[1][0]:cut_ranges[1][1],
            cut_ranges[2][0]:cut_ranges[2][1]
        ]
        return sliced_matrix

    def _find_valid_positions(self, best_set, compatibility_matrix):
        valid_positions = []
        sliced_matrix = self._get_sliced_matrix(best_set['cut_ranges'], compatibility_matrix)
        
        for i in range(sliced_matrix.shape[0]):
            for j in range(sliced_matrix.shape[1]):
                for k in range(sliced_matrix.shape[2]):
                    if sliced_matrix[i, j, k] == 1:
                        relative_positions = [pos - 24 for pos in [i, j, k]]
                        valid_positions.append({
                            'matrix_indices': (i, j, k),
                            'relative_positions': relative_positions,
                            'mutations': best_set['mutations']
                        })
                        print(f"Valid position found: matrix{(i, j, k)} -> relative{relative_positions}")
        
        return valid_positions

    def create_compatibility_tensor(
        self,
        sequence: str,
        mutation_options: List[List[Dict]]
    ) -> Tuple[np.ndarray, List[Tuple[int, int]]]:
        """
        Given a sequence and mutation options (a list of mutation lists per site),
        create an M-dimensional compatibility tensor and corresponding cut ranges.

        For each restriction site, the fragment reassembly range is computed (using the first
        mutation option's stored 'fragment_reassembly_range'). That range defines the valid reassembly
        positions for that site. For each candidate position, the sticky end is extracted using 
        self.get_sticky_end, converted to an index via self._seq_to_index, and then used to check
        pairwise compatibility with candidate sticky ends from other sites (via self.compatibility_matrix).

        The resulting tensor has shape [m1, m2, m3, ...] where each mᵢ is the number of valid 
        reassembly positions (i.e. range size) for that site. A given cell is set to 1 if, for that 
        combination of reassembly positions, all pairs of sticky ends are compatible; otherwise 0.

        Args:
            sequence: The full DNA sequence.
            mutation_options: A list (per site) of lists of mutation dictionaries.
                            (Each mutation dictionary is assumed to have a 'fragment_reassembly_range'
                            key containing a FragmentReassemblyRange instance.)

        Returns:
            A tuple (compatibility_tensor, cut_ranges) where:
            - compatibility_tensor is an M-dimensional numpy array of 0s and 1s.
            - cut_ranges is a list of (start, end_exclusive) tuples (one per site).
        """
        cut_ranges = []
        # For each site, derive the reassembly range from the first mutation option.
        for site_muts in mutation_options:
            if site_muts:
                # Use the stored fragment reassembly range.
                frag_range = site_muts[0]['fragment_reassembly_range']
                # For slicing, we use (start, end+1) since the end is inclusive in our definition.
                cut_ranges.append((frag_range.start, frag_range.end + 1))
            else:
                cut_ranges.append((0, 0))
        
        # For each site, compute the candidate reassembly positions and their sticky end indices.
        dims = []           # Will hold the number of candidate positions per site.
        sticky_indices = [] # List (per site) of sticky end indices for each candidate position.
        for (start, end) in cut_ranges:
            size = end - start  # since 'end' is exclusive
            dims.append(size)
            site_sticky_indices = []
            for pos in range(start, end):
                # Extract the sticky end from the sequence at the candidate reassembly position.
                sticky_end = self.get_sticky_end(sequence, pos)
                # Convert the 4-nt sticky end to its corresponding index (0-255).
                idx = self._seq_to_index(sticky_end)
                site_sticky_indices.append(idx)
            sticky_indices.append(site_sticky_indices)
        
        tensor_shape = tuple(dims)
        # Create an empty tensor of the required shape.
        comp_tensor = np.zeros(tensor_shape, dtype=np.uint8)
        
        # Iterate over every combination of candidate positions using np.ndindex.
        # Each combination is a tuple (i1, i2, ..., i_M) where i_k is the candidate index for site k.
        for combination in np.ndindex(tensor_shape):
            # For each site, get the sticky index for the candidate position.
            candidate_sticky_idxs = [
                sticky_indices[site_idx][combination[site_idx]]
                for site_idx in range(len(combination))
            ]
            # Check pairwise compatibility among all candidate sticky ends.
            all_compatible = True
            num_sites = len(candidate_sticky_idxs)
            for i in range(num_sites):
                for j in range(i + 1, num_sites):
                    idx_i = candidate_sticky_idxs[i]
                    idx_j = candidate_sticky_idxs[j]
                    # If the preloaded compatibility matrix indicates incompatibility, mark as not compatible.
                    if self.compatibility_matrix[idx_i, idx_j] != 1:
                        all_compatible = False
                        break
                if not all_compatible:
                    break
            comp_tensor[combination] = 1 if all_compatible else 0
        
        if self.verbose:
            print(f"Created compatibility tensor with shape: {tensor_shape}")
            print(f"Cut ranges (per site): {cut_ranges}")
        
        return comp_tensor, cut_ranges
