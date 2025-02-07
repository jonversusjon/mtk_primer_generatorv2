from typing import Dict, List, Optional, Tuple, Set
import numpy as np
from dataclasses import dataclass
from services.base import GoldenGateDesigner

@dataclass
class CutSiteRange:
    """Represents the valid cut range for a restriction site"""
    start_pos: int  # 5'-most valid cut position
    end_pos: int    # 3'-most valid cut position
    mutable_positions: List[int]  # Positions of mutable bases in restriction site
    
    @property
    def range_size(self) -> int:
        """Number of possible cut positions in range"""
        return self.end_pos - self.start_pos + 1

def calculate_cut_range(site_mutations: List[Dict]) -> CutSiteRange:
    """Calculate valid cut range based on mutable positions in a restriction site."""

    # Get all mutable positions - CHANGE THIS LINE
    mutable_positions = [mut['rs_index'] for mut in site_mutations]
    
    # Find extremes
    five_prime_pos = min(mutable_positions)
    three_prime_pos = max(mutable_positions)
    
    # Calculate range as before
    start_pos = five_prime_pos - 24
    end_pos = three_prime_pos - 15
    
    return CutSiteRange(
        start_pos=start_pos,
        end_pos=end_pos,
        mutable_positions=mutable_positions
    )
    
class MutationOptimizer(GoldenGateDesigner):
    """
    Optimizes mutation selection based on compatibility and codon usage frequency.
    Uses a precomputed binary compatibility table for fast lookups.
    """
    
    def __init__(self, verbose: bool = False):
        """
        Initialize the optimizer with a precomputed compatibility table.
        
        Args:
            compatibility_table_path: Path to binary compatibility lookup table
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
            # Read binary file
            with open(path, 'rb') as f:
                binary_data = f.read()
            
            # Verify file size is correct (256×256 = 65,536 bits = 8,192 bytes)
            expected_size = 256 * 256 // 8  # Convert bits to bytes
            if len(binary_data) != expected_size:
                raise ValueError(
                    f"Invalid compatibility table size. Expected {expected_size} bytes, "
                    f"got {len(binary_data)} bytes"
                )
                
            # Convert binary data to boolean array
            # Each byte contains 8 compatibility values
            compatibility_bits = np.unpackbits(
                np.frombuffer(binary_data, dtype=np.uint8)
            )
            
            # Reshape into 256×256 matrix
            compatibility_matrix = compatibility_bits.reshape(256, 256)
            
            if self.verbose:
                # Print some basic validation info
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


    def _seq_to_index(self, seq: str) -> int:
        """
        Converts a 4-nucleotide sequence to its corresponding matrix index.
        Sequences are indexed in alphabetical order (A=0, C=1, G=2, T=3).
        Each nucleotide position contributes 4^(3-pos) to the final index.
        Example: 'ACTG' -> A=0, C=1, T=3, G=2 -> 0×64 + 1×16 + 3×4 + 2×1 = 30
        
        Args:
            seq: String containing exactly 4 nucleotides (ACTG only)
            
        Returns:
            Integer index from 0 to 255 corresponding to sequence position in matrix
            
        Raises:
            ValueError: If sequence is not exactly 4 nucleotides or contains invalid characters
        """
        # Input validation
        if len(seq) != 4:
            raise ValueError(f"Sequence must be exactly 4 nucleotides, got {len(seq)}: {seq}")
        
        seq = seq.upper()
        if not all(nt in 'ACTG' for nt in seq):
            raise ValueError(f"Sequence must contain only A, C, G, T, got: {seq}")
        
        # Nucleotide to value mapping (alphabetical order)
        NT_VALUES = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        
        # Calculate index using positional values
        # Each position contributes 4^(3-pos) to final index
        index = 0
        for pos, nt in enumerate(seq):
            power = 3 - pos  # 3,2,1,0 for positions 0,1,2,3
            index += NT_VALUES[nt] * (4 ** power)
        
        if self.verbose:
            print(f"Converted sequence {seq} to index {index}")
        
        return index

    def get_sticky_end(self, sequence: str, cut_pos: int) -> str:
        """Get 4bp sticky end that would result from cutting at position"""
        return sequence[cut_pos:cut_pos + 4]

    def create_compatibility_tensor(self, sequence: str, mutation_options: List[List[Dict]]) -> Tuple[np.ndarray, List[CutSiteRange]]:
        print("\n=== Starting create_compatibility_tensor ===")
        print(f"Input sequence length: {len(sequence)}")
        print(f"Number of mutation options: {len(mutation_options)}")
        
        # Calculate cut ranges for each site
        cut_ranges = [calculate_cut_range(site_muts) 
                    for site_muts in mutation_options 
                    if site_muts]  # Skip empty mutation lists
        
        print(f"\nCalculated cut ranges: {len(cut_ranges)}")
        for i, r in enumerate(cut_ranges):
            print(f"Range {i}: start={r.start_pos}, end={r.end_pos}, size={r.range_size}")
            print(f"Mutable positions: {r.mutable_positions}")
        
        if not cut_ranges:
            print("No cut ranges found, returning empty arrays")
            return np.array([]), []
                
        # Create tensor with dimensions based on range sizes
        tensor_shape = [r.range_size for r in cut_ranges]
        print(f"\nTensor shape: {tensor_shape}")
        compatibility_tensor = np.ones(tensor_shape, dtype=np.int8)
        print(f"Created initial compatibility tensor with shape: {compatibility_tensor.shape}")
        
        # Iterate through all positions in the tensor
        print("\nStarting tensor population loop...")
        for indices in np.ndindex(*tensor_shape):
            
            # Get corresponding sticky ends for each position
            sticky_ends = []
            for idx, r in zip(indices, cut_ranges):
                pos = r.start_pos + idx
                sticky_end = self.get_sticky_end(sequence, pos)
                sticky_ends.append(sticky_end)
                        
            # Check compatibility between all pairs
            for i in range(len(sticky_ends)):
                for j in range(i + 1, len(sticky_ends)):
                    idx1_matrix = self._seq_to_index(sticky_ends[i])
                    idx2_matrix = self._seq_to_index(sticky_ends[j])
                    
                    # Check compatibility matrix access
                    try:
                        is_compatible = self.compatibility_matrix[idx1_matrix, idx2_matrix]
                    except Exception as e:
                        print(f"    ERROR accessing compatibility matrix: {str(e)}")
                        print(f"    Matrix type: {type(self.compatibility_matrix)}")
                        print(f"    Matrix dtype: {self.compatibility_matrix.dtype}")
                        raise
                    
                    if not is_compatible:
                        compatibility_tensor[indices] = 0
                        break
                
                if compatibility_tensor[indices] == 0:
                    break
                    
        print("\n=== Tensor creation completed ===")
        print(f"Final tensor shape: {compatibility_tensor.shape}")
        print(f"Compatible combinations: {np.sum(compatibility_tensor)}")
        print(f"Tensor density: {np.mean(compatibility_tensor):.2%}")
        
        return compatibility_tensor, cut_ranges