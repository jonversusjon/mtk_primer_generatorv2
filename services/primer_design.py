from typing import Dict, List, Optional, Any, Tuple, Union
import primer3
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from .base import PrimerDesignLogger
from .utils import GoldenGateUtils
import numpy as np
from contextlib import contextmanager

class PrimerDesigner(PrimerDesignLogger):
    def __init__(
        self, 
        part_num_left: List[str],
        part_num_right: List[str],
        kozak: str = "MTK",
        verbose: bool = False
    ):
        super().__init__(verbose=verbose)
        self.state = {
            'current_operation': '',
            'primers_designed': 0,
            'current_mutation': None
        }
        
        self.part_num_left = part_num_left
        self.part_num_right = part_num_right
        self.kozak = kozak
        
        # Default parameters for primer design
        self.default_params = {
            'tm_threshold': 45.0,
            'min_3p_match': 10,
            'max_mismatches': 1,
            'mv_conc': 50.0,
            'dv_conc': 1.5,
            'dntp_conc': 0.2,
            'dna_conc': 250.0,
            'min_tm': 57
        }
    
        self.bsmbi_site = "CGTCTC"
        self.spacer = "GAA"
        self.utils = GoldenGateUtils()

    @contextmanager
    def debug_context(self, operation: str):
        """Context manager for debugging operations."""
        prev_operation = self.state['current_operation']
        self.state['current_operation'] = operation
        try:
            yield
        finally:
            self.state['current_operation'] = prev_operation

    def design_primers(
        self,
        sequence: str,
        seq_index: int,
        mutations: Optional[List[Dict]] = None,
        compatibility_matrices: List[np.ndarray] = None,
        template_seq: Optional[str] = None,
        primer_names: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """
        Main orchestrator function for primer design.
        """
        with self.debug_context("design_primers"):
            results = {}
            
            # 1. Design mutation primers if mutations are provided
            if mutations and compatibility_matrices:
                mutation_primers = self._design_mutation_primers(
                    sequence, 
                    mutations, 
                    compatibility_matrices,
                    primer_names[seq_index] if primer_names else None
                )
                if mutation_primers:
                    results['mutation_primers'] = mutation_primers

            # 2. Generate edge primers
            mtk_partend_sequences = self.utils.get_mtk_partend_sequences()
            edge_primers = self._generate_edge_primers(
                sequence,
                seq_index,
                mtk_partend_sequences,
                primer_names[seq_index] if primer_names else None
            )
            if edge_primers:
                results['edge_primers'] = edge_primers

            # 3. Perform off-target analysis if template sequence is provided
            if template_seq and results:
                self._analyze_off_targets_for_all_primers(results, template_seq)

            return results

    def _design_mutation_primers(
        self, 
        target_seq: str, 
        mutation_sets: List[Dict], 
        comp_matrices: List[np.ndarray],
        primer_name: Optional[str] = None
    ) -> Optional[Dict[str, Dict[str, str]]]:
        """
        Design primers for mutation sets based on compatibility matrices.
        """
        with self.debug_context("design_mutation_primers"):
            for set_index, mutation_set in enumerate(mutation_sets):
                self.state['current_mutation'] = set_index
                
                # Find valid overhang combinations
                valid_coords = np.argwhere(comp_matrices[set_index] == 1)
                if valid_coords.size == 0:
                    self.logger.debug(f"No valid overhangs for mutation set {set_index}")
                    continue

                # Apply mutations to sequence
                mutated_seq = self._apply_mutations(target_seq, mutation_set)
                
                # Design primers using first valid combination
                primers = self._construct_mutation_primers(
                    mutated_seq,
                    mutation_set,
                    valid_coords[0],
                    primer_name
                )
                
                if primers and self._validate_all_primers(primers):
                    self.state['primers_designed'] += len(primers)
                    return primers
            
            return None

    def _apply_mutations(self, target_seq: str, mutation_set: Dict) -> str:
        """Apply a set of mutations to the target sequence."""
        with self.debug_context("apply_mutations"):
            mutated_seq = list(target_seq)
            
            for mut in mutation_set.values():
                pos = mut["position"] - 1
                alt_seq = mut["alternative_sequence"]
                orig_seq = mut["original_sequence"]
                
                mutated_seq[pos:pos + len(orig_seq)] = alt_seq
                
                if self.verbose:
                    self._log_mutation(target_seq, ''.join(mutated_seq), pos, len(orig_seq))
            
            return ''.join(mutated_seq)

    def _construct_mutation_primers(
        self,
        mutated_seq: str,
        mutation_set: Dict,
        selected_coord: np.ndarray,
        primer_name: Optional[str] = None,
        annealing_length: int = 18
    ) -> Optional[Dict[str, Dict[str, str]]]:
        """Construct primers for a mutation set using selected coordinates."""
        with self.debug_context("construct_mutation_primers"):
            primers = {}
            
            for mut in mutation_set.values():
                position = mut["position"] - 1
                
                # Design forward and reverse primers
                forward = self._construct_primer(
                    mutated_seq,
                    position,
                    mut,
                    selected_coord[0],
                    annealing_length,
                    is_reverse=False,
                    primer_name=primer_name
                )
                
                reverse = self._construct_primer(
                    mutated_seq,
                    position,
                    mut,
                    selected_coord[0],
                    annealing_length,
                    is_reverse=True,
                    primer_name=primer_name
                )
                
                if not forward or not reverse:
                    return None
                    
                primers[mut["site"]] = {
                    "forward": forward,
                    "reverse": reverse,
                    "mutation_info": mut
                }

            return primers

    def _construct_primer(
        self,
        sequence: str,
        position: int,
        mut: Dict,
        selected_coord: int,
        annealing_length: int,
        is_reverse: bool,
        primer_name: Optional[str] = None
    ) -> Optional[str]:
        """Construct a single primer with proper components."""
        # Get binding sequence
        if is_reverse:
            start_pos = position - annealing_length + 1
            if start_pos < 0:
                return None
            binding = sequence[start_pos:position + 1]
            binding = str(Seq(binding).reverse_complement())
            # Use bottom extended sequence for reverse primer
            extended_seq = str(mut['overhangs']['bottom_extended_sequences'][selected_coord])
        else:
            if position + annealing_length > len(sequence):
                return None
            binding = sequence[position:position + annealing_length]
            # Use top extended sequence for forward primer
            extended_seq = str(mut['overhangs']['top_extended_sequences'][selected_coord])

        # Construct complete primer
        primer_seq = self.spacer + self.bsmbi_site + extended_seq + binding
        print(f"primer_seq: {primer_seq}")
        # Add primer name if provided
        if primer_name:
            suffix = "_RV" if is_reverse else "_FW"
            return (f"{primer_name}{suffix}", primer_seq)
        return primer_seq

    def _validate_all_primers(self, primers: Dict[str, Dict[str, Union[str, Dict]]]) -> bool:
        """Validate all primers in a set."""
        with self.debug_context("validate_primers"):
            for site_primers in primers.values():
                if not self._validate_primer_pair(site_primers["forward"], site_primers["reverse"]):
                    return False
            return True

    def _validate_primer_pair(self, forward: str, reverse: str) -> bool:
        """Validate a pair of primers."""
        # Extract actual sequences if primers are named tuples
        forward_seq = forward[1] if isinstance(forward, tuple) else forward
        reverse_seq = reverse[1] if isinstance(reverse, tuple) else reverse
        
        # Check BsmBI sites
        forward_count = forward_seq.count(self.bsmbi_site)
        reverse_count = reverse_seq.count(self.bsmbi_site)
        if forward_count != 1 or reverse_count != 1:
            return False
            
        # Additional validations could be added here
        return True

    def _generate_edge_primers(
        self,
        sequence: str,
        seq_index: int,
        part_end_dict: Dict[str, str],
        primer_name: Optional[str] = None
    ) -> Dict[str, Tuple[str, str]]:
        """Generate edge primers for assembly."""
        with self.debug_context("generate_edge_primers"):
            primers = {}
            
            # Forward primer
            forward_name, forward_seq = self._construct_edge_primer(
                sequence,
                part_end_dict,
                self.part_num_left[seq_index],
                "forward",
                primer_name
            )
            
            # Reverse primer
            reverse_name, reverse_seq = self._construct_edge_primer(
                sequence,
                part_end_dict,
                self.part_num_right[seq_index],
                "reverse",
                primer_name
            )
            
            if forward_name and reverse_name:
                primers['forward'] = (forward_name, forward_seq)
                primers['reverse'] = (reverse_name, reverse_seq)
                
            return primers

    def _analyze_off_targets_for_all_primers(
        self,
        primer_results: Dict[str, Any],
        template_seq: str
    ) -> None:
        """Perform off-target analysis for all primers."""
        with self.debug_context("analyze_off_targets"):
            try:
                for primer_type, primers in primer_results.items():
                    if primer_type == 'mutation_primers':
                        for site_primers in primers.values():
                            self._analyze_primer_pair(site_primers, template_seq)
                    elif primer_type == 'edge_primers':
                        self._analyze_primer_pair(primers, template_seq)
            except Exception as e:
                self.logger.error(f"Error in off-target analysis: {str(e)}")

    def _analyze_primer_pair(
        self,
        primer_pair: Dict[str, Any],
        template_seq: str
    ) -> None:
        """Analyze a pair of primers for off-target binding."""
        for orientation in ['forward', 'reverse']:
            try:
                if isinstance(primer_pair[orientation], tuple):
                    binding_seq = primer_pair[orientation][1]
                else:
                    binding_seq = primer_pair[orientation]
                    
                off_targets = self._find_off_targets(binding_seq, template_seq)
                primer_pair[f'{orientation}_off_targets'] = off_targets
            except Exception as e:
                self.logger.error(f"Error analyzing {orientation} primer: {str(e)}")
                primer_pair[f'{orientation}_off_targets'] = []
    
# # services/primer_design.py
# from typing import Dict, List, Optional, Any, Tuple
# import primer3
# from Bio.Seq import Seq
# from Bio.SeqUtils import MeltingTemp as mt
# from .base import PrimerDesignLogger
# from .utils import GoldenGateUtils
# import numpy as np

# class PrimerDesigner(PrimerDesignLogger):
#     def __init__(
#         self, 
#         part_num_left: List[str],
#         part_num_right: List[str],
#         kozak: str = "MTK",
#         verbose: bool = False):
#         super().__init__(verbose=verbose)
#         self.state = {
#             'current_operation': '',
#             'primers_designed': 0,
#             'current_mutation': None
#         }
        
#         self.part_num_left = part_num_left
#         self.part_num_right = part_num_right
#         self.kozak = kozak
        
#         # Default parameters
#         self.default_params = {
#             'tm_threshold': 45.0,
#             'min_3p_match': 10,
#             'max_mismatches': 1,
#             'mv_conc': 50.0,
#             'dv_conc': 1.5,
#             'dntp_conc': 0.2,
#             'dna_conc': 250.0,
#             'min_tm': 57
#         }
    
#         self.bsmbi_site = "CGTCTC"
#         self.spacer = "GAA"
#         self.utils = GoldenGateUtils()
        
#     def design_primers(
#         self,
#         sequence: str,
#         seq_index: int,
#         mutations: Optional[List[Dict]] = None,
#         compatibility_matrices: List[np.ndarray] = None,
#         template_seq: Optional[str] = None,
#     ):
#         if mutations:
#             # print(f"mutations: {mutations}")
#             internal_primers = self.design_mutation_primers(sequence, mutations, compatibility_matrices)
        
#         return
 

#     def design_mutation_primers(self, target_seq, mutation_sets, comp_matrices, annealing_length=18):
#         """Main entry point for designing mutation primers."""
#         self.logger.debug(f"Starting primer design with {len(mutation_sets)} mutation sets")
#         for set_index, mutation_set in enumerate(mutation_sets):
#             # Find valid overhang combinations
#             valid_coords = self._find_valid_overhangs(comp_matrices[set_index])
#             if valid_coords is None:
#                 self.logger.debug(f"No valid overhangs found for mutation set {set_index}")
#                 continue

#             # Design primers for this mutation set
#             primers = self._design_primers_for_set(
#                 target_seq, 
#                 mutation_set, 
#                 valid_coords[0],  # Use first valid combination
#                 annealing_length
#             )
            
#             if primers:
#                 return primers
        
#         return None

#     def _find_valid_overhangs(self, comp_matrix):
#         """Find valid overhang combinations from compatibility matrix."""
#         valid_coords = np.argwhere(comp_matrix == 1)
#         return valid_coords if valid_coords.size > 0 else None

#     def _apply_mutations(self, target_seq, mutation_set):
#         """Apply mutations to target sequence and log changes."""
#         mutated_seq = list(target_seq)
        
#         for mut in mutation_set.values():
#             pos = mut["position"] - 1
#             alt_seq = mut["alternative_sequence"]
#             orig_seq = mut["original_sequence"]
            
#             # Apply mutation
#             mutated_seq[pos:pos + len(orig_seq)] = alt_seq
            
#             if self.verbose:
#                 self._log_mutation(target_seq, mutated_seq, pos, len(orig_seq))
        
#         return "".join(mutated_seq)

    def _log_mutation(self, original_seq, mutated_seq, position, length):
        """Log details of each mutation."""
        for i in range(position, position + length):
            old_base = original_seq[i] if i < len(original_seq) else "N/A"
            new_base = mutated_seq[i] if i < len(mutated_seq) else "N/A"
            if old_base != new_base:
                self.logger.info(f"Mutated base at position {i+1}: {old_base} -> {new_base}")

#     def _get_primer_components(self, mut, selected_coord_idx):
#         """Extract overhang sequences for primer design."""
#         print(f"mut: {mut}")
#         overhangs = mut["overhangs"]
#         return {
#             "forward_overhang": str(overhangs["top_strand_overhangs"][selected_coord_idx]),
#             "reverse_overhang": str(overhangs["bottom_strand_overhangs"][selected_coord_idx]),
#             "extra_forward": str(overhangs["top_extra_forward"][selected_coord_idx]),
#             "extra_reverse": str(overhangs["bottom_extra_reverse"][selected_coord_idx])
#         }

#     def _design_forward_primer(self, mutated_seq, position, components, annealing_length):
#         """Design forward primer with given components."""
#         if position + annealing_length > len(mutated_seq):
#             return None

#         binding = mutated_seq[position:position + annealing_length]
#         complete_overhang = components["extra_forward"] + components["forward_overhang"]
        
#         primer = self.spacer + self.bsmbi_site + complete_overhang + binding
#         self.logger.debug(f"Forward primer: {primer}")
#         return primer

#     def _design_reverse_primer(self, mutated_seq, position, components, annealing_length):
#         """Design reverse primer with given components."""
#         start_index = position - annealing_length + 1
#         if start_index < 0:
#             return None

#         binding = mutated_seq[start_index:position + 1]
#         rev_binding = self.utils.reverse_complement(binding)
#         complete_overhang = components["extra_reverse"] + components["reverse_overhang"]
        
#         primer = self.spacer + self.bsmbi_site + complete_overhang + rev_binding
#         self.logger.debug(f"Reverse primer: {primer}")
#         return primer


#     def _design_primers_for_set(self, target_seq, mutation_set, selected_coord, annealing_length):
#         """Design primers for a complete mutation set."""
#         # First apply all mutations to get final sequence
#         mutated_seq = self._apply_mutations(target_seq, mutation_set)
#         primers = {}
        
#         print(f"******** MUTATION SET *********")
#         print(f"mutation_set: {mutation_set}")
        
#         # Iterate through each mutation site
#         for j, mut in enumerate(mutation_set.values()):
#             position = mut["position"] - 1
            
#             # Get the six_nuc_seq from top_extended_sequences
#             six_nuc_seq = str(mut['overhangs']['top_extended_sequences'][selected_coord[j]])
#             rev_six_nuc_seq = str(Seq(six_nuc_seq).reverse_complement())
            
#             # Construct primers
#             forward = self.spacer + self.bsmbi_site + six_nuc_seq + self._get_binding_seq_for(position, annealing_length, mutated_seq)
#             reverse = self.spacer + self.bsmbi_site + rev_six_nuc_seq + self._get_binding_seq_rev(position, annealing_length, mutated_seq)
            
#             if not forward or not reverse:
#                 return None
                
#             primers[mut["site"]] = {
#                 "forward": forward,
#                 "reverse": reverse
#             }

#         return primers

#     def _get_binding_seq_for(self, position, annealing_length, sequence):
#         """
#         Get the forward binding sequence starting at a given position.
        
#         Args:
#             position (int): Starting position in the sequence (0-based)
#             annealing_length (int): Length of the binding sequence to extract
#             sequence (str): Full DNA sequence to extract from
            
#         Returns:
#             str: The binding sequence, or None if invalid
#         """
#         # Validate inputs
#         if not isinstance(position, int) or position < 0:
#             return None
            
#         if not isinstance(annealing_length, int) or annealing_length < 1:
#             return None
            
#         if not sequence or not isinstance(sequence, str):
#             return None
            
#         # Check if we have enough sequence length
#         if position + annealing_length > len(sequence):
#             return None
            
#         # Extract the binding sequence starting at position
#         binding_seq = sequence[position:position + annealing_length]
        
#         # Validate the binding sequence only contains valid DNA bases
#         valid_bases = set('ATCG')
#         if not all(base.upper() in valid_bases for base in binding_seq):
#             return None
            
#         return binding_seq

#     def _get_binding_seq_rev(self, position, annealing_length, sequence):
#         """
#         Get the reverse binding sequence ending at a given position.
        
#         Args:
#             position (int): End position in the sequence (0-based)
#             annealing_length (int): Length of the binding sequence to extract
#             sequence (str): Full DNA sequence to extract from
            
#         Returns:
#             str: The reverse binding sequence, or None if invalid
#         """
#         # Validate inputs
#         if not isinstance(position, int) or position < 0:
#             return None
            
#         if not isinstance(annealing_length, int) or annealing_length < 1:
#             return None
            
#         if not sequence or not isinstance(sequence, str):
#             return None
            
#         # Check if we have enough sequence length before position
#         start_pos = position - annealing_length + 1
#         if start_pos < 0:
#             return None
            
#         # Extract the binding sequence ending at position
#         binding_seq = sequence[start_pos:position + 1]
        
#         # Validate the binding sequence only contains valid DNA bases
#         valid_bases = set('ATCG')
#         if not all(base.upper() in valid_bases for base in binding_seq):
#             return None
            
#         # Get the reverse complement
#         complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
#         rev_complement = ''.join(complement[base.upper()] for base in reversed(binding_seq))
            
#         return rev_complement

#     def find_off_targets_for_primer(
#         self,
#         binding_seq: str,
#         template_seq: Optional[str] = None,
#         **kwargs
#     ) -> List[Dict[str, float]]:
#         """Find potential off-target binding sites using primer3."""
#         with self.debug_context("find_off_targets"):
#             if not isinstance(binding_seq, str):
#                 raise TypeError(f"Expected string for binding_seq, got {type(binding_seq)}")

#             if not template_seq:
#                 self.logger.warning(f"Skipping off-target analysis for primer: {binding_seq}")
#                 return []

#             # Merge default parameters with provided kwargs
#             params = {**self.default_params, **kwargs}
            
#             return self._analyze_off_targets(
#                 binding_seq.upper(),
#                 template_seq.upper(),
#                 params
#             )

#     def _analyze_off_targets(
#         self,
#         binding_seq: str,
#         template_seq: str,
#         params: Dict
#     ) -> List[Dict[str, float]]:
#         """Internal method to analyze off-target binding."""
#         off_targets = []
#         primer_len = len(binding_seq)
        
#         for i in range(len(template_seq) - primer_len + 1):
#             window = template_seq[i:i + primer_len]
            
#             # Check 3' end match
#             match_3p = sum(1 for a, b in zip(
#                 binding_seq[-params['min_3p_match']:],
#                 window[-params['min_3p_match']:]
#             ) if a == b)
            
#             if match_3p < params['min_3p_match']:
#                 continue

#             # Calculate heterodimer properties
#             try:
#                 result = primer3.bindings.calc_heterodimer(
#                     seq1=binding_seq,
#                     seq2=window,
#                     mv_conc=params['mv_conc'],
#                     dv_conc=params['dv_conc'],
#                     dntp_conc=params['dntp_conc'],
#                     dna_conc=params['dna_conc']
#                 )

#                 if result.tm >= params['tm_threshold']:
#                     mismatches = sum(1 for a, b in zip(binding_seq, window) if a != b)
#                     if mismatches <= params['max_mismatches']:
#                         off_targets.append({
#                             "index": i,
#                             "tm": round(result.tm, 2),
#                             "mismatches": mismatches
#                         })
#             except Exception as e:
#                 self.logger.error(f"Error in heterodimer calculation: {str(e)}")

#         if off_targets:
#             self._log_off_target_results(binding_seq, off_targets)
            
#         return off_targets



#     def _design_primer_pair(
#         self,
#         seq: Seq,
#         nucleotide_index: int,
#         new_nucleotide: str,
#         spacer: str,
#         shift: int,
#         min_tm: float
#     ) -> Optional[Tuple[Dict, Dict]]:
#         """Design a single pair of primers."""
#         left = nucleotide_index - (6 - 1) + shift
#         right = nucleotide_index + 1 + shift

#         if not self._validate_primer_positions(left, right, len(seq)):
#             return None

#         six_nuc_seq = self._get_modified_sequence(
#             seq, nucleotide_index, new_nucleotide, left, right
#         )

#         binding_seqs = self._find_binding_sequences(seq, right, left, min_tm)
#         if not all(binding_seqs):
#             return None

#         primers = self._construct_primer_pair(
#             seq, spacer, self.bsmbi_site, six_nuc_seq, *binding_seqs
#         )
        
#         if self._validate_primer_pair(primers):
#             return self._format_primer_pair(primers)
            
#         return None

#     def generate_GG_edge_primers(
#         self,
#         seq: Seq,
#         part_end_dict: Dict[str, str],
#         part_num: str,
#         primer_orientation: str,
#         kozak: str,
#         primer_name: Optional[str] = None
#     ) -> Tuple[Optional[str], Optional[str]]:
#         """Generate Golden Gate edge primers."""
#         with self.debug_context("generate_edge_primers"):
#             try:
#                 seq = Seq(seq)
#                 part_info = self._get_part_info(part_num, primer_orientation, kozak)
                
#                 part_specific_seq = part_end_dict.get(part_info['key'])
#                 if not part_specific_seq:
#                     raise ValueError(f"No sequence found for {part_info['key']}")

#                 primer = self._construct_edge_primer(
#                     seq, part_specific_seq, primer_orientation, 
#                     part_info['suffix'], primer_name
#                 )
                
#                 return primer
                
#             except Exception as e:
#                 self.logger.error(f"Error generating edge primer: {str(e)}")
#                 return None, None

#     def _log_off_target_results(self, binding_seq: str, off_targets: List[Dict]):
#         """Log off-target analysis results."""
#         self.logger.info("\n🔍 Off-Target Primer Found")
#         self.logger.info(f"Primer Binding Sequence: {binding_seq}")
#         self.logger.info(f"Total Off-Targets: {len(off_targets)}")
        
#         for off_target in off_targets:
#             self.logger.info(f"📌 Off-Target at Index: {off_target['index']}")
#             self.logger.info(f"🔥 Tm: {off_target['tm']}°C")
#             self.logger.info(f"❌ Mismatches: {off_target['mismatches']}")

#     def _validate_inputs(self, seq: Seq, nucleotide_index: int, new_nucleotide: str):
#         """Validate input parameters."""
#         if not isinstance(seq, (str, Seq)):
#             raise TypeError(f"seq should be string or Seq, got {type(seq)}")
#         if not (0 <= nucleotide_index < len(seq)):
#             raise IndexError(f"Index {nucleotide_index} is out of bounds")
#         if not isinstance(new_nucleotide, str) or len(new_nucleotide) != 1:
#             raise TypeError(f"new_nucleotide should be single character")

#     def _get_part_info(self, part_num: str, primer_orientation: str, kozak: str) -> Dict[str, str]:
#         """
#         Determines part-specific key and primer name suffix based on part number and orientation.
        
#         Args:
#             part_num: The part number ('2', '3', '3a', etc.)
#             primer_orientation: Direction of the primer ('forward' or 'reverse')
#             kozak: Type of kozak sequence ('canonical' or 'mtk')
        
#         Returns:
#             Dict containing:
#                 - 'key': The key to look up in part_end_dict
#                 - 'suffix': The suffix to use in the primer name
#         """
#         # Special case for part 2 reverse with canonical kozak
#         if part_num == '2' and primer_orientation == 'reverse' and kozak == 'canonical':
#             return {
#                 'key': '2starreverse',
#                 'suffix': '2*'
#             }
        
#         # Special case for parts 3 and 3a forward with canonical kozak
#         elif part_num in ['3', '3a'] and primer_orientation == 'forward' and kozak == 'canonical':
#             return {
#                 'key': f'{part_num}starforward',
#                 'suffix': f'{part_num}*'
#             }
        
#         # Default case
#         else:
#             return {
#                 'key': f'{part_num}{primer_orientation}',
#                 'suffix': part_num
#             }

#     def _validate_primer_positions(self, left: int, right: int, seq_length: int) -> bool:
#         """Validate primer positions are within sequence bounds."""
#         return 0 <= left < seq_length and 0 <= right <= seq_length

#     def _get_modified_sequence(
#         self,
#         seq: Seq,
#         nucleotide_index: int,
#         new_nucleotide: str,
#         left: int,
#         right: int
#     ) -> str:
#         """Get the modified sequence with the new nucleotide."""
#         modified_seq = seq[:nucleotide_index] + new_nucleotide + seq[nucleotide_index + 1:]
#         return str(modified_seq[left:right])

#     def _find_binding_sequences(
#         self,
#         seq: Seq,
#         right: int,
#         left: int,
#         min_tm: float
#     ) -> Tuple[str, str]:
#         """Find binding sequences for forward and reverse primers."""
#         forward_binding = self._find_binding_sequence(seq, right, min_tm)
#         reverse_binding = self._find_binding_sequence(
#             seq.reverse_complement(),
#             left - 1,
#             min_tm
#         )
#         return forward_binding, reverse_binding

#     def _find_binding_sequence(self, seq: Seq, start: int, min_tm: float) -> str:
#         """Find the shortest binding sequence with Tm >= min_tm."""
#         from Bio.SeqUtils import MeltingTemp as mt
        
#         for i in range(15, 36):  # Typical primer length range
#             binding_seq = seq[start:start + i]
#             if mt.Tm_NN(binding_seq) >= min_tm:
#                 return str(binding_seq)
#         return ""

#     def _construct_primer_pair(
#         self,
#         seq: Seq,
#         spacer: str,
#         bsmbi_site: str,
#         six_nuc_seq: str,
#         forward_binding: str,
#         reverse_binding: str
#     ) -> Dict[str, Dict]:
#         """Construct forward and reverse primers."""
#         primers = {
#             "forward": self._construct_primer(
#                 spacer, bsmbi_site, six_nuc_seq, forward_binding, False
#             ),
#             "reverse": self._construct_primer(
#                 spacer, bsmbi_site, six_nuc_seq, reverse_binding, True
#             )
#         }
#         return primers

#     def _construct_primer(
#         self,
#         spacer: str,
#         six_nuc_seq: str,
#         binding_seq: str,
#         reverse: bool = False
#     ) -> Dict[str, str]:
#         """Construct a single primer."""
#         if reverse:
#             from Bio.Seq import Seq
#             overhang_sequence = self.bsmbi_site + str(Seq(six_nuc_seq).reverse_complement())
#         else:
#             overhang_sequence = self.bsmbi_site + six_nuc_seq

#         primer_sequence = spacer + overhang_sequence + binding_seq
        
#         return {
#             "binding_sequence": str(binding_seq),
#             "overhang_sequence": str(overhang_sequence),
#             "primer_sequence": str(primer_sequence)
#         }

#     def _validate_primer_pair(self, primers: Dict[str, Dict]) -> bool:
#         """
#         Validate both primers in a pair.
        
#         Args:
#             primers: Dictionary containing forward and reverse primer information
#             bsmbi_site: The BsmBI recognition sequence to check for
            
#         Returns:
#             bool: True if both primers have exactly one BsmBI site
#         """
#         forward_site_count = primers["forward"]["overhang_sequence"].count(self.bsmbi_site)
#         reverse_site_count = primers["reverse"]["overhang_sequence"].count(self.bsmbi_site)
#         return forward_site_count == 1 and reverse_site_count == 1

#     def _format_primer_pair(self, primers: Dict[str, Dict]) -> Tuple[Dict, Dict]:
#         """Format primer pair with additional information."""
#         from Bio.SeqUtils import MeltingTemp as mt
        
#         formatted_primers = []
#         for primer_dict in [primers["forward"], primers["reverse"]]:
#             formatted = {
#                 "primer_sequence": primer_dict["primer_sequence"],
#                 "binding_sequence": primer_dict["binding_sequence"],
#                 "overhang_sequence": primer_dict["overhang_sequence"],
#                 "tm": round(mt.Tm_NN(primer_dict["primer_sequence"]), 1)
#             }
#             formatted_primers.append(formatted)
        
#         return tuple(formatted_primers)

#     def _perform_off_target_analysis(
#         self,
#         primers: Dict[str, List[Dict]],
#         template_seq: str
#     ) -> None:
#         """Perform off-target analysis for all primers."""
#         for primer_list in ["forward", "reverse"]:
#             for primer_dict in primers[primer_list]:
#                 try:
#                     primer_dict["off_target_analysis"] = self.find_off_targets_for_primer(
#                         primer_dict["binding_sequence"],
#                         template_seq=template_seq
#                     )
#                 except Exception as e:
#                     self.logger.error(f"Error in off-target analysis: {str(e)}")
#                     primer_dict["off_target_analysis"] = []

#     def _construct_edge_primer(
#         self,
#         seq: Seq,
#         part_specific_seq: str,
#         primer_orientation: str,
#         suffix: str,
#         primer_name: Optional[str] = None
#     ) -> Tuple[str, str]:
#         """Construct an edge primer with proper naming."""
#         n_length = 20  # Default primer binding length
        
#         if primer_orientation == 'forward':
#             edge_primer = part_specific_seq + str(seq[:n_length])
#             formatted_name = f"MTK{suffix}_{primer_name}_FW" if primer_name else f"MTK{suffix}_FW"
#         elif primer_orientation == 'reverse':
#             edge_primer = part_specific_seq + str(seq[-n_length:].reverse_complement())
#             formatted_name = f"MTK{suffix}_{primer_name}_RV" if primer_name else f"MTK{suffix}_RV"
#         else:
#             raise ValueError("primer_orientation must be 'forward' or 'reverse'")
            
#         return formatted_name, edge_primer
    
#     def _design_internal_mutation_primers(
#         self,
#         sequence: Seq,
#         seq_index: int,
#         mutation_options: List,
#         part_num_left: List[str],
#         part_num_right: List[str],
#         primer_name: Optional[List[str]],
#         template_seq: Optional[str],
#         kozak: str
#     ) -> List[List[str]]:
#         """Design primers based on prioritized mutations, with assembly site placement."""
#         with self.debug_context("design_internal_mutation_primers"):
#             primer_data = []

#             # Iterate over the mutation options (prioritized mutation list)
#             for mutation in mutation_options:
#                 # Use the mutation dictionary directly (it already contains the necessary info)
#                 mutation_data = mutation  
                
#                 # Extract site index and other mutation details
#                 site_index = mutation_data.get("site_index")
#                 mutated_base = mutation_data.get("new_nt")
#                 mutation_seq = mutation_data.get("mutated_seq")
#                 # Provide a default if 'mutable_positions' is missing (adjust if needed)
#                 mutable_positions = mutation_data.get("mutable_positions", [])
                
#                 # Design primers directly based on mutation data
#                 primer_sets = self._design_primers_for_mutation(
#                     seq=sequence,
#                     mutation=mutation,
#                     part_num_left=part_num_left[seq_index],
#                     part_num_right=part_num_right[seq_index],
#                     primer_name=primer_name[seq_index] if primer_name else None,
#                     kozak=kozak,
#                     verbose=self.verbose
#                 )
                
#                 # Format and append the primer sets
#                 primer_data.extend(
#                     self.primer_selector.format_primers_for_output(primer_sets)
#                 )

#             # Add edge primers for assembly site placement
#             mtk_partend_sequences = self.utils.get_mtk_partend_sequences()
#             primer_data.extend(self._get_edge_primers(
#                 sequence, seq_index, mtk_partend_sequences,
#                 part_num_left, part_num_right, kozak, primer_name
#             ))

#             return primer_data



#     def _design_primers_for_mutation(
#         self,
#         seq: Seq,
#         mutation: Dict,
#         part_num_left: str,
#         part_num_right: str,
#         primer_name: Optional[str],
#         kozak: str,
#         verbose: bool
#     ) -> List[Dict]:
#         """Design primers based on a single mutation."""
#         # This function will handle the design of primers based on the mutation.

#         mutation_data = mutation['mutation']
#         mutated_seq = mutation_data['mutated_seq']
#         mutated_base = mutation_data['new_nt']
#         rs_index = mutation_data['rs_index']
        
#         # Design primers for the mutation
#         primer_data = []

#         # Example primer design steps (this should be replaced with your actual logic)
#         # We generate the forward and reverse primers based on mutation data:
#         forward_primer = self._generate_forward_primer(seq, mutation, part_num_left, part_num_right, kozak)
#         reverse_primer = self._generate_reverse_primer(seq, mutation, part_num_left, part_num_right, kozak)

#         # Store primers with associated metadata (e.g., primer name)
#         primer_data.append({
#             "forward": forward_primer,
#             "reverse": reverse_primer,
#             "site_index": mutation['site_index'],
#             "mutation_seq": mutated_seq,
#             "primer_name": primer_name
#         })

#         return primer_data

#     def _reverse_complement(self, seq: str) -> str:
#         """Return the reverse complement of a DNA sequence."""
#         complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
#         return "".join(complement.get(nt.upper(), nt) for nt in reversed(seq))


#     def _generate_forward_primer(
#         self,
#         seq: Seq,
#         mutation: Dict,
#         part_num_left: str,
#         part_num_right: str,
#         kozak: str,
#         primer_type: str
#     ) -> str:
#         """
#         Generate the forward primer based on mutation details and primer type.

#         For internal primers:
#         - Extract the binding sequence from the region defined by the mutation's
#             fragment reassembly range (set by the compatibility matrix).
#         - Append " GAA" and the BsmBI recognition site.
        
#         For external primers:
#         - Extract the binding sequence (as above).
#         - Prepend the overhang looked up using the mtk left part number.
        
#         The optional kozak sequence is prepended in either case.
#         """
#         # Determine the binding region.
#         if "fragment_reassembly_range" in mutation and mutation["fragment_reassembly_range"]:
#             # Expecting an object (or dict) with 'start' and 'end'
#             frag_range = mutation["fragment_reassembly_range"]
#             binding_seq = str(seq)[frag_range.start : frag_range.end + 1]
#         else:
#             # Fallback: use site_index and a default length (this should rarely happen)
#             site_index = mutation.get("site_index", 0)
#             default_length = 20  # Not used if comp matrix determined binding region
#             binding_seq = str(seq)[site_index : site_index + default_length]

#         if primer_type.lower() == "internal":
#             # For internal primers, add spacer " GAA" and BsmBI recognition site.
#             bsmbi_site = "CGTCTC"  # Example BsmBI recognition sequence
#             primer_body = binding_seq + " GAA" + bsmbi_site
#         elif primer_type.lower() == "external":
#             # For external primers, lookup the overhang based on the left part number.
#             overhang = self.utils.get_overhang(part_num_left)
#             primer_body = overhang + binding_seq
#         else:
#             raise ValueError(f"Unknown primer type: {primer_type}")

#         # Prepend the kozak sequence if provided.
#         if kozak:
#             return kozak + primer_body
#         else:
#             return primer_body


#     def _generate_reverse_primer(
#         self,
#         seq: Seq,
#         mutation: Dict,
#         part_num_left: str,
#         part_num_right: str,
#         kozak: str,
#         primer_type: str
#     ) -> str:
#         """
#         Generate the reverse primer based on mutation details and primer type.

#         For internal primers:
#         - Extract the binding sequence from the region defined by the mutation's
#             fragment reassembly range.
#         - Compute its reverse complement.
#         - Append " GAA" and the reverse complement of the BsmBI recognition site.
        
#         For external primers:
#         - Extract the binding sequence (as above) and compute its reverse complement.
#         - Append the overhang looked up using the mtk right part number.
        
#         The optional kozak sequence is prepended in either case.
#         """
#         if "fragment_reassembly_range" in mutation and mutation["fragment_reassembly_range"]:
#             frag_range = mutation["fragment_reassembly_range"]
#             binding_seq = str(seq)[frag_range.start : frag_range.end + 1]
#         else:
#             site_index = mutation.get("site_index", 0)
#             default_length = 20
#             binding_seq = str(seq)[site_index : site_index + default_length]

#         # Get reverse complement of the binding sequence.
#         binding_seq_rc = self._reverse_complement(binding_seq)

#         if primer_type.lower() == "internal":
            
#             # Compute the reverse complement of the BsmBI recognition site.
#             bsmbi_site_rc = self._reverse_complement(self.bsmbi_site)
#             primer_body = binding_seq_rc + " GAA" + bsmbi_site_rc
#         elif primer_type.lower() == "external":
#             overhang = self.utils.get_overhang(part_num_right)
#             primer_body = binding_seq_rc + overhang
#         else:
#             raise ValueError(f"Unknown primer type: {primer_type}")

#         if kozak:
#             return kozak + primer_body
#         else:
#             return primer_body


#     def _get_edge_primers(
#         self,
#         sequence: Seq,
#         index: int,
#         part_end_dict: Dict[str, str],
#         part_num_left: List[str],
#         part_num_right: List[str],
#         kozak: str,
#         primer_name: Optional[List[str]]
#     ) -> List[List[str]]:
#         """Generate edge primers for a sequence."""
#         with self.debug_context("get_edge_primers"):
#             primer_data = []
#             left_name, left_primer = self.primer_designer.generate_GG_edge_primers(
#                 sequence, part_end_dict, part_num_left[index], "forward", kozak, primer_name
#             )
#             right_name, right_primer = self.primer_designer.generate_GG_edge_primers(
#                 sequence, part_end_dict, part_num_right[index], "reverse", kozak, primer_name
#             )
#             primer_data.append([left_name, left_primer, f"Amplicon_{index + 1}"])
#             primer_data.append([right_name, right_primer, f"Amplicon_{index + 1}"])
#             return primer_data