# services/primer_design.py
from typing import Dict, List, Optional, Any, Tuple
import primer3
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from .base import GoldenGateDesigner


class PrimerDesigner(GoldenGateDesigner):
    def __init__(self, verbose: bool = False):
        super().__init__(verbose=verbose)
        self.state = {
            'current_operation': '',
            'primers_designed': 0,
            'current_mutation': None
        }
        # Default parameters
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

    def find_off_targets_for_primer(
        self,
        binding_seq: str,
        template_seq: Optional[str] = None,
        **kwargs
    ) -> List[Dict[str, float]]:
        """Find potential off-target binding sites using primer3."""
        with self.debug_context("find_off_targets"):
            if not isinstance(binding_seq, str):
                raise TypeError(f"Expected string for binding_seq, got {type(binding_seq)}")

            if not template_seq:
                self.logger.warning(f"Skipping off-target analysis for primer: {binding_seq}")
                return []

            # Merge default parameters with provided kwargs
            params = {**self.default_params, **kwargs}
            
            return self._analyze_off_targets(
                binding_seq.upper(),
                template_seq.upper(),
                params
            )

    def _analyze_off_targets(
        self,
        binding_seq: str,
        template_seq: str,
        params: Dict
    ) -> List[Dict[str, float]]:
        """Internal method to analyze off-target binding."""
        off_targets = []
        primer_len = len(binding_seq)
        
        for i in range(len(template_seq) - primer_len + 1):
            window = template_seq[i:i + primer_len]
            
            # Check 3' end match
            match_3p = sum(1 for a, b in zip(
                binding_seq[-params['min_3p_match']:],
                window[-params['min_3p_match']:]
            ) if a == b)
            
            if match_3p < params['min_3p_match']:
                continue

            # Calculate heterodimer properties
            try:
                result = primer3.bindings.calc_heterodimer(
                    seq1=binding_seq,
                    seq2=window,
                    mv_conc=params['mv_conc'],
                    dv_conc=params['dv_conc'],
                    dntp_conc=params['dntp_conc'],
                    dna_conc=params['dna_conc']
                )

                if result.tm >= params['tm_threshold']:
                    mismatches = sum(1 for a, b in zip(binding_seq, window) if a != b)
                    if mismatches <= params['max_mismatches']:
                        off_targets.append({
                            "index": i,
                            "tm": round(result.tm, 2),
                            "mismatches": mismatches
                        })
            except Exception as e:
                self.logger.error(f"Error in heterodimer calculation: {str(e)}")

        if off_targets:
            self._log_off_target_results(binding_seq, off_targets)
            
        return off_targets

    def design_primers_for_mutation(
        self,
        seq: Seq,
        nucleotide_index: int,
        new_nucleotide: str,
        spacer: str,
        bsmbi_site: str,
        min_tm: float = 57,
        template_seq: Optional[str] = None,
    ) -> Dict[str, List[Dict]]:
        """Design primers for a specific mutation."""
        with self.debug_context("design_primers_for_mutation"):
            self._validate_inputs(seq, nucleotide_index, new_nucleotide)
            
            self.state['current_operation'] = 'primer_design'
            self.state['current_mutation'] = {
                'position': nucleotide_index,
                'new_base': new_nucleotide
            }
            
            primers = {"forward": [], "reverse": []}
            
            for shift in range(6):
                primer_pair = self._design_primer_pair(
                    seq, nucleotide_index, new_nucleotide,
                    spacer, bsmbi_site, shift, min_tm
                )
                if primer_pair:
                    primers["forward"].append(primer_pair[0])
                    primers["reverse"].append(primer_pair[1])
                    self.state['primers_designed'] += 2

            if template_seq:
                self._perform_off_target_analysis(primers, template_seq)

            return primers

    def _design_primer_pair(
        self,
        seq: Seq,
        nucleotide_index: int,
        new_nucleotide: str,
        spacer: str,
        bsmbi_site: str,
        shift: int,
        min_tm: float
    ) -> Optional[Tuple[Dict, Dict]]:
        """Design a single pair of primers."""
        left = nucleotide_index - (6 - 1) + shift
        right = nucleotide_index + 1 + shift

        if not self._validate_primer_positions(left, right, len(seq)):
            return None

        six_nuc_seq = self._get_modified_sequence(
            seq, nucleotide_index, new_nucleotide, left, right
        )

        binding_seqs = self._find_binding_sequences(seq, right, left, min_tm)
        if not all(binding_seqs):
            return None

        primers = self._construct_primer_pair(
            seq, spacer, bsmbi_site, six_nuc_seq, *binding_seqs
        )
        
        if self._validate_primer_pair(primers):
            return self._format_primer_pair(primers)
            
        return None

    def generate_GG_edge_primers(
        self,
        seq: Seq,
        part_end_dict: Dict[str, str],
        part_num: str,
        primer_orientation: str,
        kozak: str,
        primer_name: Optional[str] = None
    ) -> Tuple[Optional[str], Optional[str]]:
        """Generate Golden Gate edge primers."""
        with self.debug_context("generate_edge_primers"):
            try:
                seq = Seq(seq)
                part_info = self._get_part_info(part_num, primer_orientation, kozak)
                
                part_specific_seq = part_end_dict.get(part_info['key'])
                if not part_specific_seq:
                    raise ValueError(f"No sequence found for {part_info['key']}")

                primer = self._construct_edge_primer(
                    seq, part_specific_seq, primer_orientation, 
                    part_info['suffix'], primer_name
                )
                
                return primer
                
            except Exception as e:
                self.logger.error(f"Error generating edge primer: {str(e)}")
                return None, None

    def _log_off_target_results(self, binding_seq: str, off_targets: List[Dict]):
        """Log off-target analysis results."""
        self.logger.info("\n🔍 Off-Target Primer Found")
        self.logger.info(f"Primer Binding Sequence: {binding_seq}")
        self.logger.info(f"Total Off-Targets: {len(off_targets)}")
        
        for off_target in off_targets:
            self.logger.info(f"📌 Off-Target at Index: {off_target['index']}")
            self.logger.info(f"🔥 Tm: {off_target['tm']}°C")
            self.logger.info(f"❌ Mismatches: {off_target['mismatches']}")

    def _validate_inputs(self, seq: Seq, nucleotide_index: int, new_nucleotide: str):
        """Validate input parameters."""
        if not isinstance(seq, (str, Seq)):
            raise TypeError(f"seq should be string or Seq, got {type(seq)}")
        if not (0 <= nucleotide_index < len(seq)):
            raise IndexError(f"Index {nucleotide_index} is out of bounds")
        if not isinstance(new_nucleotide, str) or len(new_nucleotide) != 1:
            raise TypeError(f"new_nucleotide should be single character")

    def _get_part_info(self, part_num: str, primer_orientation: str, kozak: str) -> Dict[str, str]:
        """
        Determines part-specific key and primer name suffix based on part number and orientation.
        
        Args:
            part_num: The part number ('2', '3', '3a', etc.)
            primer_orientation: Direction of the primer ('forward' or 'reverse')
            kozak: Type of kozak sequence ('canonical' or 'mtk')
        
        Returns:
            Dict containing:
                - 'key': The key to look up in part_end_dict
                - 'suffix': The suffix to use in the primer name
        """
        # Special case for part 2 reverse with canonical kozak
        if part_num == '2' and primer_orientation == 'reverse' and kozak == 'canonical':
            return {
                'key': '2starreverse',
                'suffix': '2*'
            }
        
        # Special case for parts 3 and 3a forward with canonical kozak
        elif part_num in ['3', '3a'] and primer_orientation == 'forward' and kozak == 'canonical':
            return {
                'key': f'{part_num}starforward',
                'suffix': f'{part_num}*'
            }
        
        # Default case
        else:
            return {
                'key': f'{part_num}{primer_orientation}',
                'suffix': part_num
            }

    def _validate_primer_positions(self, left: int, right: int, seq_length: int) -> bool:
        """Validate primer positions are within sequence bounds."""
        return 0 <= left < seq_length and 0 <= right <= seq_length

    def _get_modified_sequence(
        self,
        seq: Seq,
        nucleotide_index: int,
        new_nucleotide: str,
        left: int,
        right: int
    ) -> str:
        """Get the modified sequence with the new nucleotide."""
        modified_seq = seq[:nucleotide_index] + new_nucleotide + seq[nucleotide_index + 1:]
        return str(modified_seq[left:right])

    def _find_binding_sequences(
        self,
        seq: Seq,
        right: int,
        left: int,
        min_tm: float
    ) -> Tuple[str, str]:
        """Find binding sequences for forward and reverse primers."""
        forward_binding = self._find_binding_sequence(seq, right, min_tm)
        reverse_binding = self._find_binding_sequence(
            seq.reverse_complement(),
            left - 1,
            min_tm
        )
        return forward_binding, reverse_binding

    def _find_binding_sequence(self, seq: Seq, start: int, min_tm: float) -> str:
        """Find the shortest binding sequence with Tm >= min_tm."""
        from Bio.SeqUtils import MeltingTemp as mt
        
        for i in range(15, 36):  # Typical primer length range
            binding_seq = seq[start:start + i]
            if mt.Tm_NN(binding_seq) >= min_tm:
                return str(binding_seq)
        return ""

    def _construct_primer_pair(
        self,
        seq: Seq,
        spacer: str,
        bsmbi_site: str,
        six_nuc_seq: str,
        forward_binding: str,
        reverse_binding: str
    ) -> Dict[str, Dict]:
        """Construct forward and reverse primers."""
        primers = {
            "forward": self._construct_primer(
                spacer, bsmbi_site, six_nuc_seq, forward_binding, False
            ),
            "reverse": self._construct_primer(
                spacer, bsmbi_site, six_nuc_seq, reverse_binding, True
            )
        }
        return primers

    def _construct_primer(
        self,
        spacer: str,
        bsmbi_site: str,
        six_nuc_seq: str,
        binding_seq: str,
        reverse: bool = False
    ) -> Dict[str, str]:
        """Construct a single primer."""
        if reverse:
            from Bio.Seq import Seq
            overhang_sequence = bsmbi_site + str(Seq(six_nuc_seq).reverse_complement())
        else:
            overhang_sequence = bsmbi_site + six_nuc_seq

        primer_sequence = spacer + overhang_sequence + binding_seq
        
        return {
            "binding_sequence": str(binding_seq),
            "overhang_sequence": str(overhang_sequence),
            "primer_sequence": str(primer_sequence)
        }

    def _validate_primer_pair(self, primers: Dict[str, Dict], bsmbi_site: str) -> bool:
        """
        Validate both primers in a pair.
        
        Args:
            primers: Dictionary containing forward and reverse primer information
            bsmbi_site: The BsmBI recognition sequence to check for
            
        Returns:
            bool: True if both primers have exactly one BsmBI site
        """
        forward_site_count = primers["forward"]["overhang_sequence"].count(bsmbi_site)
        reverse_site_count = primers["reverse"]["overhang_sequence"].count(bsmbi_site)
        return forward_site_count == 1 and reverse_site_count == 1

    def _format_primer_pair(self, primers: Dict[str, Dict]) -> Tuple[Dict, Dict]:
        """Format primer pair with additional information."""
        from Bio.SeqUtils import MeltingTemp as mt
        
        formatted_primers = []
        for primer_dict in [primers["forward"], primers["reverse"]]:
            formatted = {
                "primer_sequence": primer_dict["primer_sequence"],
                "binding_sequence": primer_dict["binding_sequence"],
                "overhang_sequence": primer_dict["overhang_sequence"],
                "tm": round(mt.Tm_NN(primer_dict["primer_sequence"]), 1)
            }
            formatted_primers.append(formatted)
        
        return tuple(formatted_primers)

    def _perform_off_target_analysis(
        self,
        primers: Dict[str, List[Dict]],
        template_seq: str
    ) -> None:
        """Perform off-target analysis for all primers."""
        for primer_list in ["forward", "reverse"]:
            for primer_dict in primers[primer_list]:
                try:
                    primer_dict["off_target_analysis"] = self.find_off_targets_for_primer(
                        primer_dict["binding_sequence"],
                        template_seq=template_seq
                    )
                except Exception as e:
                    self.logger.error(f"Error in off-target analysis: {str(e)}")
                    primer_dict["off_target_analysis"] = []

    def _construct_edge_primer(
        self,
        seq: Seq,
        part_specific_seq: str,
        primer_orientation: str,
        suffix: str,
        primer_name: Optional[str] = None
    ) -> Tuple[str, str]:
        """Construct an edge primer with proper naming."""
        n_length = 20  # Default primer binding length
        
        if primer_orientation == 'forward':
            edge_primer = part_specific_seq + str(seq[:n_length])
            formatted_name = f"MTK{suffix}_{primer_name}_FW" if primer_name else f"MTK{suffix}_FW"
        elif primer_orientation == 'reverse':
            edge_primer = part_specific_seq + str(seq[-n_length:].reverse_complement())
            formatted_name = f"MTK{suffix}_{primer_name}_RV" if primer_name else f"MTK{suffix}_RV"
        else:
            raise ValueError("primer_orientation must be 'forward' or 'reverse'")
            
        return formatted_name, edge_primer