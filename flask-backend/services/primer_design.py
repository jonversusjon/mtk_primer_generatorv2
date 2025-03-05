from typing import Dict, List, Optional, Any, Tuple, Union
import primer3
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from .utils import GoldenGateUtils
import numpy as np
from config.logging_config import logger
from services.base import debug_context


class PrimerDesigner:
    """
    Handles primer design for Golden Gate assembly.
    """

    def __init__(
        self,
        kozak: str = "MTK",
        verbose: bool = False,
    ):
        self.logger = logger.getChild("PrimerDesigner")
        self.utils = GoldenGateUtils()

        self.verbose = verbose
        self.state = {
            'current_operation': '',
            'primers_designed': 0,
            'current_mutation': None
        }

        self.kozak = kozak
        self.part_end_dict = self.utils.get_mtk_partend_sequences()

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

    def design_mutation_primers(
        self,
        full_sequence: str,
        mutation_sets: List[Dict],
        comp_matrices: List[np.ndarray],
        primer_name: Optional[str] = None
    ) -> Optional[Dict[str, Dict[str, str]]]:
        """
        Design primers for mutation sets based on compatibility matrices.
        """
        with debug_context("design_mutation_primers"):
            for set_index, mutation_set in enumerate(mutation_sets):
                self.state['current_mutation'] = set_index

                # Find valid overhang combinations
                valid_coords = np.argwhere(comp_matrices[set_index] == 1)
                if valid_coords.size == 0:
                    continue

                # Apply mutations to sequence
                mutated_seq = self._apply_mutations(
                    full_sequence, mutation_set)

                # Design primers using first valid combination
                primers = self._construct_mutation_primers(
                    mutated_seq,
                    mutation_set,
                    valid_coords[0],
                    primer_name
                )

                if primers:
                    self.state['primers_designed'] += len(primers)
                    return primers

            return None

    def _apply_mutations(self, target_seq: str, mutation_set: Dict) -> str:
        """Apply a set of mutations to the target sequence."""
        with debug_context("apply_mutations"):
            mutated_seq = list(target_seq)

            for mut in mutation_set.values():
                pos = mut["position"] - 1
                alt_seq = mut["alternative_sequence"]
                orig_seq = mut["original_sequence"]

                mutated_seq[pos:pos + len(orig_seq)] = alt_seq

                if self.verbose:
                    self._log_mutation(target_seq, ''.join(
                        mutated_seq), pos, len(orig_seq))

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
        with debug_context("construct_mutation_primers"):
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

                if forward and reverse:
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
        mutation_length = len(mut["original_sequence"])

        if is_reverse:
            # For the reverse primer, anneal immediately downstream of the mutation.
            start_pos = position + mutation_length
            end_pos = start_pos + annealing_length
            if end_pos > len(sequence):
                return None
            binding = sequence[start_pos:end_pos]
            binding = str(Seq(binding).reverse_complement())
            # Use bottom extended sequence for reverse primer
            extended_seq = str(
                mut['overhangs']['bottom_extended_sequences'][selected_coord])
        else:
            # For the forward primer, anneal immediately upstream of the mutation.
            end_pos = position
            start_pos = end_pos - annealing_length
            if start_pos < 0:
                return None
            binding = sequence[start_pos:end_pos]
            # Use top extended sequence for forward primer
            extended_seq = str(
                mut['overhangs']['top_extended_sequences'][selected_coord])

        # Construct the complete primer: spacer + BsmbI site + extended (overhang) + binding region
        primer_seq = self.spacer + self.bsmbi_site + extended_seq + binding

        if primer_name:
            suffix = "_RV" if is_reverse else "_FW"
            return (f"{primer_name}{suffix}", primer_seq)

        return primer_seq

    def _analyze_off_targets_for_all_primers(
        self,
        primer_results: Dict[str, Any],
        template_seq: str
    ) -> None:
        """Perform off-target analysis for all primers."""
        with debug_context("analyze_off_targets"):
            for primer_type, primers in primer_results.items():
                if primer_type == 'mutation_primers':
                    for site_primers in primers.values():
                        self._analyze_primer_pair(
                            site_primers, template_seq)
                elif primer_type == 'edge_primers':
                    self._analyze_primer_pair(primers, template_seq)

    def _analyze_primer_pair(
        self,
        primer_pair: Dict[str, Any],
        template_seq: str
    ) -> None:
        """Analyze a pair of primers for off-target binding."""
        for orientation in ['forward', 'reverse']:
            if orientation in primer_pair:
                if isinstance(primer_pair[orientation], tuple):
                    binding_seq = primer_pair[orientation][1]
                else:
                    binding_seq = primer_pair[orientation]

                off_targets = self._find_off_targets(binding_seq, template_seq)
                primer_pair[f'{orientation}_off_targets'] = off_targets

    def _log_mutation(self, original_seq, mutated_seq, position, length):
        """Log details of each mutation."""
        for i in range(position, position + length):
            old_base = original_seq[i] if i < len(original_seq) else "N/A"
            new_base = mutated_seq[i] if i < len(mutated_seq) else "N/A"
            if old_base != new_base:
                logger.info(
                    f"Mutated base at position {i+1}: {old_base} -> {new_base}")

    def _construct_edge_primer(
            self,
            sequence: str,
            part_num: str,
            forward_or_reverse: str,
            primer_name: str):
        part_end_dict = self.part_end_dict

        name = f"{primer_name}_{forward_or_reverse}"
        seq = f"{primer_name} placeholder primer sequence"

        return name, seq

    def _find_off_targets(self, binding_seq: str, template_seq: str):
        return []

    def _calculate_gc_content(self, sequence: str) -> float:
        """
        Calculates the GC content percentage of a DNA sequence.
        """
        if not sequence:
            return 0.0
        sequence = sequence.upper()
        gc_count = sequence.count("G") + sequence.count("C")
        return (gc_count / len(sequence)) * 100.0

    def _calculate_tm(self, sequence: str) -> float:
        """
        Calculates the approximate melting temperature (Tm) of a DNA sequence.
        """
        if not sequence:
            return 0.0
        sequence = sequence.upper()
        length = len(sequence)
        a_count = sequence.count("A")
        t_count = sequence.count("T")
        g_count = sequence.count("G")
        c_count = sequence.count("C")

        if length < 14:
            tm = (a_count + t_count) * 2 + (g_count + c_count) * 4
        else:
            tm = 64.9 + (41 * (g_count + c_count - 16.4)) / length
        return round(tm, 2)

    def generate_GG_edge_primers(self, idx, sequence, overhang_5_prime, overhang_3_prime, primer_name):
        """
        Generate Golden Gate edge primers for a sequence.
        """
        # Convert to string for consistent handling
        seq_str = str(sequence)
        seq_length = len(seq_str)

        # Default primer length
        primer_length = 20  # Excluding overhangs and recognition site

        # Define enzyme recognition sequences
        enzyme_rec_seq = "GGTCTC"  # BsaI

        # 5' Primer (forward)
        forward_primer_binding = seq_str[:primer_length]
        forward_primer = overhang_5_prime + enzyme_rec_seq + "A" + forward_primer_binding

        # 3' Primer (reverse)
        reverse_primer_binding = str(
            Seq(seq_str[-primer_length:]).reverse_complement())
        reverse_primer = overhang_3_prime + enzyme_rec_seq + "A" + reverse_primer_binding

        # Calculate properties
        forward_tm = self._calculate_tm(forward_primer_binding)
        reverse_tm = self._calculate_tm(reverse_primer_binding)

        forward_gc = self._calculate_gc_content(forward_primer_binding)
        reverse_gc = self._calculate_gc_content(reverse_primer_binding)

        # Create result dictionary
        primers = {
            "forward_primer": {
                "sequence": forward_primer,
                "binding_region": forward_primer_binding,
                "tm": forward_tm,
                "gc_content": forward_gc,
                "length": len(forward_primer)
            },
            "reverse_primer": {
                "sequence": reverse_primer,
                "binding_region": reverse_primer_binding,
                "tm": reverse_tm,
                "gc_content": reverse_gc,
                "length": len(reverse_primer)
            },
            "product_size": seq_length
        }

        return primers
