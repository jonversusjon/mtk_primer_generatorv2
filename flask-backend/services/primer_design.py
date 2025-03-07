from typing import Dict, List, Optional, Any, Tuple, Union
import primer3
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from .utils import GoldenGateUtils
import numpy as np
from config.logging_config import logger
from services.base import debug_context
from services.debug.debug_utils import MutationDebugger, debug_function, visualize_matrix

from dataclasses import dataclass, field
from typing import Optional, List, Dict
from models.primer import Primer, MutationPrimer
from functools import wraps
import logging


class PrimerDesigner:
    """
    Handles primer design for Golden Gate assembly.
    """

    def __init__(
        self,
        kozak: str = "MTK",
        verbose: bool = False,
        debug: bool = False,
    ):
        self.logger = logger.getChild("PrimerDesigner")
        self.utils = GoldenGateUtils()

        self.verbose = verbose
        self.debug = debug

        # Initialize debugger if debug mode is enabled
        if self.debug:
            self.debugger = MutationDebugger(
                parent_logger=logger,
                use_custom_format=True
            )
            # Manually set propagate to False to avoid duplicate logs
            if hasattr(self.debugger.logger, 'propagate'):
                self.debugger.logger.propagate = False

            self.logger.info("Debug mode enabled for PrimerDesigner")

            # Replace our debug_function decorator with an actual implementation
            self._setup_debug_wrappers()

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

        # Debug validation if debug mode is on
        if self.debug:
            self.debugger.validate(
                self.part_end_dict is not None,
                "MTK part end sequences loaded successfully",
                {"kozak": self.kozak}
            )

    def _setup_debug_wrappers(self):
        """Set up debug function wrappers dynamically"""
        # Replace our debug_function decorator with an actual implementation
        def apply_debug_wrapper(method):
            @wraps(method)
            def wrapped(self, *args, **kwargs):
                # Log function start
                func_name = method.__name__
                self.debugger.log_function_start(func_name, args, kwargs)

                # Run the function
                try:
                    result = method(self, *args, **kwargs)

                    # Basic validation
                    if result is not None:
                        if isinstance(result, tuple):
                            for i, item in enumerate(result):
                                self.debugger.validate(
                                    item is not None,
                                    f"Return value {i+1} is not None"
                                )
                        else:
                            self.debugger.validate(
                                result is not None, "Return value is not None")

                    # Log function end
                    self.debugger.log_function_end(func_name, result)
                    return result

                except Exception as e:
                    self.debugger.log_error(
                        f"Exception in {func_name}: {str(e)}")
                    raise

            return wrapped

        # Apply debug wrapper to all methods decorated with @debug_function
        for attr_name in dir(self.__class__):
            if attr_name.startswith('__'):
                continue

            attr = getattr(self.__class__, attr_name)
            if callable(attr) and hasattr(attr, '__wrapped__'):
                setattr(self, attr_name, apply_debug_wrapper(
                    attr).__get__(self, self.__class__))

    @debug_function
    def design_mutation_primers(
        self,
        full_sequence: str,
        mutation_sets: List[Dict],
        comp_matrices: List[np.ndarray],
        primer_name: Optional[str] = None
    ) -> Optional[List[MutationPrimer]]:
        """
        Design primers for mutation sets based on compatibility matrices.
        Returns a list of MutationPrimer objects.
        """
        if self.debug:
            self.debugger.validate(
                isinstance(full_sequence, str) and len(full_sequence) > 0,
                "Input sequence is valid",
                {"sequence_length": len(full_sequence)}
            )

            self.debugger.validate(
                isinstance(mutation_sets, list) and len(mutation_sets) > 0,
                f"Received {len(mutation_sets)} mutation sets",
                {"first_set_sites": list(
                    mutation_sets[0].keys()) if mutation_sets else None}
            )

            self.debugger.validate(
                isinstance(comp_matrices, list) and len(
                    comp_matrices) == len(mutation_sets),
                f"Received {len(comp_matrices)} compatibility matrices",
                {"first_matrix_shape":
                    comp_matrices[0].shape if comp_matrices else None}
            )

        with debug_context("design_mutation_primers"):
            for set_index, mutation_set in enumerate(mutation_sets):
                self.state['current_mutation'] = set_index

                if self.debug:
                    self.debugger.log_step(
                        "Process Mutation Set",
                        f"Processing mutation set {set_index+1}/{len(mutation_sets)}",
                        {"sites": list(mutation_set.keys())}
                    )

                # Find valid overhang combinations
                valid_coords = np.argwhere(comp_matrices[set_index] == 1)

                if self.debug:
                    valid_combinations = np.count_nonzero(
                        comp_matrices[set_index])
                    self.debugger.validate(
                        valid_coords.size > 0,
                        f"Found {valid_combinations} valid overhang combinations",
                        {"matrix_size": comp_matrices[set_index].size}
                    )

                    if valid_coords.size > 0:
                        self.debugger.log_step(
                            "Matrix Visualization",
                            f"Compatibility matrix for set {set_index+1}",
                            visualize_matrix(comp_matrices[set_index])
                        )

                if valid_coords.size == 0:
                    if self.debug:
                        self.debugger.log_warning(
                            f"No valid overhang combinations for set {set_index+1}")
                    continue

                # Apply mutations to the sequence
                mutated_seq = self._apply_mutations(
                    full_sequence, mutation_set)

                if self.debug:
                    self.debugger.validate(
                        mutated_seq is not None and len(
                            mutated_seq) == len(full_sequence),
                        f"Successfully applied mutations to sequence",
                        {"original_length": len(
                            full_sequence), "mutated_length": len(mutated_seq)}
                    )

                # Construct MutationPrimer objects using the first valid combination
                mutation_primers = self._construct_mutation_primers(
                    mutated_seq,
                    mutation_set,
                    valid_coords[0],
                    primer_name
                )

                if self.debug:
                    self.debugger.validate(
                        mutation_primers is not None,
                        f"Successfully constructed mutation primers",
                        {"primer_count": len(mutation_primers)
                         if mutation_primers else 0}
                    )

                if mutation_primers:
                    return mutation_primers

            if self.debug:
                self.debugger.log_warning(
                    "Failed to design primers for any mutation set")
            return None

    @debug_function
    def _apply_mutations(self, target_seq: str, mutation_set: Dict) -> str:
        """
        Apply a set of mutations to the target sequence.
        """
        if self.debug:
            self.debugger.log_step(
                "Apply Mutations",
                f"Applying {len(mutation_set)} mutations to sequence",
                {"sequence_length": len(target_seq)}
            )

        with debug_context("apply_mutations"):
            mutated_seq = list(target_seq)

            for site_key, mut in mutation_set.items():
                pos = mut["position"] - 1
                alt_seq = mut["alternative_sequence"]
                orig_seq = mut["original_sequence"]

                if self.debug:
                    self.debugger.log_step(
                        "Mutation Details",
                        f"Applying mutation at site {site_key}, position {pos+1}",
                        {
                            "original": orig_seq,
                            "alternative": alt_seq,
                            "context": target_seq[max(0, pos-5):min(len(target_seq), pos+len(orig_seq)+5)]
                        }
                    )

                mutated_seq[pos:pos + len(orig_seq)] = alt_seq

                if self.verbose or self.debug:
                    self._log_mutation(target_seq, ''.join(
                        mutated_seq), pos, len(orig_seq))

            result = ''.join(mutated_seq)

            if self.debug:
                self.debugger.validate(
                    len(result) == len(target_seq),
                    f"Mutation applied successfully, sequence length preserved",
                    {"original_length": len(target_seq),
                     "mutated_length": len(result)}
                )

            return result

    @debug_function
    def _construct_mutation_primers(
        self,
        mutated_seq: str,
        mutation_set: Dict,
        selected_coord: np.ndarray,
        primer_name: Optional[str] = None,
        annealing_length: int = 18
    ) -> List[MutationPrimer]:
        """
        Construct MutationPrimer objects for a mutation set using selected coordinates.
        """
        if self.debug:
            self.debugger.log_step(
                "Construct Primers",
                f"Constructing primers with annealing length {annealing_length}",
                {"selected_coord": selected_coord}
            )

        with debug_context("construct_mutation_primers"):
            mutation_primers = []

            for site_key, mut in mutation_set.items():
                position = mut["position"] - 1

                if self.debug:
                    self.debugger.log_step(
                        "Design Primer Pair",
                        f"Designing primers for site {site_key} at position {position+1}",
                        {"mutation": mut["alternative_sequence"]}
                    )

                # Design forward and reverse primers using the new helper
                forward_tuple = self._construct_primer(
                    mutated_seq,
                    position,
                    mut,
                    selected_coord[0],
                    annealing_length,
                    is_reverse=False,
                    primer_name=primer_name
                )

                if self.debug and forward_tuple:
                    self.debugger.log_step(
                        "Forward Primer",
                        f"Created forward primer for {site_key}",
                        {"name": forward_tuple[0],
                            "sequence": forward_tuple[1]}
                    )

                reverse_tuple = self._construct_primer(
                    mutated_seq,
                    position,
                    mut,
                    selected_coord[0],
                    annealing_length,
                    is_reverse=True,
                    primer_name=primer_name
                )

                if self.debug and reverse_tuple:
                    self.debugger.log_step(
                        "Reverse Primer",
                        f"Created reverse primer for {site_key}",
                        {"name": reverse_tuple[0],
                            "sequence": reverse_tuple[1]}
                    )

                if forward_tuple and reverse_tuple:
                    forward_primer = Primer(
                        name=forward_tuple[0], sequence=forward_tuple[1])
                    reverse_primer = Primer(
                        name=reverse_tuple[0], sequence=reverse_tuple[1])
                    mutation_primer = MutationPrimer(
                        site=mut["site"],
                        position=mut["position"],
                        forward=forward_primer,
                        reverse=reverse_primer,
                        mutation_info=mut
                    )
                    mutation_primers.append(mutation_primer)

                    if self.debug:
                        # Calculate and log primer properties
                        f_binding = forward_tuple[1].split(
                            self.bsmbi_site + mut['overhangs']['top_extended_sequences'][selected_coord[0]])[1]
                        r_binding = reverse_tuple[1].split(
                            self.bsmbi_site + mut['overhangs']['bottom_extended_sequences'][selected_coord[0]])[1]

                        f_tm = self._calculate_tm(f_binding)
                        r_tm = self._calculate_tm(r_binding)
                        f_gc = self._calculate_gc_content(f_binding)
                        r_gc = self._calculate_gc_content(r_binding)

                        self.debugger.log_step(
                            "Primer Properties",
                            f"Primer properties for {site_key}",
                            {
                                "forward_tm": f_tm,
                                "reverse_tm": r_tm,
                                "forward_gc": f_gc,
                                "reverse_gc": r_gc,
                                "forward_length": len(f_binding),
                                "reverse_length": len(r_binding)
                            }
                        )
                else:
                    if self.debug:
                        self.debugger.log_warning(
                            f"Failed to create primer pair for site {site_key}")

            if self.debug:
                self.debugger.validate(
                    len(mutation_primers) > 0,
                    f"Successfully created {len(mutation_primers)} primer pairs",
                    {"expected": len(mutation_set)}
                )

            return mutation_primers

    @debug_function
    def _construct_primer(
        self,
        sequence: str,
        position: int,
        mut: Dict,
        selected_coord: int,
        annealing_length: int,
        is_reverse: bool,
        primer_name: Optional[str] = None
    ) -> Optional[Tuple[str, str]]:
        """
        Construct a primer that both incorporates the mutation directly and ensures 
        the mutation is within the BsmBI-generated overhang.
        """
        direction = "reverse" if is_reverse else "forward"

        if self.debug:
            self.debugger.log_step(
                f"Construct {direction.capitalize()} Primer",
                f"Designing {direction} primer at position {position+1}",
                {
                    "annealing_length": annealing_length,
                    "selected_coord": selected_coord
                }
            )

        original_seq = mut["original_sequence"]
        alternative_seq = mut["alternative_sequence"]
        mutation_length = len(original_seq)

        # Set default primer name if none provided
        if not primer_name:
            primer_name = f"Mut_{mut['site']}"
        primer_suffix = "_RV" if is_reverse else "_FW"
        primer_name = f"{primer_name}{primer_suffix}"

        # Get the appropriate extended sequence from overhangs
        if is_reverse:
            extended_seq = str(
                mut['overhangs']['bottom_extended_sequences'][selected_coord])

            # Determine optimal binding site that crosses the mutation
            # For reverse primers, we want to start binding from a position that includes the mutation
            binding_start = position - annealing_length // 2  # Start before the mutation
            binding_end = position + mutation_length + \
                annealing_length // 2  # End after the mutation

            if binding_start < 0 or binding_end > len(sequence):
                if self.debug:
                    self.debugger.log_warning(
                        f"Invalid binding region for reverse primer: [{binding_start}, {binding_end}]",
                        {"sequence_length": len(sequence)}
                    )
                return None

            # Create a binding region that spans the mutation
            # But incorporate the mutated sequence instead of the original
            binding_before_mutation = sequence[binding_start:position]
            binding_after_mutation = sequence[position +
                                              mutation_length:binding_end]

            # Construct binding sequence with the mutation incorporated
            binding = binding_before_mutation + alternative_seq + binding_after_mutation
            binding = str(Seq(binding).reverse_complement())

        else:
            extended_seq = str(
                mut['overhangs']['top_extended_sequences'][selected_coord])

            # Determine optimal binding site that crosses the mutation
            binding_start = position - annealing_length // 2  # Start before the mutation
            binding_end = position + mutation_length + \
                annealing_length // 2  # End after the mutation

            if binding_start < 0 or binding_end > len(sequence):
                if self.debug:
                    self.debugger.log_warning(
                        f"Invalid binding region for forward primer: [{binding_start}, {binding_end}]",
                        {"sequence_length": len(sequence)}
                    )
                return None

            # Create a binding region that spans the mutation
            # But incorporate the mutated sequence instead of the original
            binding_before_mutation = sequence[binding_start:position]
            binding_after_mutation = sequence[position +
                                              mutation_length:binding_end]

            # Construct binding sequence with the mutation incorporated
            binding = binding_before_mutation + alternative_seq + binding_after_mutation

        # Construct the complete primer: spacer + BsmbI site + extended sequence (which includes mutation) + binding region
        primer_seq = self.spacer + self.bsmbi_site + extended_seq + binding

        if self.debug:
            # Validate the constructed primer
            self.debugger.validate(
                len(primer_seq) > len(self.spacer +
                                      self.bsmbi_site + extended_seq),
                f"Primer contains valid binding region",
                {
                    "spacer": self.spacer,
                    "enzyme_site": self.bsmbi_site,
                    "extended_seq": extended_seq,
                    "binding_region": binding,
                    "total_length": len(primer_seq)
                }
            )

            # Check if the mutation is incorporated in the primer
            mutation_in_primer = alternative_seq in binding if not is_reverse else alternative_seq in str(
                Seq(binding).reverse_complement())
            self.debugger.validate(
                mutation_in_primer,
                f"Mutation is incorporated in the primer binding region",
                {"mutation": alternative_seq}
            )

        return (primer_name, primer_seq)

    @debug_function
    def _analyze_off_targets_for_all_primers(
        self,
        primer_results: Dict[str, Any],
        template_seq: str
    ) -> None:
        """Perform off-target analysis for all primers."""
        if self.debug:
            self.debugger.log_step(
                "Off-Target Analysis",
                f"Analyzing off-targets for {len(primer_results)} primer types",
                {"template_length": len(template_seq)}
            )

        with debug_context("analyze_off_targets"):
            for primer_type, primers in primer_results.items():
                if self.debug:
                    self.debugger.log_step(
                        "Analyze Primer Type",
                        f"Processing {primer_type}"
                    )

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

                if self.debug:
                    self.debugger.log_step(
                        f"{orientation.capitalize()} Off-Targets",
                        f"Found {len(off_targets)} off-targets",
                        {"first_offsets": off_targets[:3]
                            if off_targets else []}
                    )

    def _log_mutation(self, original_seq, mutated_seq, position, length):
        """Log details of each mutation."""
        for i in range(position, position + length):
            old_base = original_seq[i] if i < len(original_seq) else "N/A"
            new_base = mutated_seq[i] if i < len(mutated_seq) else "N/A"
            if old_base != new_base:
                msg = f"Mutated base at position {i+1}: {old_base} -> {new_base}"
                if self.debug:
                    self.debugger.log_step("Mutation Detail", msg)
                else:
                    logger.info(msg)

    @debug_function
    def _construct_edge_primer(
            self,
            sequence: str,
            part_num: str,
            forward_or_reverse: str,
            primer_name: str):
        part_end_dict = self.part_end_dict

        if self.debug:
            self.debugger.log_step(
                "Edge Primer Construction",
                f"Constructing {forward_or_reverse} edge primer for part {part_num}",
                {"sequence_length": len(sequence)}
            )

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

    @debug_function
    def generate_GG_edge_primers(self, idx, sequence, mtk_part_left, mtk_part_right, primer_name):
        """
        Generate Golden Gate edge primers for a sequence with optimal primer length calculation.

        Args:
            idx: Index or identifier for the sequence
            sequence: DNA sequence for which primers are to be designed
            overhang_5_prime: 5' overhang sequence
            overhang_3_prime: 3' overhang sequence
            primer_name: Base name for the primers

        Returns:
            Dictionary containing primer information
        """
        if self.debug:
            self.debugger.log_step(
                "Generate Edge Primers",
                f"Generating primers for sequence {idx}",
                {
                    "sequence_length": len(sequence),
                    "left_part": mtk_part_left,
                    "right_part": mtk_part_right
                }
            )

        # Convert to string for consistent handling
        seq_str = str(sequence)
        seq_length = len(seq_str)

        # Calculate optimal primer lengths instead of using fixed length
        forward_length = self._calculate_optimal_primer_length(
            seq_str, 0, 'forward')
        reverse_length = self._calculate_optimal_primer_length(
            seq_str, len(seq_str), 'reverse')

        if self.debug:
            self.debugger.log_step(
                "Calculated Lengths",
                f"Determined optimal primer lengths",
                {
                    "forward_length": forward_length,
                    "reverse_length": reverse_length
                }
            )

        overhang_5_prime = self.utils.get_mtk_partend_sequence(
            mtk_part_left, "forward", kozak=self.kozak)
        overhang_3_prime = self.utils.get_mtk_partend_sequence(
            mtk_part_right, "reverse", kozak=self.kozak)

        if self.debug:
            self.debugger.validate(
                overhang_5_prime is not None and overhang_3_prime is not None,
                f"Successfully retrieved MTK part end overhangs",
                {
                    "5_prime_overhang": overhang_5_prime,
                    "3_prime_overhang": overhang_3_prime
                }
            )

        # 5' Primer (forward)
        forward_primer_binding = seq_str[:forward_length]
        forward_primer = overhang_5_prime + forward_primer_binding

        # 3' Primer (reverse)
        reverse_primer_binding = str(
            Seq(seq_str[-reverse_length:]).reverse_complement())
        reverse_primer = overhang_3_prime + reverse_primer_binding

        if self.debug:
            self.debugger.log_step(
                "Constructed Primers",
                f"Created forward and reverse primers",
                {
                    "forward_primer": forward_primer,
                    "forward_binding": forward_primer_binding,
                    "reverse_primer": reverse_primer,
                    "reverse_binding": reverse_primer_binding
                }
            )

        # Calculate properties
        forward_tm = self._calculate_tm(forward_primer_binding)
        reverse_tm = self._calculate_tm(reverse_primer_binding)

        forward_gc = self._calculate_gc_content(forward_primer_binding)
        reverse_gc = self._calculate_gc_content(reverse_primer_binding)

        if self.debug:
            self.debugger.log_step(
                "Primer Properties",
                f"Calculated primer properties",
                {
                    "forward_tm": forward_tm,
                    "reverse_tm": reverse_tm,
                    "forward_gc": forward_gc,
                    "reverse_gc": reverse_gc
                }
            )

        # Create result dictionary
        primers = {
            "forward_primer": {
                "name": f"{primer_name}_F",
                "sequence": forward_primer,
                "binding_region": forward_primer_binding,
                "tm": forward_tm,
                "gc_content": forward_gc,
                "length": len(forward_primer)
            },
            "reverse_primer": {
                "name": f"{primer_name}_R",
                "sequence": reverse_primer,
                "binding_region": reverse_primer_binding,
                "tm": reverse_tm,
                "gc_content": reverse_gc,
                "length": len(reverse_primer)
            },
            "product_size": seq_length
        }

        if self.debug:
            self.debugger.validate(
                'forward_primer' in primers and 'reverse_primer' in primers,
                f"Successfully created edge primers",
                {"product_size": primers["product_size"]}
            )

        return primers

    @debug_function
    def _calculate_optimal_primer_length(self, sequence, position, direction='forward'):
        """
        Calculate the optimal primer length based on sequence properties.

        Args:
            sequence: The DNA sequence
            position: Start position for calculation (0 for forward, len(seq) for reverse)
            direction: 'forward' or 'reverse'

        Returns:
            Optimal primer length
        """
        if self.debug:
            self.debugger.log_step(
                "Calculate Primer Length",
                f"Determining optimal {direction} primer length from position {position}",
                {"sequence_length": len(sequence)}
            )

        # Set minimum and maximum primer lengths
        min_length = 18
        max_length = 30
        target_tm = 60  # Target melting temperature in °C

        # Initialize with minimum length
        optimal_length = min_length

        if direction == 'forward':
            # For forward primers, start from position and extend right
            for length in range(min_length, min(max_length + 1, len(sequence) - position)):
                primer_seq = sequence[position:position + length]
                tm = self._calculate_tm(primer_seq)

                if self.debug:
                    self.debugger.log_step(
                        "Length Iteration",
                        f"Testing length {length}",
                        {"sequence": primer_seq, "tm": tm, "target": target_tm},
                        level=logging.DEBUG
                    )

                if tm >= target_tm:
                    optimal_length = length
                    break
        else:  # reverse
            # For reverse primers, start from position and extend left
            for length in range(min_length, min(max_length + 1, position + 1)):
                if position - length < 0:
                    break
                primer_seq = sequence[position - length:position]
                tm = self._calculate_tm(primer_seq)

                if self.debug:
                    self.debugger.log_step(
                        "Length Iteration",
                        f"Testing length {length}",
                        {"sequence": primer_seq, "tm": tm, "target": target_tm},
                        level=logging.DEBUG
                    )

                if tm >= target_tm:
                    optimal_length = length
                    break

        if self.debug:
            self.debugger.validate(
                optimal_length >= min_length,
                f"Calculated optimal primer length: {optimal_length}",
                {"direction": direction, "min_length": min_length,
                    "max_length": max_length}
            )

        return optimal_length
