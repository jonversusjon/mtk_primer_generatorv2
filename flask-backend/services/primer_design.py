import numpy as np
from Bio.Seq import Seq
from .utils import GoldenGateUtils
from config.logging_config import logger
from services.debug.debug_utils import MutationDebugger, visualize_matrix
from models.primer import Primer, MutationPrimer
import logging
from services.debug.debug_mixin import DebugMixin


class PrimerDesigner(DebugMixin):
    """
    Handles primer design for Golden Gate assembly.
    """

    def __init__(self, kozak: str = "MTK", verbose: bool = False, debug: bool = False):
        self.logger = logger.getChild("PrimerDesigner")
        self.utils = GoldenGateUtils()
        self.verbose = verbose
        self.debug = debug
        self.debugger = None

        if self.debug:
            self.debugger = MutationDebugger(
                parent_logger=logger, use_custom_format=True)
            if hasattr(self.debugger.logger, 'propagate'):
                self.debugger.logger.propagate = False
            self.logger.info("Debug mode enabled for PrimerDesigner")

        self.state = {'current_operation': '',
                      'primers_designed': 0, 'current_mutation': None}
        self.kozak = kozak
        self.part_end_dict = self.utils.get_mtk_partend_sequences()
        self.default_params = {
            'tm_threshold': 45.0, 'min_3p_match': 10, 'max_mismatches': 1,
            'mv_conc': 50.0, 'dv_conc': 1.5, 'dntp_conc': 0.2,
            'dna_conc': 250.0, 'min_tm': 57
        }
        self.bsmbi_site = "CGTCTC"
        self.spacer = "GAA"

        # Validate that the MTK part end sequences have been loaded.
        self.validate(self.part_end_dict is not None,
                      "MTK part end sequences loaded successfully", {"kozak": self.kozak})

    @DebugMixin.debug_wrapper
    def design_mutation_primers(self, full_sequence: str, mutation_sets: list, comp_matrices: list, primer_name: str = None):
        """
        Designs mutation primers for the provided full sequence given one or more mutation sets.
        Uses the compatibility matrix to pick valid overhang combinations.
        """
        # Validate inputs.
        self.validate(isinstance(full_sequence, str) and full_sequence,
                      "Input sequence is valid", {"sequence_length": len(full_sequence)})
        self.validate(isinstance(mutation_sets, list) and mutation_sets,
                      f"Received {len(mutation_sets)} mutation sets",
                      {"first_set_sites": list(mutation_sets[0].keys()) if mutation_sets else None})
        self.validate(isinstance(comp_matrices, list) and len(comp_matrices) == len(mutation_sets),
                      f"Received {len(comp_matrices)} compatibility matrices",
                      {"first_matrix_shape": comp_matrices[0].shape if comp_matrices else None})

        for set_index, mutation_set in enumerate(mutation_sets):
            self.state['current_mutation'] = set_index
            self.log_step("Process Mutation Set", f"Processing mutation set {set_index+1}/{len(mutation_sets)}",
                          {"sites": list(mutation_set.keys())})

            # Find valid overhang combinations using the compatibility matrix.
            valid_coords = np.argwhere(comp_matrices[set_index] == 1)
            valid_combinations = np.count_nonzero(comp_matrices[set_index])
            self.validate(valid_coords.size > 0,
                          f"Found {valid_combinations} valid overhang combinations",
                          {"matrix_size": comp_matrices[set_index].size})

            if valid_coords.size > 0:
                self.log_step("Matrix Visualization", f"Compatibility matrix for set {set_index+1}",
                              visualize_matrix(comp_matrices[set_index]))
            else:
                self.debugger and self.debugger.log_warning(
                    f"No valid overhang combinations for set {set_index+1}")
                continue

            # Apply mutations to create the mutated sequence.
            mutated_seq = self._apply_mutations(full_sequence, mutation_set)
            self.validate(mutated_seq is not None and len(mutated_seq) == len(full_sequence),
                          "Successfully applied mutations to sequence",
                          {"original_length": len(full_sequence), "mutated_length": len(mutated_seq)})

            # Choose the first valid coordinate row from the compatibility matrix.
            chosen_coords = valid_coords[0]
            position_keys = list(mutation_set.keys())

            for i, site_key in enumerate(position_keys):
                site_info = mutation_set[site_key]
                overhang_index = int(chosen_coords[i].item())
                try:
                    overhang_data = site_info["overhangs"]["overhang_options"][overhang_index]
                except (KeyError, IndexError) as e:
                    self.log_step("Overhang Error", f"Could not retrieve overhang for site {site_key}",
                                  {"error": str(e)}, level=logging.ERROR)
                    continue

                self.log_step(
                    "Chosen Overhang",
                    f"For site {site_key} (position {site_info['position']}), using overhang index {overhang_index}",
                    {
                        "top_overhang": overhang_data.get("top_overhang"),
                        "bottom_overhang": overhang_data.get("bottom_overhang"),
                        "top_extended": overhang_data.get("top_extended"),
                        "bottom_extended": overhang_data.get("bottom_extended"),
                    }
                )

            # Construct MutationPrimer objects using the chosen overhang row.
            mutation_primers = self._construct_mutation_primers(
                original_seq=full_sequence,
                mutated_seq=mutated_seq,
                mutation_set=mutation_set,
                valid_coords=chosen_coords,
                primer_name=primer_name
            )
            self.validate(mutation_primers is not None, "Successfully constructed mutation primers",
                          {"primer_count": len(mutation_primers) if mutation_primers else 0})

            if mutation_primers:
                return mutation_primers

        self.debugger and self.debugger.log_warning(
            "Failed to design primers for any mutation set")
        return None

    @DebugMixin.debug_wrapper
    def _apply_mutations(self, target_seq: str, mutation_set: dict) -> str:
        """
        Applies the mutations specified in mutation_set to target_seq.
        """
        self.log_step("Apply Mutations", f"Applying {len(mutation_set)} mutations to sequence", {
            "sequence_length": len(target_seq)})
        mutated_seq = list(target_seq)
        for site_key, mut in mutation_set.items():
            pos = mut["position"] - 1  # converting to 0-indexed.
            alt_seq = mut["alternative_sequence"]
            orig_seq = mut["original_sequence"]
            self.log_step("Mutation Details", f"Applying mutation at site {site_key}, position {pos+1}",
                          {"original": orig_seq, "alternative": alt_seq,
                           "context": target_seq[max(0, pos-5):min(len(target_seq), pos+len(orig_seq)+5)]})
            mutated_seq[pos:pos + len(orig_seq)] = alt_seq
            if self.verbose or self.debug:
                self._log_mutation(target_seq, ''.join(
                    mutated_seq), pos, len(orig_seq))
        result = ''.join(mutated_seq)
        self.validate(len(result) == len(target_seq),
                      "Mutation applied successfully, sequence length preserved",
                      {"original_length": len(target_seq), "mutated_length": len(result)})
        return result

    def _log_mutation(self, original_seq, mutated_seq, position, length):
        for i in range(position, position + length):
            old_base = original_seq[i] if i < len(original_seq) else "N/A"
            new_base = mutated_seq[i] if i < len(mutated_seq) else "N/A"
            if old_base != new_base:
                msg = f"Mutated base at position {i+1}: {old_base} -> {new_base}"
                if self.debugger:
                    self.debugger.log_step("Mutation Detail", msg)
                else:
                    self.logger.info(msg)

    @DebugMixin.debug_wrapper
    def _construct_mutation_primers(self, original_seq: str, mutated_seq: str, mutation_set: dict,
                                    valid_coords: np.ndarray, primer_name: str = None, flank: int = 10) -> list:
        self.log_step("Construct Primers", f"Constructing primers with flank length {flank}",
                      {"selected_coord": valid_coords})
        mutation_primers = []

        for i, (site_key, mut) in enumerate(mutation_set.items()):
            selected_coord = int(valid_coords[i].item())
            position = mut["position"] - 1
            self.log_step("Design Primer Pair", f"Designing primers for site {site_key} at position {position+1}",
                          {"mutation": mut["alternative_sequence"], "selected_coord": selected_coord})

            forward_tuple = self._construct_primer(mutated_seq, position, mut, selected_coord,
                                                   flank, is_reverse=False, primer_name=primer_name)
            reverse_tuple = self._construct_primer(mutated_seq, position, mut, selected_coord,
                                                   flank, is_reverse=True, primer_name=primer_name)

            if forward_tuple and reverse_tuple:
                forward_primer = Primer(
                    name=forward_tuple[0], sequence=forward_tuple[1])
                reverse_primer = Primer(
                    name=reverse_tuple[0], sequence=reverse_tuple[1])
                mutation_primer = MutationPrimer(site=mut["site"], position=mut["position"],
                                                 forward=forward_primer, reverse=reverse_primer, mutation_info=mut)
                mutation_primers.append(mutation_primer)
            else:
                self.debugger and self.debugger.log_warning(
                    f"Failed to create primer pair for site {site_key}")

        self.validate(len(mutation_primers) > 0, f"Successfully created {len(mutation_primers)} primer pairs",
                      {"expected": len(mutation_set)})
        return mutation_primers

    @DebugMixin.debug_wrapper
    def _construct_primer(self, mutated_seq: str, position: int, mut: dict, selected_coord: int,
                          flank: int, is_reverse: bool, primer_name: str = None):
        """
        Constructs a mutagenic primer using the mutated sequence as the binding region.
        The final primer consists of:
          - A non-annealing tail: spacer + BsmBI recognition site + an extended overhang (from the compatibility matrix)
          - An annealing region: a segment of the mutated sequence that spans from (position - flank) to (position + mutation_length + flank)
        For reverse primers, the annealing region is reverse-complemented.
        """
        direction = "reverse" if is_reverse else "forward"
        self.log_step(f"Construct {direction.capitalize()} Primer", f"Designing {direction} primer at position {position+1}",
                      {"flank": flank, "selected_coord": selected_coord})

        if not primer_name:
            primer_name = f"Mut_{mut['site']}"
        primer_suffix = "_RV" if is_reverse else "_FW"
        primer_name = f"{primer_name}{primer_suffix}"

        mutation_length = len(mut["original_sequence"])
        binding_start = max(0, position - flank)
        binding_end = min(len(mutated_seq), position + mutation_length + flank)

        if binding_end - binding_start < (mutation_length + 2 * flank):
            self.debugger and self.debugger.log_warning(f"Binding region length insufficient for {direction} primer",
                                                        {"binding_start": binding_start, "binding_end": binding_end,
                                                         "expected_length": mutation_length + 2 * flank})
            return None

        binding_region = mutated_seq[binding_start:binding_end]
        if is_reverse:
            binding_region = str(Seq(binding_region).reverse_complement())

        extended_seq = str(mut['overhangs']['overhang_options'][selected_coord][
            'bottom_extended' if is_reverse else 'top_extended'])

        primer_seq = self.spacer + self.bsmbi_site + extended_seq + binding_region
        self.validate(len(primer_seq) > len(self.spacer + self.bsmbi_site + extended_seq),
                      "Primer contains valid binding region",
                      {"spacer": self.spacer,
                       "enzyme_site": self.bsmbi_site,
                       "extended_seq": extended_seq,
                       "binding_region": binding_region,
                       "total_length": len(primer_seq),
                       "binding_length": len(binding_region)})
        return (primer_name, primer_seq)

    @DebugMixin.debug_wrapper
    def generate_GG_edge_primers(self, idx, sequence, mtk_part_left, mtk_part_right, primer_name):
        self.log_step("Generate Edge Primers", f"Generating primers for sequence {idx}",
                      {"sequence_length": len(sequence), "left_part": mtk_part_left, "right_part": mtk_part_right})

        seq_str = str(sequence)
        seq_length = len(seq_str)
        forward_length = self._calculate_optimal_primer_length(
            seq_str, 0, 'forward')
        reverse_length = self._calculate_optimal_primer_length(
            seq_str, len(seq_str), 'reverse')

        self.log_step("Calculated Lengths", "Determined optimal primer lengths",
                      {"forward_length": forward_length, "reverse_length": reverse_length})

        overhang_5_prime = self.utils.get_mtk_partend_sequence(
            mtk_part_left, "forward", kozak=self.kozak)
        overhang_3_prime = self.utils.get_mtk_partend_sequence(
            mtk_part_right, "reverse", kozak=self.kozak)
        self.validate(overhang_5_prime is not None and overhang_3_prime is not None,
                      "Successfully retrieved MTK part end overhangs",
                      {"5_prime_overhang": overhang_5_prime, "3_prime_overhang": overhang_3_prime})

        forward_primer_binding = seq_str[:forward_length]
        forward_primer = overhang_5_prime + forward_primer_binding
        reverse_primer_binding = str(
            Seq(seq_str[-reverse_length:]).reverse_complement())
        reverse_primer = overhang_3_prime + reverse_primer_binding

        primers = {
            "forward_primer": {
                "name": f"{primer_name}_F",
                "sequence": forward_primer,
                "binding_region": forward_primer_binding,
                "tm": self._calculate_tm(forward_primer_binding),
                "gc_content": self._calculate_gc_content(forward_primer_binding),
                "length": len(forward_primer)
            },
            "reverse_primer": {
                "name": f"{primer_name}_R",
                "sequence": reverse_primer,
                "binding_region": reverse_primer_binding,
                "tm": self._calculate_tm(reverse_primer_binding),
                "gc_content": self._calculate_gc_content(reverse_primer_binding),
                "length": len(reverse_primer)
            },
            "product_size": seq_length
        }

        self.validate('forward_primer' in primers and 'reverse_primer' in primers,
                      "Successfully created edge primers",
                      {"product_size": primers["product_size"]})
        return primers

    def _calculate_optimal_primer_length(self, sequence, position, direction='forward'):
        self.log_step("Calculate Primer Length", f"Determining optimal {direction} primer length from position {position}",
                      {"sequence_length": len(sequence)})

        min_length = 18
        max_length = 30
        target_tm = 60
        optimal_length = min_length

        if direction == 'forward':
            for length in range(min_length, min(max_length + 1, len(sequence) - position)):
                primer_seq = sequence[position:position + length]
                tm = self._calculate_tm(primer_seq)
                self.log_step("Length Iteration", f"Testing length {length}",
                              {"sequence": primer_seq, "tm": tm, "target": target_tm})
                if tm >= target_tm:
                    optimal_length = length
                    break
        else:
            for length in range(min_length, min(max_length + 1, position + 1)):
                if position - length < 0:
                    break
                primer_seq = sequence[position - length:position]
                tm = self._calculate_tm(primer_seq)
                self.log_step("Length Iteration", f"Testing length {length}",
                              {"sequence": primer_seq, "tm": tm, "target": target_tm})
                if tm >= target_tm:
                    optimal_length = length
                    break

        self.validate(optimal_length >= min_length, f"Calculated optimal primer length: {optimal_length}",
                      {"direction": direction, "min_length": min_length, "max_length": max_length})
        return optimal_length

    def _calculate_gc_content(self, sequence: str) -> float:
        if not sequence:
            return 0.0
        sequence = sequence.upper()
        gc_count = sequence.count("G") + sequence.count("C")
        return (gc_count / len(sequence)) * 100.0

    def _calculate_tm(self, sequence: str) -> float:
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
