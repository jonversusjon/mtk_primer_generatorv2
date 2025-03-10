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
        Uses the compatibility matrix to pick valid overhang combinations and constructs MutationPrimer objects.
        """
        # (Input validation code omitted for brevity)

        for set_index, mutation_set in enumerate(mutation_sets):
            self.state['current_mutation'] = set_index
            self.log_step("Process Mutation Set", f"Processing mutation set {set_index+1}/{len(mutation_sets)}",
                          {"sites": list(mutation_set.keys())})

            # Use the compatibility matrix to obtain valid overhang combinations.
            valid_coords = np.argwhere(comp_matrices[set_index] == 1)
            self.log_step("Compatibility Matrix", f"Valid coordinates for mutation set {set_index+1}",
                          {"valid_coords": valid_coords.tolist()})
            self.validate(valid_coords.size > 0,
                          f"Found {np.count_nonzero(comp_matrices[set_index])} valid overhang combinations",
                          {"matrix_size": comp_matrices[set_index].size})
            self.log_step("Matrix Visualization", f"Compatibility matrix for set {set_index+1}",
                          visualize_matrix(comp_matrices[set_index]))

            # we now know, for the current mutation set, which overhangs are compatible with each other
            # amongst the N sites to mutate
            for i, site_to_mutate_data in enumerate(mutation_set.values()):
                print(f"site_to_mutate_data: {site_to_mutate_data}")
                # design primer for the ith mutation
                # the mutated context for designing primers is contained in the alternative codon dict
                site, position, original_codon_sequence, alternative_codon_sequence, mutated_base_index, overhangs, primer_context = site_to_mutate_data.values()

                mutated_seq = self._apply_mutations(
                    full_sequence, mutation_set)
                self.validate(mutated_seq is not None and len(mutated_seq) == len(full_sequence),
                              "Successfully applied mutations to sequence",
                              {"original_length": len(full_sequence), "mutated_length": len(mutated_seq)})
                # annealing_region = mutated_seq[sticky_end_placement[0]:sticky_end_placement[0] + self.binding_length]
                # the forward primer is self.spacer + self.bsmbi_site + annealing_region

            # # Log chosen sticky end details (drilling down into the codon structure).
            # for i, site_key in enumerate(position_keys):
            #     site_info = mutation_set[site_key]
            #     overhang_index = int(chosen_coords[i].item())
            #     try:
            #         if "overhangs" in site_info:
            #             overhang_options = site_info["overhangs"].get(
            #                 "overhang_options", [])
            #             if overhang_options and overhang_index < len(overhang_options):
            #                 chosen_option = overhang_options[overhang_index]
            #                 top_info = chosen_option.get("top_overhang")
            #                 bottom_info = chosen_option.get("bottom_overhang")
            #                 self.log_step(
            #                     "Chosen Overhang",
            #                     f"For site {site_key} (position {site_info['position']}), using overhang option index {overhang_index}",
            #                     {
            #                         "top_overhang": top_info,
            #                         "bottom_overhang": bottom_info,
            #                     }
            #                 )
            #             else:
            #                 raise KeyError(
            #                     "No valid overhang options available.")
            #         else:
            #             # Fallback logging for legacy structure.
            #             if "alternative_codons" in site_info:
            #                 alternative = site_info["alternative_codons"][overhang_index]
            #             elif "codons" in site_info:
            #                 alternative = site_info["codons"][0]["alternative_codons"][overhang_index]
            #             else:
            #                 raise KeyError(
            #                     "No alternative codon information found.")
            #             sticky_ends = alternative.get("sticky_ends", {})
            #             sticky_key = list(sticky_ends.keys())[
            #                 0] if sticky_ends else None
            #             top_option = sticky_ends.get(sticky_key, {}).get(
            #                 "top_strand", []) if sticky_key else None
            #             bottom_option = sticky_ends.get(sticky_key, {}).get(
            #                 "bottom_strand", []) if sticky_key else None
            #             self.log_step(
            #                 "Chosen Sticky Ends",
            #                 f"For site {site_key} (position {site_info['position']}), using alternative codon index {overhang_index}",
            #                 {
            #                     "alternative_codon": alternative["seq"],
            #                     "top_sticky": top_option[overhang_index] if top_option and len(top_option) > overhang_index else None,
            #                     "bottom_sticky": bottom_option[overhang_index] if bottom_option and len(bottom_option) > overhang_index else None,
            #                 }
            #             )
                # except (KeyError, IndexError) as e:
                #     self.log_step("Sticky End Error", f"Could not retrieve sticky end for site {site_key}",
                #                   {"error": str(e)}, level=logging.ERROR)
                #     continue

            # Construct MutationPrimer objects.
            mutation_primers = self._construct_mutation_primers(
                original_seq=full_sequence,
                mutated_seq=mutated_seq,
                mutation_set=mutation_set,
                valid_coords=[],
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
        self.log_step("Apply Mutations", f"Applying {len(mutation_set)} mutations to sequence",
                      {"sequence_length": len(target_seq)})
        mutated_seq = list(target_seq)
        for site_key, mut in mutation_set.items():
            pos = mut["position"] - 1  # Convert to 0-indexed.
            alt_codon_seq = mut["alternative_codon_sequence"]
            orig_codon_seq = mut["original_codon_sequence"]
            self.log_step("Mutation Details", f"Applying mutation at site {site_key}, position {pos+1}",
                          {"original": orig_codon_seq, "alternative": alt_codon_seq,
                           "context_sequence": target_seq[max(0, pos-5):min(len(target_seq), pos+len(orig_codon_seq)+5)]})
            mutated_seq[pos:pos + len(orig_codon_seq)] = alt_codon_seq
            if self.verbose or self.debug:
                self._log_mutation(target_seq, ''.join(
                    mutated_seq), pos, len(orig_codon_seq))
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
                                    valid_coords: np.ndarray, primer_name: str = None, binding_length: int = 10) -> list:
        self.log_step("Construct Primers", f"Constructing primers with binding length {binding_length}",
                      {"selected_coord": valid_coords})
        mutation_primers = []

        for i, (site_key, mut) in enumerate(mutation_set.items()):
            selected_coord = int(valid_coords[i].item())
            position = mut["position"] - 1
            self.log_step("Design Primer Pair", f"Designing primers for site {site_key} at position {position+1}",
                          {"mutation": mut["alternative_codon_sequence"], "selected_coord": selected_coord})

            forward_tuple = self._construct_primer(mutated_seq, position, mut, selected_coord,
                                                   binding_length, is_reverse=False, primer_name=primer_name)
            reverse_tuple = self._construct_primer(mutated_seq, position, mut, selected_coord,
                                                   binding_length, is_reverse=True, primer_name=primer_name)

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
                          binding_length: int, is_reverse: bool, primer_name: str = None):
        """
        Constructs a mutagenic primer using the mutated sequence as the binding region.
        The final primer consists of:
        - A non-annealing tail: spacer + BsmBI recognition site
        - An annealing region: a segment of the mutated sequence starting at the sticky end start index.
        For reverse primers, the annealing region is reverse-complemented.
        """
        direction = "reverse" if is_reverse else "forward"
        self.log_step(f"Construct {direction.capitalize()} Primer",
                      f"Designing {direction} primer at site position {position+1}",
                      {"binding_length": binding_length, "selected_coord": selected_coord})

        if not primer_name:
            primer_name = f"Mut_{mut.get('site', 'unknown')}"
        primer_suffix = "_RV" if is_reverse else "_FW"
        primer_name = f"{primer_name}{primer_suffix}"

        # Use the new structure if available: expect an "overhangs" key with "overhang_options"
        if "overhangs" in mut:
            overhang_options = mut["overhangs"].get("overhang_options", [])
            if not overhang_options or selected_coord >= len(overhang_options):
                self.debugger and self.debugger.log_warning(
                    f"No valid overhang option found for site {mut.get('site', 'unknown')}")
                return None
            overhang_option = overhang_options[selected_coord]
            # Choose the strand info based on primer orientation.
            strand_info = overhang_option.get(
                "bottom_overhang") if is_reverse else overhang_option.get("top_overhang")
            if not strand_info:
                self.debugger and self.debugger.log_warning(
                    f"No sticky end information found for site {mut.get('site', 'unknown')}")
                return None
            overhang_start_index = strand_info.get("start_index")
            if overhang_start_index is None:
                self.debugger and self.debugger.log_warning(
                    f"Sticky end start index missing for site {mut.get('site', 'unknown')}")
                return None
        else:
            # Fallback to legacy structure if needed (e.g., using alternative_codons)
            if "alternative_codons" in mut:
                alternative_codons = mut["alternative_codons"]
            elif "codons" in mut and len(mut["codons"]) > 0:
                alternative_codons = mut["codons"][0].get(
                    "alternative_codons", [])
            else:
                alternative_codons = []
            if not alternative_codons or selected_coord >= len(alternative_codons):
                self.debugger and self.debugger.log_warning(
                    f"No valid alternative codon found for site {mut.get('site', 'unknown')}")
                return None
            alternative = alternative_codons[selected_coord]
            sticky_ends = alternative.get("sticky_ends", {})
            # Assume a key (like "position_2") exists; adjust if necessary.
            sticky_key = list(sticky_ends.keys())[0] if sticky_ends else None
            if not sticky_key:
                self.debugger and self.debugger.log_warning(
                    f"No sticky ends found for alternative codon {alternative['seq']} at site {mut.get('site', 'unknown')}")
                return None
            if is_reverse:
                sticky_list = sticky_ends.get(
                    sticky_key, {}).get("bottom_strand", [])
            else:
                sticky_list = sticky_ends.get(
                    sticky_key, {}).get("top_strand", [])
            if not sticky_list or selected_coord >= len(sticky_list):
                self.debugger and self.debugger.log_warning(
                    f"No valid sticky end option found for site {mut.get('site', 'unknown')}")
                return None
            overhang_option = sticky_list[selected_coord]
            overhang_start_index = overhang_option.get("start_index")
            if overhang_start_index is None:
                self.debugger and self.debugger.log_warning(
                    f"Sticky end start index missing for site {mut.get('site', 'unknown')}")
                return None

        # Calculate the binding (annealing) region from the mutated sequence.
        binding_region = mutated_seq[overhang_start_index:
                                     overhang_start_index + binding_length]
        if len(binding_region) < binding_length:
            self.debugger and self.debugger.log_warning(
                f"Binding region length insufficient for {direction} primer",
                {"start_index": overhang_start_index, "requested_length": binding_length,
                 "actual_length": len(binding_region)}
            )
            return None

        if is_reverse:
            binding_region = str(Seq(binding_region).reverse_complement())

        # The final primer: non-annealing tail (spacer + BsmBI site) + binding region.
        primer_seq = self.spacer + self.bsmbi_site + binding_region

        self.validate(len(primer_seq) > len(self.spacer + self.bsmbi_site),
                      "Primer contains valid binding region",
                      {"spacer": self.spacer,
                       "enzyme_site": self.bsmbi_site,
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
