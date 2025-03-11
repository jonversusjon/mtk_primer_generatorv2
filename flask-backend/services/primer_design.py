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
            'tm_threshold': 55.0, 'min_3p_match': 10, 'max_mismatches': 1,
            'mv_conc': 50.0, 'dv_conc': 1.5, 'dntp_conc': 0.2,
            'dna_conc': 250.0, 'min_tm': 57
        }
        self.bsmbi_site = "CGTCTC"
        self.spacer = "GAA"
        self.max_binding_length = 30

        # Validate that the MTK part end sequences have been loaded.
        self.validate(self.part_end_dict is not None,
                      "MTK part end sequences loaded successfully", {"kozak": self.kozak})

    @DebugMixin.debug_wrapper
    def design_mutation_primers(self, mutation_sets: list, comp_matrices: list, primer_name: str = None):
        """
        Designs mutation primers for the provided full sequence given one or more mutation sets.
        Uses the compatibility matrix to pick valid overhang combinations and constructs MutationPrimer objects.
        """

        for set_index, mutation_set in enumerate(mutation_sets):
            self.state['current_mutation'] = set_index

            # Use the compatibility matrix to obtain valid overhang combinations.
            valid_coords = np.argwhere(comp_matrices[set_index] == 1)
            selected_coords = valid_coords[np.random.choice(
                valid_coords.shape[0])]
            selected_coords = selected_coords.tolist()

            # Use log_step to output debug info:
            self.log_step(
                "Debug Info", f"Valid coordinates for mutation set {set_index+1}: {valid_coords.tolist()}")
            self.log_step(
                "Debug Info", f"Selected coordinates for mutation set {set_index+1}: {selected_coords}")

            self.log_step("Compatibility Matrix", f"Valid coordinates for mutation set {set_index+1}",
                          {"valid_coords": valid_coords.tolist()})
            self.validate(valid_coords.size > 0,
                          f"Found {np.count_nonzero(comp_matrices[set_index])} valid overhang combinations",
                          {"matrix_size": comp_matrices[set_index].size})
            self.log_step("Matrix Visualization", f"Compatibility matrix for set {set_index+1}",
                          visualize_matrix(comp_matrices[set_index]))

            # we now know, for the current mutation set, which overhangs are compatible with each other
            # amongst the N sites to mutate

            for i, site_to_mutate_data in enumerate(mutation_set):
                overhang_options = site_to_mutate_data["overhangs"].get(
                    "overhang_options", [])

                # Ensure valid_coords has an entry for this index
                if i >= len(valid_coords):
                    raise IndexError(f"valid_coords index {i} out of range")

                selected_overhang = selected_coords[i]
                overhang_data = overhang_options[selected_overhang]
                print(f"overhang_data: {overhang_data}")
                # Extract overhang sequences safely
                try:
                    top_overhang_seq = overhang_data["top_overhang"]["seq"]
                    bottom_overhang_seq = overhang_data["bottom_overhang"]["seq"]
                    self.log_step(
                        "Overhang Sequences", f"Top overhang: {top_overhang_seq}, Bottom overhang: {bottom_overhang_seq}",)
                except KeyError as e:
                    raise KeyError(
                        f"Missing expected key in overhang_data: {e}")

                # Use extracted sequences as needed

                print(f"site_to_mutate_data: {site_to_mutate_data}")

                # design primer for the ith mutation
                # the mutated context for designing primers is contained in the alternative codon dict

            # Construct MutationPrimer objects.
            mutation_primers = self._construct_mutation_primers(
                mutation_set=mutation_set,
                selected_coords=selected_coords,
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
    def _construct_mutation_primers(self, mutation_set: list,
                                    selected_coords: list,
                                    primer_name: str = None,
                                    min_binding_length: int = 10) -> list:
        """
        For each mutation in mutation_set, design forward and reverse primers such that:
        - The forward primer is on the left (5') side of the overhang/mutation site
            (on the sense strand).
        - The reverse primer is on the right (3') side of the overhang/mutation site
            (on the sense strand), but we reverse-complement that region to get the actual
            primer that anneals to the antisense strand.

        'overhang_start_index' is the 0-based index where the 4-nt sticky end begins (e.g. 32).
        """

        self.log_step("Construct Primers",
                      f"Constructing primers with binding length {min_binding_length}",
                      {"selected_coord": selected_coords})
        mutation_primers = []

        def reverse_complement(seq: str) -> str:
            complement = str.maketrans('ATCGatcg', 'TAGCtagc')
            return seq.translate(complement)[::-1]

        def calculate_gc(seq: str) -> float:
            seq = seq.upper()
            if len(seq) == 0:
                return 0.0
            gc_count = seq.count("G") + seq.count("C")
            return round((gc_count / len(seq)) * 100, 2)

        for i, mutation_set_obj in enumerate(mutation_set):
            selected_overhang = selected_coords[i]
            overhang_data = mutation_set_obj["overhangs"]["overhang_options"][selected_overhang]
            print(f"overhang_data: {overhang_data}")
            mutated_context = mutation_set_obj["mutated_context"]
            overhang_start_index = overhang_data["overhang_start_index"]

            tm_threshold = self.default_params["tm_threshold"]

            # ========== FORWARD PRIMER ==========
            f_5prime_cap = overhang_start_index - 1

            # Expand the annealing region if Tm is too low.
            f_seq_length = min_binding_length
            f_anneal = mutated_context[f_5prime_cap:f_5prime_cap + f_seq_length]
            while self._calculate_tm(f_anneal) < tm_threshold and f_seq_length < self.max_binding_length:
                f_seq_length += 1
                f_anneal = mutated_context[f_5prime_cap:f_5prime_cap + f_seq_length]

            self.log_step(
                "Forward Primer Sequence",
                f"Forward primer annealing region[1:5]: {f_anneal[1:5]}, "
                f"overhang_seq: {overhang_data['top_overhang']['seq']}"
            )

            self.validate(
                str(f_anneal[1:5]).strip().upper() == str(
                    overhang_data["top_overhang"]["seq"]).strip().upper(),
                f"Forward primer annealing region[1:5] ({f_anneal[1:5].strip().upper()}) "
                f"matches expected overhang sequence ({overhang_data['top_overhang']['seq'].strip().upper()})"
            )

            f_primer_sequence = self.spacer + self.bsmbi_site + f_anneal

            # ========== REVERSE PRIMER ==========
            r_5prime_cap = overhang_start_index + 5

            # Expand the annealing region if Tm is too low.
            r_seq_length = min_binding_length
            r_anneal = mutated_context[r_5prime_cap:r_5prime_cap + r_seq_length]
            while self._calculate_tm(r_anneal) < tm_threshold and r_seq_length < self.max_binding_length:
                r_seq_length += 1
                r_anneal = mutated_context[r_5prime_cap -
                                           r_seq_length:r_5prime_cap]

            r_anneal = reverse_complement(r_anneal)
            r_primer_sequence = self.spacer + self.bsmbi_site + r_anneal
            self.log_step(
                "Reverse Primer Sequence",
                f"Reverse primer annealing region[1:5]: {r_anneal[1:5]}, "
                f"overhang_seq: {overhang_data['bottom_overhang']['seq']}"
            )
            self.log_step(
                "r_anneal", f"Reverse primer annealing region: {r_anneal}")
            self.validate(
                str(r_anneal[1:5]).strip().upper() == str(
                    overhang_data["bottom_overhang"]["seq"]).strip().upper(),
                f"Reverse primer annealing region[1:5] ({repr(str(r_anneal[1:5]).strip().upper())}) "
                f"matches expected overhang sequence ({repr(str(overhang_data['bottom_overhang']['seq']).strip().upper())})"
            )

            # ========== CREATE PRIMER OBJECTS ==========
            f_primer = Primer(
                name=primer_name +
                "_forward" if primer_name else f"primer_{i}_forward",
                sequence=f_primer_sequence,
                binding_region=f_anneal,
                tm=self._calculate_tm(f_anneal),
                gc_content=calculate_gc(f_anneal),
                length=len(f_primer_sequence)
            )

            r_primer = Primer(
                name=primer_name +
                "_reverse" if primer_name else f"primer_{i}_reverse",
                sequence=r_primer_sequence,
                binding_region=r_anneal,
                tm=self._calculate_tm(r_anneal),
                gc_content=calculate_gc(r_anneal),
                length=len(r_primer_sequence)
            )

            mutation_primer = MutationPrimer(
                site=mutation_set_obj["site"],
                position=mutation_set_obj["position"],
                forward=f_primer,
                reverse=r_primer,
                mutation_info=mutation_set_obj
            )

            self.log_step("Primer Design Result",
                          f"Designed primer pair for mutation {mutation_set_obj['site']} at position {mutation_set_obj['position']}")
            self.log_step("Mutated Context", mutated_context)
            self.log_step("Overhang", overhang_data)

            self.log_step("Forward Primer Name", f_primer.name)
            self.log_step("Forward Primer Sequence", f_primer.sequence)
            self.log_step("Forward Primer Binding Region",
                          f_primer.binding_region)
            self.log_step("Forward Primer Tm", f_primer.tm)
            self.log_step("Forward Primer GC Content", f_primer.gc_content)
            self.log_step("Forward Primer Length", f_primer.length)

            self.log_step("Reverse Primer Name", r_primer.name)
            self.log_step("Reverse Primer Sequence", r_primer.sequence)
            self.log_step("Reverse Primer Binding Region",
                          r_primer.binding_region)
            self.log_step("Reverse Primer Tm", r_primer.tm)
            self.log_step("Reverse Primer GC Content", r_primer.gc_content)
            self.log_step("Reverse Primer Length", r_primer.length)

            mutation_primers.append(mutation_primer)

        return mutation_primers

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
