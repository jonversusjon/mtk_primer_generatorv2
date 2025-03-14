import numpy as np
from Bio.Seq import Seq
from .utils import GoldenGateUtils
from config.logging_config import logger
from services.debug.debug_utils import MutationDebugger, visualize_matrix
from models.primer import Primer, MutationPrimer
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

        self.state = {
            'current_operation': '',
            'primers_designed': 0,
            'current_mutation': None
        }
        self.kozak = kozak
        self.part_end_dict = self.utils.get_mtk_partend_sequences()
        self.default_params = {
            'tm_threshold': 55.0,
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
        self.max_binding_length = 30

        self.validate(
            self.part_end_dict is not None,
            "MTK part end sequences loaded successfully",
            {"kozak": self.kozak}
        )

    @DebugMixin.debug_wrapper
    def design_mutation_primers(self, mutation_sets: list, comp_matrices: list, primer_name: str = None):
        """
        Designs mutation primers for the provided mutation sets using compatibility matrices.
        """
        for set_index, mutation_set in enumerate(mutation_sets):
            self.state['current_mutation'] = set_index
            comp_matrix = comp_matrices[set_index]
            valid_coords = np.argwhere(comp_matrix == 1)
            self.validate(
                valid_coords.size > 0,
                f"Found {np.count_nonzero(comp_matrix)} valid overhang combinations",
                {"matrix_size": comp_matrix.size}
            )

            # Randomly select one valid coordinate set
            selected_coords = valid_coords[np.random.choice(
                valid_coords.shape[0])].tolist()

            self.log_step(
                "Debug Info", f"Mutation set {set_index+1}: Valid coordinates: {valid_coords.tolist()}")
            self.log_step(
                "Debug Info", f"Mutation set {set_index+1}: Selected coordinates: {selected_coords}")
            self.log_step("Matrix Visualization", f"Compatibility matrix for set {set_index+1}",
                          visualize_matrix(comp_matrix))

            # Process each mutation site in the set
            for i, site_data in enumerate(mutation_set):
                overhang_options = site_data["overhangs"].get(
                    "overhang_options", [])
                if i >= len(valid_coords):
                    raise IndexError(f"valid_coords index {i} out of range")
                selected_overhang = selected_coords[i]
                overhang_data = overhang_options[selected_overhang]

                # Log overhang sequences
                try:
                    top_seq = overhang_data["top_overhang"]["seq"]
                    bottom_seq = overhang_data["bottom_overhang"]["seq"]
                    self.log_step("Overhang Sequences",
                                  f"Top: {top_seq}, Bottom: {bottom_seq}")
                except KeyError as e:
                    raise KeyError(
                        f"Missing expected key in overhang_data: {e}")

                self.log_step(
                    "Site Data", f"Processing mutation site: {site_data}")

            # Construct MutationPrimer objects for the current mutation set
            mutation_primers = self._construct_mutation_primers(
                mutation_set=mutation_set,
                selected_coords=selected_coords,
                primer_name=primer_name
            )
            self.validate(
                mutation_primers is not None,
                "Successfully constructed mutation primers",
                {"primer_count": len(mutation_primers)
                 if mutation_primers else 0}
            )
            if mutation_primers:
                return mutation_primers

        if self.debugger:
            self.debugger.log_warning(
                "Failed to design primers for any mutation set")
        return None

    @DebugMixin.debug_wrapper
    def _construct_mutation_primers(self, mutation_set: list, selected_coords: list,
                                    primer_name: str = None, min_binding_length: int = 10) -> list:
        """
        Constructs forward and reverse mutation primers for each mutation in mutation_set.
        """
        self.log_step("Construct Primers", f"Binding length {min_binding_length}", {
                      "selected_coords": selected_coords})
        mutation_primers = []

        for i, mut_obj in enumerate(mutation_set):
            selected_overhang = selected_coords[i]
            overhang_data = mut_obj["overhangs"]["overhang_options"][selected_overhang]
            mutated_context = mut_obj["mutated_context"]
            overhang_start = overhang_data["overhang_start_index"]

            tm_threshold = self.default_params["tm_threshold"]

            # Design forward primer
            f_5prime = overhang_start - 1
            f_seq_length = min_binding_length
            f_anneal = mutated_context[f_5prime:f_5prime + f_seq_length]
            while self.utils.calculate_tm(f_anneal) < tm_threshold and f_seq_length < self.max_binding_length:
                f_seq_length += 1
                f_anneal = mutated_context[f_5prime:f_5prime + f_seq_length]

            self.log_step("Forward Primer Sequence",
                          f"Annealing region: {f_anneal[1:5]} vs expected: {overhang_data['top_overhang']['seq']}")
            self.validate(
                f_anneal[1:5].strip().upper(
                ) == overhang_data["top_overhang"]["seq"].strip().upper(),
                f"Forward annealing region mismatch: got {f_anneal[1:5]}",
            )
            f_primer_seq = self.spacer + self.bsmbi_site + f_anneal

            # Design reverse primer
            r_5prime = overhang_start + 5
            r_seq_length = min_binding_length
            r_anneal = mutated_context[r_5prime:r_5prime + r_seq_length]
            while self.utils.calculate_tm(r_anneal) < tm_threshold and r_seq_length < self.max_binding_length:
                r_seq_length += 1
                r_anneal = mutated_context[r_5prime - r_seq_length:r_5prime]
            r_anneal = self.utils.reverse_complement(r_anneal)
            r_primer_seq = self.spacer + self.bsmbi_site + r_anneal

            self.log_step("Reverse Primer Sequence",
                          f"Annealing region: {r_anneal[1:5]} vs expected: {overhang_data['bottom_overhang']['seq']}")
            self.validate(
                r_anneal[1:5].strip().upper(
                ) == overhang_data["bottom_overhang"]["seq"].strip().upper(),
                f"Reverse annealing region mismatch: got {r_anneal[1:5]}",
            )

            # Create primer objects
            f_primer = Primer(
                name=(primer_name +
                      "_forward") if primer_name else f"primer_{i}_forward",
                sequence=f_primer_seq,
                binding_region=f_anneal,
                tm=self.utils.calculate_tm(f_anneal),
                gc_content=self.utils.gc_content(f_anneal),
                length=len(f_primer_seq)
            )
            r_primer = Primer(
                name=(primer_name +
                      "_reverse") if primer_name else f"primer_{i}_reverse",
                sequence=r_primer_seq,
                binding_region=r_anneal,
                tm=self.utils.calculate_tm(r_anneal),
                gc_content=self.utils.gc_content(r_anneal),
                length=len(r_primer_seq)
            )

            mutation_primer = MutationPrimer(
                site=mut_obj["site"],
                position=mut_obj["position"],
                forward=f_primer,
                reverse=r_primer,
                mutation_info=mut_obj
            )

            self.log_step("Primer Design Result",
                          f"Designed primer pair for site {mut_obj['site']} at position {mut_obj['position']}")
            mutation_primers.append(mutation_primer)

        return mutation_primers

    @DebugMixin.debug_wrapper
    def generate_GG_edge_primers(self, idx, sequence, mtk_part_left, mtk_part_right, primer_name):
        self.log_step("Generate Edge Primers", f"Sequence {idx}",
                      {"length": len(sequence), "left_part": mtk_part_left, "right_part": mtk_part_right})
        seq_str = str(sequence)
        seq_length = len(seq_str)
        f_length = self.utils.calculate_optimal_primer_length(
            seq_str, 0, 'forward')
        r_length = self.utils.calculate_optimal_primer_length(
            seq_str, len(seq_str), 'reverse')

        self.log_step("Calculated Lengths", "Optimal primer lengths",
                      {"forward_length": f_length, "reverse_length": r_length})

        overhang_5p = self.utils.get_mtk_partend_sequence(
            mtk_part_left, "forward", kozak=self.kozak)
        overhang_3p = self.utils.get_mtk_partend_sequence(
            mtk_part_right, "reverse", kozak=self.kozak)
        self.validate(
            overhang_5p is not None and overhang_3p is not None,
            "Retrieved MTK part end overhangs",
            {"5_prime": overhang_5p, "3_prime": overhang_3p}
        )

        f_binding = seq_str[:f_length]
        r_binding = str(Seq(seq_str[-r_length:]).reverse_complement())
        forward_primer = overhang_5p + f_binding
        reverse_primer = overhang_3p + r_binding

        primers = {
            "forward_primer": {
                "name": f"{primer_name}_F",
                "sequence": forward_primer,
                "binding_region": f_binding,
                "tm": self.utils.calculate_tm(f_binding),
                "gc_content": self.utils.gc_content(f_binding),
                "length": len(forward_primer)
            },
            "reverse_primer": {
                "name": f"{primer_name}_R",
                "sequence": reverse_primer,
                "binding_region": r_binding,
                "tm": self.utils.calculate_tm(r_binding),
                "gc_content": self.utils.gc_content(r_binding),
                "length": len(reverse_primer)
            },
            "product_size": seq_length
        }

        self.validate(
            'forward_primer' in primers and 'reverse_primer' in primers,
            "Edge primers created successfully",
            {"product_size": seq_length}
        )
        return primers
