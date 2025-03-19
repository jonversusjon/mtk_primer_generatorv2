import numpy as np
from Bio.Seq import Seq
from .utils import GoldenGateUtils
from log_utils import logger
from models import Primer, MutationPrimerPair, EdgePrimerPair, MutationSet, MutationSetCollection, OverhangOption
import logging
from typing import List

class PrimerDesigner():
    """
    Handles primer design for Golden Gate assembly.
    """

    def __init__(self, kozak: str = "MTK", verbose: bool = False, debug: bool = False):
        self.utils = GoldenGateUtils()
        self.verbose = verbose
        self.debug = debug

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

        logger.validate(
            self.part_end_dict is not None,
            "MTK part end sequences loaded successfully",
            {"kozak": self.kozak}
        )

    def design_mutation_primers(self, mutation_sets: MutationSetCollection,
                                primer_name: str = None, max_results: int = 1):
        """
        Designs mutation primers for the provided mutation sets using compatibility matrices.
        Returns a dict mapping each mutation set index to a list of MutationPrimer objects.
        """
        all_primers = {}
        mut_site_keys = mutation_sets.sites_to_mutate
        
        for idx, mut_set in enumerate(mutation_sets.sets):
            comp_matrix = mut_set.compatibility

            # Validate and log valid overhang combinations.
            valid_coords = np.argwhere(comp_matrix == 1)
            logger.validate(
                valid_coords.size > 0,
                f"Found {np.count_nonzero(comp_matrix)} valid overhang combination(s)",
                {"matrix_size": comp_matrix.size}
            )
            logger.log_step("Valid Coordinates",
                            f"Mutation set {idx+1}: Valid coordinates: {valid_coords.tolist()}")
            logger.log_step("Matrix Visualization",
                            f"Compatibility matrix for set {idx+1}",
                            logger.visualize_matrix(comp_matrix))

            # Determine which coordinate combinations to process.
            if max_results == 0:
                coords_to_process = valid_coords.tolist()
                logger.log_step("Processing All Coordinates",
                                f"Processing all {len(coords_to_process)} valid coordinate combination(s).")
            else:
                sample_size = min(max_results, valid_coords.shape[0])
                selected_indices = np.random.choice(valid_coords.shape[0], size=sample_size, replace=False)
                coords_to_process = [valid_coords[i].tolist() for i in selected_indices]
                logger.log_step("Random Coordinate Selection",
                                f"Randomly selected {sample_size} coordinate combination(s): {coords_to_process}")

            mut_set_primer_sets = []
            for coords in coords_to_process:
                # each set of coords represent a valid set of mut_primer pairs
                logger.log_step("Coordinate Processing", f"Processing coordinate combination: {coords}")
                # Ensure the length of the coordinate combination matches the number of mutation sites.
                if len(coords) != len(mut_set.mutations):
                    raise ValueError(
                        f"Coordinate length {len(coords)} does not match number of mutation sites {len(mut_set)}"
                    )

                # Construct MutationPrimer objects for the current combination.
                mut_set_primer_set: List[MutationPrimerPair] = self._construct_mutation_primer_set(
                    mut_site_keys=mut_site_keys,
                    mutation_set=mut_set,
                    selected_coords=coords,
                    primer_name=primer_name
                )
                
                logger.validate(
                    mut_set_primer_set is not None,
                    "Successfully constructed mutation primers",
                    {"primer_count": len(mut_set_primer_set) if mut_set_primer_set else 0}
                )
                
                if mut_set_primer_set:
                    logger.log_step("Constructed Primers",
                                    f"Constructed mutation primer pairs for combination {coords}: {mut_set_primer_set}")
                    mut_set_primer_sets.extend(mut_set_primer_set)

            logger.log_step("Set Summary",
                            f"Total primer pairs constructed for mutation set {idx+1}: {len(mut_set_primer_set)}")
            all_primers[idx] = mut_set_primer_set

        if not all_primers:
            if self.debug:
                logger.log_step("Design Failure", "Failed to design primers for any mutation set", level=logging.WARNING)
            return None

        return all_primers


    def _construct_mutation_primer_set(
        self,
        mut_site_keys: List[str],
        mutation_set: MutationSet,
        selected_coords: list,
        primer_name: str = None,
        min_binding_length: int = 10) -> list[MutationPrimerPair]:
        """
        Constructs forward and reverse mutation primers for each mutation in the mutation_set.
        Expects mutation_set to be a list of Mutation objects (Pydantic models) and uses their
        overhang_options (a list of OverhangOption objects) via attribute access.
        """
        logger.log_step("Construct Primers", f"Binding length: {min_binding_length}", {"selected_coords": selected_coords})
        mutation_primers = []
        # TODO: Make sure that mutation_sets is made up of pydantic objects
        for i, mutation in enumerate(mutation_set.mutations):
            selected_overhang = selected_coords[i]
            # Access overhang_options as a list of OverhangOption objects.
            overhang_options = mutation.overhang_options
            if selected_overhang >= len(overhang_options):
                raise IndexError(f"Selected overhang index {selected_overhang} out of range for mutation site {i}")
            overhang_data: OverhangOption = overhang_options[selected_overhang]
            mutated_context = mutation.mut_context
            overhang_start = overhang_data.overhang_start_index

            tm_threshold = self.default_params["tm_threshold"]

            # Design forward primer
            f_5prime = overhang_start - 1
            f_seq_length = min_binding_length
            f_anneal = mutated_context[f_5prime:f_5prime + f_seq_length]
            while self.utils.calculate_tm(f_anneal) < tm_threshold and f_seq_length < self.max_binding_length:
                f_seq_length += 1
                f_anneal = mutated_context[f_5prime:f_5prime + f_seq_length]

            logger.log_step("Forward Primer",
                            f"Annealing region: {f_anneal[1:5]} vs expected: {overhang_data.top_overhang}")
            logger.validate(
                f_anneal[1:5].strip().upper() == overhang_data.top_overhang.strip().upper(),
                f"Forward annealing region mismatch: got {f_anneal[1:5]}"
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

            logger.log_step("Reverse Primer",
                            f"Annealing region: {r_anneal[1:5]} vs expected: {overhang_data.bottom_overhang}")
            logger.validate(
                r_anneal[1:5].strip().upper() == overhang_data.bottom_overhang.strip().upper(),
                f"Reverse annealing region mismatch: got {r_anneal[1:5]}"
            )

            # Create Primer objects.
            f_primer = Primer(
                name=(primer_name + "_forward") if primer_name else f"primer_{i}_forward",
                sequence=f_primer_seq,
                binding_region=f_anneal,
                tm=self.utils.calculate_tm(f_anneal),
                gc_content=self.utils.gc_content(f_anneal),
                length=len(f_primer_seq)
            )
            r_primer = Primer(
                name=(primer_name + "_reverse") if primer_name else f"primer_{i}_reverse",
                sequence=r_primer_seq,
                binding_region=r_anneal,
                tm=self.utils.calculate_tm(r_anneal),
                gc_content=self.utils.gc_content(r_anneal),
                length=len(r_primer_seq)
            )

            # For MutationPrimerPair, we need site and position.
            # If mut_obj does not have these, use fallbacks.
            site_val = mut_site_keys[i]
            position_val = mutation.first_mut_idx

            # Construct the MutationPrimerPair using the Mutation object as mutation_info.
            mutation_primer_pair = MutationPrimerPair(
                site=site_val,
                position=position_val,
                forward=f_primer,
                reverse=r_primer,
                mutation=mutation
            )
            logger.log_step("Primer Pair Constructed",
                            f"Designed primer pair for site {site_val} at position {position_val}")
            mutation_primers.append(mutation_primer_pair)

        return mutation_primers



    def generate_GG_edge_primers(self, idx, sequence, mtk_part_left, mtk_part_right, primer_name) -> EdgePrimerPair:
        logger.log_step("Generate Edge Primers",
                        f"Sequence {idx}",
                        {"length": len(sequence), "left_part": mtk_part_left, "right_part": mtk_part_right})
        seq_str = str(sequence)
        seq_length = len(seq_str)
        f_length = self.utils.calculate_optimal_primer_length(seq_str, 0, 'forward')
        r_length = self.utils.calculate_optimal_primer_length(seq_str, len(seq_str), 'reverse')
        logger.log_step("Optimal Primer Lengths",
                        "Calculated optimal primer lengths",
                        {"forward_length": f_length, "reverse_length": r_length})
        
        overhang_5p = self.utils.get_mtk_partend_sequence(mtk_part_left, "forward", kozak=self.kozak)
        overhang_3p = self.utils.get_mtk_partend_sequence(mtk_part_right, "reverse", kozak=self.kozak)
        logger.validate(
            overhang_5p is not None and overhang_3p is not None,
            "Retrieved MTK part end overhangs",
            {"5_prime": overhang_5p, "3_prime": overhang_3p}
        )
        
        f_binding = seq_str[:f_length]
        r_binding = str(Seq(seq_str[-r_length:]).reverse_complement())
        
        forward_seq = overhang_5p + f_binding
        reverse_seq = overhang_3p + r_binding

        f_primer = Primer(
            name=f"{primer_name}_F",
            sequence=forward_seq,
            binding_region=f_binding,
            tm=self.utils.calculate_tm(f_binding),
            gc_content=self.utils.gc_content(f_binding),
            length=len(forward_seq)
        )
        r_primer = Primer(
            name=f"{primer_name}_R",
            sequence=reverse_seq,
            binding_region=r_binding,
            tm=self.utils.calculate_tm(r_binding),
            gc_content=self.utils.gc_content(r_binding),
            length=len(reverse_seq)
        )
        
        logger.validate(
            f_primer is not None and r_primer is not None,
            "Edge primers created successfully",
            {"product_size": seq_length}
        )
        
        edge_primers = EdgePrimerPair(
            forward=f_primer,
            reverse=r_primer
        )
        return edge_primers

