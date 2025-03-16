# services/protocol.py

from typing import List, Dict, Optional, Tuple
from .sequence_prep import SequencePreparator
from .rs_detector import RestrictionSiteDetector
from .mutation_analyzer import MutationAnalyzer
from .mutation_optimizer import MutationOptimizer
from .primer_design import PrimerDesigner
from .reactions import ReactionOrganizer
from .utils import GoldenGateUtils
from config.logging_config import logger
from .base import debug_context
from services.debug.debug_utils import MutationDebugger, visualize_matrix
from services.debug.debug_mixin import DebugMixin
from models.sequences import SequenceToDomesticate
from pydantic import validate_call


class GoldenGateProtocol(DebugMixin):
    """
    Orchestrates the Golden Gate protocol by managing sequence preparation,
    primer design, mutation analysis, and optimization.
    class SequenceToDomesticate(BaseModel):
        sequence_index: int
        primer_name: Optional[str] = None
        sequence: str
        mtk_part_left: str
        mtk_part_right: str
        restriction_sites: List[RestrictionSite]
    """

    @validate_call
    def __init__(
        self,
        sequencesToDomesticate: List[SequenceToDomesticate],
        codon_usage_dict: Dict[str, Dict[str, float]],
        max_mutations: int,
        template_seq: Optional[str] = None,
        kozak: str = "MTK",
        output_tsv_path: str = "designed_primers.tsv",
        max_results: int = 1,
        verbose: bool = False,
        debug: bool = False,
    ):
        self.logger = logger.getChild("GoldenGateProtocol")
        self.debug = debug
        self.debugger = None

        if self.debug:
            self.debugger = MutationDebugger(
                parent_logger=logger, use_custom_format=True)
            if hasattr(self.debugger.logger, 'propagate'):
                self.debugger.logger.propagate = False
            self.logger.info("Debug mode enabled for PrimerDesigner")

        self.utils = GoldenGateUtils()
        self.sequence_preparator = SequencePreparator()
        self.rs_analyzer = RestrictionSiteDetector()
        self.mutation_analyzer = MutationAnalyzer(
            codon_usage_dict=codon_usage_dict,
            max_mutations=max_mutations,
            verbose=verbose,
            debug=True,
        )
        self.mutation_optimizer = MutationOptimizer(
            verbose=verbose, debug=True)
        self.primer_designer = PrimerDesigner(
            kozak=kozak, verbose=verbose, debug=True)
        self.reaction_organizer = ReactionOrganizer()
        self.logger.debug(
            f"GoldenGateProtocol initialized with codon_usage_dict: {codon_usage_dict}")
        if verbose:
            self.logger.info("GoldenGateProtocol is running in verbose mode.")

        self.sequencesToDomesticate = sequencesToDomesticate
        self.template_seq = template_seq
        self.kozak = kozak
        self.verbose = verbose
        self.codon_usage_dict = codon_usage_dict
        self.max_mutations = max_mutations
        self.output_tsv_path = output_tsv_path

        self.mtk_partend_sequences = self.utils.get_mtk_partend_sequences()
        self.state = {
            'current_sequence_index': 0,
            'current_step': '',
            'mutations_found': [],
            'primers_designed': []
        }

        self.max_results = max_results

    @DebugMixin.debug_wrapper
    def create_gg_protocol(self, progress_callback) -> dict:
        """
        Main function to orchestrate the Golden Gate protocol creation.
        Returns:
            dict: A dictionary containing protocol details.
        """
        print(f"Starting Golden Gate protocol creation...")
        logger.info("Starting Golden Gate protocol creation...")

        result_data = {}

        for idx, seq_object in enumerate(self.sequencesToDomesticate):
            single_seq, mtk_part_left, mtk_part_right, primer_name = _getSequenceData(
                seq_object, idx+1)

            sequence_data = {
                "sequence_index": idx,
                "processed_sequence": None,
                "mtk_part_left": mtk_part_left,
                "mtk_part_right": mtk_part_right,
                "restriction_sites": [],
                "mutations": None,
                "PCR_reactions": {},
                "messages": [],
                "errors": None
            }

            print(
                f"Processing sequence {idx+1}/{len(self.sequencesToDomesticate)}")

            # 1. Preprocess sequence (remove start/stop codons, etc.)
            with debug_context("Preprocessing sequence"):
                processed_seq, message, _ = self.sequence_preparator.preprocess_sequence(
                    single_seq, mtk_part_left)

                if message:
                    sequence_data["messages"].append(message)

                sequence_data["processed_sequence"] = str(
                    processed_seq) if processed_seq else str(single_seq)

            # 2. Find restriction sites
            with debug_context("Finding restriction sites"):
                sites_to_mutate = self.rs_analyzer.find_sites_to_mutate(
                    processed_seq, idx)

                sequence_data["restriction_sites"] = sites_to_mutate

            # 3. Mutation analysis and mutation primer design
            mutation_primers = {}

            if sites_to_mutate:
                with debug_context("Mutation analysis"):
                    mutation_options = self.mutation_analyzer.get_all_mutations(
                        sites_to_mutate)
                    optimized_mutations, compatibility_matrices = None, None

                    if mutation_options:
                        optimized_mutations, compatibility_matrices = self.mutation_optimizer.optimize_mutations(
                            mutation_options=mutation_options
                        )

                        sequence_data["mutations"] = {
                            "all_mutation_options": optimized_mutations,
                            "compatibility": compatibility_matrices
                        }

                with debug_context("Mutation primer design"):
                    # TODO: add valid solution tracking to enable returning multiple solutions up to
                    # max_results number of solutions to the user
                    mutation_primers = self.primer_designer.design_mutation_primers(
                        mutation_sets=optimized_mutations,
                        comp_matrices=compatibility_matrices,
                        primer_name=primer_name,
                        max_results=self.max_results,
                    )
                    sequence_data["mutation_primers"] = mutation_primers

            print(f"Mutation primers designed: {mutation_primers}")
            print(f"designing edge primers...")
            # 4. Generate edge primers
            edge_primers = self.primer_designer.generate_GG_edge_primers(
                idx, processed_seq, mtk_part_left, mtk_part_right, primer_name
            )
            sequence_data["edge_primers"] = edge_primers
            self.log_step("Sequence Data",
                          "Edge primers generated.",
                          {"edge_forward": edge_primers["forward_primer"],
                           "edge_reverse": edge_primers["reverse_primer"]})
            # 5. Group primers into PCR reactions
            print(f"Grouping primers into PCR reactions...")
            self.log_step("Group PCR Reactions",
                          "sequence_data: {sequence_data}")
            sequence_data["PCR_reactions"] = self.reaction_organizer.group_primers_into_pcr_reactions(
                sequence_data)
            print(f"Finished grouping primers into PCR reactions...")
            # Store the processed sequence data for this sequence number
            result_data[idx] = sequence_data
            self.utils.print_object_schema(
                result_data, indent=0, name="result_data")

        # Pydantic v2 validation of result_data
        from models.protocols import MTKDomesticationProtocol
        try:
            validated_protocol = MTKDomesticationProtocol.model_validate(
                {"result_data": result_data})
            return validated_protocol.model_dump()
        except Exception as e:
            self.logger.error(f"Validation error in protocol creation: {e}")
            raise e


@DebugMixin.debug_wrapper
def _getSequenceData(seq_object: SequenceToDomesticate, i: int) -> Tuple[str, str, str, str]:
    """Extracts sequence data from the provided SequenceToDomesticate model."""
    single_seq = seq_object.sequence
    mtk_part_left = seq_object.mtk_part_left
    mtk_part_right = seq_object.mtk_part_right
    primer_name = seq_object.primer_name or f"Sequence_{i}"
    return single_seq, mtk_part_left, mtk_part_right, primer_name
