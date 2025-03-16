# services/protocol.py

from typing import List, Dict, Optional

from models import DomesticationResult, SequenceToDomesticate
from services import (
    SequencePreparator,
    RestrictionSiteDetector,
    MutationAnalyzer,
    MutationOptimizer,
    PrimerDesigner,
    ReactionOrganizer,
)
from .utils import GoldenGateUtils
from services.debug import MutationDebugger, DebugMixin, debug_context
from config.logging_config import logger

class GoldenGateProtocol(DebugMixin):
    """
    Orchestrates the Golden Gate protocol by managing sequence preparation,
    primer design, mutation analysis, and optimization.
    """

    def __init__(
        self,
        sequences_to_domesticate: List[SequenceToDomesticate],
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

        self.sequences_to_domesticate = sequences_to_domesticate
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
        print("Starting Golden Gate protocol creation...")
        logger.info("Starting Golden Gate protocol creation...")

        result_data = {}

        for idx, seq_to_dom in enumerate(self.sequences_to_domesticate):

            dom_result = DomesticationResult(
                sequence_index=idx,
                mtk_part_left = seq_to_dom.mtk_part_left,
                mtk_part_right = seq_to_dom.mtk_part_right,
            )

            print(
                f"Processing sequence {idx+1}/{len(self.sequences_to_domesticate)}")

            # 1. Preprocess sequence (remove start/stop codons, etc.)
            with debug_context("Preprocessing sequence"):
                processed_seq, message, _ = self.sequence_preparator.preprocess_sequence(
                    seq_to_dom.sequence, seq_to_dom.mtk_part_left)

                if message:
                    dom_result.messages.append(message)

                dom_result.processed_sequence = str(
                    processed_seq) if processed_seq else str(seq_to_dom.sequence)

            # 2. Find restriction sites
            with debug_context("Finding restriction sites"):
                sites_to_mutate = self.rs_analyzer.find_sites_to_mutate(
                    processed_seq, idx)

                dom_result.restriction_sites = sites_to_mutate

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

                        dom_result.mutation_options = {
                            "all_mutation_options": optimized_mutations,
                            "compatibility": compatibility_matrices
                        }

                with debug_context("Mutation primer design"):
                    mutation_primers = self.primer_designer.design_mutation_primers(
                        mutation_sets=optimized_mutations,
                        comp_matrices=compatibility_matrices,
                        primer_name=seq_to_dom.primer_name,
                        max_results=self.max_results,
                    )
                    dom_result.mutation_primers = mutation_primers

            self.log_step(f"Mutation primers designed: {mutation_primers}")
            self.log_step("designing edge primers...")
            
            # 4. Generate edge primers
            dom_result.edge_primers = self.primer_designer.generate_GG_edge_primers(
                idx, processed_seq, seq_to_dom.mtk_part_left, seq_to_dom.mtk_part_right, seq_to_dom.primer_name
            )
            
            self.log_step("Sequence Data",  
                          "Edge primers generated.",
                          f"edge_forward: {dom_result.edge_primers.forward}",
                          f"edge_reverse: {dom_result.edge_primers.reverse}")
            
            # 5. Group primers into PCR reactions
            print("Grouping primers into PCR reactions...")
            self.log_step("Group PCR Reactions",
                          "sequence_data: {sequence_data}")
            dom_result.PCR_reactions = self.reaction_organizer.group_primers_into_pcr_reactions(
                dom_result)
            print("Finished grouping primers into PCR reactions...")
            
            result_data[idx] = dom_result
            
        return result_data