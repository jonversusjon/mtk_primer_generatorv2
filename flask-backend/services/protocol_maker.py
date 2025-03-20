# services/protocol.py

from typing import List, Dict, Optional

from models import DomesticationResult, SequenceToDomesticate, MutationSetCollection
from services import (
    SequencePreparator,
    RestrictionSiteDetector,
    MutationAnalyzer,
    MutationOptimizer,
    PrimerDesigner,
    ReactionOrganizer,
)
from .utils import GoldenGateUtils
from log_utils import logger

class ProtocolMaker():
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
        self.debug = debug

        self.utils = GoldenGateUtils()
        self.sequence_preparator = SequencePreparator()
        self.rs_analyzer = RestrictionSiteDetector(codon_dict=codon_usage_dict)
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
        logger.debug(f"Protocol maker initialized with codon_usage_dict: {codon_usage_dict}")
        if verbose:
            logger.log_step("Verbose Mode", "Protocol maker is running in verbose mode.")

        self.sequences_to_domesticate: List[SequenceToDomesticate] = sequences_to_domesticate
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

    def create_gg_protocol(self, progress_callback) -> dict:
        """
        Main function to orchestrate the Golden Gate protocol creation.
        Returns:
            dict: A dictionary containing protocol details.
        """
        print("Starting Golden Gate protocol creation...")
        logger.log_step("Protocol Start", "Starting Golden Gate protocol creation...")

        result_data = {}

        for idx, seq_to_dom in enumerate(self.sequences_to_domesticate):
            logger.log_step("Process Sequence", f"Processing sequence {idx+1} of {len(self.sequences_to_domesticate)}",
                            {"sequence_index": idx})
            dom_result = DomesticationResult(
                sequence_index=idx,
                mtk_part_left=seq_to_dom.mtk_part_left,
                mtk_part_right=seq_to_dom.mtk_part_right,
            )

            # 1. Preprocess sequence (remove start/stop codons, etc.)
            logger.log_step("Preprocessing", f"Preprocessing sequence at index {idx}")
            processed_seq, message, _ = self.sequence_preparator.preprocess_sequence(
                seq_to_dom.sequence, seq_to_dom.mtk_part_left)
            if message:
                dom_result.messages.append(message)
            dom_result.processed_sequence = str(processed_seq) if processed_seq else str(seq_to_dom.sequence)

            # 2. Find restriction sites
            logger.log_step("Restriction Site Detection", f"Detecting restriction sites for sequence index {idx}")
            sites_to_mutate = self.rs_analyzer.find_sites_to_mutate(processed_seq, idx)
            dom_result.restriction_sites = sites_to_mutate

            # 3. Mutation analysis and mutation primer design
            mutation_primers = {}
            if sites_to_mutate:
                logger.log_step("Mutation Analysis", f"Analyzing mutations for sequence index {idx}")
                mutation_options = self.mutation_analyzer.get_all_mutations(sites_to_mutate)
                optimized_mutations: MutationSetCollection = None
                if mutation_options:
                    logger.log_step("Mutation Optimization", f"Optimizing mutations for sequence index {idx}")
                    optimized_mutations = self.mutation_optimizer.optimize_mutations(
                        mutation_options=mutation_options
                    )

                    logger.log_step("Primer Design", f"Designing mutation primers for sequence index {idx}")
                    mutation_primers = self.primer_designer.design_mutation_primers(
                        mutation_sets=optimized_mutations,
                        primer_name=seq_to_dom.primer_name,
                        max_results=self.max_results,
                    )
                    dom_result.mut_primers = mutation_primers
                logger.log_step("Mutation Primers", f"Mutation primers designed: {mutation_primers}")

            # 4. Generate edge primers
            logger.log_step("Edge Primer Design", f"Designing edge primers for sequence index {idx}")
            dom_result.edge_primers = self.primer_designer.generate_GG_edge_primers(
                idx, processed_seq, seq_to_dom.mtk_part_left, seq_to_dom.mtk_part_right, seq_to_dom.primer_name
            )
            logger.log_step("Edge Primer Result", "Edge primers generated.",
                            {"edge_forward": dom_result.edge_primers.forward,
                             "edge_reverse": dom_result.edge_primers.reverse})

            # 5. Group primers into PCR reactions
            print("Grouping primers into PCR reactions...")
            logger.log_step("PCR Reaction Grouping", "Grouping primers into PCR reactions using designed primers.")
            dom_result.PCR_reactions = self.reaction_organizer.group_primers_into_pcr_reactions(dom_result)
            print("Finished grouping primers into PCR reactions...")

            result_data[idx] = dom_result

        return result_data
