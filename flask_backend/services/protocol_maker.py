# services/protocol.py

from typing import Dict, Optional, List
from functools import partial

from flask_backend.models import (
    DomesticationResult,
    SequenceToDomesticate,
    MutationSetCollection,
    RestrictionSite
)
from flask_backend.services import (
    SequencePreparator,
    RestrictionSiteDetector,
    MutationAnalyzer,
    MutationOptimizer,
    PrimerDesigner,
    ReactionOrganizer,
)
from flask_backend.services.utils import GoldenGateUtils
from flask_backend.logging import logger

class ProtocolMaker():
    """
    Orchestrates the Golden Gate protocol by managing sequence preparation,
    primer design, mutation analysis, and optimization.
    """

    def __init__(
        self,
        request_idx: int,
        sequence_to_domesticate: SequenceToDomesticate,
        codon_usage_dict: Dict[str, Dict[str, float]],
        max_mutations: int,
        template_seq: Optional[str] = None,
        kozak: str = "MTK",
        output_tsv_path: str = "designed_primers.tsv",
        max_results: str = "one",
        verbose: bool = False,
        debug: bool = False,
        job_id: Optional[str] = None,
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
        logger.debug(f"Protocol maker for sequence {request_idx+1} initialized with codon_usage_dict: {codon_usage_dict}")
        if verbose:
            logger.log_step("Verbose Mode", "Protocol maker is running in verbose mode.")

        self.request_idx = request_idx
        self.seq_to_dom: SequenceToDomesticate = sequence_to_domesticate
        self.template_seq = template_seq
        self.kozak = kozak
        self.verbose = verbose
        self.codon_usage_dict = codon_usage_dict
        self.max_mutations = max_mutations
        self.output_tsv_path = output_tsv_path
        self.max_results = max_results
        self.job_id = job_id

    def create_gg_protocol(self, send_update) -> dict:
        """
        Main function to orchestrate the Golden Gate protocol creation.
        Returns:
            dict: A dictionary containing protocol details.
        """
        # Initialize the result object
        logger.log_step("Protocol Start - Process Sequence", f"Processing sequence {self.request_idx+1}")        
        dom_result = DomesticationResult(
            sequence_index=self.request_idx,
            mtk_part_left=self.seq_to_dom.mtk_part_left,
            mtk_part_right=self.seq_to_dom.mtk_part_right,
        )
        
        # 1 - Preprocess the sequence
        logger.log_step("Preprocessing", f"Preprocessing sequence at index {self.request_idx+1}")
        processed_seq, valid_seq = self.sequence_preparator.preprocess_sequence(
            self.seq_to_dom.sequence, 
            self.seq_to_dom.mtk_part_left,
            partial(send_update, step="Preprocessing")
        )
        # TODO: need a way to gracefully handle invalid sequences that snuck past
        # the frontend validation.
        if not valid_seq:
            logger.log_step("Preprocessing Error", "Invalid sequence detected during preprocessing.")
            return dom_result
        dom_result.processed_sequence = str(processed_seq) if processed_seq else str(self.seq_to_dom.sequence)
    
        # 2 - Restriction Site Detection
        logger.log_step("Restriction Site Detection", f"Detecting restriction sites for sequence {self.request_idx+1}")
        restriction_sites: List[RestrictionSite] = self.rs_analyzer.find_restriction_sites(
            processed_seq,
            partial(send_update, step="Restriction Sites")
        )
        dom_result.restriction_sites = restriction_sites

        mutation_primers = {}
        if restriction_sites:
            # 3 - Mutation Analysis
            logger.log_step("Mutation Analysis", f"Analyzing mutations for sequence {self.request_idx+1}")
            mutation_options = self.mutation_analyzer.get_all_mutations(
                restriction_sites,
                partial(send_update, step="Mutation Analysis")
            )
            
            mutation_sets: MutationSetCollection = None 
            if mutation_options:
                logger.log_step("Mutation Optimization", f"Optimizing mutations for sequence {self.request_idx+1}")
                mutation_sets: MutationSetCollection = self.mutation_optimizer.optimize_mutations(
                    mutation_options,
                    partial(send_update, step="Mutation Analysis")
                )

                # 4 - Primer Design
                logger.log_step("Primer Design", f"Designing mutation primers for sequence {self.request_idx+1}")
                mutation_primers = self.primer_designer.design_mutation_primers(
                    mutation_sets,
                    self.seq_to_dom.primer_name if self.seq_to_dom.primer_name else f"Primer{self.request_idx+1}",
                    self.max_results if self.max_results else "one",
                    partial(send_update, step="Primer Design")
                )
                dom_result.mut_primers = mutation_primers
        
        # Generate edge primers
        logger.log_step("Edge Primer Design", f"Designing edge primers for sequence {self.request_idx+1}")
        dom_result.edge_primers = self.primer_designer.generate_GG_edge_primers(
            self.request_idx,
            processed_seq,
            self.seq_to_dom.mtk_part_left,
            self.seq_to_dom.mtk_part_right,
            self.seq_to_dom.primer_name,
            partial(send_update, step="Primer Design")
        )
        logger.log_step("Edge Primer Result", "Edge primers generated.",
                        {"edge_forward": dom_result.edge_primers.forward,
                        "edge_reverse": dom_result.edge_primers.reverse})
        
        logger.log_step("PCR Reaction Grouping", "Grouping primers into PCR reactions using designed primers.")
        
        nested_reactions = self.reaction_organizer.group_primers_into_pcr_reactions(
            dom_result,
            partial(send_update, step="PCR Reaction Grouping")
        )

        # Flatten the reactions into a list of PCRReaction
        dom_result.PCR_reactions = [
            reaction
            for mutation_set in nested_reactions["mutation_sets"]
            for solution in mutation_set["solutions"]
            for reaction in solution["reactions"]
        ]

        print("Finished grouping primers into PCR reactions...")

        return dom_result