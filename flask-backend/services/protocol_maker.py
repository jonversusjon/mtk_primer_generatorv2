# services/protocol.py

from typing import Dict, Optional

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
        
        self.mtk_partend_sequences = self.utils.get_mtk_partend_sequences()
        self.state = {
            'current_sequence_index': 0,
            'current_step': '',
            'mutations_found': [],
            'primers_designed': []
        }


    def create_gg_protocol(self, progress_callback) -> dict:
        """
        Main function to orchestrate the Golden Gate protocol creation.
        Returns:
            dict: A dictionary containing protocol details.
        """
        print("Starting Golden Gate protocol creation...")
        
        progress_callback(
            step="Protocol Start",
            message=f"Starting Golden Gate protocol creation for sequence {self.request_idx+1}..."
        )
        
        logger.log_step("Protocol Start", "Starting Golden Gate protocol creation...")
        logger.log_step("Process Sequence", f"Processing sequence {self.request_idx+1}")
        
        dom_result = DomesticationResult(
            sequence_index=self.request_idx,
            mtk_part_left=self.seq_to_dom.mtk_part_left,
            mtk_part_right=self.seq_to_dom.mtk_part_right,
        )

        # 1. Preprocess sequence (remove start/stop codons, etc.)
        progress_callback(
            step="Preprocessing",
            message=f"Preprocessing sequence at index {self.request_idx+1}..."
        )
        logger.log_step("Preprocessing",
                        f"Preprocessing sequence at index {self.request_idx+1}")
        
        processed_seq, message, _ = self.sequence_preparator.preprocess_sequence(
            self.seq_to_dom.sequence, self.seq_to_dom.mtk_part_left)
        if message:
            dom_result.messages.append(message)
            
        dom_result.processed_sequence = str(processed_seq) if processed_seq else str(self.seq_to_dom.sequence)
        progress_callback(
            step="Preprocessing",
            message=f"Finished preprocessing sequence at index {self.request_idx+1}... {message}"
        )
        
        # 2. Find restriction sites
        logger.log_step("Restriction Site Detection",
                        f"Detecting restriction sites for sequence {self.request_idx+1}")
        sites_to_mutate = self.rs_analyzer.find_sites_to_mutate(processed_seq, self.request_idx)
        dom_result.restriction_sites = sites_to_mutate

        # 3. Mutation analysis and mutation primer design
        mutation_primers = {}
        if sites_to_mutate:
            logger.log_step("Mutation Analysis",
                            f"Analyzing mutations for sequence {self.request_idx+1}")
            mutation_options = self.mutation_analyzer.get_all_mutations(sites_to_mutate)
            optimized_mutations: MutationSetCollection = None
            if mutation_options:
                logger.log_step("Mutation Optimization",
                                f"Optimizing mutations for sequence {self.request_idx+1}")
                optimized_mutations = self.mutation_optimizer.optimize_mutations(
                    mutation_options=mutation_options
                )

                logger.log_step("Primer Design",
                                f"Designing mutation primers for sequence {self.request_idx+1}")
                mutation_primers = self.primer_designer.design_mutation_primers(
                    mutation_sets=optimized_mutations,
                    primer_name=self.seq_to_dom.primer_name,
                    max_results_str=self.max_results,
                )
                dom_result.mut_primers = mutation_primers
            logger.log_step("Mutation Primers",
                            f"Mutation primers designed: {mutation_primers}")

        # 4. Generate edge primers
        logger.log_step("Edge Primer Design",
                        f"Designing edge primers for sequence {self.request_idx+1}")
        dom_result.edge_primers = self.primer_designer.generate_GG_edge_primers(
            self.request_idx,
            processed_seq,
            self.seq_to_dom.mtk_part_left,
            self.seq_to_dom.mtk_part_right,
            self.seq_to_dom.primer_name
        )
        logger.log_step("Edge Primer Result", "Edge primers generated.",
                        {"edge_forward": dom_result.edge_primers.forward,
                        "edge_reverse": dom_result.edge_primers.reverse})

        # 5. Group primers into PCR reactions
        print("Grouping primers into PCR reactions...")
        logger.log_step("PCR Reaction Grouping",
                        "Grouping primers into PCR reactions using designed primers.")
        
        dom_result.PCR_reactions = self.reaction_organizer.group_primers_into_pcr_reactions(dom_result)
        print("Finished grouping primers into PCR reactions...")

        return dom_result

