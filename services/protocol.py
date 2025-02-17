# services/protocol.py
from typing import List, Dict, Optional, Union
from .base import PrimerDesignLogger
from .sequence_prep import SequencePreparator
from .primer_design import PrimerDesigner
from .primer_select import PrimerSelector
from .mutation_analyzer import MutationAnalyzer
from .mutation_optimizer import MutationOptimizer
from .utils import GoldenGateUtils
from itertools import product
from .townsend import generate_GG_protocol

class GoldenGateProtocol(PrimerDesignLogger):
    def __init__(
        self,
        seq: List[str],
        codon_usage_dict: Dict[str, Dict[str, float]],
        part_num_left: List[str],
        part_num_right: List[str],
        max_mutations: int,
        primer_name: Optional[List[str]] = None,
        template_seq: Optional[str] = None,
        kozak: str = "MTK",
        output_tsv_path: str = "designed_primers.tsv",
        verbose: bool = False):
        super().__init__(verbose=verbose)
        self.utils = GoldenGateUtils()
        self.sequence_preparator = SequencePreparator()
        self.primer_designer = PrimerDesigner(
            part_num_left=part_num_left,
            part_num_right=part_num_right,
            kozak=kozak,
            verbose=verbose
        )

        self.primer_selector = PrimerSelector()
        print(f"initialize GoldenGateProtocol - verbose: {verbose}")

        self.mutation_optimizer = MutationOptimizer(verbose=verbose)
        self.mutation_analyzer = MutationAnalyzer(
            sequence=seq, 
            codon_usage_dict=codon_usage_dict,
            max_mutations=max_mutations, 
            verbose=verbose)
        
        self.seq=seq
        primer_name=primer_name
        self.template_seq=template_seq
        self.verbose=verbose
        self.codon_usage_dict=codon_usage_dict
        self.max_mutations=max_mutations
        self.output_tsv_path=output_tsv_path
        self.part_num_left=part_num_left
        
        self.state = {
            'current_sequence_index': 0,
            'current_step': '',
            'mutations_found': [],
            'primers_designed': []
        }

    def create_gg_protocol(
        self,
    ) -> List[List[str]]:
        """
        Main function to orchestrate the Golden Gate protocol creation.
        Each step is clearly defined and its output feeds into the next step.
        """
        self.logger.info("Starting Golden Gate protocol creation...")
        
        # Process each sequence and collect all primer data
        all_primer_data = []
        for i, single_seq in enumerate(self.seq):
            self.logger.info("Running Townsend primer generator")
            protocol = generate_GG_protocol(single_seq, self.part_num_left[0], True)
            # print(f"intended error: {intended_erro}")
            try:
                # 1. Remove start/stop codons and find restriction sites
                with self.debug_context("Preprocessing sequence"):
                    processed_seq = self.sequence_preparator.preprocess_sequence(single_seq)
                    sites_to_mutate = self.sequence_preparator.find_and_summarize_sites(processed_seq, i)
                
                if sites_to_mutate:
                    with self.debug_context("Mutation analysis"):
                        # 2. Generate all possible silent mutations using MutationAnalyzer
                        mutation_options = self.mutation_analyzer.get_all_mutations(
                            sites_to_mutate=sites_to_mutate)
                        
                        # 3. Optimize and prioritize mutations if there are valid options
                        if mutation_options:
                            optimized_mutations, compatibility_matrices = self.mutation_optimizer.optimize_mutations(
                                sequence=processed_seq, mutation_options=mutation_options
                            )
                        else:
                            optimized_mutations, compatibility_matrices = None, None

                        
                    with self.debug_context("Primer design"):
                        # 4. Design primers using PrimerDesigner
                        primers = self.primer_designer.design_primers(
                            sequence=processed_seq,
                            seq_index=i,
                            mutations=optimized_mutations,
                            compatibility_matrices=compatibility_matrices,
                            template_seq=self.template_seq
                        )
                else:
                    
                    self.logger.info(f"No restriction sites found in sequence {i}")
                    primers = self.primer_designer.generate_GG_edge_primers(
                        processed_seq, i, 
                    )
                    all_primer_data.extend(primers)
                    continue
                
                all_primer_data.extend(primers)
                
            except Exception as e:
                self.logger.error(f"Error processing sequence {i}: {str(e)}")
                continue
            
        # Step 3: Save results and return
        self._save_primers_to_tsv(all_primer_data, self.output_tsv_path)
        return all_primer_data
    
    def _save_primers_to_tsv(self, primer_data: List[List[str]], output_tsv_path: str) -> None:
        """Saves primer data to a TSV file."""
        with self.debug_context("save_primers_to_tsv"):
            if not primer_data:
                self.logger.warning("No primer data to save.")
                return

            try:
                with open(output_tsv_path, "w") as tsv_file:
                    tsv_file.write("Primer Name\tSequence\tAmplicon\n")
                    for row in primer_data:
                        tsv_file.write("\t".join(map(str, row)) + "\n")
            except IOError as e:
                self.logger.error(f"Error writing to file {output_tsv_path}: {e}")
                raise