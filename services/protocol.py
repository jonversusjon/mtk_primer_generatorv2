# services/protocol.py
from Bio.Seq import Seq
from typing import List, Dict, Optional, Union
from .base import GoldenGateDesigner
from .sequence_prep import SequencePreparator
from .primer_design import PrimerDesigner
from .primer_select import PrimerSelector
from .mutation_analyzer import MutationAnalyzer
from .mutation_optimizer import MutationOptimizer
from .utils import GoldenGateUtils
from itertools import product

class GoldenGateProtocol(GoldenGateDesigner):
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
        self.mutation_optimizer = MutationOptimizer(
            verbose=verbose)
        self.mutation_analyzer = MutationAnalyzer(
            sequence=seq, 
            codon_usage_dict=codon_usage_dict,
            max_mutations=max_mutations, 
            verbose=verbose)
        
        self.seq=seq
        primer_name=primer_name
        template_seq=template_seq
        self.verbose=verbose
        self.codon_usage_dict=codon_usage_dict
        self.max_mutations=max_mutations
        self.output_tsv_path=output_tsv_path
        
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
            try:
                # 1. Remove start/stop codons and find restriction sites
                with self.debug_context("Preprocessing sequence"):
                    processed_seq = self.sequence_preparator.preprocess_sequence(single_seq)
                    sites_to_mutate = self.sequence_preparator.find_and_summarize_sites(processed_seq, i)
                
                if sites_to_mutate:
                    with self.debug_context("Mutation analysis"):
                        # 2. Generate all possible silent mutations using MutationAnalyzer
                        mutation_options = self.mutation_analyzer.get_all_mutations(sites_to_mutate=sites_to_mutate)
                        
                        print(f"mutation_options: {mutation_options}")
                        if not mutation_options:
                            self.logger.warning(f"No valid mutations found for sequence {i}")
                            primers = self.primer_designer.generate_GG_edge_primers(
                                processed_seq, i,
                            )
                            all_primer_data.extend(primers)
                            continue
                        
                        # 3. Optimize and prioritize mutations using MutationOptimizer
                        optimized_mutations = self.mutation_optimizer.optimize_mutations(
                            mutation_options=mutation_options,
                        )
                    with self.debug_context("Primer design"):
                        # 4. Design primers using PrimerDesigner
                        primers = self.primer_designer.design_primers_for_mutation(
                            sequence=processed_seq,
                            mutations=optimized_mutations,
                            seq_index=i,
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
    
    
    
    
    
    
    
    

    def _process_single_sequence(
        self,
        single_seq: str,
        seq_index: int,
        part_num_left: List[str],
        part_num_right: List[str],
        codon_usage_dict: Dict[str, Dict[str, float]],
        max_mutations: int,
        primer_name: Optional[List[str]],
        template_seq: Optional[str],
        kozak: str
    ) -> List[List[str]]:
        """Process a single sequence and return its primer data."""
        
        sequence = self.sequence_preparator._preprocess_sequence(single_seq)
        sites_to_mutate = self.sequence_preparator._find_and_summarize_sites(sequence, seq_index)
        
        self.state['current_step'] = 'mutation_analysis'
        if self.verbose: print("=== Starting mutation_analysis ===")
        if self.mutation_analyzer is None:
            self.mutation_analyzer = MutationAnalyzer(
                sequence=sequence,
                codon_usage_dict=codon_usage_dict,
                sites_to_mutate=sites_to_mutate,
                max_mutations=max_mutations,
                verbose=self.verbose
            )
            
        all_mutations = self.mutation_analyzer._get_all_mutations()
        print(f"all_mutations: {all_mutations}")
        # print(f"                best_mutations: {best_mutations}")
        # self.state['current_step'] = 'primer_design'
        # print(f"                seq_index: {seq_index}")
        # print(f"                best_mutations: {best_mutations}")
        # print(f"                part_num_left: {part_num_left}")
        # print(f"                part_num_right: {part_num_right}")
        # print(f"                primer_name: {primer_name}")
        # print(f"                kozak: {kozak}")
        
        # # Calculate cut ranges and create compatibility tensor
        # compatibility_tensor = self.mutation_optimizer.create_compatibility_tensor(
        #     sequence=sequence_str,
        #     mutation_options=site_mutations
        # )

        # best_mutations, valid_positions = self.mutation_optimizer.find_optimal_mutations(mutation_result, compatibility_tensor)
        # internal_mutation_primers = self._design_internal_mutation_primers(
        #     sequence=sequence,
        #     seq_index=seq_index,
        #     mutation_options=best_mutations,
        #     part_num_left=part_num_left,
        #     part_num_right=part_num_right,
        #     primer_name=primer_name,
        #     template_seq=template_seq,
        #     kozak=kozak
        # )

        # print(f'internal_mutation_primers: {internal_mutation_primers}')
        # return internal_mutation_primers

            
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