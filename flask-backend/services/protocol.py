# services/protocol.py
from typing import List, Dict, Optional, Union
from .sequence_prep import SequencePreparator
from .primer_design import PrimerDesigner
from .primer_select import PrimerSelector
from .mutation_analyzer import MutationAnalyzer
from .mutation_optimizer import MutationOptimizer
from .utils import GoldenGateUtils
from itertools import product
from config.logging_config import logger
from .base import debug_context

class GoldenGateProtocol:
    """
    Orchestrates the Golden Gate protocol by managing sequence preparation,
    primer design, mutation analysis, and optimization.
    """

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
        verbose: bool = False
    ):
        self.logger = logger.getChild("GoldenGateProtocol")
        self.utils = GoldenGateUtils()
        self.sequence_preparator = SequencePreparator()
        self.primer_designer = PrimerDesigner(
            part_num_left=part_num_left,
            part_num_right=part_num_right,
            kozak=kozak,
            verbose=verbose
        )
        self.primer_selector = PrimerSelector()
        self.mutation_optimizer = MutationOptimizer(verbose=verbose)
        print(f"GoldenGateProtocol initialized with codon_usage_dict: {codon_usage_dict}")
        self.mutation_analyzer = MutationAnalyzer(
            sequence=seq,
            codon_usage_dict=codon_usage_dict,
            max_mutations=max_mutations,
            verbose=verbose
        )

        """
        Initialize GoldenGateProtocol with provided sequences and parameters.
        
        Args:
            seq (List[str]): List of sequences to process.
            codon_usage_dict (Dict[str, Dict[str, float]]): Codon usage table.
            part_num_left (List[str]): Left part numbers for assembly.
            part_num_right (List[str]): Right part numbers for assembly.
            max_mutations (int): Maximum mutations allowed per site.
            primer_name (Optional[List[str]]): Names of primers.
            template_seq (Optional[str]): Template sequence for primer design.
            kozak (str): Kozak sequence type, default is "MTK".
            output_tsv_path (str): Path to output TSV file.
            verbose (bool): If True, provides user-facing logs in production.
        """
        self.seq = seq
        self.primer_name = primer_name
        self.template_seq = template_seq
        self.verbose = verbose
        self.codon_usage_dict = codon_usage_dict
        self.max_mutations = max_mutations
        self.output_tsv_path = output_tsv_path
        self.part_num_left = part_num_left
        self.part_num_right = part_num_right

        self.state = {
            'current_sequence_index': 0,
            'current_step': '',
            'mutations_found': [],
            'primers_designed': []
        }

        self.logger.debug(f"GoldenGateProtocol initialized with verbose={verbose}")
        if self.verbose:
            self.logger.info("GoldenGateProtocol is running in verbose mode.")

    def create_gg_protocol(self) -> Dict:
        """
        Main function to orchestrate the Golden Gate protocol creation.
        """
        logger.info("Starting Golden Gate protocol creation...")
        result_data = {
            'primers': [],
            'restriction_sites': [],
            'mutations': [],
            'sequence_analysis': [],
            'has_errors': False,
            'sequence_errors': {}
        }

        for i, single_seq in enumerate(self.seq):
            try:
                print(f"Processing sequence {i+1}/{len(self.seq)}")
                sequence_data = {
                    'sequence_index': i,
                    'processed_sequence': None,
                    'restriction_sites': [],
                    'mutations': None,
                    'primers': None
                }

                # 1. Remove start/stop codons and find restriction sites
                with debug_context("Preprocessing sequence"):
                    processed_seq, message, success = self.sequence_preparator.preprocess_sequence(
                        single_seq)
                    
                    if not success:
                        result_data['has_errors'] = True
                        result_data['sequence_errors'][i] = message
                        logger.error(f"Preprocessing failed for sequence {i}: {message}")
                        continue
                    
                with debug_context("Finding restriction sites"):  
                    sites_to_mutate = self.sequence_preparator.find_and_summarize_sites(
                        processed_seq, i)
                    sequence_data['processed_sequence'] = processed_seq
                    sequence_data['restriction_sites'] = sites_to_mutate
                        
                if sites_to_mutate:
                    with debug_context("Mutation analysis"):
                        # 2. Generate all possible silent mutations using MutationAnalyzer
                        mutation_options = self.mutation_analyzer.get_all_mutations(
                            sites_to_mutate=sites_to_mutate)

                        # 3. Optimize and prioritize mutations if there are valid options
                        if mutation_options:
                            optimized_mutations, compatibility_matrices = self.mutation_optimizer.optimize_mutations(
                                sequence=processed_seq, mutation_options=mutation_options
                            )
                            sequence_data['mutations'] = {
                                'all_mutation_options': optimized_mutations,
                                'compatibility': compatibility_matrices
                            }
                            
                        else:
                            optimized_mutations, compatibility_matrices = None, None

                    with debug_context("Primer design"):
                        # 4. Design primers using PrimerDesigner
                        primers = self.primer_designer.design_primers(
                            sequence=processed_seq,
                            seq_index=i,
                            mutations=optimized_mutations,
                            compatibility_matrices=compatibility_matrices,
                            template_seq=self.template_seq
                        )
                else:

                    logger.info(
                        f"No restriction sites found in sequence {i}")
                    primers = self.primer_designer.generate_GG_edge_primers(
                        processed_seq, i,
                    )
                    continue

                sequence_data['primers'] = primers
                result_data['sequence_analysis'].append(sequence_data)
                result_data['primers'].extend(primers)
            
            except Exception as e:
                logger.error(f"Unhandled error processing sequence {i}: {str(e)}")
                sequence_data['error'] = str(e)
                result_data['sequence_analysis'].append(sequence_data)
                result_data['has_errors'] = True
                result_data['sequence_errors'][i] = str(e)
                continue

        # Step 3: Save results and return
        if not result_data['has_errors']:
            self._save_primers_to_tsv(result_data['primers'], self.output_tsv_path)
                  
        return result_data

    def _save_primers_to_tsv(self, primer_data: List[List[str]], output_tsv_path: str) -> None:
        """Saves primer data to a TSV file."""
        with debug_context("save_primers_to_tsv"):
            if not primer_data:
                logger.warning("No primer data to save.")
                return

            try:
                with open(output_tsv_path, "w") as tsv_file:
                    tsv_file.write("Primer Name\tSequence\tAmplicon\n")
                    for row in primer_data:
                        tsv_file.write("\t".join(map(str, row)) + "\n")
            except IOError as e:
                logger.error(
                    f"Error writing to file {output_tsv_path}: {e}")
                raise
            
    def generate_GG_edge_primers(self, processed_seq: str, i: int,):
        return []
