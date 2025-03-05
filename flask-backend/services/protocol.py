# services/protocol.py
import json
from typing import List, Dict, Optional
from .sequence_prep import SequencePreparator
from .primer_design import PrimerDesigner
from .primer_select import PrimerSelector
from .mutation_analyzer import MutationAnalyzer
from .mutation_optimizer import MutationOptimizer
from .utils import GoldenGateUtils
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
        self.logger.debug(
            f"GoldenGateProtocol initialized with codon_usage_dict: {codon_usage_dict}")
        if verbose:
            self.logger.info("GoldenGateProtocol is running in verbose mode.")
        self.mutation_analyzer = MutationAnalyzer(
            sequence=seq,
            codon_usage_dict=codon_usage_dict,
            max_mutations=max_mutations,
            verbose=verbose
        )

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

    def create_gg_protocol(self) -> dict:
        """
        Main function to orchestrate the Golden Gate protocol creation.
        Returns:
            dict: A dictionary containing protocol details.
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
            sequence_data = {
                'sequence_index': i,
                'processed_sequence': None,
                'restriction_sites': [],
                'mutations': None,
                'primers': None,
                'preprocessing_message': None
            }

            try:
                print(f"Processing sequence {i+1}/{len(self.seq)}")

                # 1. Preprocess sequence (remove start/stop codons, etc.)
                with debug_context("Preprocessing sequence"):
                    processed_seq, message, success = self.sequence_preparator.preprocess_sequence(
                        single_seq)
                    if message:
                        sequence_data['preprocessing_message'] = message
                    # If preprocessing returns None, fallback to the original sequence.
                    if processed_seq is None:
                        processed_seq = single_seq
                    # Always convert Seq object to string when storing in result
                    sequence_data['processed_sequence'] = str(
                        processed_seq) if processed_seq is not None else None

                # 2. Find restriction sites
                with debug_context("Finding restriction sites"):
                    sites_to_mutate = self.sequence_preparator.find_and_summarize_sites(
                        processed_seq, i)
                    sequence_data['restriction_sites'] = sites_to_mutate

                # 3. Mutation analysis and primer design
                if sites_to_mutate:
                    with debug_context("Mutation analysis"):
                        mutation_options = self.mutation_analyzer.get_all_mutations(
                            sites_to_mutate=sites_to_mutate
                        )
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
                        primers = self.primer_designer.design_primers(
                            sequence=processed_seq,
                            seq_index=i,
                            mutations=optimized_mutations,
                            compatibility_matrices=compatibility_matrices,
                            template_seq=self.template_seq
                        )
                else:
                    logger.info(f"No restriction sites found in sequence {i}")
                    try:
                        overhang_5_prime = self.part_num_left[i]
                        overhang_3_prime = self.part_num_right[i]
                        primers = self.primer_designer.generate_GG_edge_primers(
                            i, processed_seq, overhang_5_prime, overhang_3_prime
                        )
                    except Exception as e:
                        logger.error(
                            f"Error generating edge primers for sequence {i}: {e}")
                        primers = []

                sequence_data['primers'] = primers
                result_data['sequence_analysis'].append(sequence_data)
                result_data['primers'].extend(primers)

            except Exception as e:
                logger.exception(
                    f"Unhandled error processing sequence {i}: {str(e)}")
                sequence_data['error'] = str(e)
                result_data['sequence_analysis'].append(sequence_data)
                result_data['has_errors'] = True
                result_data['sequence_errors'][i] = str(e)
                continue

        # Save primers to TSV if no critical errors occurred
        if not result_data['has_errors']:
            self._save_primers_to_tsv(
                result_data['primers'], self.output_tsv_path)

        # Convert any remaining non-serializable objects and return the dictionary
        return self.utils.convert_non_serializable(result_data)

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
                logger.error(f"Error writing to file {output_tsv_path}: {e}")
                raise

    def generate_GG_edge_primers(self, processed_seq: str, i: int):
        return []
