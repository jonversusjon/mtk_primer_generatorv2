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
    def __init__(self, verbose: bool = False):
        super().__init__(verbose=verbose)
        self.utils = GoldenGateUtils()
        self.sequence_preparator = SequencePreparator()
        self.primer_designer = PrimerDesigner()
        self.primer_selector = PrimerSelector()
        self.mutation_analyzer = MutationAnalyzer()
        self.mutation_optimizer = MutationOptimizer(verbose=verbose)
        
        self.state = {
            'current_sequence_index': 0,
            'current_step': '',
            'mutations_found': [],
            'primers_designed': []
        }

    def create_gg_protocol(
        self,
        seq: List[str],
        part_num_left: List[str],
        part_num_right: List[str],
        codon_usage_dict: Dict[str, Dict[str, float]],
        max_mutations: int,
        primer_name: Optional[List[str]] = None,
        template_seq: Optional[str] = None,
        kozak: str = "MTK",
        output_tsv_path: str = "designed_primers.tsv",
        verbose: bool = False,
    ) -> List[List[str]]:
        """
        Main function to orchestrate the Golden Gate protocol creation.
        
        Args:
            seq: List of DNA sequences to process
            part_num_left: List of left part numbers
            part_num_right: List of right part numbers
            codon_usage_dict: Dictionary mapping codons to their usage frequencies
            max_mutations: Maximum number of mutations allowed
            primer_name: Optional list of primer names
            template_seq: Optional template sequence for comparison
            kozak: Kozak sequence (default: "MTK")
            output_tsv_path: Path to save output TSV file
            verbose: Whether to print detailed progress
            
        Returns:
            List of lists containing primer data [name, sequence, amplicon]
        """
        self._validate_inputs(seq, part_num_left, part_num_right, primer_name)
        
        with self.debug_context("create_gg_protocol"):
            self.logger.info("Starting Golden Gate protocol creation...")
            primer_data = []
            
            for i, single_seq in enumerate(seq):
                self.state['current_sequence_index'] = i
                self.state['current_step'] = 'preprocessing'
                
                try:
                    primer_data.extend(
                        self._process_single_sequence(
                            single_seq=single_seq,
                            seq_index=i,
                            part_num_left=part_num_left,
                            part_num_right=part_num_right,
                            codon_usage_dict=codon_usage_dict,
                            max_mutations=max_mutations,
                            primer_name=primer_name,
                            template_seq=template_seq,
                            kozak=kozak
                        )
                    )
                except Exception as e:
                    self.logger.error(f"Error processing sequence {i}: {str(e)}")
                    continue

            self._save_primers_to_tsv(primer_data, output_tsv_path)
            return primer_data

    def _validate_inputs(self, seq: List[str], part_num_left: List[str], 
                        part_num_right: List[str], primer_name: Optional[List[str]]) -> None:
        """Validate input parameters for consistency."""
        if len(part_num_left) != len(seq):
            raise ValueError("Length of part_num_left must match number of sequences")
        if len(part_num_right) != len(seq):
            raise ValueError("Length of part_num_right must match number of sequences")
        if primer_name and len(primer_name) != len(seq):
            raise ValueError("Length of primer_name must match number of sequences")

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
        
        sequence = self._preprocess_sequence(single_seq)
        sites_to_mutate = self._find_and_summarize_sites(sequence, seq_index)

        self.state['current_step'] = 'mutation_analysis'
        mutation_result = self._analyze_mutations(
            sequence=sequence,
            sites_to_mutate=sites_to_mutate,
            codon_usage_dict=codon_usage_dict,
        )
        print(f"mutation_result final: {mutation_result}")
        
        
        if mutation_result['status'] == 'error':
            self.logger.warning(f"Mutation analysis failed: {mutation_result['message']}")
            return []
            
        selected_mutations = mutation_result['selected_mutations']
        
        self.state['current_step'] = 'primer_design'
        primer_data = self._design_primers(
            sequence=sequence,
            seq_index=seq_index,
            mutation_options=selected_mutations,
            part_num_left=part_num_left,
            part_num_right=part_num_right,
            primer_name=primer_name,
            template_seq=template_seq,
            kozak=kozak
        )

        return primer_data

    def _preprocess_sequence(self, sequence: str) -> Seq:
        """Converts sequence to uppercase and adjusts for frame and codons."""
        with self.debug_context("preprocess_sequence"):
            sequence = Seq(sequence.upper())
            return self.sequence_preparator.adjust_sequence_for_frame_and_codons(sequence)

    def _find_and_summarize_sites(self, sequence: Seq, index: int) -> List[Dict]:
        """Finds and summarizes restriction sites needing mutation."""
        with self.debug_context("find_and_summarize_sites"):
            sites_to_mutate = self.sequence_preparator.find_bsmbi_bsai_sites(index, sequence, verbose=self.verbose)
            if sites_to_mutate:
                self.sequence_preparator.summarize_bsmbi_bsai_sites(sites_to_mutate)
            return sites_to_mutate

    def _analyze_mutations(
        self,
        sequence: Union[Seq, str],
        sites_to_mutate: List[Dict],
        codon_usage_dict: Dict[str, Dict[str, float]],
        max_mutations: int = 1,
    ) -> Dict:
        """
        Analyze mutations for the sequence and evaluate their compatibility.
        
        Args:
            sequence: Input DNA sequence
            sites_to_mutate: List of restriction sites to mutate
            codon_usage_dict: Dictionary of codon usage frequencies
            compatibility_table_path: Path to binary compatibility lookup table
        
        Returns:
            Dictionary containing:
            - status: 'success' or 'error'
            - message: Description of result
            - mutation_options: All possible mutations
            - selected_mutations: Optimal compatible mutations
            - compatibility_score: Score indicating overall compatibility
            - cut_positions: List of selected cut positions
        """
        with self.debug_context("analyze_mutations"):
            try:
                # Convert Seq to str if necessary
                sequence_str = str(sequence) if isinstance(sequence, Seq) else sequence
                
                # Get all possible mutation options
                mutation_result = self.mutation_analyzer.gather_mutation_options(
                    seq=sequence_str,
                    sites_to_mutate=sites_to_mutate,
                    codon_usage_dict=codon_usage_dict,
                    max_mutations=max_mutations,
                    verbose=self.verbose
                )
                
                if not mutation_result or 'restriction_sites' not in mutation_result:
                    return {
                        'status': 'error',
                        'message': 'No mutation options generated',
                        'mutation_options': [],
                        'selected_mutations': None,
                    }

                restriction_sites = mutation_result['restriction_sites']
                
                # Process mutation options for each site to get valid ranges
                site_mutations = []
                for site in restriction_sites:
                    site_muts = []
                    for codon in site['codons']:
                        if 'mutations' in codon:
                            for mutation in codon['mutations']:
                                mut_copy = mutation.copy()
                                mut_copy.update({
                                    'site_start': site['start_index'],
                                    'enzyme': site['enzyme'],
                                    'strand': site['strand']
                                })
                                site_muts.append(mut_copy)
                    if site_muts:
                        site_mutations.append(site_muts)
                
                if not site_mutations:
                    return {
                        'status': 'error',
                        'message': 'No valid mutations found',
                        'mutation_options': mutation_result,
                    }
                
                # Calculate cut ranges and create compatibility tensor
                compatibility_tensor, cut_ranges = self.mutation_optimizer.create_compatibility_tensor(
                    sequence=sequence_str,
                    mutation_options=site_mutations
                )
                
                print(f'cut_ranges: {cut_ranges}')
                # if not selected_indices:
                #     return {
                #         'status': 'error',
                #         'message': 'No compatible mutation combination found',
                #         'mutation_options': mutation_result,
                #         'selected_mutations': None,
                #         'compatibility_score': 0
                #     }
                
                # return {
                #     'status': 'success',
                #     'message': 'Compatible mutations found',
                #     'mutation_options': mutation_result,
                #     'selected_mutations': selected_mutations,
                #     'compatibility_score': 1,
                #     'cut_positions': list(selected_indices)
                # }
                   
                return {}
             
            except Exception as e:
                if self.verbose:
                    print(f"Error in mutation analysis: {str(e)}")
                return {
                    'status': 'error',
                    'message': f'Error during mutation analysis: {str(e)}',
                    'mutation_options': mutation_result if 'mutation_result' in locals() else {},
                    'selected_mutations': None,
                    'compatibility_score': 0
                }
            
    def _design_primers(
        self,
        sequence: Seq,
        seq_index: int,
        mutation_options: List,
        part_num_left: List[str],
        part_num_right: List[str],
        primer_name: Optional[List[str]],
        template_seq: Optional[str],
        kozak: str
    ) -> List[List[str]]:
        """Design primers based on selected mutations."""
        with self.debug_context("design_primers"):
            primer_data = []
            
            if mutation_options:
                primer_sets = self.primer_designer.get_all_possible_internal_primers(
                    seq=sequence,
                    seq_index=seq_index,
                    proposed_mutations=mutation_options,
                    part_num_left=part_num_left[seq_index],
                    part_num_right=part_num_right[seq_index],
                    primer_name=primer_name[seq_index] if primer_name else None,
                    verbose=self.verbose
                )
                internal_primers = self.primer_selector.select_best_internal_primers(primer_sets, template_seq)
                primer_data.extend(self.primer_selector.format_primers_for_output(
                    internal_primers,
                    primer_name[seq_index] if primer_name else None
                ))

            # Add edge primers
            mtk_partend_sequences = self.utils.get_mtk_partend_sequences()
            primer_data.extend(self._get_edge_primers(
                sequence, seq_index, mtk_partend_sequences,
                part_num_left, part_num_right, kozak, primer_name
            ))

            return primer_data

    def _get_edge_primers(
        self,
        sequence: Seq,
        index: int,
        part_end_dict: Dict[str, str],
        part_num_left: List[str],
        part_num_right: List[str],
        kozak: str,
        primer_name: Optional[List[str]]
    ) -> List[List[str]]:
        """Generate edge primers for a sequence."""
        with self.debug_context("get_edge_primers"):
            primer_data = []
            left_name, left_primer = self.primer_designer.generate_GG_edge_primers(
                sequence, part_end_dict, part_num_left[index], "forward", kozak, primer_name
            )
            right_name, right_primer = self.primer_designer.generate_GG_edge_primers(
                sequence, part_end_dict, part_num_right[index], "reverse", kozak, primer_name
            )
            primer_data.append([left_name, left_primer, f"Amplicon_{index + 1}"])
            primer_data.append([right_name, right_primer, f"Amplicon_{index + 1}"])
            return primer_data

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