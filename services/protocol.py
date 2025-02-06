# services/protocol.py
from Bio.Seq import Seq
from typing import List, Dict, Optional
from .base import GoldenGateDesigner
from .sequence_prep import SequencePreparator
from .primer_design import PrimerDesigner
from .primer_select import PrimerSelector
from .mutation import MutationAnalyzer
from .utils import GoldenGateUtils


class GoldenGateProtocol(GoldenGateDesigner):
    def __init__(self, verbose: bool = False):
        super().__init__(verbose=verbose)
        self.utils = GoldenGateUtils()
        self.primer_designer = PrimerDesigner()
        self.primer_selector = PrimerSelector()
        self.mutation_analyzer = MutationAnalyzer()
        self.sequence_preparator = SequencePreparator()
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
        """Main function to orchestrate the Golden Gate protocol creation."""
        
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
        mutation_options = self._analyze_mutations(
            sequence, sites_to_mutate, codon_usage_dict, template_seq
        )

        self.state['current_step'] = 'primer_design'
        primer_data = self._design_primers(
            sequence=sequence,
            seq_index=seq_index,
            mutation_options=mutation_options,
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

    def _find_and_summarize_sites(self, sequence: Seq, index: int) -> Dict:
        """Finds and summarizes restriction sites needing mutation."""
        with self.debug_context("find_and_summarize_sites"):
            sites_to_mutate = self.sequence_preparator.find_bsmbi_bsai_sites(index, sequence, verbose=self.verbose)
            if sites_to_mutate:
                self.sequence_preparator.summarize_bsmbi_bsai_sites(sites_to_mutate)
            return sites_to_mutate

    def _analyze_mutations(
        self,
        sequence: Seq,
        sites_to_mutate: Dict,
        codon_usage_dict: Dict,
        template_seq: Optional[str]
    ) -> List:
        """Analyze and select mutations for the sequence."""
        with self.debug_context("analyze_mutations"):
            mutation_options = self.mutation_analyzer.gather_mutation_options(
                seq=sequence,
                sites_to_mutate=sites_to_mutate,
                codon_usage_dict=codon_usage_dict,
                spacer="GAA",
                bsmbi_site="CGTCTC",
                min_tm=57,
                template_seq=template_seq,
                verbose=self.verbose
            )
            
            selected_mutations = []
            while mutation_options and not selected_mutations:
                selected_mutations = self.mutation_analyzer.find_best_mutation_set(mutation_options, self.verbose)
                
            return selected_mutations

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

    def _save_primers_to_tsv(self, primer_data: List[List[str]], output_tsv_path: str):
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

# from Bio.Seq import Seq
# from typing import List, Dict, Optional
# from .sequence_prep import (
#     adjust_sequence_for_frame_and_codons,
#     find_bsmbi_bsai_sites,
#     summarize_bsmbi_bsai_sites
# )
# from .primer_design import generate_GG_edge_primers, get_all_possible_internal_primers
# from .primer_select import select_best_internal_primers, format_primers_for_output
# from .mutation import *
# from .utils import get_mtk_partend_sequences, get_codon_usage_dict


# def create_gg_protocol(
#     seq: List[str],
#     part_num_left: List[str],
#     part_num_right: List[str],
#     codon_usage_dict: Dict[str, Dict[str, float]],
#     max_mutations: int,
#     primer_name: Optional[List[str]] = None,
#     template_seq: Optional[str] = None,
#     kozak: str = "MTK",
#     verbose: bool = False,
#     output_tsv_path: str = "designed_primers.tsv",
# ) -> List[List[str]]:
#     """Main function to orchestrate the Golden Gate protocol creation."""

#     print("Starting Golden Gate protocol creation...")
#     primer_data = []
    
#     for i, single_seq in enumerate(seq):
#         single_seq = preprocess_sequence(single_seq)
#         sites_to_mutate = find_and_summarize_sites(single_seq, i, verbose)

#         num_sites_to_mutate = sum(len(sites)
#                                   for sites in sites_to_mutate.values())
#         print(f"Number of sites to mutate: {num_sites_to_mutate}")

#         # Step 1: Gather mutation options sorted by codon usage frequency
#         mutation_options = gather_mutation_options(
#             seq=single_seq,
#             sites_to_mutate=sites_to_mutate,
#             codon_usage_dict=codon_usage_dict,
#             spacer="GAA",
#             bsmbi_site="CGTCTC",
#             min_tm=57,
#             template_seq=template_seq,
#             verbose=verbose
#         )
#         # rank_and_print_mutation_sets(mutation_options, verbose)
#         print(f"\n🔍 mutation_options (before selection): {mutation_options}")
        
#         # Step 2: Iterate while mutation options exist
#         selected_mutations = []
#         while mutation_options and not selected_mutations:
#             selected_mutations = find_best_mutation_set(
#                 mutation_options, verbose)
#             print(f"✅ selected_mutations: {selected_mutations}")  # Debugging output

#         if not selected_mutations:
#             print("❌ No valid mutation subset found! Skipping primer design.")
            
#         # Step 3: Proceed with primer design
#         if selected_mutations:
#             primer_sets = get_all_possible_internal_primers(
#                 seq=single_seq,
#                 seq_index=i,
#                 proposed_mutations=selected_mutations,
#                 part_num_left=part_num_left[i],
#                 part_num_right=part_num_right[i],
#                 primer_name=primer_name[i],
#                 verbose=verbose
#             )
#             internal_primers = select_best_internal_primers(
#                 primer_sets, template_seq)
#             primer_data += format_primers_for_output(
#                 internal_primers, primer_name[i])

#         else:
#             print(
#                 f"❌ No valid mutation subset found for sequence {i + 1}. Skipping primer design.")

#        # Step 4: Design edge primers
#         mtk_partend_sequences = get_mtk_partend_sequences()
#         primer_data += get_edge_primers(
#             single_seq, i, mtk_partend_sequences, part_num_left, part_num_right, kozak, primer_name
#         )

#         # Save and return results
#         save_primers_to_tsv(primer_data, output_tsv_path)

#         return primer_data


# def preprocess_sequence(sequence: str) -> Seq:
#     """Converts sequence to uppercase and adjusts for frame and codons."""
#     sequence = Seq(sequence.upper())
#     return adjust_sequence_for_frame_and_codons(sequence)


# def find_and_summarize_sites(sequence: Seq, index: int, verbose: bool) -> Dict:
#     """Finds and summarizes restriction sites needing mutation."""
#     sites_to_mutate = find_bsmbi_bsai_sites(index, sequence, verbose=verbose)
#     if sites_to_mutate:
#         summarize_bsmbi_bsai_sites(sites_to_mutate)
#     return sites_to_mutate


# def get_edge_primers(
#     sequence: Seq,
#     index: int,
#     part_end_dict: Dict[str, str],
#     part_num_left: List[str],
#     part_num_right: List[str],
#     kozak: str,
#     primer_name: Optional[List[str]]
# ) -> List[List[str]]:
#     """Fetches edge primers for a given sequence and returns primer data."""
#     primer_data = []
#     left_name, left_primer = generate_GG_edge_primers(
#         sequence, part_end_dict, part_num_left[index], "forward", kozak, primer_name
#     )
#     right_name, right_primer = generate_GG_edge_primers(
#         sequence, part_end_dict, part_num_right[index], "reverse", kozak, primer_name
#     )
#     primer_data.append([left_name, left_primer, f"Amplicon_{index + 1}"])
#     primer_data.append([right_name, right_primer, f"Amplicon_{index + 1}"])

#     return primer_data


# def save_primers_to_tsv(primer_data: List[List[str]], output_tsv_path: str):
#     """Saves primer data to a TSV file."""
#     if not primer_data:
#         print("Warning: No primer data to save.")
#         return  # Exit early if no data

#     try:
#         with open(output_tsv_path, "w") as tsv_file:
#             tsv_file.write("Primer Name\tSequence\tAmplicon\n")
#             for row in primer_data:
#                 tsv_file.write("\t".join(map(str, row)) + "\n")
#     except IOError as e:
#         print(f"Error writing to file {output_tsv_path}: {e}")
