# services/sequence_prep.py
from Bio.Seq import Seq
from prettytable import PrettyTable
import re
from typing import Dict, Optional, Union
from .base import GoldenGateDesigner


class SequencePreparator(GoldenGateDesigner):
    def __init__(self, verbose: bool = False):
        super().__init__(verbose=verbose)
        self.state = {
            'current_sequence': None,
            'adjustments_made': [],
            'restriction_sites': {}
        }

    def adjust_sequence_for_frame_and_codons(
        self, 
        seq: Union[str, Seq],
    ) -> Seq:
        """Adjusts sequence to be in frame and removes start/stop codons."""
        with self.debug_context("adjust_sequence_for_frame_and_codons"):
            if isinstance(seq, str):
                seq = Seq(seq)
            
            self.state['current_sequence'] = seq
            self.state['adjustments_made'] = []

            # Adjust for frame
            original_length = len(seq)
            seq = seq[:len(seq) - (len(seq) % 3)]
            if len(seq) != original_length:
                self.state['adjustments_made'].append('frame_adjusted')
                self.logger.info(f"Trimmed {original_length - len(seq)} nucleotides for frame adjustment")

            # Check codons
            translated_seq = seq.translate()

            # Handle stop codon
            if translated_seq[-1] == '*':
                self.state['adjustments_made'].append('stop_codon_removed')
                self.logger.info('Stop codon removed from the end of the sequence')
                seq = seq[:-3]

            # Handle start codon
            if translated_seq[0] == 'M':
                self.state['adjustments_made'].append('start_codon_removed')
                self.logger.info('Start codon removed from the beginning of the sequence')
                seq = seq[3:]

            return seq

    def find_bsmbi_bsai_sites(
        self,
        sequence_index: int,
        seq: Union[str, Seq],
        verbose: bool = False,
    ) -> Dict:
        """Identifies restriction sites within a DNA sequence."""
        with self.debug_context("find_bsmbi_bsai_sites"):
            recognition_sequences = {
                'BsmBI': 'CGTCTC',
                'BsaI': 'GGTCTC'
            }

            found_sites = {}
            seq_obj = Seq(str(seq).upper())
            
            for enzyme, recognition_seq in recognition_sequences.items():
                site_details = self._find_sites_for_enzyme(
                    seq_obj, 
                    recognition_seq
                )
                found_sites[enzyme] = site_details

            self.state['restriction_sites'] = found_sites
            total_sites = sum(len(sites) for sites in found_sites.values())
            
            self.logger.info(
                f"Found {total_sites} restriction sites in sequence {sequence_index + 1}"
            )
            
            return found_sites

    def _find_sites_for_enzyme(
        self,
        seq: Seq,
        recognition_seq: str
    ) -> list:
        """Helper method to find sites for a specific enzyme."""
        site_details = []
        
        # Forward strand matches
        forward_matches = re.finditer(f"(?={recognition_seq})", str(seq))
        for match in forward_matches:
            site_details.append({
                'position': match.start(),
                'sequence': recognition_seq,
                'strand': '+'
            })

        # Reverse strand matches
        rev_comp = str(Seq(recognition_seq).reverse_complement())
        reverse_matches = re.finditer(f"(?={rev_comp})", str(seq))
        for match in reverse_matches:
            site_details.append({
                'position': match.start(),
                'sequence': rev_comp,
                'strand': '-'
            })

        return site_details

    def summarize_bsmbi_bsai_sites(self, site_details: Dict):
        """Creates a formatted summary of restriction sites."""
        with self.debug_context("summarize_restriction_sites"):
            site_type_descriptions = {
                'BsmBI': 'BsmBI Restriction Site',
                'BsaI': 'BsaI Restriction Site',
            }

            table = PrettyTable()
            table.field_names = ["Site Type", "Number of Instances", "Position(s)"]

            for site_type, details_list in site_details.items():
                if not details_list:
                    continue
                    
                site_type_description = site_type_descriptions.get(site_type, site_type)
                positions = ', '.join(str(details['position']) for details in details_list)
                number_of_instances = len(details_list)
                
                table.add_row([
                    site_type_description,
                    number_of_instances,
                    positions
                ])

            self.logger.info("\nRestriction Site Analysis Summary:")
            self.logger.info(str(table))

# import re
# from Bio.Seq import Seq
# from prettytable import PrettyTable


# def adjust_sequence_for_frame_and_codons(seq, output_widget=None):
#     """
#     Adjusts the sequence to be in frame by trimming excess nucleotides and removes start/stop codons if present.

#     Args:
#         seq (Bio.Seq.Seq): The DNA sequence to adjust.
#         output_widget (ipywidgets.Output, optional): Output widget for displaying messages.

#     Returns:
#         Bio.Seq.Seq: The adjusted DNA sequence.
#     """
#     seq = seq[:len(seq) - (len(seq) % 3)]  # Adjust for frame

#     # Translate the sequence for codon checks
#     translated_seq = seq.translate()

#     # Remove stop codon at the end if present
#     if translated_seq[-1] == '*':
#         message = 'Stop codon removed from the end of the sequence.'
#         if output_widget:
#             with output_widget:
#                 print(message)
#         else:
#             print(message)
#         seq = seq[:-3]

#     # Remove start codon at the beginning if present
#     if translated_seq[0] == 'M':
#         message = 'Start codon removed from the beginning of the sequence.'
#         if output_widget:
#             with output_widget:
#                 print(message)
#         else:
#             print(message)
#         seq = seq[3:]

#     return seq



# def find_bsmbi_bsai_sites(i, seq, verbose=False):
#     """
#     Identifies the indices of specified restriction sites within a given DNA sequence.

#     Args:
#         i (int): The sequence number from the user input (used for logging/output).
#         seq (str or Bio.Seq.Seq): The DNA sequence to search.
#         verbose (bool): Whether to print verbose output.

#     Returns:
#         dict: A dictionary mapping restriction site names to lists of dictionaries,
#               where each dictionary represents a unique site and contains its position,
#               sequence, and strand.
#     """
#     recognition_sequences = {
#         'BsmBI': 'CGTCTC',
#         'BsaI': 'GGTCTC'
#     }

#     found_sites = {}
#     seq_obj = Seq(seq.upper())  # Ensure seq is a Bio.Seq.Seq object

#     for enzyme, recognition_seq in recognition_sequences.items():
#         forward_matches = re.finditer(f"(?={recognition_seq})", str(seq_obj))
#         reverse_matches = re.finditer(f"(?={str(Seq(recognition_seq).reverse_complement())})", str(seq_obj))

#         site_details = []
#         for match in forward_matches:
#             site_details.append({
#                 'position': match.start()+1,
#                 'sequence': recognition_seq,
#                 'strand': '+'
#             })
#         for match in reverse_matches:
#             site_details.append({
#                 'position': match.start()+1,
#                 'sequence': str(Seq(recognition_seq).reverse_complement()),
#                 'strand': '-'
#             })

#         found_sites[enzyme] = site_details

#     # Summarize and log
#     found_sites_count = sum(len(lst) for lst in found_sites.values())
#     if verbose:
#         if found_sites_count > 0:
#             print(f'Sequence {i+1} BsmBI/BsaI sites found: {found_sites}')
#         else:
#             print(f'Found no BsmBI/BsaI sites for sequence {i}')

#     return found_sites


# def summarize_bsmbi_bsai_sites(site_details):
#     """
#     Summarizes restriction site details in a table format.

#     Args:
#         site_details (dict): Dictionary containing restriction site details.

#     Returns:
#         None
#     """
#     print('\n=====================================================')
#     print('Restriction Site Analysis Summary:')
#     print('=====================================================\n')

#     table = PrettyTable()
#     table.field_names = ["Site Type", "Number of Instances", "Position(s)"]

#     # Define a mapping from site type keys to descriptive names
#     site_type_descriptions = {
#         'BsmBI': 'BsmBI Restriction Site',
#         'BsaI': 'BsaI Restriction Site',
#     }

#     for site_type, details_list in site_details.items():
#         if not details_list:  # Skip if the list is empty
#             continue
#         # Use the mapping to convert site_type to a more meaningful description
#         site_type_description = site_type_descriptions.get(site_type, site_type)
#         positions = ', '.join(str(details['position']) for details in details_list)
#         number_of_instances = len(details_list)
#         table.add_row([site_type_description, number_of_instances, positions])

#     print(table)
