# services/sequence_prep.py
from Bio.Seq import Seq
from prettytable import PrettyTable
import re
from typing import Dict, Optional, Union, List
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
             
            self.logger.info("\n")    
            print("\nRestriction Site Analysis Summary:")
            print(str(table))
        
    def preprocess_sequence(self, sequence: str) -> Seq:
        """Converts sequence to uppercase and adjusts for frame and codons."""
        with self.debug_context("preprocess_sequence"):
            sequence = Seq(sequence.upper())
            return self.adjust_sequence_for_frame_and_codons(sequence)
        
    def find_and_summarize_sites(self, sequence: Seq, index: int) -> List[Dict]:
        """Finds and summarizes restriction sites needing mutation."""
        with self.debug_context("find_and_summarize_sites"):
            sites_to_mutate = self.find_bsmbi_bsai_sites(index, sequence, verbose=self.verbose)
            if sites_to_mutate:
                self.summarize_bsmbi_bsai_sites(sites_to_mutate)
            return sites_to_mutate