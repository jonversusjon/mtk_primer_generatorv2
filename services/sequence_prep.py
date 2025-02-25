# services/sequence_prep.py
from Bio.Seq import Seq, CodonTable
from prettytable import PrettyTable
import re
from typing import Dict, Optional, Union, List
from collections import defaultdict
import logging
from config.logging_config import logger
from .base import debug_context

class SequencePreparator:
    """
    Handles sequence preprocessing for Golden Gate assembly,
    including frame adjustments and restriction site analysis.
    """

    def __init__(self, verbose: bool = False):
        self.logger = logger.getChild("SequencePreparator")
        
        """
        Initialize sequence preparation parameters.
        
        Args:
            verbose (bool): If True, provide user-facing logs in production.
        """
        self.verbose = verbose

        self.state = {
            'current_sequence': None,
            'adjustments_made': [],
            'restriction_sites': {}
        }

        self.logger.debug("SequencePreparator initialized.")

        if self.verbose:
            self.logger.info("SequencePreparator is running in verbose mode.")

    def adjust_sequence_for_frame_and_codons(
        self, 
        seq: Union[str, Seq],
    ) -> Seq:
        """Adjusts sequence to be in frame and removes start/stop codons."""
        with debug_context("adjust_sequence_for_frame_and_codons"):
            if isinstance(seq, str):
                seq = Seq(seq)
            print(f"seq: {seq}")
            self.state['current_sequence'] = seq
            self.state['adjustments_made'] = []

            # Adjust for frame
            original_length = len(seq)
            seq = seq[:len(seq) - (len(seq) % 3)]
            if len(seq) != original_length:
                self.state['adjustments_made'].append('frame_adjusted')
                logger.info(f"Trimmed {original_length - len(seq)} nucleotides for frame adjustment")

            # Check codons
            translated_seq = seq.translate()

            # Handle stop codon
            if translated_seq[-1] == '*':
                self.state['adjustments_made'].append('stop_codon_removed')
                logger.info('Stop codon removed from the end of the sequence')
                seq = seq[:-3]

            # Handle start codon
            if translated_seq[0] == 'M':
                self.state['adjustments_made'].append('start_codon_removed')
                logger.info('Start codon removed from the beginning of the sequence')
                seq = seq[3:]

            return seq
        

    def find_bsmbi_bsai_sites(self, sequence, verbose):
        """
        Finds restriction enzyme recognition sites within a DNA sequence,
        and records the codons spanned by the site.

        Args:
            sequence (str or Bio.Seq): The DNA sequence.
            recognition_seq (str): The recognition sequence of the restriction enzyme.

        Returns:
            A list of dictionaries, where each dictionary represents a
            recognition site and contains:
                'position':  The 1-based start position of the site in the sequence.
                'sequence':  The recognized sequence (may be reverse complement).
                'frame':     The reading frame (0, 1, or 2).
                'strand':    '+' for forward strand, '-' for reverse strand.
                'codons':    A list of codons spanned by the recognition site.
                            Each codon is a tuple: (codon_sequence, codon_position, amino_acid).
                'enzyme':    The name of the enzyme (e.g., 'BsmBI' or 'BsaI').
        """
        
        seq = str(sequence)  # Ensure seq is a string for consistent indexing

        with debug_context("find_bsmbi_bsai_sites"):
            recognition_sequences = {
                'BsmBI': 'CGTCTC',
                'BsaI': 'GGTCTC'
            }
            seq_obj = Seq(seq.upper())
            
            sites_to_mutate = []
            for enzyme, recognition_seq in recognition_sequences.items():
                site_details = self._find_sites_for_enzyme(seq_obj, recognition_seq)
                for site in site_details:
                    site['enzyme'] = enzyme
                    sites_to_mutate.append(site)
                    
            sites_to_mutate.sort(key=lambda site: site['position'])

            return sites_to_mutate


    def get_codons(self, seq, start_index, length, frame):
        """
        Extracts codons spanned by the recognition site, handling edge cases
        and strand orientation.

        Args:
            seq (Seq): The full template sequence being domesticated.
            start_index (int): The 0-based start index of the recognition site.
            length (int): The length of the recognition site sequence.
            frame (int): The frame offset (0, 1, or 2).

        Returns:
            List of tuples (codon_sequence, codon_position, amino_acid)
        """
        with debug_context("find_codons"):
            codons = []
            translation_table = CodonTable.unambiguous_dna_by_id[1]  # Standard genetic code

            if frame == 0:
                codon_positions = [start_index, start_index + 3]
            elif frame == 1:
                codon_positions = [start_index - 1, start_index + 2, start_index + 5]
            elif frame == 2:
                codon_positions = [start_index - 2, start_index + 1, start_index + 4]
            else:
                raise ValueError("Frame must be 0, 1, or 2")

            # Extract codons, ensuring we stay within bounds
            for codon_pos in codon_positions:
                if 0 <= codon_pos <= len(seq) - 3:  # Ensure valid triplet extraction
                    codon_dict = {}
                    codon_seq = seq[codon_pos:codon_pos + 3]
                    codon_dict["codon_seq"] = str(codon_seq)
                    codon_dict["amino_acid"] = translation_table.forward_table.get(str(codon_seq), 'X')  # 'X' for unknown/stop
                    codon_dict["position"] = codon_pos + 1
                    codons.append(codon_dict)  # Store 1-based codon position
            
            if self.verbose:
                logger.info(
                    f"Extracted codons for sequence of length {len(seq)} with start_index {start_index}, "
                    f"length {length}, and frame {frame}: {codons}"
                )
            return codons


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
            index = match.start()
            frame = (index) % 3  # Calculate frame (0, 1, or 2)
            codons = self.get_codons(seq, index, len(recognition_seq), frame)
            site_details.append({
                'position': index+1,
                'sequence': recognition_seq,
                'frame': frame,
                'codons': codons,
                'strand': '+'
            })

        # Reverse strand matches
        rev_comp = str(Seq(recognition_seq).reverse_complement())
        reverse_matches = re.finditer(f"(?={rev_comp})", str(seq))
        for match in reverse_matches:
            index = match.start()
            frame = (index) % 3  # Calculate frame (0, 1, or 2)
            codons = self.get_codons(seq, index, len(recognition_seq), frame)
            site_details.append({
                'position': index+1,
                'sequence': rev_comp,
                'frame': frame,
                'codons': codons,
                'strand': '-'
            })

        return site_details


    def summarize_bsmbi_bsai_sites(self, site_details):
        """Creates a formatted summary of restriction sites."""
        with debug_context("summarize_restriction_sites"):
            site_type_descriptions = {
                'BsmBI': 'BsmBI Restriction Site',
                'BsaI': 'BsaI Restriction Site',
            }

            # Ensure site_details is structured as a dictionary
            if isinstance(site_details, list):
                grouped_sites = defaultdict(list)
                for site in site_details:
                    # Infer site type from sequence (if applicable) or provide a key in `site`
                    site_type = "BsmBI" if site["sequence"] == "GAGACC" else "BsaI"
                    grouped_sites[site_type].append(site)
                site_details = grouped_sites

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

            logger.info("\n")
            logger.info("\nRestriction Site Analysis Summary:")
            logger.info(f"\n{str(table)}")

    def preprocess_sequence(self, sequence: str) -> Seq:
        """Converts sequence to uppercase and adjusts for frame and codons."""
        with debug_context("preprocess_sequence"):
            sequence = Seq(sequence.upper())
            return self.adjust_sequence_for_frame_and_codons(sequence)
        
    def find_and_summarize_sites(self, sequence: Seq, index: int) -> List[Dict]:
        """Finds and summarizes restriction sites needing mutation."""
        with debug_context("find_and_summarize_sites"):
            sites_to_mutate = self.find_bsmbi_bsai_sites(sequence, self.verbose)
            if sites_to_mutate:
                self.summarize_bsmbi_bsai_sites(sites_to_mutate)
            return sites_to_mutate