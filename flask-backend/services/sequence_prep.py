# services/sequence_prep.py
from Bio.Seq import Seq, CodonTable
from prettytable import PrettyTable
import re
from typing import Dict, Optional, Union, List, Tuple
from collections import defaultdict
from config.logging_config import logger
from .base import debug_context


class SequencePreparator:
    """
    Sequence Preparation Module

    This module handles preprocessing of DNA sequences for Golden Gate assembly,
    including the removal of start/stop codons, ensuring sequences are in the correct reading frame,
    and identifying restriction enzyme recognition sites (specifically BsmBI and BsaI) that require mutation.

    Input Parameters:
        - sequence (str or Bio.Seq.Seq): The DNA sequence to be prepared.
        - verbose (bool): Controls detailed logging during site identification.

    Output Data Structures:
        The methods provide structured outputs:
        - preprocess_sequence: Returns a tuple (cleaned_sequence (str), message (str), in_frame (bool)).
        - find_bsmbi_bsai_sites: Returns a sorted list of dictionaries, each containing:
            • 'position' (int): Start position of the recognition site.
            • 'sequence' (str): Recognition site sequence.
            • 'enzyme' (str): Enzyme name ('BsmBI' or 'BsaI').
            • 'frame' (int): Reading frame (0, 1, or 2).
            • 'strand' (str): Strand orientation ('+' or '-').
            • 'codons' (List[Dict]): Codons overlapping the recognition site with details:
                - 'codon_seq' (str): Codon sequence.
                - 'position' (int): Codon position within the sequence.
                - 'amino_acid' (str): Amino acid encoded by the codon.
            • 'context_sequence' (str): Sequence context surrounding the recognition site.
            • 'context_recognition_site_indices' (List[int]): Positions of the mutated bases within the context.

        - get_codons: Returns a list of codon dictionaries as described above.
    """

    def __init__(self, verbose: bool = False):
        self.logger = logger.getChild("SequencePreparator")
        self.verbose = verbose
        self.state = {
            'current_sequence': None,
            'adjustments_made': [],
            'restriction_sites': {}
        }

    def preprocess_sequence(self, sequence: Union[str, Seq]) -> Tuple[Optional[Seq], str, bool]:
        """
        Processes a DNA sequence by removing start/stop codons and ensuring proper frame.
        """
        with debug_context("pre_process_sequence"):
            # Convert to Seq object and uppercase
            if isinstance(sequence, str):
                sequence = Seq(sequence.upper())
            else:
                sequence = sequence.upper()

            cleaned_sequence = sequence
            sequence_length = len(sequence)

            # Initialize tracking variables
            trim_start_codon = False
            trim_stop_codon = False
            in_frame = sequence_length % 3 == 0

            # Check for start codon
            if len(cleaned_sequence) >= 3 and cleaned_sequence[0:3] == "ATG":
                cleaned_sequence = cleaned_sequence[3:]
                trim_start_codon = True
                logger.info(
                    'Start codon removed from the beginning of the sequence')

            # Check for stop codons
            stop_codons = ["TAA", "TAG", "TGA"]
            if len(cleaned_sequence) >= 3 and cleaned_sequence[-3:] in stop_codons:
                cleaned_sequence = cleaned_sequence[:-3]
                trim_stop_codon = True
                logger.info('Stop codon removed from the end of the sequence')

            # Handle frame adjustment
            final_length = len(cleaned_sequence)
            remainder = final_length % 3

            if remainder != 0:
                if trim_start_codon:
                    # If we removed start codon, trim from the end
                    cleaned_sequence = cleaned_sequence[:-remainder]
                elif trim_stop_codon:
                    # If we removed stop codon, trim from the beginning
                    cleaned_sequence = cleaned_sequence[remainder:]

            # Create message
            if not in_frame:
                if trim_start_codon and trim_stop_codon:
                    message = "Provided sequence does not appear to be in frame, using provided start codon to infer translation frame. Stop and start codons detected and have been removed."
                elif trim_start_codon:
                    message = "Provided sequence does not appear to be in frame, using provided start codon to infer translation frame. Start codon has been removed."
                elif trim_stop_codon:
                    message = "Provided sequence does not appear to be in frame, using provided stop codon to infer frame. Stop codon has been removed."
                else:
                    message = "Provided sequence does not appear to be in frame. If this is not intended, please check the sequence."
                    return str(sequence), message, False
            else:
                if trim_start_codon and trim_stop_codon:
                    message = "Start and stop codons detected and removed."
                elif trim_start_codon:
                    message = "Start codon detected and removed."
                elif trim_stop_codon:
                    message = "Stop codon detected and removed."
                else:
                    message = "Sequence is in frame, no codon adjustments needed."

        return str(cleaned_sequence), message, True

    def find_bsmbi_bsai_sites(self, sequence, verbose):
        """
        Finds BsmBI and BsaI restriction enzyme recognition sites within a DNA sequence.
        """
        seq = str(sequence)

        with debug_context("find_bsmbi_bsai_sites"):
            recognition_sequences = {
                'BsmBI': 'CGTCTC',
                'BsaI': 'GGTCTC'
            }
            seq_obj = Seq(seq.upper())

            sites_to_mutate = []
            for enzyme, recognition_seq in recognition_sequences.items():
                site_details = self._find_sites_for_enzyme(
                    seq_obj, recognition_seq)
                for site in site_details:
                    site['enzyme'] = enzyme
                    sites_to_mutate.append(site)

            sites_to_mutate.sort(key=lambda site: site['position'])
            return sites_to_mutate

    def get_codons(self, seq, start_index, length, frame):
        """
        Extracts codons spanned by the recognition site.
        """
        with debug_context("find_codons"):
            codons = []
            translation_table = CodonTable.unambiguous_dna_by_id[1]

            if frame == 0:
                codon_positions = [start_index, start_index + 3]
            elif frame == 1:
                codon_positions = [start_index - 1,
                                   start_index + 2, start_index + 5]
            elif frame == 2:
                codon_positions = [start_index - 2,
                                   start_index + 1, start_index + 4]
            else:
                return []

            # Extract codons
            for codon_pos in codon_positions:
                if 0 <= codon_pos <= len(seq) - 3:
                    codon_dict = {}
                    codon_seq = seq[codon_pos:codon_pos + 3]
                    codon_dict["codon_seq"] = str(codon_seq)
                    codon_dict["amino_acid"] = translation_table.forward_table.get(
                        str(codon_seq), 'X')
                    # Removed +1 for 0-indexing
                    codon_dict["position"] = codon_pos
                    codons.append(codon_dict)

            return codons

    def _find_sites_for_enzyme(self, seq: Seq, recognition_seq: str) -> list:
        """Helper method to find sites for a specific enzyme."""
        site_details = []
        seq_str = str(seq)

        forward_matches = list(re.finditer(
            re.escape(recognition_seq), seq_str))

        for match in forward_matches:
            index = match.start()
            frame = index % 3
            codons = self.get_codons(seq, index, len(recognition_seq), frame)

            # Get context sequence (30bp upstream and 30bp downstream)
            start_context = max(0, index - 30)
            end_context = min(len(seq_str), index + len(recognition_seq) + 30)
            context_seq = seq_str[start_context:end_context]

            # Calculate mutated base indices relative to the context sequence
            context_recognition_site_indices = [
                i - start_context for i in range(index, index + len(recognition_seq))]

            site_details.append({
                'position': index,  # Removed +1 for 0-indexing
                'recognition_sequence': recognition_seq,
                'frame': frame,
                'codons': codons,
                'strand': '+',
                'context_sequence': context_seq,
                'context_recognition_site_indices': context_recognition_site_indices
            })

        # Reverse strand matches
        rev_comp = str(Seq(recognition_seq).reverse_complement())
        reverse_matches = list(re.finditer(re.escape(rev_comp), seq_str))

        for match in reverse_matches:
            index = match.start()
            frame = index % 3
            codons = self.get_codons(seq, index, len(recognition_seq), frame)

            # Get context sequence (30bp upstream and 30bp downstream)
            start_context = max(0, index - 30)
            end_context = min(len(seq_str), index + len(recognition_seq) + 30)
            context_seq = seq_str[start_context:end_context]

            # Calculate mutated base indices relative to the context sequence
            context_recognition_site_indices = [
                i - start_context for i in range(index, index + len(rev_comp))]

            site_details.append({
                'position': index,  # Removed +1 for 0-indexing
                'sequence': rev_comp,
                'frame': frame,
                'codons': codons,
                'strand': '-',
                'context_sequence': context_seq,
                'context_recognition_site_indices': context_recognition_site_indices
            })

        return site_details

    def summarize_bsmbi_bsai_sites(self, site_details):
        """Creates a formatted summary of restriction sites."""
        with debug_context("summarize_restriction_sites"):
            site_type_descriptions = {
                'BsmBI': 'BsmBI Restriction Site',
                'BsaI': 'BsaI Restriction Site',
            }

            # Convert list to dictionary if needed
            if isinstance(site_details, list):
                grouped_sites = defaultdict(list)
                for site in site_details:
                    site_type = site.get('enzyme')
                    grouped_sites[site_type].append(site)
                site_details = grouped_sites

            table = PrettyTable()
            table.field_names = ["Site Type",
                                 "Number of Instances", "Position(s)"]

            for site_type, details_list in site_details.items():
                if not details_list:
                    continue

                site_type_description = site_type_descriptions.get(
                    site_type, site_type)
                # Add +1 to position for display only
                positions = ', '.join(
                    str(details['position'] + 1) for details in details_list)
                number_of_instances = len(details_list)

                table.add_row([
                    site_type_description,
                    number_of_instances,
                    positions
                ])

            logger.info("\nRestriction Site Analysis Summary:")
            logger.info(f"\n{str(table)}")

    def find_and_summarize_sites(self, sequence: Seq, index: int) -> List[Dict]:
        """Finds and summarizes restriction sites needing mutation."""
        with debug_context("find_and_summarize_sites"):
            sites_to_mutate = self.find_bsmbi_bsai_sites(
                sequence, self.verbose)
            if sites_to_mutate:
                self.summarize_bsmbi_bsai_sites(sites_to_mutate)
            return sites_to_mutate
