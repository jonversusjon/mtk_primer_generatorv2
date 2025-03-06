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
    Handles sequence preprocessing for Golden Gate assembly,
    including frame adjustments and restriction site analysis.
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
                    return sequence, message, False
            else:
                if trim_start_codon and trim_stop_codon:
                    message = "Start and stop codons detected and removed."
                elif trim_start_codon:
                    message = "Start codon detected and removed."
                elif trim_stop_codon:
                    message = "Stop codon detected and removed."
                else:
                    message = "Sequence is in frame, no codon adjustments needed."

        return cleaned_sequence, message, True

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
                    codon_dict["position"] = codon_pos + 1
                    codons.append(codon_dict)

            return codons

    def _find_sites_for_enzyme(self, seq: Seq, recognition_seq: str) -> list:
        """Helper method to find sites for a specific enzyme."""
        site_details = []
        seq_str = str(seq)
        print(f"seq_str: {seq_str}")
        print(f"recognition_seq: {recognition_seq}")
        # Forward strand matches
        forward_matches = list(re.finditer(
            re.escape(recognition_seq), seq_str))

        print(f"forward_matches: {forward_matches}")
        for match in forward_matches:
            index = match.start()
            frame = index % 3
            codons = self.get_codons(seq, index, len(recognition_seq), frame)
            site_details.append({
                'position': index + 1,
                'sequence': recognition_seq,
                'frame': frame,
                'codons': codons,
                'strand': '+'
            })

        # Reverse strand matches
        rev_comp = str(Seq(recognition_seq).reverse_complement())
        reverse_matches = list(re.finditer(re.escape(rev_comp), seq_str))
        print(f"reverse_matches: {reverse_matches}")
        for match in reverse_matches:
            index = match.start()
            frame = index % 3
            codons = self.get_codons(seq, index, len(recognition_seq), frame)
            site_details.append({
                'position': index + 1,
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
                positions = ', '.join(
                    str(details['position']) for details in details_list)
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
