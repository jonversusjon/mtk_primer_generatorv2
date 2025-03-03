# services/sequence_prep.py
from Bio.Seq import Seq, CodonTable
from prettytable import PrettyTable
import re
from typing import Dict, Optional, Union, List, Tuple

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
        

    def preprocess_sequence(self, sequence: Union[str, Seq]) -> Tuple[Optional[Seq], str, bool]:
        """
        Processes a DNA sequence by:
        1. Converting to uppercase
        2. Removing start codon (ATG) if present
        3. Removing stop codons (TAA, TAG, TGA) if present
        4. Ensuring the sequence is in frame (divisible by 3)
        
        Args:
            sequence: DNA sequence as string or Seq object
            
        Returns:
            Tuple of (processed_sequence, message, success)
            Where success is False if processing failed
        """
        with debug_context("pre_process_sequence"):
            # Convert to Seq object and uppercase if needed
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
            message_parts = []
            
            # Check for start codon
            if len(cleaned_sequence) >= 3 and cleaned_sequence[0:3] == "ATG":
                cleaned_sequence = cleaned_sequence[3:]
                trim_start_codon = True
                logger.info('Start codon removed from the beginning of the sequence')
            
            # Check for stop codons (TAA, TAG, TGA)
            stop_codons = ["TAA", "TAG", "TGA"]
            if len(cleaned_sequence) >= 3 and cleaned_sequence[-3:] in stop_codons:
                cleaned_sequence = cleaned_sequence[:-3]
                trim_stop_codon = True
                logger.info('Stop codon removed from the end of the sequence')
            
            # Handle frame adjustment based on codon presence
            final_length = len(cleaned_sequence)
            remainder = final_length % 3
            
            if remainder != 0:
                if trim_start_codon:
                    # If we removed start codon, trim from the end to make divisible by 3
                    bases_to_trim = remainder
                    cleaned_sequence = cleaned_sequence[:-bases_to_trim]
                    logger.info(f"Trimmed {bases_to_trim} nucleotides from end for frame adjustment")
                elif trim_stop_codon:
                    # If we removed stop codon, trim from the beginning to make divisible by 3
                    bases_to_trim = remainder
                    cleaned_sequence = cleaned_sequence[bases_to_trim:]
                    logger.info(f"Trimmed {bases_to_trim} nucleotides from start for frame adjustment")
                else:
                    # If no start/stop codon was found, we don't know how to infer frame
                    # This will be handled in the message logic
                    pass
            
            # Create human readable message
            if not in_frame:
                if trim_start_codon and trim_stop_codon:
                    message = "Provided sequence does not appear to be in frame, using provided start codon to infer translation frame. Stop and start codons detected and have been removed."
                elif trim_start_codon:
                    message = "Provided sequence does not appear to be in frame, using provided start codon to infer translation frame. Start codon has been removed."
                elif trim_stop_codon:
                    message = "Provided sequence does not appear to be in frame, using provided stop codon to infer frame. Stop codon has been removed."
                else:
                    # If no codons were found and sequence is not in frame, return error message
                    message = "Provided sequence does not appear to be in frame, please check your sequence and try again"
                    # Return None for sequence, the error message, and success=False
                    return None, message, False
            else:
                # Sequence is already in frame, we just report codon trimming
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
        """
        TODO: MAKE SURE THIS RETURNS THE DICTIONARY WITH ENZYME AS PRIMARY KEY
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
            seq_str = str(seq)
            
            # Simple string search first to confirm presence
            print(f"Sequence length: {len(seq_str)}")
            print(f"First 30 bases: {seq_str[:30]}")
            print(f"Last 30 bases: {seq_str[-30:]}")

            # Check if substring exists directly
            print(f"Manual check for reverse site: {'gagacc' in seq_str.lower()}")
            print(f"Position: {seq_str.lower().find('gagacc')}")
            # Direct string search instead of regex first
            fwd_positions = []
            pos = 0
            while True:
                pos = seq_str.find(recognition_seq, pos)
                if pos == -1:
                    break
                fwd_positions.append(pos)
                pos += 1
                
            rev_comp = str(Seq(recognition_seq).reverse_complement())
            print(f"Recognition sequence: '{recognition_seq}'")
            print(f"Reverse complement: '{rev_comp}'")
            rev_positions = []
            pos = 0
            while True:
                pos = seq_str.find(rev_comp, pos)
                if pos == -1:
                    break
                rev_positions.append(pos)
                pos += 1
                
            print(f"Direct string search for {recognition_seq}: Found at positions {fwd_positions}")
            print(f"Direct string search for {rev_comp}: Found at positions {rev_positions}")
            
            # Forward strand matches using regex
            forward_matches = list(re.finditer(re.escape(recognition_seq), seq_str))
            print(f"Regex search for {recognition_seq}: Found {len(forward_matches)} matches")
            
            for match in forward_matches:
                index = match.start()
                frame = (index) % 3  # Calculate frame (0, 1, or 2)
                codons = self.get_codons(seq, index, len(recognition_seq), frame)
                site_details.append({
                    'position': index+1,  # Convert to 1-based indexing
                    'sequence': recognition_seq,
                    'frame': frame,
                    'codons': codons,
                    'strand': '+'
                })

            # Reverse strand matches
            reverse_matches = list(re.finditer(re.escape(rev_comp), seq_str))
            print(f"Regex search for {rev_comp}: Found {len(reverse_matches)} matches")
            
            for match in reverse_matches:
                index = match.start()
                frame = (index) % 3  # Calculate frame (0, 1, or 2)
                codons = self.get_codons(seq, index, len(recognition_seq), frame)
                site_details.append({
                    'position': index+1,  # Convert to 1-based indexing
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
        
    def find_and_summarize_sites(self, sequence: Seq, index: int) -> List[Dict]:
        """Finds and summarizes restriction sites needing mutation."""
        with debug_context("find_and_summarize_sites"):
            sites_to_mutate = self.find_bsmbi_bsai_sites(sequence, self.verbose)
            if sites_to_mutate:
                self.summarize_bsmbi_bsai_sites(sites_to_mutate)
            return sites_to_mutate