# services/sequence_prep.py
from Bio.Seq import Seq
from typing import Optional, Tuple
from log_utils import logger
import logging

class SequencePreparator:
    """
    Sequence Preparation Module

    This module handles preprocessing of DNA sequences for Golden Gate assembly,
    including the removal of start/stop codons and ensuring sequences are in the correct reading frame.
    """

    def __init__(self, verbose: bool = False, debug: bool = False):
        self.verbose = verbose
        self.debug = debug
        self.state = {
            'current_sequence': None,
            'adjustments_made': [],
            'restriction_sites': {}
        }
        if self.debug:
            logger.log_step("Initialization", "Debug mode enabled for SequencePreparator")

    @logger.log_function
    def preprocess_sequence(
        self,
        sequence: str,
        matk_part_left: str
    ) -> Tuple[Optional[Seq], str, bool]:
        """
        Processes a DNA sequence by removing start/stop codons and ensuring proper frame.
        """
        with logger.debug_context("pre_process_sequence"):
            logger.log_step("Conversion", "Converting input sequence to uppercase and Seq object")
            if isinstance(sequence, str):
                sequence = Seq(sequence.upper())
            else:
                sequence = sequence.upper()

            cleaned_sequence = sequence
            sequence_length = len(sequence)
            logger.log_step("Conversion", f"Sequence converted. Original length: {sequence_length}")


            # Initialize tracking variables
            trim_start_codon = False
            trim_stop_codon = False
            in_frame = sequence_length % 3 == 0
            logger.log_step("Validation", f"Initial frame check: in_frame = {in_frame}")

            # Check for start codon only if matk_part_left is "3" or "3a"
            if matk_part_left in {"3", "3a"} and len(cleaned_sequence) >= 3 and cleaned_sequence[:3] == "ATG":
                cleaned_sequence = cleaned_sequence[3:]
                trim_start_codon = True
                logger.log_step("Start Codon Removal", "Start codon removed from the beginning of the sequence")
            else:
                logger.log_step("Start Codon Check", "No start codon removal performed")

            # Check for stop codons
            stop_codons = {"TAA", "TAG", "TGA"}
            if len(cleaned_sequence) >= 3 and cleaned_sequence[-3:] in stop_codons:
                cleaned_sequence = cleaned_sequence[:-3]
                trim_stop_codon = True
                logger.log_step("Stop Codon Removal", "Stop codon removed from the end of the sequence")
            else:
                logger.log_step("Stop Codon Check", "No stop codon removal performed")

            # Handle frame adjustment
            final_length = len(cleaned_sequence)
            remainder = final_length % 3

            if remainder != 0:
                if trim_start_codon:
                    logger.log_step("Frame Adjustment", f"Trimming {remainder} bases from the end due to frame remainder after start codon removal")
                    cleaned_sequence = cleaned_sequence[:-remainder]
                elif trim_stop_codon:
                    logger.log_step("Frame Adjustment", f"Trimming {remainder} bases from the beginning due to frame remainder after stop codon removal")
                    cleaned_sequence = cleaned_sequence[remainder:]
                else:
                    logger.log_step("Frame Adjustment", f"Sequence length {final_length} is not a multiple of 3 and no codon removal detected", level=logging.WARNING)

            logger.log_step("Sequence Length", f"Original length: {sequence_length}, Cleaned length: {len(cleaned_sequence)}")

            # Create message based on what was adjusted
            if not in_frame:
                if trim_start_codon and trim_stop_codon:
                    message = ("Provided sequence does not appear to be in frame, using provided start codon "
                               "to infer translation frame. Stop and start codons detected and have been removed.")
                elif trim_start_codon:
                    message = ("Provided sequence does not appear to be in frame, using provided start codon "
                               "to infer translation frame. Start codon has been removed.")
                elif trim_stop_codon:
                    message = ("Provided sequence does not appear to be in frame, using provided stop codon "
                               "to infer frame. Stop codon has been removed.")
                else:
                    message = ("Provided sequence does not appear to be in frame. If this is not intended, please check the sequence.")
                    logger.log_step("Frame Warning", "Sequence not in frame and no codon trimming performed", level=logging.ERROR)
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

            logger.log_step("Preprocessing Complete", f"Final cleaned sequence: {str(cleaned_sequence)}")
        return str(cleaned_sequence), message, True
