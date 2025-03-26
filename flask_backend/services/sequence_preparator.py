# services/sequence_prep.py
import logging
from typing import Optional, Tuple
from Bio.Seq import Seq

from flask_backend.logging import logger

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

    def preprocess_sequence(
        self,
        sequence: str,
        matk_part_left: str,
        progress_callback: Optional[callable] = None
    ) -> Tuple[Optional[Seq], str, bool]:
        """
        Processes a DNA sequence by removing start/stop codons and ensuring proper frame.
        
        Args:
            sequence: Input DNA sequence
            matk_part_left: MTK part number on the left side
            progress_callback: Optional callback for progress reporting
            
        Returns:
            Tuple containing (processed sequence, message, success flag)
        """
        # Report initial progress
        if progress_callback:
            progress_callback(0, "Starting sequence preprocessing", "progress")
            
        with logger.debug_context("pre_process_sequence"):
            # Step 1: Initial conversion and setup (0-20%)
            logger.log_step("Conversion", "Converting input sequence to uppercase and Seq object")
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
            logger.log_step("Validation", f"Initial frame check: in_frame = {in_frame}")
            
            if progress_callback:
                progress_callback(20, "Frame checked, looking for start/stop codons", "progress")

            # Step 2: Start codon check (20-40%)
            # Check for start codon only if matk_part_left is "3" or "3a"
            if matk_part_left in {"3", "3a"} and len(cleaned_sequence) >= 3 and cleaned_sequence[:3] == "ATG":
                cleaned_sequence = cleaned_sequence[3:]
                trim_start_codon = True
                logger.log_step("Start Codon Removal", "Start codon removed from the beginning of the sequence")
                if progress_callback:
                    progress_callback(40, "Start codon detected and removed", "progress")
            else:
                logger.log_step("Start Codon Check", "No start codon removal performed")
                if progress_callback:
                    progress_callback(40, "No start codon detected or removal needed", "progress")

            # Step 3: Stop codon check (40-60%)
            # Check for stop codons
            stop_codons = {"TAA", "TAG", "TGA"}
            if len(cleaned_sequence) >= 3 and cleaned_sequence[-3:] in stop_codons:
                cleaned_sequence = cleaned_sequence[:-3]
                trim_stop_codon = True
                logger.log_step("Stop Codon Removal", "Stop codon removed from the end of the sequence")
                if progress_callback:
                    progress_callback(60, "Stop codon detected and removed", "progress")
            else:
                logger.log_step("Stop Codon Check", "No stop codon removal performed")
                if progress_callback:
                    progress_callback(60, "No stop codon detected or removal needed", "progress")

            # Step 4: Frame adjustment (60-80%)
            # Handle frame adjustment
            final_length = len(cleaned_sequence)
            remainder = final_length % 3

            if remainder != 0:
                if trim_start_codon:
                    logger.log_step("Frame Adjustment", f"Trimming {remainder} bases from the end due to frame remainder after start codon removal")
                    cleaned_sequence = cleaned_sequence[:-remainder]
                    if progress_callback:
                        progress_callback(80, f"Adjusting frame by trimming {remainder} bases from the end", "progress")
                elif trim_stop_codon:
                    logger.log_step("Frame Adjustment", f"Trimming {remainder} bases from the beginning due to frame remainder after stop codon removal")
                    cleaned_sequence = cleaned_sequence[remainder:]
                    if progress_callback:
                        progress_callback(80, f"Adjusting frame by trimming {remainder} bases from the beginning", "progress")
                else:
                    logger.log_step("Frame Adjustment", f"Sequence length {final_length} is not a multiple of 3 and no codon removal detected", level=logging.WARNING)
                    if progress_callback:
                        progress_callback(80, "Warning: Sequence not in frame and no codons to trim", "progress")
            else:
                if progress_callback:
                    progress_callback(80, "Sequence already in frame, no adjustment needed", "progress")

            logger.log_step("Sequence Length", f"Original length: {sequence_length}, Cleaned length: {len(cleaned_sequence)}")

            # Step 5: Final message creation (80-100%)
            # Create message based on what was adjusted
            if not in_frame:
                if trim_start_codon and trim_stop_codon:
                    message = ("Provided sequence does not appear to be in frame, using provided start codon "
                               "to infer translation frame. Stop and start codons detected and have been removed.")
                    notification_count = 3
                elif trim_start_codon:
                    message = ("Provided sequence does not appear to be in frame, using provided start codon "
                               "to infer translation frame. Start codon has been removed.")
                    notification_count = 2
                elif trim_stop_codon:
                    message = ("Provided sequence does not appear to be in frame, using provided stop codon "
                               "to infer frame. Stop codon has been removed.")
                    notification_count = 2
                else:
                    message = ("Provided sequence does not appear to be in frame. If this is not intended, please check the sequence.")
                    notification_count = 1
                    logger.log_step("Frame Warning", "Sequence not in frame and no codon trimming performed", level=logging.ERROR)
                    if progress_callback:
                        progress_callback(100, "Error: Sequence not in frame and cannot be corrected", "progress", notification_count=notification_count)
                    return str(sequence), message, False
            else:
                if trim_start_codon and trim_stop_codon:
                    message = "Start and stop codons detected and removed."
                    notification_count = 2
                elif trim_start_codon:
                    message = "Start codon detected and removed."
                    notification_count = 1
                elif trim_stop_codon:
                    message = "Stop codon detected and removed."
                    notification_count = 1
                else:
                    message = "Sequence is in frame, no codon adjustments needed."
                    notification_count = 0

            logger.log_step("Preprocessing Complete", f"Final cleaned sequence: {str(cleaned_sequence)}")
            
            if progress_callback:
                progress_callback(100, f"Preprocessing complete: {message}", "progress", notification_count=notification_count, display_message=message)
                
        return str(cleaned_sequence), True