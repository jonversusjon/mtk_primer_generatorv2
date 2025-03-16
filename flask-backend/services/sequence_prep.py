# services/sequence_prep.py
from Bio.Seq import Seq
from typing import Optional, Tuple
from config.logging_config import logger
from services.debug import DebugMixin, MutationDebugger, debug_context


class SequencePreparator:
    """
    Sequence Preparation Module

    This module handles preprocessing of DNA sequences for Golden Gate assembly,
    including the removal of start/stop codons, ensuring sequences are in the correct reading frame,

    """

    def __init__(self, verbose: bool = False, debug: bool = False):
        self.logger = logger.getChild("SequencePreparator")
        self.debug = debug
        self.debugger = None

        if self.debug:
            self.debugger = MutationDebugger(
                parent_logger=logger, use_custom_format=True)
            if hasattr(self.debugger.logger, 'propagate'):
                self.debugger.logger.propagate = False
            self.logger.info("Debug mode enabled for PrimerDesigner")

        self.verbose = verbose
        self.state = {
            'current_sequence': None,
            'adjustments_made': [],
            'restriction_sites': {}
        }

    @DebugMixin.debug_wrapper
    def preprocess_sequence(
        self,
        sequence: str,
        matk_part_left: str
            ) -> Tuple[Optional[Seq], str, bool]:
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

            # Check for start codon only if matk_part_left is "3" or "3a"
            if matk_part_left in {"3", "3a"} and len(cleaned_sequence) >= 3 and cleaned_sequence[:3] == "ATG":
                cleaned_sequence = cleaned_sequence[3:]
                trim_start_codon = True
                logger.info(
                    'Start codon removed from the beginning of the sequence')

            # Check for stop codons
            stop_codons = {"TAA", "TAG", "TGA"}
            if len(cleaned_sequence) >= 3 and cleaned_sequence[-3:] in stop_codons:
                cleaned_sequence = cleaned_sequence[:-3]
                trim_stop_codon = True
                logger.info('Stop codon removed from the end of the sequence')

            # Handle frame adjustment
            final_length = len(cleaned_sequence)
            remainder = final_length % 3

            if remainder != 0:
                if trim_start_codon:
                    # Trim from the end
                    cleaned_sequence = cleaned_sequence[:-remainder]
                elif trim_stop_codon:
                    # Trim from the beginning
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
