# services/sequence_prep.py
from Bio.Seq import Seq, CodonTable
from prettytable import PrettyTable
import re
from typing import Dict, Optional, Union, List, Tuple
from collections import defaultdict
from config.logging_config import logger
from .base import debug_context
from models.sequences import SequenceToDomesticate
from services.debug.debug_mixin import DebugMixin
from services.debug.debug_utils import MutationDebugger


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
    def preprocess_sequence(self, sequence: SequenceToDomesticate, matk_part_left: str) -> Tuple[Optional[Seq], str, bool]:
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
