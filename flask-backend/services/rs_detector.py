from Bio.Seq import Seq, CodonTable
import re
from typing import List
from config.logging_config import logger
from models import RestrictionSite, Codon
from debug import MutationDebugger, DebugMixin, debug_context
from .utils import GoldenGateUtils

class RestrictionSiteDetector:
    """
    Restriction Site Detector Module

    This module identifies restriction enzyme recognition sites (BsmBI and BsaI)
    within a DNA sequence, extracting context information, codons overlapping the site,
    and providing a summary of the sites found.
    """

    def __init__(self, verbose: bool = False, debug: bool = False):
        self.logger = logger.getChild("RestrictionSiteDetector")
        self.verbose = verbose
        self.debug = debug
        self.debugger = None

        if self.debug:
            self.debugger = MutationDebugger(
                parent_logger=logger, use_custom_format=True)
            if hasattr(self.debugger.logger, 'propagate'):
                self.debugger.logger.propagate = False
            self.logger.info("Debug mode enabled for RestrictionSiteDetector")
        self.utils = GoldenGateUtils()
        
    @DebugMixin.debug_wrapper
    def find_sites_to_mutate(self, sequence: str, index: int) -> List[RestrictionSite]:
        """
        Finds both BsmBI and BsaI restriction enzyme recognition sites on both strands of a DNA sequence.
        """
        seq_str = str(sequence).upper()
        with debug_context("find_bsmbi_bsai_sites"):
            recognition_sequences = {
                'BsmBI': 'CGTCTC',
                'BsaI': 'GGTCTC'
            }
            sites_to_mutate = []

            # For each enzyme, search for both forward and reverse complement matches
            for enzyme, rec_seq in recognition_sequences.items():
                patterns = [
                    (rec_seq, '+', 'recognition_sequence'),
                    (str(Seq(rec_seq).reverse_complement()), '-', 'sequence')
                ]
                for pattern, strand, pattern_key in patterns:
                    for match in re.finditer(re.escape(pattern), seq_str):
                        index = match.start()
                        frame = index % 3

                        # Extract context (30bp upstream and downstream)
                        start_context = max(0, index - 30)
                        end_context = min(
                            len(seq_str), index + len(pattern) + 30)
                        context_seq = seq_str[start_context:end_context]

                        relative_index = index - start_context

                        codons = self.get_codons(
                            context_seq, relative_index, frame)
                        context_recognition_site_indices = [
                            i - start_context for i in range(index, index + len(pattern))]

                        site_to_mutate = RestrictionSite(
                            position = index,
                            frame = frame,
                            codons = codons,
                            strand = strand,
                            context_seq = context_seq,
                            context_rs_indices = context_recognition_site_indices,
                            context_first_base = start_context,
                            context_last_base = end_context,
                            enzyme = enzyme
                        )
                        
                        sites_to_mutate.append(site_to_mutate)

            sites_to_mutate.sort(key=lambda site: site.position)

            self.log_step(
                "Restriction sites found: \n"
                f"{self.utils.summarize_restriction_sites(sites_to_mutate)}"
            )
            
            return sites_to_mutate

    @DebugMixin.debug_wrapper
    def get_codons(self, context_seq: str, recognition_start_index: int, frame: int) -> List[Codon]:
        """
        Extracts codons spanned by the recognition site from a context sequence and 
        immediately stores a tuple (of 0s and 1s) for each codon indicating which 
        bases are within the restriction enzyme recognition site.
        """
        with debug_context("find_codons"):
            codons = []
            translation_table = CodonTable.unambiguous_dna_by_id[1]

            # Define codon positions and the corresponding overlap tuples based on the reading frame.
            if frame == 0:
                codon_positions = [
                    recognition_start_index,
                    recognition_start_index + 3
                ]
                # Each tuple represents the overlap (three bases) for the respective codon.
                codon_bases_in_rs = [
                    (1, 1, 1),  # first codon: all bases are in the recognition site
                    (1, 1, 1)   # second codon: all bases are in the recognition site
                ]
            elif frame == 1:
                codon_positions = [
                    recognition_start_index - 1,
                    recognition_start_index + 2,
                    recognition_start_index + 5
                ]
                codon_bases_in_rs = [
                    (0, 1, 1),  # first codon: only positions 1 and 2 are in the site
                    (1, 1, 1),  # second codon: all bases are in the site
                    (1, 0, 0)   # third codon: only position 0 is in the site
                ]
            elif frame == 2:
                codon_positions = [
                    recognition_start_index - 2,
                    recognition_start_index + 1,
                    recognition_start_index + 4
                ]
                codon_bases_in_rs = [
                    (0, 0, 1),  # first codon: only position 2 is in the site
                    (1, 1, 1),  # second codon: all bases are in the site
                    (1, 1, 0)   # third codon: only positions 0 and 1 are in the site
                ]
            else:
                return []

            # Iterate through codon positions and create Codon objects with the appropriate overlap tuple.
            for codon_index, pos in enumerate(codon_positions):
                if 0 <= pos <= len(context_seq) - 3:
                    codon_seq = context_seq[pos: pos + 3]
                    # The overlap tuple is taken directly from our list.
                    overlap = codon_bases_in_rs[codon_index]

                    codon = Codon(
                        amino_acid=translation_table.forward_table.get(str(codon_seq), 'X'),
                        context_position=pos,
                        codon_sequence=str(codon_seq),
                        recognition_overlap=overlap  # Store the 3-dimensional tuple.
                    )

                    codons.append(codon)

            return codons



