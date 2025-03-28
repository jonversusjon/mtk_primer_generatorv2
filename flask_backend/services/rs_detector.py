import logging
import re
from typing import List, Dict
from Bio.Seq import Seq, CodonTable

from flask_backend.models import RestrictionSite, Codon
from flask_backend.logging import logger
from flask_backend.services.utils import GoldenGateUtils

class RestrictionSiteDetector():
    """
    Restriction Site Detector Module

    This module identifies restriction enzyme recognition sites (BsmBI and BsaI)
    within a DNA sequence, extracting context information, codons overlapping the site,
    and providing a summary of the sites found.
    """

    def __init__(
        self,
        codon_dict: Dict,
        verbose: bool = False, debug: bool = False):
        self.verbose = verbose
        self.debug = debug
        logger.log_step("Initialization", f"Initializing RestrictionSiteDetector with verbose={verbose} and debug={debug}")
        self.utils = GoldenGateUtils()

        self.codon_dict = codon_dict
        
    def find_restriction_sites(
        self,
        sequence: str,
        send_update: callable,
        ) -> List[RestrictionSite]:
        """
        Finds both BsmBI and BsaI restriction enzyme recognition sites on both strands of a DNA sequence.
        """
        seq_str = str(sequence).upper()

        recognition_sequences = {
            'BsmBI': 'CGTCTC',
            'BsaI': 'GGTCTC'
        }
        restriction_sites = []

        # For each enzyme, search for both forward and reverse complement matches
        for enzyme, rec_seq in recognition_sequences.items():
            logger.log_step("Enzyme Processing", f"Processing enzyme {enzyme} with recognition sequence {rec_seq}")
            patterns = [
                (rec_seq, '+', 'recognition_sequence'),
                (str(Seq(rec_seq).reverse_complement()), '-', 'sequence')
            ]
            for pattern, strand, pattern_key in patterns:
                for match in re.finditer(re.escape(pattern), seq_str):
                    found_index = match.start()
                    frame = found_index % 3
                    logger.log_step("Match Found", f"Found {enzyme} match at index {found_index} on strand {strand} with frame {frame}")
                    
                    # Extract context (30bp upstream and downstream)
                    start_context = max(0, found_index - 30)
                    end_context = min(len(seq_str), found_index + len(pattern) + 30)
                    context_seq = seq_str[start_context:end_context]
                    logger.log_step("Context Extraction", f"Extracted context from {start_context} to {end_context}")

                    relative_index = found_index - start_context
                    codons = self.get_codons(context_seq, relative_index, frame)
                    context_recognition_site_indices = [i - start_context for i in range(found_index, found_index + len(pattern))]

                    restriction_site = RestrictionSite(
                        position=found_index,
                        frame=frame,
                        codons=codons,
                        strand=strand,
                        recognition_seq=pattern,
                        context_seq=context_seq,
                        context_rs_indices=context_recognition_site_indices,
                        
                        context_first_base=start_context,
                        context_last_base=end_context,
                        enzyme=enzyme
                    )
                    
                    restriction_sites.append(restriction_site)
                    logger.log_step("Site Added", f"Added site at position {found_index} for enzyme {enzyme}")

        restriction_sites.sort(key=lambda site: site.position)
        logger.log_step("Result", f"Total sites found: {len(restriction_sites)}")
                
        if restriction_sites:
            send_update(message=f"Found {len(restriction_sites)} restriction sites", prog=100, restriction_sites=restriction_sites, notification_count=len(restriction_sites))
        else:
            send_update(message="No restriction sites found", prog=100, restriction_sites=restriction_sites, callout="No site mutations needed", notification_count=0)
                    
        return restriction_sites

    def get_codons(self, context_seq: str, recognition_start_index: int, frame: int) -> List[Codon]:
        """
        Extracts codons spanned by the recognition site from a context sequence and 
        immediately stores a tuple (of 0s and 1s) for each codon indicating which 
        bases are within the restriction enzyme recognition site.
        """
        logger.log_step("Codon Extraction", f"Extracting codons from context sequence with recognition_start_index={recognition_start_index} and frame={frame}")
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
                [0, 1, 2], # first codon: all bases are in the recognition site
                [0, 1, 2]  # second codon: all bases are in the recognition site
            ]

        elif frame == 1:
            codon_positions = [
                recognition_start_index - 1,
                recognition_start_index + 2,
                recognition_start_index + 5
            ]
            codon_bases_in_rs = [
                [1, 2],     # first codon: only positions 1 and 2 are in the site
                [0, 1, 2],  # second codon: all bases are in the site
                [0]         # third codon: only position 0 is in the site
            ]
        elif frame == 2:
            codon_positions = [
                recognition_start_index - 2,
                recognition_start_index + 1,
                recognition_start_index + 4
            ]
            codon_bases_in_rs = [
                [0],        # first codon: only position 2 is in the site
                [0, 1, 2],  # second codon: all bases are in the site
                [0, 1]      # third codon: only positions 0 and 1 are in the site
            ]
        else:
            logger.log_step("Codon Extraction", f"Invalid frame {frame} encountered, returning empty codon list", level=logging.WARNING)
            return []

        # Iterate through codon positions and create Codon objects with the appropriate overlap tuple.
        for codon_index, pos in enumerate(codon_positions):
            print(f"Codon index: {codon_index}, Position: {pos}")
            if 0 <= pos <= len(context_seq) - 3:
                codon_seq = context_seq[pos: pos + 3]
                # The overlap list is taken directly from our list.
                overlap = codon_bases_in_rs[codon_index]
                amino_acid = self.utils.get_amino_acid(codon_seq)
                usage = self.utils.get_codon_usage(codon_seq, amino_acid, self.codon_dict)

                codon = Codon(
                    amino_acid=translation_table.forward_table.get(str(codon_seq), 'X'),
                    context_position=pos,
                    codon_sequence=str(codon_seq),
                    rs_overlap=overlap,
                    usage=usage,
                )
                codons.append(codon)
                logger.log_step("Codon Processed", f"Processed codon {codon_seq} at position {pos} with overlap {overlap}")
            else:
                logger.log_step("Codon Skipped", f"Skipping codon at position {pos} due to insufficient length", level=logging.WARNING)

        logger.log_step("Codon Extraction", f"Extracted {len(codons)} codon(s) from context sequence")
        return codons
