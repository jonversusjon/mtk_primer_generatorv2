from Bio.Seq import Seq, CodonTable
import re
from prettytable import PrettyTable
from typing import Dict, List, Union
from collections import defaultdict
from config.logging_config import logger
from .base import debug_context


class RestrictionSiteDetector:
    """
    Restriction Site Detector Module

    This module identifies restriction enzyme recognition sites (BsmBI and BsaI)
    within a DNA sequence, extracting context information, codons overlapping the site,
    and providing a summary of the sites found.
    """

    def __init__(self, verbose: bool = False):
        self.logger = logger.getChild("RestrictionSiteDetector")
        self.verbose = verbose

    def get_codons(self, context_seq: str, recognition_start_index: int, length: int, frame: int) -> List[Dict]:
        """
        Extracts codons spanned by the recognition site from a context sequence.
        """
        with debug_context("find_codons"):
            codons = []
            translation_table = CodonTable.unambiguous_dna_by_id[1]

            # Determine codon positions based on reading frame
            if frame == 0:
                codon_positions = [recognition_start_index,
                                   recognition_start_index + 3]
            elif frame == 1:
                codon_positions = [recognition_start_index - 1,
                                   recognition_start_index + 2, recognition_start_index + 5]
            elif frame == 2:
                codon_positions = [recognition_start_index - 2,
                                   recognition_start_index + 1, recognition_start_index + 4]
            else:
                return []

            for pos in codon_positions:
                if 0 <= pos <= len(context_seq) - 3:
                    codon_seq = context_seq[pos: pos + 3]
                    codons.append({
                        "codon_sequence": str(codon_seq),
                        "amino_acid": translation_table.forward_table.get(str(codon_seq), 'X'),
                        "context_position": pos  # relative position in context_seq
                    })

            return codons

    def find_bsmbi_bsai_sites(self, sequence: Union[str, Seq]) -> List[Dict]:
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
                            context_seq, relative_index, len(pattern), frame)
                        context_recognition_site_indices = [
                            i - start_context for i in range(index, index + len(pattern))]

                        site_details = {
                            'position': index,
                            'frame': frame,
                            'codons': codons,
                            'strand': strand,
                            'context_sequence': context_seq,
                            'context_recognition_site_indices': context_recognition_site_indices,
                            'context_first_base': start_context,
                            'context_last_base': end_context,
                            'enzyme': enzyme
                        }
                        site_details[pattern_key] = pattern
                        sites_to_mutate.append(site_details)

            sites_to_mutate.sort(key=lambda site: site['position'])
            
            # Pydantic validation of the sites
            from models.mtk import RestrictionSite
            validated_sites = []
            for site in sites_to_mutate:
                try:
                    validated_site = RestrictionSite.parse_obj(site)
                    validated_sites.append(validated_site.dict())
                except Exception as e:
                    self.logger.error(f"Validation error in restriction site: {e}")
                    raise e
            return validated_sites

    def summarize_bsmbi_bsai_sites(self, site_details: Union[List[Dict], Dict[str, List[Dict]]]) -> None:
        """
        Creates a formatted summary of restriction sites.
        """
        with debug_context("summarize_restriction_sites"):
            site_type_descriptions = {
                'BsmBI': 'BsmBI Restriction Site',
                'BsaI': 'BsaI Restriction Site',
            }

            if isinstance(site_details, list):
                grouped_sites = defaultdict(list)
                for site in site_details:
                    enzyme = site.get('enzyme')
                    grouped_sites[enzyme].append(site)
                site_details = grouped_sites

            table = PrettyTable()
            table.field_names = ["Site Type",
                                 "Number of Instances", "Position(s)"]

            for enzyme, details_list in site_details.items():
                if not details_list:
                    continue

                enzyme_desc = site_type_descriptions.get(enzyme, enzyme)
                positions = ', '.join(
                    str(detail['position'] + 1) for detail in details_list)
                table.add_row([enzyme_desc, len(details_list), positions])

            logger.info("\nRestriction Site Analysis Summary:")
            logger.info(f"\n{table}")

    def find_sites_to_mutate(self, sequence: Seq, index: int) -> List[Dict]:
        """
        Finds and summarizes restriction sites needing mutation.
        """
        with debug_context("find_sites_to_mutate"):
            sites_to_mutate = self.find_bsmbi_bsai_sites(sequence)
            if sites_to_mutate:
                self.summarize_bsmbi_bsai_sites(sites_to_mutate)
            return sites_to_mutate
