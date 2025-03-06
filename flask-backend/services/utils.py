# services/utils.py
import json
import os
import csv
from functools import lru_cache
from Bio.Seq import Seq
from Bio.Data import CodonTable
from typing import Dict, List, Optional, Any
import numpy as np
from config.logging_config import logger
from .base import debug_context


class GoldenGateUtils:
    def __init__(self, verbose: bool = False):
        self.verbose = verbose
        self.logger = logger.getChild("GoldenGateUtils")
        self.data_dir = os.path.join(
            os.path.dirname(__file__), "../static/data")
        self.codon_tables_dir = os.path.join(
            self.data_dir, "codon_usage_tables")

    def load_json_file(self, filename: str) -> Optional[Dict]:
        """Loads a JSON file from the static/data directory."""
        with debug_context("load_json_file"):
            filepath = os.path.join(self.data_dir, filename)

            try:
                with open(filepath, "r") as file:
                    return json.load(file)
            except FileNotFoundError:
                logger.error(f"File not found: {filepath}")
                return None
            except json.JSONDecodeError:
                logger.error(f"Invalid JSON format in: {filepath}")
                return None

    @lru_cache(maxsize=10)
    def get_codon_usage_dict(self, species: str) -> Optional[Dict]:
        """Loads a codon usage table for a species."""
        with debug_context("get_codon_usage_dict"):
            filename = os.path.join(self.codon_tables_dir, f"{species}.json")
            try:
                with open(filename, "r") as f:
                    return json.load(f)
            except FileNotFoundError:
                logger.error(
                    f"Codon usage table not found for species: {species}")
                return None

    @lru_cache(maxsize=1)
    def get_mtk_partend_sequences(self) -> Optional[Dict]:
        """Loads MTK part-end sequences."""
        return self.load_json_file("mtk_partend_sequences.json")

    def get_mtk_partend_sequence(self, mtk_part_num: str, primer_direction: str, kozak: str = "MTK") -> Optional[str]:
        """
        Retrieve the correct overhang sequence based on the part number, direction, and kozak preference.

        Args:
            mtk_part_num (str): The part number (e.g., "2", "3", "4a").
            primer_direction (str): Either "forward" or "reverse".
            kozak (str): Either "MTK" (default) or "Canonical".

        Returns:
            Optional[str]: The corresponding primer sequence if found, otherwise None.
        """
        mtk_sequences = self.get_mtk_partend_sequences()
        if not mtk_sequences:
            return None

        # Construct the standard key
        base_key = f"{mtk_part_num}{primer_direction}"

        # Check if a star version should be used for Canonical Kozak cases
        if kozak == "Canonical" and mtk_part_num in {"2", "3", "3a"} and primer_direction == "forward":
            star_key = f"{mtk_part_num}star{primer_direction}"
            return mtk_sequences.get(star_key, mtk_sequences.get(base_key))

        # Default behavior (MTK)
        return mtk_sequences.get(base_key, None)

    def get_available_species(self) -> List[str]:
        """Gets list of available species from codon usage tables."""
        with debug_context("get_available_species"):
            if os.path.exists(self.codon_tables_dir):
                return [f[:-5] for f in os.listdir(self.codon_tables_dir)
                        if f.endswith('.json')]
            logger.warning("Codon usage tables directory not found")
            return []

    def reverse_complement(self, seq: str) -> str:
        """Returns the reverse complement of a DNA sequence."""
        return str(Seq(seq).reverse_complement())

    def get_amino_acid(self, codon: str) -> str:
        """Translates a codon to its amino acid."""
        return str(Seq(codon).translate())

    def get_codons_for_amino_acid(self, amino_acid: str):
        """Returns a list of codons that encode the given amino acid."""
        amino_acid = amino_acid.upper()
        table = CodonTable.unambiguous_dna_by_id[1]

        if amino_acid == '*':
            return table.stop_codons

        return [codon for codon, aa in table.forward_table.items() if aa == amino_acid]

    def export_primers_to_tsv(
        self,
        forward_primers: List[tuple],
        reverse_primers: List[tuple],
        filename: str = "primers.tsv"
    ) -> None:
        """Exports primers to TSV file."""
        with debug_context("export_primers"):
            filepath = os.path.join(self.data_dir, filename)

            try:
                with open(filepath, mode='w', newline='') as file:
                    writer = csv.writer(file, delimiter='\t')

                    for name, sequence in forward_primers:
                        writer.writerow([
                            name,
                            sequence,
                            "Generated for Golden Gate Assembly"
                        ])

                    for name, sequence in reverse_primers:
                        writer.writerow([
                            name,
                            sequence,
                            "Generated for Golden Gate Assembly"
                        ])

                logger.info(f"Primers exported to {filepath}")

            except Exception as e:
                logger.error(f"Error exporting primers: {str(e)}")
                raise

    def gc_content(self, seq: str) -> float:
        """Computes GC content of a DNA sequence."""
        if not seq:
            return 0
        gc_count = sum(1 for nt in seq.upper() if nt in "GC")
        return gc_count / len(seq)

    def translate_codon(self, codon: str) -> str:
        """Translates a codon into its corresponding amino acid."""
        with debug_context("translate_codon"):
            try:
                codon = codon.upper().replace("U", "T")
                standard_table = CodonTable.unambiguous_dna_by_id[1]
                return standard_table.forward_table.get(codon, "?")
            except Exception as e:
                logger.error(f"Error translating codon {codon}: {str(e)}")
                return "?"

    def seq_to_index(self, seq: str) -> int:
        """Converts a 4-nucleotide sequence to its corresponding matrix index."""
        seq = seq.upper()
        NT_VALUES = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        index = 0
        for pos, nt in enumerate(seq):
            power = 3 - pos
            index += NT_VALUES[nt] * (4 ** power)

        return index

    def get_recognition_site_bases(self, frame, codon_index):
        """
        Returns a list of indices indicating which bases in the codon
        are within the restriction enzyme recognition site.
        """
        if frame == 0:
            return [0, 1, 2]
        elif frame == 1:
            if codon_index == 0:
                return [1, 2]
            elif codon_index == 1:
                return [0, 1, 2]
            elif codon_index == 2:
                return [0]
        elif frame == 2:
            if codon_index == 0:
                return [2]
            elif codon_index == 1:
                return [0, 1, 2]
            elif codon_index == 2:
                return [0, 1]
        else:
            return []

    def _load_compatibility_table(self, path: str) -> np.ndarray:
        """Loads the binary compatibility table into a numpy array."""
        with open(path, 'rb') as f:
            binary_data = f.read()

        compatibility_bits = np.unpackbits(
            np.frombuffer(binary_data, dtype=np.uint8)
        )
        compatibility_matrix = compatibility_bits.reshape(256, 256)

        if self.verbose:
            logger.info(
                f"Loaded compatibility matrix with shape: {compatibility_matrix.shape}")
            logger.info(
                f"Number of compatible pairs: {np.sum(compatibility_matrix)}")
            logger.info(f"Matrix density: {np.mean(compatibility_matrix):.2%}")

        return compatibility_matrix

    def get_codon_usage(self, codon: str, amino_acid: str, codon_usage_dict: dict, default_usage: float = 0.0) -> float:
        """Retrieves codon usage frequency."""
        # Convert DNA (T) to RNA (U) for lookup
        codon_rna = codon.replace("T", "U")

        # Retrieve codon usage or default
        if amino_acid in codon_usage_dict:
            return codon_usage_dict[amino_acid].get(codon_rna, default_usage)
        return default_usage

    def convert_non_serializable(self, obj):
        """Convert non-serializable objects (Seq, ndarray, tuple) to JSON-compatible types."""
        if isinstance(obj, Seq):
            return str(obj)  # Convert BioPython Seq to string
        elif isinstance(obj, np.ndarray):
            return obj.tolist()  # Convert NumPy arrays to lists
        elif isinstance(obj, dict):
            return {k: self.convert_non_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self.convert_non_serializable(v) for v in obj]
        elif isinstance(obj, tuple):
            # Convert tuple to list for JSON
            return [self.convert_non_serializable(v) for v in obj]
        return obj  # Return as-is if it's already JSON serializable
