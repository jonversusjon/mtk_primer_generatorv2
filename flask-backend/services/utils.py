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

    def package_form_data(self, request: Any) -> tuple:
        """Parses and restructures form data."""
        with debug_context("package_form_data"):
            try:
                # First check if it's JSON data
                if request.is_json:
                    data = request.json
                    return data, None

                # Then try form data
                elif request.form:
                    # Get data from request
                    data = self._get_request_data(request)

                    if data is None:
                        return None, "No data received"

                    # Extract fields assuming valid data
                    packaged_data = {
                        'templateSequence': data.get('templateSequence', ''),
                        'numSequences': int(data.get('numSequences', 1)),
                        'species': data.get('species', ''),
                        'kozak': data.get('kozak', 'MTK'),
                        'sequences': data.get('sequences', []),
                        'verboseMode': data.get('verboseMode', False)
                    }

                    return packaged_data, None

                else:
                    return None, "No data received"

            except Exception as e:
                logger.error(f"Error packaging form data: {str(e)}")
                return None, str(e)

    def _get_request_data(self, request: Any) -> dict:
        """Get data from request object, handling different content types."""
        try:
            if request.is_json:
                return request.get_json(force=True)
            elif request.form:
                return dict(request.form)
            else:
                raw_data = request.get_data()
                if raw_data:
                    return json.loads(raw_data)
                return None

        except Exception as e:
            logger.error(f"Error getting request data: {str(e)}")
            return None

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
        """Convert non-serializable objects (Seq, ndarray) to JSON-compatible types."""
        if isinstance(obj, Seq):
            return str(obj)  # Convert BioPython Seq to string
        elif isinstance(obj, np.ndarray):
            return obj.tolist()  # Convert NumPy arrays to lists
        elif isinstance(obj, dict):
            return {k: self.convert_non_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self.convert_non_serializable(v) for v in obj]
        return obj  # Return as-is if it's already JSON serializable
