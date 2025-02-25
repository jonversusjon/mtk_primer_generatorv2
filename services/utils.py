# services/utils.py
import json
import os
import csv
import logging
from functools import lru_cache
from Bio.Seq import Seq
from Bio.Data import CodonTable
from typing import Dict, List, Optional, Any
from Bio.Data import CodonTable
import numpy as np
import logging
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
        """
        Returns the reverse complement of a DNA sequence.

        Args:
            seq (str): The DNA sequence.

        Returns:
            str: The reverse complement of the input sequence.
        """
        return str(Seq(seq).reverse_complement())

    def get_amino_acid(self, codon: str) -> str:
        """Translates a codon to its amino acid."""
        return str(Seq(codon).translate())

    def get_codons_for_amino_acid(self, amino_acid: str):
        """
        Returns a list of codons (as strings) that encode the given amino acid,
        based on the standard genetic code (NCBI table 1).

        Args:
            amino_acid (str): One-letter amino acid code (e.g., 'E' for Glutamic acid).

        Returns:
            list[str]: A list of codon strings.

        Raises:
            TypeError: If amino_acid is not a string.
            ValueError: If amino_acid is not a single valid character.
        """

        if not isinstance(amino_acid, str):
            raise TypeError("amino_acid must be a string")

        if len(amino_acid) != 1:
            raise ValueError(
                "amino_acid must be a single character representing the amino acid")

        valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY*")
        amino_acid = amino_acid.upper()

        if amino_acid not in valid_amino_acids:
            raise ValueError(
                f"'{amino_acid}' is not a valid amino acid one-letter code")

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
        """Parses and restructures form data from HTMX."""
        with debug_context("package_form_data"):
            print(f"Got request data: {request}")

            try:
                # First check if it's JSON data
                if request.is_json:
                    print("Received JSON data")
                    data = request.json
                    print(f"JSON data: {data}")
                    return data, None

                # Then try form data
                elif request.form:
                    print("Received form data")
                    # Determine data format and get data
                    data = self._get_request_data(request)

                    if data is None:
                        return None, "No data received"

                    if not isinstance(data, dict):
                        return None, f"Invalid data format: expected dict, got {type(data)}"

                    # Extract and validate main fields
                    try:
                        packaged_data = self._extract_form_fields(data)
                        logger.debug(f"Extracted form fields: {packaged_data}")

                        if packaged_data is None:
                            return None, "Failed to extract form fields"

                        # Ensure verboseMode exists with a default value
                        packaged_data['verboseMode'] = packaged_data.get('verboseMode', False)
                        logger.debug(f"package_data verbose: {packaged_data['verboseMode']}")

                        return packaged_data, None

                    except Exception as e:
                        logger.error(f"Error in _extract_form_fields: {str(e)}")
                        return None, str(e)

                # If neither JSON nor form data
                else:
                    print("No data received")
                    print(f"Request data: {request.get_data()}")
                    return None, "No data received"

            except Exception as e:
                print(f"Error packaging form data: {str(e)}")
                return None, str(e)

    def _get_request_data(self, request: Any) -> dict:
        """Get data from request object, handling different content types."""
        try:
            raw_data = request.get_data()

            if request.is_json:
                # Handle JSON data
                try:
                    json_data = request.get_json(force=True)
                    logger.debug(f"Parsed JSON data: {json_data}")
                    return json_data
                except Exception as e:
                    logger.error(f"Error parsing JSON: {str(e)}")
                    return None
            elif request.form:
                # Handle form data
                form_data = dict(request.form)
                logger.debug(f"Form data: {form_data}")
                return form_data
            elif raw_data:
                # Try to parse raw data as JSON
                try:
                    json_data = json.loads(raw_data)
                    logger.debug(f"Parsed raw data as JSON: {json_data}")
                    return json_data
                except json.JSONDecodeError as e:
                    logger.error(f"Failed to parse raw data as JSON: {e}")
                    return None
            else:
                logger.error("No data found in request")
                return None
                
        except Exception as e:
            logger.error(f"Error getting request data: {str(e)}")
            return None

    def _extract_form_fields(self, data: dict) -> dict:
        """Extract and validate form fields from the request data."""
        try:
            logger.debug("Extracting form fields from data: %s", data)
            
            # Initialize the packaged data with defaults
            packaged_data = {
                'templateSequence': data.get('templateSequence', ''),
                'numSequences': int(data.get('numSequences', 1)),
                'species': data.get('species', ''),
                'kozak': data.get('kozak', 'MTK'),
                'sequences': [],
                'verboseMode': data.get('verboseMode', False)
            }
            
            # Extract sequence information
            sequences = data.get('sequences', [])
            logger.debug("Found sequences: %s", sequences)
            
            # Validate each sequence
            for seq in sequences:
                if not isinstance(seq, dict):
                    logger.error(f"Invalid sequence format: {seq}")
                    continue
                    
                sequence_data = {
                    'primerName': seq.get('primerName', ''),
                    'mtkPart': seq.get('mtkPart', ''),
                    'sequence': seq.get('sequence', '')
                }
                
                # Log the sequence data being added
                logger.debug("Adding sequence data: %s", sequence_data)
                packaged_data['sequences'].append(sequence_data)
                
            # Check if we have at least one sequence
            if not packaged_data['sequences']:
                logger.error("No valid sequences found in form data")
                packaged_data['sequences'] = [{
                    'primerName': '',
                    'mtkPart': '',
                    'sequence': ''
                }]
                
            return packaged_data
            
        except Exception as e:
            logger.error(f"Error in _extract_form_fields: {str(e)}")
            raise

    def _get_first_value(
        self,
        data: Dict,
        key: str,
        default: Any = "",
        cast_type: type = str
    ) -> Any:
        """Safely retrieves first value from form field."""
        value = data.get(key, default)

        if isinstance(value, list):
            value = value[0] if value else default

        try:
            return cast_type(value)
        except (ValueError, TypeError):
            logger.error(
                f"Failed to cast {key}={value} to {cast_type.__name__}")
            return default

    def validate_form_data(
        self,
        data: Dict,
        required_fields: Optional[List[str]] = None,
        field_validators: Optional[Dict] = None
    ) -> List[str]:
        """Validates form data against requirements and rules."""
        with debug_context("validate_form_data"):
            error_messages = []

            if required_fields:
                missing_fields = [
                    field for field in required_fields
                    if field not in data
                ]
                if missing_fields:
                    error_messages.append(
                        f"Missing required fields: {', '.join(missing_fields)}"
                    )

            if field_validators:
                for field, validator in field_validators.items():
                    if field in data:
                        validation_error = validator(data[field])
                        if validation_error:
                            error_messages.append(validation_error)

            return error_messages

    def seq_to_index(self, seq: str) -> int:
        """
        Converts a 4-nucleotide sequence to its corresponding matrix index.
        Sequences are indexed in alphabetical order (A=0, C=1, G=2, T=3).

        Args:
            seq: String containing exactly 4 nucleotides (ACTG only)

        Returns:
            Integer index from 0 to 255 corresponding to sequence position in matrix

        Raises:
            ValueError: If sequence is not exactly 4 nucleotides or contains invalid characters
        """
        if len(seq) != 4:
            raise ValueError(
                f"Sequence must be exactly 4 nucleotides, got {len(seq)}: {seq}")

        seq = seq.upper()
        if not all(nt in 'ACTG' for nt in seq):
            raise ValueError(
                f"Sequence must contain only A, C, G, T, got: {seq}")

        NT_VALUES = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        index = 0
        for pos, nt in enumerate(seq):
            power = 3 - pos
            index += NT_VALUES[nt] * (4 ** power)

        return index

    def get_recognition_site_bases(self, frame, codon_index):
        """
        Given the reading frame and the codon index (0 for first codon, etc.),
        returns a list of indices (0, 1, 2) indicating which bases in the codon
        are within the 6-base restriction enzyme recognition site.

        Args:
            frame (int): The reading frame (0, 1, or 2).
            codon_index (int): The index of the codon in the list (0 for first, etc.)

        Returns:
            List[int]: The indices of the codon that are part of the recognition site.
        """
        if frame == 0:
            # With frame 0, the recognition site aligns with codon boundaries:
            return [0, 1, 2]
        elif frame == 1:
            # For frame=1, we assume three codons:
            # First codon: only its last two bases (indices 1 and 2) are within the site.
            # Second codon: all three bases are in.
            # Third codon: only its first base (index 0) is in.
            if codon_index == 0:
                return [1, 2]
            elif codon_index == 1:
                return [0, 1, 2]
            elif codon_index == 2:
                return [0]
        elif frame == 2:
            # For frame=2, we assume three codons:
            # First codon: only its third base (index 2) is within the site.
            # Second codon: all three bases are in.
            # Third codon: only its first two bases (indices 0 and 1) are in.
            if codon_index == 0:
                return [2]
            elif codon_index == 1:
                return [0, 1, 2]
            elif codon_index == 2:
                return [0, 1]
        else:
            raise ValueError("Frame must be 0, 1, or 2")

    def _load_compatibility_table(self, path: str) -> np.ndarray:
        """
        Loads the binary compatibility table into a numpy array.

        The binary file contains a 256×256 compatibility matrix where each element
        represents whether two 4-nucleotide sequences are compatible. The sequences
        are ordered alphabetically, so AA[AA] is at [0,0] and TT[TT] is at [255,255].

        Args:
            path: Path to the binary compatibility table file

        Returns:
            256x256 numpy array where element [i,j] indicates if sequence i is compatible with j

        Raises:
            FileNotFoundError: If compatibility table file cannot be found
            ValueError: If table dimensions or content are invalid
        """
        try:
            with open(path, 'rb') as f:
                binary_data = f.read()

            expected_size = 256 * 256 // 8  # 65,536 bits = 8,192 bytes
            if len(binary_data) != expected_size:
                raise ValueError(
                    f"Invalid compatibility table size. Expected {expected_size} bytes, "
                    f"got {len(binary_data)} bytes"
                )

            compatibility_bits = np.unpackbits(
                np.frombuffer(binary_data, dtype=np.uint8)
            )
            compatibility_matrix = compatibility_bits.reshape(256, 256)

            if self.verbose:
                logger.info(
                    f"Loaded compatibility matrix with shape: {compatibility_matrix.shape}")
                logger.info(
                    f"Number of compatible pairs: {np.sum(compatibility_matrix)}")
                logger.info(
                    f"Matrix density: {np.mean(compatibility_matrix):.2%}")

            return compatibility_matrix

        except FileNotFoundError:
            raise FileNotFoundError(
                f"Compatibility table not found at {path}. "
                "Please ensure the binary file is in the correct location."
            )
        except Exception as e:
            raise ValueError(f"Error loading compatibility table: {str(e)}")

    def get_codon_usage(self, codon: str, amino_acid: str, codon_usage_dict: dict, default_usage: float = 0.0) -> float:
        """
        Converts a DNA codon (T → U) to RNA and retrieves its codon usage frequency.

        Args:
            codon (str): A three-letter DNA codon (e.g., "GAA").
            amino_acid (str): The corresponding amino acid (single-letter code).
            codon_usage_dict (dict): Dictionary mapping amino acids to codon usage frequencies.
            default_usage (float, optional): Default usage value if the codon is not found. Defaults to 1.0.

        Returns:
            float: The codon usage frequency.
        """
        if not isinstance(codon, str) or len(codon) != 3:
            raise ValueError(
                f"Invalid codon: {codon}. Must be a three-letter DNA string.")

        if amino_acid not in codon_usage_dict:
            print(
                f"Warning: Amino acid {amino_acid} not found in codon usage dictionary.")
            return default_usage

        # Convert DNA (T) to RNA (U) for lookup
        codon_rna = codon.replace("T", "U")

        # Retrieve codon usage, return default if not found
        usage_value = codon_usage_dict[amino_acid].get(
            codon_rna, default_usage)

        return usage_value
