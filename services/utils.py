# services/utils.py
import json
import os
import csv
import logging
from functools import lru_cache
from Bio.Seq import Seq
from Bio.Data import CodonTable
from typing import Dict, List, Optional, Any
from .base import GoldenGateDesigner


class GoldenGateUtils(GoldenGateDesigner):
    def __init__(self, verbose: bool = False):
        super().__init__(verbose=verbose)
        self.data_dir = os.path.join(os.path.dirname(__file__), "../static/data")
        self.codon_tables_dir = os.path.join(self.data_dir, "codon_usage_tables")

    def load_json_file(self, filename: str) -> Optional[Dict]:
        """Loads a JSON file from the static/data directory."""
        with self.debug_context("load_json_file"):
            filepath = os.path.join(self.data_dir, filename)
            
            try:
                with open(filepath, "r") as file:
                    return json.load(file)
            except FileNotFoundError:
                self.logger.error(f"File not found: {filepath}")
                return None
            except json.JSONDecodeError:
                self.logger.error(f"Invalid JSON format in: {filepath}")
                return None

    @lru_cache(maxsize=10)
    def get_codon_usage_dict(self, species: str) -> Optional[Dict]:
        """Loads a codon usage table for a species."""
        with self.debug_context("get_codon_usage_dict"):
            filename = os.path.join(self.codon_tables_dir, f"{species}.json")
            try:
                with open(filename, "r") as f:
                    return json.load(f)
            except FileNotFoundError:
                self.logger.error(f"Codon usage table not found for species: {species}")
                return None

    @lru_cache(maxsize=1)
    def get_mtk_partend_sequences(self) -> Optional[Dict]:
        """Loads MTK part-end sequences."""
        return self.load_json_file("mtk_partend_sequences.json")

    def get_available_species(self) -> List[str]:
        """Gets list of available species from codon usage tables."""
        with self.debug_context("get_available_species"):
            if os.path.exists(self.codon_tables_dir):
                return [f[:-5] for f in os.listdir(self.codon_tables_dir) 
                        if f.endswith('.json')]
            self.logger.warning("Codon usage tables directory not found")
            return []

    def get_amino_acid(self, codon: str) -> str:
        """Translates a codon to its amino acid."""
        return str(Seq(codon).translate())

    def export_primers_to_tsv(
        self,
        forward_primers: List[tuple],
        reverse_primers: List[tuple],
        filename: str = "primers.tsv"
    ) -> None:
        """Exports primers to TSV file."""
        with self.debug_context("export_primers"):
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
                        
                self.logger.info(f"Primers exported to {filepath}")
                
            except Exception as e:
                self.logger.error(f"Error exporting primers: {str(e)}")
                raise

    def gc_content(self, seq: str) -> float:
        """Computes GC content of a DNA sequence."""
        if not seq:
            return 0
        gc_count = sum(1 for nt in seq.upper() if nt in "GC")
        return gc_count / len(seq)

    def translate_codon(self, codon: str) -> str:
        """Translates a codon into its corresponding amino acid."""
        with self.debug_context("translate_codon"):
            try:
                codon = codon.upper().replace("U", "T")
                standard_table = CodonTable.unambiguous_dna_by_id[1]
                return standard_table.forward_table.get(codon, "?")
            except Exception as e:
                self.logger.error(f"Error translating codon {codon}: {str(e)}")
                return "?"

    # Form handling methods
    def package_form_data(self, request: Any) -> tuple:
        """Parses and restructures form data from HTMX."""
        with self.debug_context("package_form_data"):
            try:
                # Determine data format and get data
                data = self._get_request_data(request)
                if not isinstance(data, dict):
                    return None, "Invalid data format"

                # Extract and validate main fields
                packaged_data = self._extract_form_fields(data)
                
                return packaged_data, None
                
            except Exception as e:
                self.logger.error(f"Error packaging form data: {str(e)}")
                return None, str(e)

    def _get_request_data(self, request: Any) -> Dict:
        """Extracts data from request based on content type."""
        if request.content_type == "application/json":
            return request.get_json()
        elif request.content_type in ["application/x-www-form-urlencoded", "multipart/form-data"]:
            return request.form.to_dict(flat=False)
        else:
            raise ValueError("Unsupported Media Type")

    def _extract_form_fields(self, data: Dict) -> Dict:
        """Extracts and structures form fields."""
        template_sequence = self._get_first_value(data, "templateSequence", "")
        num_sequences = self._get_first_value(data, "numSequences", 1, int)
        species = self._get_first_value(data, "species", "")
        kozak = self._get_first_value(data, "kozak", "")
        verbose_mode = self._get_first_value(data, "verbose_mode", "off").lower() == "on"

        sequences = []
        for i in range(num_sequences):
            sequences.append({
                "primerName": self._get_first_value(data, f"sequences[{i}][primerName]", ""),
                "mtkPart": self._get_first_value(data, f"sequences[{i}][mtkPart]", ""),
                "sequence": self._get_first_value(data, f"sequences[{i}][sequence]", ""),
            })

        return {
            "templateSequence": template_sequence,
            "numSequences": num_sequences,
            "species": species,
            "kozak": kozak,
            "sequences": sequences,
            "verboseMode": verbose_mode
        }

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
            self.logger.error(f"Failed to cast {key}={value} to {cast_type.__name__}")
            return default

    def validate_form_data(
        self,
        data: Dict,
        required_fields: Optional[List[str]] = None,
        field_validators: Optional[Dict] = None
    ) -> List[str]:
        """Validates form data against requirements and rules."""
        with self.debug_context("validate_form_data"):
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

# import json
# import os
# import csv
# import logging
# from functools import lru_cache
# from Bio.Seq import Seq
# from Bio.Data import CodonTable

# logging.basicConfig(level=logging.INFO)

# def load_json_file(filename):
#     """
#     Loads a JSON file from the static/data directory.

#     Args:
#         filename (str): Name of the JSON file.

#     Returns:
#         dict: Parsed JSON data, or None if an error occurs.
#     """
#     filepath = os.path.join(os.path.dirname(__file__), "../static/data", filename)

#     try:
#         with open(filepath, "r") as file:
#             return json.load(file)
#     except FileNotFoundError:
#         logging.error(f"File not found: {filepath}")
#         return None
#     except json.JSONDecodeError:
#         logging.error(f"Invalid JSON format in: {filepath}")
#         return None


# @lru_cache(maxsize=10)
# def get_codon_usage_dict(species):
#     """
#     Loads a codon usage table from a JSON file.

#     Args:
#         species (str): The name of the species (e.g., "homo_sapiens").

#     Returns:
#         dict: The codon usage table, or None if not found.
#     """
#     data_dir = "static/data/codon_usage_tables"
#     filename = os.path.join(data_dir, f"{species}.json")
#     try:
#         with open(filename, "r") as f:
#             codon_usage_table = json.load(f)
#         return codon_usage_table
#     except FileNotFoundError:
#         return None


# @lru_cache(maxsize=1)
# def get_mtk_partend_sequences():
#     """
#     Loads MTK part-end sequences from a JSON file.

#     Returns:
#         dict: The MTK part-end sequences, or None if the file is not found.
#     """
#     return load_json_file("mtk_partend_sequences.json")


# def get_available_species():
#     """
#     Gets a list of available species from JSON files in static/data/codon_usage_tables.

#     Returns:
#         list[str]: A list of species names.
#     """
#     data_dir = os.path.join(os.path.dirname(__file__), "../static/data/codon_usage_tables")

#     if os.path.exists(data_dir):
#         return [f[:-5] for f in os.listdir(data_dir) if f.endswith('.json')]
    
#     return []  # Return an empty list if the directory is missing


# def get_amino_acid(codon):
#     return Seq(codon).translate()


# def export_primers_to_tsv(forward_primers, reverse_primers, filename="primers.tsv"):
#     """
#     Exports the forward and reverse primers to a TSV file.

#     Args:
#         forward_primers (list): List of (name, sequence) tuples.
#         reverse_primers (list): List of (name, sequence) tuples.
#         filename (str): Output file name.
#     """
#     filepath = os.path.join(os.path.dirname(__file__), "../static/data", filename)

#     with open(filepath, mode='w', newline='') as file:
#         writer = csv.writer(file, delimiter='\t')

#         for name, sequence in forward_primers:
#             writer.writerow([name, sequence, "Generated for Golden Gate Assembly"])

#         for name, sequence in reverse_primers:
#             writer.writerow([name, sequence, "Generated for Golden Gate Assembly"])

#     logging.info(f"Primers exported to {filepath}")


# def gc_content(seq):
#     """
#     Computes the GC content of a DNA sequence as a fraction.
#     """
#     gc_count = sum(1 for nt in seq.upper() if nt in "GC")
#     return gc_count / len(seq) if seq else 0


# def package_form_data(request):
#     """
#     Parses and restructures form data from HTMX to create a structured dictionary.
#     """
#     try:
#         # Determine data format
#         if request.content_type == "application/json":
#             data = request.get_json()
#         elif request.content_type in ["application/x-www-form-urlencoded", "multipart/form-data"]:
#             data = request.form.to_dict(flat=False)
#         else:
#             return None, "Unsupported Media Type"

#         if not data:
#             return None, "No data received in the request."
        
#         # Extract main fields
#         template_sequence = _get_first_value(data, "templateSequence", "")
#         num_sequences = _get_first_value(data, "numSequences", 1, int)
#         species = _get_first_value(data, "species", "")
#         kozak = _get_first_value(data, "kozak", "")
#         verbose_mode = _get_first_value(data, "verbose_mode", "off").lower() == "on"

#         # 🔹 Transform sequence keys into a structured list
#         sequences = []
#         for i in range(num_sequences):
#             sequences.append({
#                 "primerName": _get_first_value(data, f"sequences[{i}][primerName]", ""),
#                 "mtkPart": _get_first_value(data, f"sequences[{i}][mtkPart]", ""),
#                 "sequence": _get_first_value(data, f"sequences[{i}][sequence]", ""),
#             })

#         # 🔹 Package data
#         packaged_data = {
#             "templateSequence": template_sequence,
#             "numSequences": num_sequences,
#             "species": species,
#             "kozak": kozak,
#             "sequences": sequences,
#             "verboseMode": verbose_mode
#         }

#         return packaged_data, None  # Success

#     except Exception as e:
#         logging.error(f"Unexpected error in package_form_data: {str(e)}", exc_info=True)
#         return None, f"An unexpected error occurred: {str(e)}"



# def _get_first_value(data, key, default="", cast_type=str):
#     """
#     Safely retrieves the first value from a form field (list or single value).
    
#     Args:
#         data (dict): Form or JSON data.
#         key (str): Key to fetch.
#         default (any): Default value if key is missing.
#         cast_type (type): Data type to cast the value to.
    
#     Returns:
#         Value of the key, properly casted.
#     """
#     value = data.get(key, default)

#     # Ensure value is a string or a single element (not a list)
#     if isinstance(value, list):
#         value = value[0] if value else default  # Use default if list is empty

#     try:
#         return cast_type(value)
#     except (ValueError, TypeError):
#         logging.error(f"Failed to cast {key}={value} to {cast_type.__name__}")
#         return default  # Return default if casting fails


# def validate_form_data(data, required_fields=None, field_validators=None):
#     """
#     Validates form data against required fields and custom validation rules.

#     Args:
#         data (dict): The parsed JSON data from the request.
#         required_fields (list, optional): List of required field names.
#         field_validators (dict, optional): Validation functions {field_name: validator_function}.

#     Returns:
#         list: List of error messages, or an empty list if no errors.
#     """
#     error_messages = []

#     # Check required fields
#     if required_fields:
#         missing_fields = [field for field in required_fields if field not in data]
#         if missing_fields:
#             error_messages.append(f"Missing required fields: {', '.join(missing_fields)}")

#     # Apply field-specific validation functions
#     if field_validators:
#         for field, validator in field_validators.items():
#             if field in data:
#                 validation_error = validator(data[field])
#                 if validation_error:
#                     error_messages.append(validation_error)

#     return error_messages  # Return all validation errors

# def translate_codon(codon):
#     """Translates a codon into its corresponding amino acid."""
#     codon = codon.upper().replace("U", "T")  # Ensure DNA format
#     standard_table = CodonTable.unambiguous_dna_by_id[1]  # Standard table
#     return standard_table.forward_table.get(codon, "?")  # Returns "?" if not found