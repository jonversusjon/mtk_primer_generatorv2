from typing import Dict, Tuple, Any, Optional, List, Union
from services.utils import GoldenGateUtils
import re
from config.logging_config import logger
from services.base import debug_context

class ProtocolValidator:
    def __init__(self):
        self.utils = GoldenGateUtils()
        self.logger = logger.getChild("ProtocolValidator")

        self.required_fields = ["sequences", "species", "max-mut-per-site"]
        self.required_sequence_fields = ["sequence", "mtkPart", "primerName"]

    def validate_protocol_inputs(self, data: Dict[str, Any]) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
        """
        Validate inputs for protocol generation.
        
        Args:
            data: Dictionary containing protocol input data
            
        Returns:
            Tuple of (validated_data, error_message). If validation fails,
            validated_data will be None and error_message will contain the error.
            If validation succeeds, error_message will be None.
        """
        # Check required top-level fields
        missing_field = self._check_required_fields(data)
        if missing_field:
            return None, f"Missing required field: {missing_field}"

        # Validate sequences
        sequences = data.get("sequences", [])
        if not sequences:
            return None, "No sequences provided"

        sequence_error = self._validate_sequences(sequences)
        if sequence_error:
            return None, sequence_error

        # Validate numerical parameters
        max_mutations = self._validate_max_mutations(data)
        if isinstance(max_mutations, str):  # Error message
            return None, max_mutations

        # Validate species
        species_error = self._validate_species(data.get("species", ""))
        if species_error:
            return None, species_error

        # Format validated data
        return self._format_validated_data(data, sequences, max_mutations), None

    def is_valid_dna_sequence(self, sequence: str) -> bool:
        """Check if sequence contains only valid DNA characters."""
        return bool(re.match(r'^[ATGCWSMKRYBDHVN]*$', sequence, re.I))

    def check_frame_and_codons(self, sequence: str) -> Tuple[bool, str, bool]:
        """
        Check if sequence is in frame and detect start/stop codons
        
        Returns:
            Tuple of (is_valid, message, is_advisory)
            Where is_valid is False only if there's a critical error
        """
        if not sequence:
            return True, "", False
            
        sequence = sequence.upper()
        sequence_length = len(sequence)
        
        # Track what we find
        has_start_codon = len(sequence) >= 3 and sequence[0:3] == "ATG"
        stop_codons = ["TAA", "TAG", "TGA"]
        has_stop_codon = len(sequence) >= 3 and sequence[-3:] in stop_codons
        in_frame = sequence_length % 3 == 0
        
        # Create appropriate message
        if not in_frame:
            if has_start_codon and has_stop_codon:
                message = "Sequence not in frame, start and stop codons detected."
                is_advisory = True
            elif has_start_codon:
                message = "Sequence not in frame, start codon detected."
                is_advisory = True
            elif has_stop_codon:
                message = "Sequence not in frame, stop codon detected."
                is_advisory = True
            else:
                # No codons found and not in frame - this is a critical error
                message = "Sequence not in frame (must be divisible by 3)."
                return False, message, False
        else:
            # Sequence is in frame
            if has_start_codon and has_stop_codon:
                message = "Start and stop codons detected."
                is_advisory = True
            elif has_start_codon:
                message = "Start codon detected."
                is_advisory = True
            elif has_stop_codon:
                message = "Stop codon detected."
                is_advisory = True
            else:
                # In frame, no codons - nothing to report
                return True, "", False
                
        # If we get here, we have an advisory message
        return True, message, True

    def _check_required_fields(self, data: Dict[str, Any]) -> Optional[str]:
        """Check if all required fields are present."""
        for field in self.required_fields:
            if field not in data:
                return field
        return None

    def _validate_sequences(self, sequences: List[Dict[str, Any]]) -> Optional[str]:
        """Validate sequence entries."""
        for i, seq_entry in enumerate(sequences):
            # Check required sequence fields
            for field in self.required_sequence_fields:
                if field not in seq_entry:
                    return f"Missing {field} in sequence {i+1}"
            
            # Validate sequence content
            sequence = seq_entry["sequence"].strip()
            if not sequence:
                return f"Empty sequence in entry {i+1}"
            
            # Validate sequence characters
            if not self.is_valid_dna_sequence(sequence):
                return f"Invalid DNA sequence in entry {i+1}"
            
            # Validate MTK part number
            mtk_part = seq_entry["mtkPart"].strip()
            if not mtk_part.isdigit() and not mtk_part in ["", "3a", "3b", "4a", "4aII", "3bI", "3bII", "8a", "8b"]:
                return f"Invalid MTK part number in sequence {i+1}"
            
            # Validate primer name
            primer_name = seq_entry["primerName"].strip()
            if not primer_name:
                return f"Empty primer name in sequence {i+1}"

        return None

    def _validate_max_mutations(self, data: Dict[str, Any]) -> Union[int, str]:
        """Validate max mutations parameter."""
        try:
            max_mutations = int(data.get("max-mut-per-site", 3))
            if max_mutations < 1:
                return "max-mut-per-site must be positive"
            return max_mutations
        except ValueError:
            return "Invalid max-mut-per-site value"

    def _validate_species(self, species: str) -> Optional[str]:
        """Validate species selection."""
        species = species.strip()
        if not species:
            return "Species not specified"
        
        available_species = self.utils.get_available_species()
        if species not in available_species:
            return f"Invalid species: {species}"
        
        return None

    def _format_validated_data(
        self,
        data: Dict[str, Any],
        sequences: List[Dict[str, Any]],
        max_mutations: int
    ) -> Dict[str, Any]:
        """Format validated data for protocol generation."""
        return {
            "seq": [entry["sequence"].strip() for entry in sequences],
            "part_num_left": [entry["mtkPart"].strip() for entry in sequences],
            "primer_name": [entry["primerName"].strip() for entry in sequences],
            "part_num_right": [entry.get("partNumRight", "").strip() for entry in sequences],
            "max_mutations": max_mutations,
            "template_seq": data.get("templateSequence", "").strip() or None,
            "kozak": data.get("kozak", "MTK").strip(),
            "verbose": data.get("verbose_mode", "off") == "on",
            "species": data.get("species", "").strip()
        }