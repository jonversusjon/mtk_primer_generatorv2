from flask import Blueprint, request, Response
from functools import wraps
from validators import ProtocolValidator
import re
from typing import Optional, Tuple, Any
from dataclasses import dataclass
from http import HTTPStatus

validation = Blueprint('validation', __name__)
validator = ProtocolValidator()

@dataclass
class ValidationResult:
    """Represents the result of a validation check"""
    is_valid: bool
    message: Optional[str] = None
    status_code: int = HTTPStatus.OK

    def to_response(self) -> Response:
        """Convert validation result to Flask response"""
        if self.is_valid:
            return Response('', status=HTTPStatus.OK)
        return Response(
            self.message or 'Validation failed',
            status=self.status_code or HTTPStatus.BAD_REQUEST
        )

def validation_handler(f):
    """Decorator to handle common validation patterns"""
    @wraps(f)
    def wrapper(*args, **kwargs):
        try:
            value = request.form.get('value', '').strip()
            result = f(value, *args, **kwargs)
            if isinstance(result, ValidationResult):
                return result.to_response()
            return Response('', HTTPStatus.OK)
        except Exception as e:
            return Response(
                str(e),
                status=HTTPStatus.INTERNAL_SERVER_ERROR
            )
    return wrapper

class SequenceValidator:
    """Handles DNA sequence validation logic"""
    
    @staticmethod
    def is_valid_dna_sequence(sequence: str) -> bool:
        """Check if sequence contains only valid DNA bases"""
        return bool(re.match(r'^[ATGCWSMKRYBDHVN]*$', sequence, re.I))

    @staticmethod
    def is_amplicon_valid(amplicon: str, template: str) -> bool:
        """Validate amplicon against template sequence"""
        if not template:
            return True
        return len(amplicon) <= len(template) and amplicon.upper() in template.upper()

    @staticmethod
    def validate_sequence_length(sequence: str, min_length: int = 80) -> bool:
        """Validate sequence meets minimum length requirement"""
        return len(sequence) >= min_length

# Validation Routes
@validation.route('/validate/sequence', methods=['POST'])
@validation_handler
def validate_sequence(sequence: str) -> ValidationResult:
    """Validate DNA sequence input"""
    template = request.form.get('template', '').strip()

    if not sequence:
        return ValidationResult(False, "Sequence cannot be empty")

    if not SequenceValidator.is_valid_dna_sequence(sequence):
        return ValidationResult(
            False,
            "Only valid DNA bases (A, T, G, C, W, S, M, K, R, Y, B, D, H, V, N) allowed"
        )

    if not SequenceValidator.validate_sequence_length(sequence):
        return ValidationResult(False, "Sequence must be at least 80 bp")

    if template and not SequenceValidator.is_amplicon_valid(sequence, template):
        return ValidationResult(False, "Sequence must be within the template")

    return ValidationResult(True)

@validation.route('/validate/mtk-part', methods=['POST'])
@validation_handler
def validate_mtk_part(part_number: str) -> ValidationResult:
    """Validate MTK part number"""
    if not part_number:
        return ValidationResult(False, "MTK Part number is required")
    
    if not part_number.isdigit():
        return ValidationResult(False, "MTK part must be a number")
    
    return ValidationResult(True)

@validation.route('/validate/primer-name', methods=['POST'])
@validation_handler
def validate_primer_name(name: str) -> ValidationResult:
    """Validate primer name"""
    # Add custom primer name validation if needed
    return ValidationResult(True)

@validation.route('/validate/species', methods=['POST'])
@validation_handler
def validate_species(species: str) -> ValidationResult:
    """Validate species input"""
    error = validator._validate_species(species)
    if error:
        return ValidationResult(False, error)
    return ValidationResult(True)

@validation.route('/validate/max-mut-per-site', methods=['POST'])
@validation_handler
def validate_max_mutations(value: str) -> ValidationResult:
    """Validate maximum mutations per site"""
    try:
        max_mut = int(value)
        if max_mut <= 0:
            return ValidationResult(False, "Must be a positive number")
        return ValidationResult(True)
    except ValueError:
        return ValidationResult(False, "Must be a number")

# Error Handlers
@validation.errorhandler(Exception)
def handle_exception(e):
    """Global error handler for validation routes"""
    return Response(
        str(e),
        status=HTTPStatus.INTERNAL_SERVER_ERROR
    )