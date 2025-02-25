from flask import Blueprint, request, Response, jsonify
from functools import wraps
import re
from typing import Optional, Tuple, Any
from dataclasses import dataclass
from http import HTTPStatus
import json

validation = Blueprint('validation', __name__)


@dataclass
class ValidationResult:
    """Represents the result of a validation check"""
    is_valid: bool
    message: Optional[str] = None
    status_code: int = HTTPStatus.OK
    field_name: Optional[str] = None
    field_index: Optional[int] = None

    def to_response(self) -> Response:
        """Convert validation result to Flask response"""
        response_data = {
            "valid": self.is_valid,
            "message": self.message or '',
            "field": self.field_name,
            "index": self.field_index
        }

        response = jsonify(response_data)
        # Always return 200 for validation responses
        response.status_code = HTTPStatus.OK

        return response


def validation_handler(f):
    """Decorator to handle common validation patterns"""
    @wraps(f)
    def wrapper(*args, **kwargs):
        try:
            # Get common parameters
            data = request.json if request.is_json else request.form
            value = data.get('value', '').strip()
            field = data.get('field', '')
            sequence_index = data.get('sequenceIndex', None)
            
            # Convert sequenceIndex to int if it exists
            if sequence_index is not None:
                try:
                    sequence_index = int(sequence_index)
                except (ValueError, TypeError):
                    sequence_index = None
            
            # Call the validation function
            result = f(value, field, sequence_index, *args, **kwargs)
            
            # Ensure we have a ValidationResult
            if not isinstance(result, ValidationResult):
                result = ValidationResult(True)
                
            # Add field information
            if not result.field_name:
                result.field_name = field
            if result.field_index is None:
                result.field_index = sequence_index
                
            return result.to_response()
            
        except Exception as e:
            # Return a formatted error response
            error_result = ValidationResult(
                is_valid=False,
                message=f"Validation error: {str(e)}",
                status_code=HTTPStatus.INTERNAL_SERVER_ERROR
            )
            return error_result.to_response()
    
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
def validate_sequence(value: str, field: str, sequence_index: Optional[int]) -> ValidationResult:
    """Validate DNA sequence input"""
    if not value:
        return ValidationResult(False, "Sequence cannot be empty")

    if not SequenceValidator.is_valid_dna_sequence(value):
        return ValidationResult(
            False,
            "Only valid DNA bases (A, T, G, C, W, S, M, K, R, Y, B, D, H, V, N) allowed"
        )

    if not SequenceValidator.validate_sequence_length(value):
        return ValidationResult(False, "Sequence must be at least 80 bp")
        
    # Template validation would go here if needed
    # We could fetch the template sequence from the form data if required

    return ValidationResult(True)

@validation.route('/validate/mtk-part', methods=['POST'])
@validation_handler
def validate_mtk_part(value: str, field: str, sequence_index: Optional[int]) -> ValidationResult:
    """Validate MTK part number"""
    if not value:
        return ValidationResult(False, "MTK Part number is required")
    
    if not value.isdigit():
        return ValidationResult(False, "MTK part must be a number")
    
    return ValidationResult(True)

@validation.route('/validate/primer-name', methods=['POST'])
@validation_handler
def validate_primer_name(value: str, field: str, sequence_index: Optional[int]) -> ValidationResult:
    """Validate primer name"""
    # Add custom primer name validation if needed
    if not value:
        return ValidationResult(False, "Primer name is required")
        
    # Add more validation logic if needed
    
    return ValidationResult(True)

@validation.route('/validate/species', methods=['POST'])
@validation_handler
def validate_species(value: str, field: str, sequence_index: Optional[int]) -> ValidationResult:
    """Validate species input"""
    if not value.strip():
        return ValidationResult(False, "Species cannot be empty")
    # Optionally, add more validation logic here if needed
    return ValidationResult(True)


@validation.route('/validate/max-mut-per-site', methods=['POST'])
@validation_handler
def validate_max_mutations(value: str, field: str, sequence_index: Optional[int]) -> ValidationResult:
    """Validate maximum mutations per site"""
    try:
        max_mut = int(value)
        if max_mut <= 0:
            return ValidationResult(False, "Must be a positive number")
        return ValidationResult(True)
    except ValueError:
        return ValidationResult(False, "Must be a number")

@validation.route('validate/form', methods=['POST'])
def validate_form():
    """Validate the entire form and return validation status"""
    try:
        # Get form data
        data = request.json if request.is_json else request.form
        
        # You could implement comprehensive form validation here
        # For now, we'll return a simple success response
        
        return jsonify({
            "valid": True,
            "message": "Form validation passed"
        })
        
    except Exception as e:
        return jsonify({
            "valid": False,
            "message": f"Form validation error: {str(e)}"
        })

# Error Handlers
@validation.errorhandler(Exception)
def handle_exception(e):
    """Global error handler for validation routes"""
    return jsonify({
        "valid": False,
        "message": f"Validation error: {str(e)}"
    }), HTTPStatus.INTERNAL_SERVER_ERROR