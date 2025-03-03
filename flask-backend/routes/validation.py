from flask import Blueprint, request, Response, render_template, make_response
from functools import wraps
import re
from typing import Optional, Tuple, Any
from dataclasses import dataclass
from http import HTTPStatus
import logging
import json

logger = logging.getLogger(__name__)
validation = Blueprint('validation', __name__)

@dataclass
class ValidationResult:
    """Represents the result of a validation check with HTMX response format"""
    is_valid: bool
    is_advisory: bool = False
    message: Optional[str] = None
    field_name: Optional[str] = None
    field_index: Optional[int] = None

    def to_response(self) -> Response:
        """Convert validation result to HTMX-friendly response"""
        # For valid results with no message, return empty response
        if self.is_valid and not self.message:
            response = make_response("")
        else:
            if self.is_valid:
                css_class = "advisory" if self.is_advisory else "valid"
            else:
                css_class = "invalid"
                
            error_html = f'<div class="validation-feedback {css_class}">{self.message}</div>'
            response = make_response(error_html)
        
        # Set response status - 200 for HTMX to process the response
        response.status_code = HTTPStatus.OK
        
        return response


def htmx_validation_handler(f):
    """Decorator to handle HTMX validation requests"""
    @wraps(f)
    def wrapper(*args, **kwargs):
        try:
            # Get form data from request
            form_data = request.form
            
            # Extract field value and name
            value = form_data.get('value', '').strip()
            field = form_data.get('name', '')  # HTMX typically sends field name as 'name'
            
            # Extract sequence index if present (for arrays of fields)
            sequence_index = form_data.get('sequenceIndex', None)
            if sequence_index is not None:
                try:
                    sequence_index = int(sequence_index)
                except (ValueError, TypeError):
                    sequence_index = None
            
            if sequence_index is None and '-' in field and field.split('-')[-1].isdigit():
                parts = field.split('-')
                sequence_index = int(parts[-1])
                field = '-'.join(parts[:-1])
                
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
            logger.exception("Validation error")
            error_html = f'<div class="validation-feedback invalid">Validation error: {str(e)}</div>'
            response = make_response(error_html)
            response.status_code = HTTPStatus.OK
            return response
    
    return wrapper


class SequenceValidator:
    """Handles DNA sequence validation logic - Unchanged from original"""
    
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
    
    @staticmethod
    def check_frame_and_codons(sequence: str) -> Tuple[bool, str, bool]:
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
    
        
    @staticmethod
    def preprocess_sequence_check(sequence: str) -> Tuple[bool, str, bool]:
        """
        Simplified check for sequence preprocessing issues
        Returns tuple of (is_valid, message, is_advisory)
        """
        if not sequence:
            return True, "", False
            
        sequence = sequence.upper()
        sequence_length = len(sequence)
        
        # Check if in frame
        in_frame = sequence_length % 3 == 0
        if not in_frame:
            return False, "Sequence length must be divisible by 3 (in frame).", False
            
        # Check for start codon
        has_start_codon = len(sequence) >= 3 and sequence[0:3] == "ATG"
        
        # Check for stop codons
        stop_codons = ["TAA", "TAG", "TGA"]
        has_stop_codon = len(sequence) >= 3 and sequence[-3:] in stop_codons
        
        # Provide appropriate advisory messages
        if has_start_codon and has_stop_codon:
            return True, "Start and stop codons detected and will be removed during processing.", True
        elif has_start_codon:
            return True, "Start codon detected and will be removed during processing.", True
        elif has_stop_codon:
            return True, "Stop codon detected and will be removed during processing.", True
            
        return True, "", False

# Validation Routes
@validation.route('/validate/sequence', methods=['POST'])
def validate_sequence():
    """Direct HTMX validation with support for multiple validation messages"""
    try:
        form_data = request.form
        # Extract the index 
        sequence_index = form_data.get('sequenceIndex')
        try:
            sequence_index = int(sequence_index) if sequence_index else 0
        except ValueError:
            sequence_index = 0
            
        key = f"sequences[{sequence_index}][sequence]"
        value = form_data.get(key, '').strip()
        
        # Collect all validation messages
        validation_messages = []
        is_valid = True
        
        # Check for empty value
        if not value:
            validation_messages.append('<span class="validation-feedback invalid">Sequence cannot be empty</span>')
            is_valid = False
        else:  
            # Check for valid DNA sequence
            if not SequenceValidator.is_valid_dna_sequence(value):
                validation_messages.append('<span class="validation-feedback invalid">Only valid DNA bases (A T G C) or extended DNA code (W S M K R Y B D H V N) allowed</span>')
                is_valid = False
                
            # Check for sequence length
            if not SequenceValidator.validate_sequence_length(value):
                validation_messages.append('<span class="validation-feedback invalid">Sequence must be at least 80 bp</span>')
                is_valid = False
                
            # Additional validation for domestication sequences
            sequence_type = form_data.get('sequence_type', 'domestication')
            if sequence_type == 'domestication':
                check_valid, message, is_advisory = SequenceValidator.preprocess_sequence_check(value)
                
                if message:
                    css_class = "advisory" if is_advisory else "invalid"
                    validation_messages.append(f'<span class="validation-feedback {css_class}">{message}</span>')
                    if not check_valid:
                        is_valid = False
        
        # If no validation issues, add success message
        if is_valid:
            validation_messages.append('<span class="validation-feedback valid">Valid sequence</span>')
        
        # Create a response that includes validation messages and sets HTMX headers for form state
        response = make_response(''.join(validation_messages))
        
        # Add HTMX headers to control the submit button state
        if is_valid:
            response.headers['HX-Trigger'] = json.dumps({
                'sequenceValidChanged': {
                    'isValid': True,
                    'index': sequence_index
                }
            })
        else:
            response.headers['HX-Trigger'] = json.dumps({
                'sequenceValidChanged': {
                    'isValid': False, 
                    'index': sequence_index
                }
            })
            
        print(f"Returning feedback: {validation_messages} with is_valid={is_valid}")
        return response
        
    except Exception as e:
        print(f"Validation error: {str(e)}")
        error_response = make_response(f'<span class="validation-feedback invalid">Validation error: {str(e)}</span>')
        error_response.headers['HX-Trigger'] = json.dumps({
            'sequenceValidChanged': {
                'isValid': False,
                'index': sequence_index if 'sequence_index' in locals() else 0
            }
        })
        return error_response

@validation.route('/validate/template', methods=['POST'])
@htmx_validation_handler
def validate_template(value: str, field: str, sequence_index: Optional[int]) -> ValidationResult:
    """Validate template sequence"""
    # Skip validation if empty (templates are optional)
    if not value:
        return ValidationResult(True, "", False, field, sequence_index)
        
    # Check for valid DNA bases
    if not SequenceValidator.is_valid_dna_sequence(value):
        return ValidationResult(
            False,
            "Only valid DNA bases (A, T, G, C, W, S, M, K, R, Y, B, D, H, V, N) allowed",
            False,
            field,
            sequence_index
        )

    # Check minimum length
    if not SequenceValidator.validate_sequence_length(value):
        return ValidationResult(False, "Template must be at least 80 bp", False, field, sequence_index)
        
    return ValidationResult(True, "Valid template sequence", False, field, sequence_index)

@validation.route('/validate/amplicon', methods=['POST'])
@htmx_validation_handler
def validate_amplicon(value: str, field: str, sequence_index: Optional[int]) -> ValidationResult:
    """Validate that the amplicon is contained within the template"""
    if not value:
        return ValidationResult(False, "Amplicon sequence cannot be empty", False, field, sequence_index)
        
    # Get the template value from the form data
    template = request.form.get('template', '')
    
    # First validate the amplicon as a DNA sequence
    if not SequenceValidator.is_valid_dna_sequence(value):
        return ValidationResult(
            False,
            "Only valid DNA bases (A, T, G, C, W, S, M, K, R, Y, B, D, H, V, N) allowed",
            False,
            field,
            sequence_index
        )
        
    # If we have a template, validate Athat the amplicon is contained within it
    if template and not SequenceValidator.is_amplicon_valid(value, template):
        return ValidationResult(
            False,
            "Amplicon sequence must be contained within the template sequence",
            False,
            field,
            sequence_index
        )
        
    return ValidationResult(True, "Valid amplicon", False, field, sequence_index)

@validation.route('/validate/mtk-part', methods=['POST'])
@htmx_validation_handler
def validate_mtk_part(value: str, field: str, sequence_index: Optional[int]) -> ValidationResult:
    """Validate MTK part number"""
    if not value:
        return ValidationResult(False, "MTK Part number is required", field, sequence_index)
    
    if not value.isdigit():
        return ValidationResult(False, "MTK part must be a number", field, sequence_index)
    
    return ValidationResult(True, "Valid MTK part number", field, sequence_index)

@validation.route('/validate/primer-name', methods=['POST'])
@htmx_validation_handler
def validate_primer_name(value: str, field: str, sequence_index: Optional[int]) -> ValidationResult:
    """Validate primer name"""
    if not value:
        return ValidationResult(False, "Primer name is required", field, sequence_index)
        
    return ValidationResult(True, "Valid primer name", field, sequence_index)

@validation.route('/validate/species', methods=['POST'])
@htmx_validation_handler
def validate_species(value: str, field: str, sequence_index: Optional[int]) -> ValidationResult:
    """Validate species input"""
    if not value.strip():
        return ValidationResult(False, "Species cannot be empty", field, sequence_index)
    
    return ValidationResult(True, "Valid species name", field, sequence_index)

@validation.route('/validate/max-mut-per-site', methods=['POST'])
@htmx_validation_handler
def validate_max_mutations(value: str, field: str, sequence_index: Optional[int]) -> ValidationResult:
    """Validate maximum mutations per site"""
    try:
        max_mut = int(value)
        if max_mut <= 0:
            return ValidationResult(False, "Must be a positive number", field, sequence_index)
        return ValidationResult(True, "Valid value", field, sequence_index)
    except ValueError:
        return ValidationResult(False, "Must be a number", field, sequence_index)

@validation.route('/validate/form', methods=['POST'])
def validate_form():
    """Validate the entire form and return HTML response"""
    try:
        # Process form data
        form_data = request.form
        
        # Validate all fields (you would implement comprehensive validation here)
        # For now, we'll return a success message
        
        # Determine if this is an HTMX request that requires partial response
        is_htmx_request = request.headers.get('HX-Request') == 'true'
        
        if is_htmx_request:
            # Return a fragment for HTMX to update a specific part of the page
            return '<div class="alert alert-success">Form validation passed!</div>'
        else:
            # Full page redirect for non-HTMX form submission
            return render_template('success.html', message="Form validation passed")
        
    except Exception as e:
        error_message = f"Form validation error: {str(e)}"
        
        # Check if HTMX request
        is_htmx_request = request.headers.get('HX-Request') == 'true'
        
        if is_htmx_request:
            return f'<div class="alert alert-danger">{error_message}</div>'
        else:
            return render_template('error.html', message=error_message)

# Error Handlers
@validation.errorhandler(Exception)
def handle_exception(e):
    """Global error handler for validation routes"""
    error_message = f"Validation error: {str(e)}"
    
    # Check if HTMX request
    is_htmx_request = request.headers.get('HX-Request') == 'true'
    
    if is_htmx_request:
        return f'<div class="alert alert-danger">{error_message}</div>', HTTPStatus.OK
    else:
        return render_template('error.html', message=error_message), HTTPStatus.INTERNAL_SERVER_ERROR