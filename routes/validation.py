# routes/validation.py
from flask import Blueprint, request, jsonify
from validators import ProtocolValidator

validation = Blueprint('validation', __name__)
validator = ProtocolValidator()

@validation.route('/validate-field', methods=['POST'])
def validate_field():
    """Validate a single form field in real-time."""
    data = request.json
    field_name = data.get('field')
    value = data.get('value')
    sequence_index = data.get('sequenceIndex')

    if not field_name or value is None:
        return jsonify({
            'valid': False,
            'message': 'Invalid request'
        })

    # Validate based on field type
    if field_name == 'sequence':
        valid = validator._is_valid_dna_sequence(value)
        message = 'Invalid DNA sequence' if not valid else ''
    
    elif field_name == 'mtkPart':
        valid = value.strip().isdigit()
        message = 'MTK part must be a number' if not valid else ''
    
    elif field_name == 'primerName':
        valid = bool(value.strip())
        message = 'Primer name is required' if not valid else ''
    
    elif field_name == 'species':
        error = validator._validate_species(value)
        valid = error is None
        message = error if error else ''
    
    elif field_name == 'max-mut-per-site':
        try:
            max_mut = int(value)
            valid = max_mut > 0
            message = 'Must be a positive number' if not valid else ''
        except ValueError:
            valid = False
            message = 'Must be a number'
    
    else:
        valid = True
        message = ''

    return jsonify({
        'valid': valid,
        'message': message,
        'field': field_name,
        'sequenceIndex': sequence_index
    })