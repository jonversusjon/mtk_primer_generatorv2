# routes/main.py
from flask import Blueprint, render_template, request, jsonify
from tests.test_data import TEST_SEQ, TEST_TEMPLATE_SEQ
from services.utils import GoldenGateUtils
from services.protocol import GoldenGateProtocol
from validators.protocol_validator import ProtocolValidator

main = Blueprint("main", __name__)

TESTING_MODE = True
utils = GoldenGateUtils()
validator = ProtocolValidator()

# protocol_maker = GoldenGateProtocol()

@main.route("/")
def home():
    return render_template(
        "index.html",
        title="Home Page",
        test_seq=TEST_SEQ,
        test_template_seq=TEST_TEMPLATE_SEQ,
        testing_mode=TESTING_MODE
    )

@main.route("/get_species")
def get_species():
    species = utils.get_available_species()
    options_html = "".join(
        f'<option value="{s}">{s}</option>' for s in species)
    return options_html

@main.route("/get_sequence_inputs", methods=["GET"])
def get_sequence_inputs():
    num_sequences = int(request.args.get("numSequences", 1))
    
    # Optionally, if in testing mode, prepare test values:
    test_primer_name = ["Test Primer 1"] if num_sequences > 0 else []
    test_part_number = ["6"] if num_sequences > 0 else []
    test_sequence = TEST_SEQ
    
    return render_template("sequence_input_tabs.html",
                           num_sequences=num_sequences,
                           test_primer_name=test_primer_name,
                           test_part_number=test_part_number,
                           test_sequence=test_sequence)

@main.route('/validate_field', methods=['POST'])
def validate_field():
    """Validate a single form field."""
    if request.content_type == "application/json":
        data = request.get_json()
    else:
        data = request.form

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


@main.route('/generate_protocol', methods=['POST'])
def generate_protocol():
    # Ensure we handle both JSON and form-encoded data
    packaged_data, error_message = utils.package_form_data(request)
    verbose = packaged_data.get("verboseMode")
    
    print(f"main route verbose: {verbose}")
    if error_message:
        return jsonify({'error': error_message}), 400

    # Extract and format required data
    seq = [entry["sequence"] for entry in packaged_data.get("sequences", [])]
    part_num_left = [entry["mtkPart"] for entry in packaged_data.get("sequences", [])]
    primer_name = [entry["primerName"] for entry in packaged_data.get("sequences", [])]

    # Use `.get()` with a fallback value
    max_mutations = int(packaged_data.get("max-mut-per-site", 3))
    template_seq = packaged_data.get("templateSequence", None)
    kozak = packaged_data.get("kozak", "MTK")
    
    species = packaged_data.get("species", "")
    
    part_num_right = [entry.get("partNumRight", "") for entry in packaged_data.get("sequences", [])]

    species_codon_usage = utils.get_codon_usage_dict(species)
    protocol_maker = GoldenGateProtocol(
        seq=seq,
        codon_usage_dict=species_codon_usage,
        part_num_left=part_num_left,
        part_num_right=part_num_right,
        max_mutations=max_mutations,
        primer_name=primer_name,
        template_seq=template_seq,
        kozak=kozak,
        verbose=verbose)
    
    protocol = protocol_maker.create_gg_protocol()

    return jsonify({'protocol': protocol})
