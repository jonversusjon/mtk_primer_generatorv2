# routes/main.py
from flask import Blueprint, render_template, request, jsonify
from tests.test_data import TEST_SEQ, TEST_TEMPLATE_SEQ
from utils import get_available_species, package_form_data
from services.protocol import create_gg_protocol

# Renamed blueprint variable for consistency:
main = Blueprint("main", __name__)

TESTING_MODE = True

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
    species = get_available_species()
    options_html = "".join(
        f'<option value="{s}">{s}</option>' for s in species)
    return options_html

@main.route("/get_sequence_inputs", methods=["GET"])
def get_sequence_inputs():
    num_sequences = int(request.args.get("numSequences", 1))
    # Optionally, if in testing mode, prepare test values:
    test_primer_name = ["Test Primer 1"] if num_sequences > 0 else []
    test_part_number = ["6"] if num_sequences > 0 else []
    test_sequence = TEST_SEQ  # Replace with your test sequence if desired
    
    return render_template("sequence_input_tabs.html",
                           num_sequences=num_sequences,
                           test_primer_name=test_primer_name,
                           test_part_number=test_part_number,
                           test_sequence=test_sequence)

@main.route('/validate_inputs', methods=['POST'])
def validate_inputs():
    data = request.get_json() or {}
    # Implement your validation logic here.
    # For example, if an error is found:
    #    return jsonify({'error': 'Validation error message'}), 400
    return jsonify({'message': 'Validation successful'})

@main.route('/generate_protocol', methods=['POST'])
def generate_protocol():
    packaged_data, error_message = package_form_data(request)
    if error_message:
        return jsonify({'error': error_message}), 400

    protocol = create_gg_protocol(packaged_data)
    return jsonify({'protocol': protocol})
