# routes/main.py
from flask import Blueprint, render_template, request, jsonify
from tests.test_data import TEST_SEQ, TEST_TEMPLATE_SEQ
from services.utils import get_available_species, package_form_data, get_codon_usage_dict
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
    test_sequence = TEST_SEQ
    
    return render_template("sequence_input_tabs.html",
                           num_sequences=num_sequences,
                           test_primer_name=test_primer_name,
                           test_part_number=test_part_number,
                           test_sequence=test_sequence)

@main.route('/validate_inputs', methods=['POST'])
def validate_inputs():
    # Check Content-Type
    if request.content_type == "application/json":
        data = request.get_json()
    elif request.content_type == "application/x-www-form-urlencoded":
        data = request.form  # HTMX sends form data like this
    else:
        return jsonify({"error": "Unsupported Media Type"}), 415

    return jsonify({"status": "success", "received": data})

@main.route('/generate_protocol', methods=['POST'])
def generate_protocol():
    print('generate_protocol')

    # Ensure we handle both JSON and form-encoded data
    packaged_data, error_message = package_form_data(request)

    print('Packaged Data Keys:', packaged_data.keys())
    
    if error_message:
        return jsonify({'error': error_message}), 400

    # Extract and format required data
    seq = [entry["sequence"] for entry in packaged_data.get("sequences", [])]
    part_num_left = [entry["mtkPart"] for entry in packaged_data.get("sequences", [])]
    primer_name = [entry["primerName"] for entry in packaged_data.get("sequences", [])]

    # Use `.get()` with a fallback value
    max_mutations = int(packaged_data.get("max-mut-per-site", 3))  # Default to 3
    template_seq = packaged_data.get("templateSequence", None)
    kozak = packaged_data.get("kozak", "MTK")
    verbose = packaged_data.get("verboseMode", "off") == "on"
    species = packaged_data.get("species", "")
    
    part_num_right = [entry.get("partNumRight", "") for entry in packaged_data.get("sequences", [])]

    species_codon_usage = get_codon_usage_dict(species)

    protocol = create_gg_protocol(
        seq=seq,
        codon_usage_dict=species_codon_usage,
        part_num_left=part_num_left,
        part_num_right=part_num_right,
        max_mutations=max_mutations,
        primer_name=primer_name,
        template_seq=template_seq,
        kozak=kozak,
        verbose=verbose
    )

    return jsonify({'protocol': protocol})
