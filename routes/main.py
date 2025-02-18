from flask import Blueprint, render_template, request, jsonify, current_app
from config.logging_config import logger  # ✅ Use centralized logger
from services.utils import GoldenGateUtils
from services.protocol import GoldenGateProtocol
from validators.protocol_validator import ProtocolValidator

main = Blueprint("main", __name__)

utils = GoldenGateUtils()
validator = ProtocolValidator()

MTK_PART_NUMS = ['', '1', '2', '3', '3a', '3b', '3c', '3d', '3e', '3f',
                 '3g', '4', '4a', '4b', '4aII', '3bI', '3bII', '5', '6', '7', '8', '8a', '8b']


@main.route("/")
def home():
    """Render the homepage with default input values from config."""
    testing_mode = current_app.config.get("TESTING", False)    
    test_template_seq = current_app.config.get("TEST_TEMPLATE_SEQ", "")
    print(f"main route => test_template_seq length: {len(test_template_seq)}")
    test_seq = current_app.config.get("TEST_SEQ", [""])
    
    return render_template(
        "index.html",
        title="Home Page",
        num_sequences=len(test_seq) if testing_mode else 1,
        mtk_part_nums=MTK_PART_NUMS,
        test_template_seq=test_template_seq if testing_mode else "",
        testing_mode=testing_mode,
    )


@main.route("/get_sequence_inputs", methods=["GET"])
def get_sequence_inputs():
    """Generate sequence input tabs based on the number of requested sequences."""

    print(f"Received args: {request.args}")

    num_sequences_raw = request.args.get("numSequences", "1")
    try:
        num_sequences = int(num_sequences_raw)
    except ValueError:
        return f"Invalid numSequences value: {num_sequences_raw}", 400

    is_testing = request.args.get("testingMode", "false").lower() == "true"

    # Get test sequences from Flask config, defaulting to an empty list if not set
    test_seq_list = current_app.config.get("TEST_SEQ", [])

    # Ensure the correct number of test sequences are preloaded
    test_sequences = [
        test_seq_list[i] if is_testing and i < len(test_seq_list) else "" 
        for i in range(num_sequences)
    ]

    test_primer_names = [f"Test Primer {i+1}" if is_testing else "" for i in range(num_sequences)]
    test_part_numbers = ["6" if is_testing else "" for _ in range(num_sequences)]
    
    return render_template(
        "sequence_input_tabs.html",
        num_sequences=num_sequences,
        test_primer_name=test_primer_names,
        test_part_number=test_part_numbers,
        test_sequences=test_sequences,
        mtk_part_nums=MTK_PART_NUMS
    )



@main.route("/get_species")
def get_species():
    """Return available species as an HTML option list."""
    try:
        species = utils.get_available_species()
        options_html = "".join(f'<option value="{s}">{s}</option>' for s in species)
        return options_html
    except Exception as e:
        logger.error(f"Error fetching species: {str(e)}", exc_info=True)
        return jsonify({'error': 'Failed to fetch species'}), 500


@main.route("/validate_field", methods=["POST"])
def validate_field():
    """Validate a single form field."""
    data = request.get_json() or request.form
    field_name = data.get("field")
    value = data.get("value", "").strip()
    sequence_index = data.get("sequenceIndex")

    if not field_name:
        return jsonify({'valid': False, 'message': 'Invalid request'})

    validation_functions = {
        "sequence": lambda v: validator._is_valid_dna_sequence(v),
        "mtkPart": lambda v: v.isdigit(),
        "primerName": lambda v: bool(v),
        "species": lambda v: validator._validate_species(v) is None,
        "max-mut-per-site": lambda v: v.isdigit() and int(v) > 0
    }

    valid = validation_functions.get(field_name, lambda v: True)(value)
    message = "Invalid value" if not valid else ""

    return jsonify({'valid': valid, 'message': message, 'field': field_name, 'sequenceIndex': sequence_index})


@main.route("/generate_protocol", methods=["POST"])
def generate_protocol():
    """Generate a Golden Gate protocol from the form data."""
    try:
        logger.info("Starting protocol generation")

        packaged_data, error_message = utils.package_form_data(request)
        if error_message or packaged_data is None:
            return jsonify({'error': error_message or 'No data received'}), 400

        sequences = packaged_data.get("sequences", [])
        if not sequences or not sequences[0].get("sequence"):
            return jsonify({'error': 'No valid sequence data provided'}), 400

        seq = [entry.get("sequence", "").strip() or "TEST_SEQUENCE" for entry in sequences]
        part_num_left = [entry.get("mtkPart", "").strip() or "6" for entry in sequences]
        primer_name = [entry.get("primerName", "").strip() or f"Test Primer {i+1}" for i, entry in enumerate(sequences)]

        max_mutations = int(packaged_data.get("max-mut-per-site", 3))
        template_seq = packaged_data.get("templateSequence", "").strip() or "TEST_TEMPLATE"
        kozak = packaged_data.get("kozak", "MTK")
        species = packaged_data.get("species", "") or "Homo sapiens"
        part_num_right = ["" for _ in sequences]  # Placeholder for right part numbers

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
            verbose=packaged_data.get("verboseMode", False)
        )

        protocol_results = protocol_maker.create_gg_protocol()

        response_data = {
            "status": "success",
            "data": {
                "primers": protocol_results,
                "summary": {
                    "total_sequences": len(seq),
                    "total_primers": len(protocol_results) if protocol_results else 0,
                }
            }
        }

        return render_template("results_partial.html", results=response_data)

    except Exception as e:
        logger.error(f"Error generating protocol: {str(e)}", exc_info=True)
        return jsonify({"error": str(e)}), 500
