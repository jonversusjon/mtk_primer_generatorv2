# main.py

from flask import Blueprint, render_template, request, jsonify, current_app
from config.logging_config import logger
from services.utils import GoldenGateUtils
from services.protocol import GoldenGateProtocol
from validators.protocol_validator import ProtocolValidator
from collections import defaultdict
from Bio.Seq import Seq
import numpy as np
import json

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
    test_seq = current_app.config.get("TEST_SEQ", [])
    available_species = utils.get_available_species()

    # When testing, show as many tabs as there are test sequences; otherwise default to 1.
    num_sequences = len(test_seq) if testing_mode and test_seq else 1

    return render_template(
        "index.html",
        title="Home Page",
        num_sequences=num_sequences,
        mtk_part_nums=MTK_PART_NUMS,
        test_template_seq=test_template_seq if testing_mode else "",
        testing_mode=testing_mode,
        test_seq=test_seq,
        available_species=available_species,
    )


@main.route('/update-sequence-count', methods=['POST'])
def update_sequence_count():
    direction = request.form.get('direction')
    current = int(request.form.get('current', 1))

    if direction == 'increase' and current < 10:
        return str(current + 1)
    elif direction == 'decrease' and current > 1:
        return str(current - 1)

    return str(current)


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
    test_seq_list = current_app.config.get("TEST_SEQ", [])

    # Ensure the correct number of test sequences are preloaded
    test_sequences = [
        test_seq_list[i] if is_testing and i < len(test_seq_list) else ""
        for i in range(num_sequences)
    ]

    test_primer_names = [
        f"Test Primer {i+1}" if is_testing else "" for i in range(num_sequences)]
    test_part_numbers = [
        "6" if is_testing else "" for _ in range(num_sequences)]

    return render_template(
        "sequence_input_tabs.html",
        num_sequences=num_sequences,
        test_primer_name=test_primer_names,
        test_part_number=test_part_numbers,
        test_sequences=test_sequences,
        mtk_part_nums=MTK_PART_NUMS
    )


def process_restriction_sites(sequence_analysis):
    """Prepare restriction site summary for Jinja template rendering."""
    grouped_sites = defaultdict(list)

    for seq in sequence_analysis:
        for site in seq.get("restriction_sites", []):
            enzyme = site["enzyme"]
            grouped_sites[enzyme].append(site["position"])

    # Convert to a format Jinja can easily loop over
    processed_sites = [
        {
            "enzyme": enzyme,
            "count": len(positions),
            "positions": sorted(positions)
        }
        for enzyme, positions in grouped_sites.items()
    ]

    return processed_sites


def convert_non_serializable(obj):
    """Convert non-serializable objects (Seq, ndarray) to JSON-compatible types."""
    if isinstance(obj, Seq):
        return str(obj)  # Convert BioPython Seq to string
    elif isinstance(obj, np.ndarray):
        return obj.tolist()  # Convert NumPy arrays to lists
    elif isinstance(obj, dict):
        return {k: convert_non_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_non_serializable(v) for v in obj]
    return obj  # Return as-is if it's already JSON serializable


@main.route("/generate_protocol", methods=["POST"])
def generate_protocol():
    print("⭐ generate_protocol route triggered")

    # Get basic form data
    species = request.form.get("species", "")
    kozak = request.form.get("kozak", "MTK")
    max_mut_per_site = int(request.form.get("max_mut_per_site", 3))
    verbose_mode = "verbose_mode" in request.form
    template_sequence = request.form.get(
        "templateSequence", "").strip() or "TEST_TEMPLATE"

    # Extract sequences data
    num_sequences = int(request.form.get("numSequences", 0))
    sequences = []
    primer_names = []
    mtk_parts = []

    for i in range(num_sequences):
        seq_key = f"sequences[{i}][sequence]"
        name_key = f"sequences[{i}][primerName]"
        part_key = f"sequences[{i}][mtkPart]"

        sequence = request.form.get(seq_key, "").strip()
        if sequence:  # Only process non-empty sequences
            sequences.append(sequence)
            primer_names.append(request.form.get(
                name_key, "").strip() or f"Test Primer {i+1}")
            mtk_parts.append(request.form.get(part_key, "").strip() or "6")

    # Check if we have valid sequences
    if not sequences:
        return jsonify({
            "error": "No valid sequences provided",
            "form_data": request.form
        }), 422

    try:
        # Create protocol
        protocol_maker = GoldenGateProtocol(
            seq=sequences,
            codon_usage_dict=utils.get_codon_usage_dict(species),
            part_num_left=mtk_parts,
            part_num_right=["" for _ in sequences],
            max_mutations=max_mut_per_site,
            primer_name=primer_names,
            template_seq=template_sequence,
            kozak=kozak,
            verbose=verbose_mode
        )

        # Generate the protocol
        result = protocol_maker.create_gg_protocol()
        print("Restriction sites in result:", result.get('restriction_sites'))

        # Check for errors in the result
        if result.get('has_errors', False):
            return jsonify({
                "error": "Protocol generation error",
                "sequence_errors": result.get('sequence_errors', {}),
                "partial_result": result,
                "form_data": request.form
            }), 422

        return jsonify({"success": True, "data": result})
    except Exception as e:
        error_response = {
            "error": "An error occurred in generate_protocol",
            "details": str(e)
        }
        return jsonify(error_response), 500
