from flask import Blueprint, request, jsonify, current_app
from config.logging_config import logger
from services.utils import GoldenGateUtils
from services.protocol import GoldenGateProtocol
from validators.protocol_validator import ProtocolValidator
from flask_cors import CORS
import json

api = Blueprint("api", __name__, url_prefix="/api")
CORS(api)  # Enable CORS for all routes in this blueprint

utils = GoldenGateUtils()
validator = ProtocolValidator()


@api.route("/species", methods=["GET"])
def get_species():
    print("get species route called")
    """Return available species as JSON"""
    try:
        species = utils.get_available_species()
        print(f"Available species: {species}")
        return jsonify({"species": species})
    except Exception as e:
        logger.error(f"Error fetching species: {str(e)}", exc_info=True)
        return jsonify({'error': 'Failed to fetch species'}), 500


@api.route("/validate/sequence", methods=["POST"])
def validate_sequence():
    """Validate a DNA sequence"""
    try:
        data = request.get_json()
        sequence = data.get("sequence", "").strip()

        if not sequence:
            return jsonify({
                "isValid": False,
                "message": "Sequence cannot be empty",
                "isAdvisory": False
            })

        # Check for valid DNA bases
        if not validator.is_valid_dna_sequence(sequence):
            return jsonify({
                "isValid": False,
                "message": "Only valid DNA bases (A, T, G, C) or extended DNA code allowed",
                "isAdvisory": False
            })

        # Check for sequence length
        if len(sequence) < 80:
            return jsonify({
                "isValid": False,
                "message": "Sequence must be at least 80 bp",
                "isAdvisory": False
            })

        # Check if sequence is in frame
        in_frame = len(sequence) % 3 == 0
        if not in_frame:
            return jsonify({
                "isValid": False,
                "message": "Sequence length must be divisible by 3 (in frame)",
                "isAdvisory": False
            })

        # Check for start/stop codons
        valid, message, is_advisory = validator.check_frame_and_codons(
            sequence)

        return jsonify({
            "isValid": valid,
            "message": message,
            "isAdvisory": is_advisory
        })

    except Exception as e:
        logger.error(f"Validation error: {str(e)}", exc_info=True)
        return jsonify({
            "isValid": False,
            "message": f"Server error: {str(e)}",
            "isAdvisory": False
        }), 500


@api.route("/protocol", methods=["POST"])
def generate_protocol():
    """Generate a Golden Gate protocol based on the provided sequences"""
    try:
        # Handle form data submitted by React
        form_data = request.form.to_dict()

        # Process sequences data - form data is flattened so we need to restructure it
        sequences = []
        num_sequences = int(form_data.get("numSequences", 0))

        for i in range(num_sequences):
            seq_key = f"sequences[{i}][sequence]"
            name_key = f"sequences[{i}][primerName]"
            part_key = f"sequences[{i}][mtkPart]"

            if seq_key in form_data:
                sequences.append({
                    "sequence": form_data.get(seq_key, "").strip(),
                    "primerName": form_data.get(name_key, "").strip(),
                    "mtkPart": form_data.get(part_key, "").strip()
                })

        # Extract other form fields
        species = form_data.get("species", "")
        kozak = form_data.get("kozak", "MTK")
        max_mut_per_site = int(form_data.get("max_mut_per_site", 3))
        verbose_mode = "verbose_mode" in form_data
        template_sequence = form_data.get("templateSequence", "").strip()

        # Validate input
        if not species:
            return jsonify({"error": "Species not specified"}), 400

        if not sequences:
            return jsonify({"error": "No sequences provided"}), 400

        # Extract list data for protocol generation
        seq_list = [seq["sequence"] for seq in sequences]
        primer_names = [seq["primerName"] for seq in sequences]
        mtk_parts = [seq["mtkPart"] for seq in sequences]

        # Create protocol
        protocol_maker = GoldenGateProtocol(
            seq=seq_list,
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

        # Process any non-JSON serializable objects in the result
        result = json.loads(json.dumps(result, default=lambda o: str(o)))

        # Check for errors
        if result.get('has_errors', False):
            return jsonify({
                "error": "Error in protocol generation",
                "sequence_errors": result.get('sequence_errors', {})
            }), 422

        # Success - return the result
        return jsonify(result)

    except Exception as e:
        logger.error(f"Error in generate_protocol: {str(e)}", exc_info=True)
        return jsonify({"error": str(e)}), 500


@api.route("/export", methods=["POST"])
def export_protocol():
    """Export protocol results as a TSV file"""
    try:
        data = request.get_json()
        primers = data.get("primers", [])

        if not primers:
            return jsonify({"error": "No primers to export"}), 400

        # Generate a unique filename
        filename = f"primers_{utils.generate_unique_id()}.tsv"
        filepath = f"static/exports/{filename}"

        # Export primers to TSV
        with open(filepath, "w") as f:
            f.write("Primer Name\tSequence\tAmplicon\n")
            for primer in primers:
                f.write(f"{primer[0]}\t{primer[1]}\t{primer[2]}\n")

        # Return the download URL
        download_url = f"/static/exports/{filename}"
        return jsonify({"download_url": download_url})

    except Exception as e:
        logger.error(f"Error exporting protocol: {str(e)}", exc_info=True)
        return jsonify({"error": str(e)}), 500
