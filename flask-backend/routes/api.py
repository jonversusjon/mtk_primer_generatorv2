from flask import Blueprint, request, jsonify, current_app
from config.logging_config import logger
from services.utils import GoldenGateUtils
from services.protocol import GoldenGateProtocol
from flask_cors import CORS
import json

api = Blueprint("api", __name__, url_prefix="/api")
CORS(api)  # Enable CORS for all routes in this blueprint

utils = GoldenGateUtils()


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


@api.route("/generate_protocol", methods=["POST"])
def generate_protocol():
    """Generate a Golden Gate protocol based on the provided sequences"""
    try:
        # Parse JSON data from request
        data = request.json
        print(f"Form data received: {data}")

        # Extract main parameters
        sequencesToDomesticate = data.get("sequencesToDomesticate", [])
        species = data.get("species", "")
        kozak = data.get("kozak", "MTK")
        max_mut_per_site = data.get("max_mut_per_site", 3)
        verbose_mode = data.get("verbose_mode", True)
        template_sequence = data.get("templateSequence", "")

        # Validate input
        if not species:
            return jsonify({"error": "Species not specified"}), 400

        if not sequencesToDomesticate:
            return jsonify({"error": "No sequences provided"}), 400

        # Create protocol with the new parameter structure
        protocol_maker = GoldenGateProtocol(
            sequencesToDomesticate=sequencesToDomesticate,
            codon_usage_dict=utils.get_codon_usage_dict(species),
            max_mutations=max_mut_per_site,
            template_seq=template_sequence,
            kozak=kozak,
            verbose=verbose_mode
        )

        # Generate the protocol
        result = protocol_maker.create_gg_protocol()
        serializable_result = utils.convert_non_serializable(result)
        print(f"Generated protocol: {serializable_result}")

        # Check for errors
        if serializable_result.get('has_errors', False):
            return jsonify({
                "error": "Error in protocol generation",
                "sequence_errors": serializable_result.get('sequence_errors', {})
            }), 422

        # Success - return the result
        return jsonify(serializable_result)

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
