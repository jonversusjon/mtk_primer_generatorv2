from flask import Blueprint, request, jsonify, current_app
from config.logging_config import logger
from services.utils import GoldenGateUtils
from services.protocol import GoldenGateProtocol
from flask_cors import CORS
import json
import logging

api = Blueprint("api", __name__, url_prefix="/api")
CORS(api)  # Enable CORS for all routes in this blueprint

utils = GoldenGateUtils()

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def get_nested_keys(item, prefix=''):
    """
    Recursively collect keys from dictionaries and lists.

    :param item: The dictionary or list to inspect.
    :param prefix: The prefix to use for nested keys.
    :return: A list of string keys, including paths.
    """
    keys = []
    if isinstance(item, dict):
        for k, v in item.items():
            full_key = f"{prefix}.{k}" if prefix else k
            keys.append(full_key)
            keys.extend(get_nested_keys(v, full_key))
    elif isinstance(item, list):
        for index, element in enumerate(item):
            indexed_prefix = f"{prefix}[{index}]"
            keys.extend(get_nested_keys(element, indexed_prefix))
    return keys


@api.route("/species", methods=["GET"])
def get_species():
    """Return available species as JSON"""
    try:
        species = utils.get_available_species()
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

        # Extract main parameters
        sequencesToDomesticate = data.get("sequencesToDomesticate", [])
        species = data.get("species", "")
        kozak = data.get("kozak", "MTK")
        max_mut_per_site = data.get("max_mut_per_site", 3)
        verbose_mode = data.get("verbose_mode", True)
        template_sequence = data.get("templateSequence", "")
        max_results = data.get("max_results", 1)

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
            max_results=max_results,
            verbose=verbose_mode
        )

        # Generate the protocol
        result = protocol_maker.create_gg_protocol()
        serializable_result = utils.convert_non_serializable(result)

        result_keys = get_nested_keys(result)
        serializable_keys = get_nested_keys(serializable_result)

        # Optionally, compare key sets to see if any keys were lost in conversion
        missing_keys = set(result_keys) - set(serializable_keys)
        if missing_keys:
            logger.warning(
                "Missing keys in serializable result: %s", missing_keys)
        else:
            logger.debug("No keys are missing in serializable result.")

        # Check for errors
        if serializable_result.get('has_errors', False):
            return jsonify({
                "error": "Error in protocol generation",
                "sequence_errors": serializable_result.get('sequence_errors', {})
            }), 422

        print(
            f"Generated protocol: {serializable_result}")
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
