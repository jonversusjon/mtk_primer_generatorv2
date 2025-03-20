# services/api.py
from flask import Blueprint, request, jsonify, current_app, Response
from config.logging_config import logger
from services.utils import GoldenGateUtils
from services.protocol_maker import GoldenGateProtocol
from flask_cors import CORS
import json
import logging
import time
from threading import Thread
from models import SequenceToDomesticate
from pydantic import BaseModel


api = Blueprint("api", __name__, url_prefix="/api")
CORS(api)  # Enable CORS for all routes

utils = GoldenGateUtils()
logger.setLevel(logging.DEBUG)

# Global dictionary to store job statuses.
job_status = {}

def convert_pydantic_to_dict(obj):
    """Recursively convert Pydantic models to dictionaries."""
    if isinstance(obj, BaseModel):  # Check if it's any Pydantic model
        # Use model_dump() for Pydantic v2 or dict() for v1
        return obj.model_dump() if hasattr(obj, 'model_dump') else obj.dict()
    elif isinstance(obj, dict):
        return {k: convert_pydantic_to_dict(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_pydantic_to_dict(item) for item in obj]
    elif isinstance(obj, tuple):
        return tuple(convert_pydantic_to_dict(item) for item in obj)
    else:
        return obj

def parse_request_data(data):
    """
    Parse and validate incoming request data for the protocol job.

    :param data: dict of request data.
    :return: tuple of:
        - sequences_to_domesticate (List[SequenceToDomesticate])
        - species (str)
        - kozak (str)
        - max_mut_per_site (int)
        - verbose_mode (bool)
        - template_sequence (str)
        - max_results (int)
    """
    sequences_json = data.get("sequencesToDomesticate", [])
    sequences_to_domesticate = [SequenceToDomesticate.model_validate(seq) for seq in sequences_json]

    species = data.get("species", "")
    kozak = data.get("kozak", "MTK")
    max_mut_per_site = data.get("max_mut_per_site", 3)
    verbose_mode = data.get("verbose_mode", True)
    template_sequence = data.get("templateSequence", "")
    max_results = data.get("max_results", 1)

    return (
        sequences_to_domesticate,
        species,
        kozak,
        max_mut_per_site,
        verbose_mode,
        template_sequence,
        max_results,
    )

def merge_primer_names(serializable_result, sequences_to_domesticate):
    """
    Merge user-submitted primer names into the final result.

    :param serializable_result: The result object (list or dict).
    :param sequences_to_domesticate: List of SequenceToDomesticate objects.
    """
    # If the result is a list, we assume one entry per sequence.
    if isinstance(serializable_result, list):
        for i, seq in enumerate(sequences_to_domesticate):
            primer_name = seq.get("primerName")
            if primer_name and i < len(serializable_result):
                serializable_result[i]["primerName"] = primer_name
    else:
        # If it's a single object, store them in another way:
        serializable_result["submittedSequences"] = sequences_to_domesticate

def update_job_status(job_status, job_id, percentage, message, step, result=None):
    """
    Update the global job_status dictionary for a given job_id.

    :param job_status: Global dict storing job statuses.
    :param job_id: Unique job identifier.
    :param percentage: Numeric progress percentage.
    :param message: Human-readable status message.
    :param step: The current step of processing.
    :param result: Optional final result data.
    """
    status_update = {
        "percentage": percentage,
        "message": message,
        "step": step
    }
    if result is not None:
        status_update["result"] = result

    job_status[job_id] = status_update

def make_progress_callback(job_status, job_id):
    """
    Create a progress_callback function scoped to a specific job_id.

    :param job_status: Global dict storing job statuses.
    :param job_id: Unique job identifier.
    :return: progress_callback function.
    """
    def progress_callback(percentage, message, step):
        update_job_status(job_status, job_id, percentage, message, step)
        # Optionally, add logging or other logic here:
        # logger.debug(f"[Job {job_id}] {percentage}% - {message} ({step})")

    return progress_callback

def run_protocol_job(job_id, data):
    """
    Runs the protocol generation job in a background thread and updates
    the job_status dictionary by passing a progress_callback into the
    protocol maker. The protocol maker is expected to call progress_callback
    at key steps in its process.

    :param job_id: Unique job identifier.
    :param data: dict containing request data.
    """
    try:
        # Initial status update
        update_job_status(job_status, job_id, 0, "Job started", "init")

        # Parse request data
        (
            sequences_to_domesticate,
            species,
            kozak,
            max_mut_per_site,
            verbose_mode,
            template_sequence,
            max_results,
        ) = parse_request_data(data)

        # Create protocol maker
        protocol_maker = GoldenGateProtocol(
            sequences_to_domesticate=sequences_to_domesticate,
            codon_usage_dict=utils.get_codon_usage_dict(species),
            max_mutations=max_mut_per_site,
            template_seq=template_sequence,
            kozak=kozak,
            max_results=max_results,
            verbose=verbose_mode,
            debug=False,
        )

        # Define progress callback
        progress_callback = make_progress_callback(job_status, job_id)

        # Run protocol generation
        result = protocol_maker.create_gg_protocol(progress_callback=progress_callback)

        # Convert to a serializable format
        serializable_result = utils.convert_non_serializable(result)

        # Merge primer names
        merge_primer_names(serializable_result, sequences_to_domesticate)

        # Check for errors
        if serializable_result.get("has_errors", False):
            update_job_status(job_status, job_id, -1, "Error in protocol generation", "error")
        else:
            update_job_status(
                job_status,
                job_id,
                100,
                "Job complete",
                "done",
                result=serializable_result
            )

    except Exception as e:
        logger.error(f"Error in protocol job [{job_id}]: {str(e)}", exc_info=True)
        update_job_status(job_status, job_id, -1, f"Error: {str(e)}", "error")



@api.route("/generate_protocol", methods=["POST"])
def generate_protocol():
    """
    Start a protocol generation job in the background and immediately
    return a jobId and the submitted data for the frontend to use.
    """
    try:
        data = request.json
        # Use jobId from data or generate one if missing
        job_id = data.get("jobId", str(int(time.time() * 1000)))
        data["jobId"] = job_id

        # Start the job in a background thread
        thread = Thread(target=run_protocol_job, args=(job_id, data))
        thread.start()

        # Immediately return the job id and the submitted data
        return jsonify({
            "jobId": job_id,
            "message": "Job started",
            "submittedData": data
        }), 202
    except Exception as e:
        logger.error(
            f"Error starting protocol generation: {str(e)}", exc_info=True)
        return jsonify({"error": str(e)}), 500


@api.route("/status/<job_id>", methods=["GET"])
def stream_status(job_id):
    """
    Stream job status updates as SSE (Server-Sent Events) so the frontend
    can receive periodic messages.
    """
    def event_stream():
        last_status = None
        while True:
            status = job_status.get(job_id)
            if status and status != last_status:
                # Convert any nested Pydantic models to dictionaries
                serializable_status = convert_pydantic_to_dict(status)
                yield f"data: {json.dumps(serializable_status)}\n\n"
                last_status = status
            # Break if job is complete (100%) or failed (-1)
            if status and (
                status.get("percentage") == 100 or
                status.get("percentage") == -1
            ):
                break
            time.sleep(1)
    return Response(event_stream(), mimetype="text/event-stream")


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


@api.route("/species", methods=["GET"])
def get_species():
    """Return available species as JSON"""
    try:
        species = utils.get_available_species()
        return jsonify({"species": species})
    except Exception as e:
        logger.error(f"Error fetching species: {str(e)}", exc_info=True)
        return jsonify({'error': 'Failed to fetch species'}), 500


@api.route('/config', methods=['GET'])
def get_config():
    """Return the currently loaded configuration from app.py."""
    config = current_app.config.get("ACTIVE_CONFIG", {})
    # print(f"Current config: {config}")
    return jsonify(config)
