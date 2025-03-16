# services/api.py
from flask import Blueprint, request, jsonify, current_app, Response
from config.logging_config import logger
from services.utils import GoldenGateUtils
from services.protocol import GoldenGateProtocol
from flask_cors import CORS
import json
import logging
import time
from threading import Thread

api = Blueprint("api", __name__, url_prefix="/api")
CORS(api)  # Enable CORS for all routes

utils = GoldenGateUtils()
logger.setLevel(logging.DEBUG)

# Global dictionary to store job statuses.
job_status = {}


def run_protocol_job(job_id, data):
    """
    Runs the protocol generation job in a background thread and updates
    the job_status dictionary by passing a progress_callback into the
    protocol maker. The protocol maker is expected to call progress_callback
    at key steps in its process.
    """
    try:
        # Initial status update
        job_status[job_id] = {"percentage": 0,
                              "message": "Job started", "step": "init"}

        # Extract parameters from request data
        sequencesToDomesticate = data.get("sequencesToDomesticate", [])
        print(f"Sequences to domesticate: {sequencesToDomesticate}")
        species = data.get("species", "")
        kozak = data.get("kozak", "MTK")
        max_mut_per_site = data.get("max_mut_per_site", 3)
        verbose_mode = data.get("verbose_mode", True)
        template_sequence = data.get("templateSequence", "")
        max_results = data.get("max_results", 1)

        # Create the protocol maker instance with your parameters.
        protocol_maker = GoldenGateProtocol(
            sequencesToDomesticate=sequencesToDomesticate,
            codon_usage_dict=utils.get_codon_usage_dict(species),
            max_mutations=max_mut_per_site,
            template_seq=template_sequence,
            kozak=kozak,
            max_results=max_results,
            verbose=verbose_mode,
            debug=False,
        )

        # Define a progress callback that updates the global job_status.
        def progress_callback(percentage, message, step):
            # Optionally, you can merge this update with previous details if
            # needed.
            job_status[job_id] = {"percentage": percentage,
                                  "message": message, "step": step}
            # For example, if you want to preserve per-sequence details, you
            # could do so here.
            # logger.debug(f"Progress update for job {job_id}: {percentage}% -
            # {message} ({step})")

        # Run the protocol generation; this function is expected to invoke
        # progress_callback
        # at each key processing step (e.g. after preprocessing, mutation
        # analysis, primer design, etc.).
        result = protocol_maker.create_gg_protocol(
            progress_callback=progress_callback)

        # Convert the result to a serializable format.
        serializable_result = utils.convert_non_serializable(result)

        # Merge primer names from the submitted data into the results.
        if isinstance(serializable_result, list):
            for i, seq in enumerate(sequencesToDomesticate):
                primer_name = seq.get("primerName")
                if primer_name and i < len(serializable_result):
                    serializable_result[i]["primerName"] = primer_name
        else:
            # If the result is a single object, you can merge data in a
            # different way.
            # For example, you could include the entire submitted sequences.
            serializable_result["submittedSequences"] = sequencesToDomesticate

        # Check for errors in the protocol generation
        if serializable_result.get("has_errors", False):
            job_status[job_id] = {
                "percentage": -1,
                "message": "Error in protocol generation",
                "step": "error"
            }
        else:
            # Final update: set status to complete and attach the result.
            job_status[job_id] = {
                "percentage": 100,
                "message": "Job complete",
                "step": "done",
                "result": serializable_result
            }
    except Exception as e:
        logger.error(f"Error in protocol job: {str(e)}", exc_info=True)
        job_status[job_id] = {
            "percentage": -1,
            "message": "Error: " + str(e),
            "step": "error"
        }


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
        # Continue streaming until the job is complete or an error occurs.
        while True:
            status = job_status.get(job_id)
            if status and status != last_status:
                # Format message as SSE: data: <json>\n\n
                yield f"data: {json.dumps(status)}\n\n"
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
    print(f"Current config: {config}")
    return jsonify(config)
