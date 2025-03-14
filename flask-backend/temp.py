# api.py
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
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Global dictionary to store job statuses.
job_status = {}


def run_protocol_job(job_id, data):
    """
    Runs the protocol generation job in a background thread and updates
    the job_status dictionary periodically.
    """
    try:
        # Initial status update
        job_status[job_id] = {"percentage": 0,
                              "message": "Job started", "step": "init"}

        # Simulate some initialization work
        time.sleep(2)
        job_status[job_id] = {
            "percentage": 25, "message": "Initializing protocol...", "step": "init"}
        time.sleep(2)
        job_status[job_id] = {
            "percentage": 50, "message": "Running protocol generation...", "step": "processing"}
        time.sleep(2)
        job_status[job_id] = {
            "percentage": 75, "message": "Finalizing protocol...", "step": "finalizing"}

        # Extract parameters from request data
        sequencesToDomesticate = data.get("sequencesToDomesticate", [])
        species = data.get("species", "")
        kozak = data.get("kozak", "MTK")
        max_mut_per_site = data.get("max_mut_per_site", 3)
        verbose_mode = data.get("verbose_mode", True)
        template_sequence = data.get("templateSequence", "")
        max_results = data.get("max_results", 1)

        # Run the actual protocol generation
        protocol_maker = GoldenGateProtocol(
            sequencesToDomesticate=sequencesToDomesticate,
            codon_usage_dict=utils.get_codon_usage_dict(species),
            max_mutations=max_mut_per_site,
            template_seq=template_sequence,
            kozak=kozak,
            max_results=max_results,
            verbose=verbose_mode
        )
        result = protocol_maker.create_gg_protocol()
        serializable_result = utils.convert_non_serializable(result)

        # Check for errors in the protocol generation
        if serializable_result.get("has_errors", False):
            job_status[job_id] = {
                "percentage": -1,
                "message": "Error in protocol generation",
                "step": "error"
            }
        else:
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
    return a jobId for the frontend to use with the SSE endpoint.
    """
    try:
        data = request.json
        # Use jobId from data or generate one if missing
        job_id = data.get("jobId", str(int(time.time() * 1000)))
        data["jobId"] = job_id

        # Start the job in a background thread
        thread = Thread(target=run_protocol_job, args=(job_id, data))
        thread.start()

        # Immediately return the job id with a 202 Accepted status.
        return jsonify({"jobId": job_id, "message": "Job started"}), 202
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
            if status and (status.get("percentage") == 100 or status.get("percentage") == -1):
                break
            time.sleep(1)
    return Response(event_stream(), mimetype="text/event-stream")
