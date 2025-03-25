from flask import Blueprint, request, jsonify, current_app, Response, stream_with_context
from functools import wraps
import time
import json

from flask_backend.services import GoldenGateUtils
from flask_backend.logging import logger
from flask_backend.celery_tasks import generate_protocol_task


api = Blueprint("api", __name__, url_prefix="/api")
utils = GoldenGateUtils()


def handle_errors(f):
    """Decorator to handle exceptions in API routes with standardized error responses."""
    @wraps(f)
    def decorated_function(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except Exception as e:
            logger.error(e, exc_info=True)
            return jsonify({"error": str(e)}), 500
    return decorated_function


def get_celery_instance():
    """Get the Celery instance from the current Flask application."""
    return current_app.extensions["celery"]


@api.route("/generate_protocol", methods=["POST"])
@handle_errors
def generate_protocol():
    """Start an asynchronous protocol generation task."""
    data = request.get_json()
    if not data:
        return jsonify({"error": "No data provided"}), 400
        
    task = generate_protocol_task.delay(data)
    return jsonify({"task_id": task.id}), 202


@api.route("/status/<job_id>")
def sse_status(job_id):
    print(f"SSE route hit for job_id={job_id}")
    def event_stream():
        while True:
            async_result = get_celery_instance().AsyncResult(job_id)
            state = async_result.state
            print(f"Current state: {state}")
            meta = async_result.info or {}
            
            # Ensure meta is a dictionary before processing.
            if not isinstance(meta, dict):
                msg = "Value of KeyError:", meta.args[0]


                logger.log_step("SSE", msg)
                # Convert non-dict meta to a dict; for example, store its string representation.
                meta = {"error": str(meta)}
            
            # print("Meta keys:", list(meta.keys()))
            
            try:
                data = json.dumps(meta)
            except TypeError as e:
                # In case any value is still not serializable, convert it.
                meta = {str(k): (v if isinstance(v, (int, float, str, bool, list, dict)) else str(v))
                        for k, v in meta.items()}
                data = json.dumps(meta)
                data['error'] = str(e)
            
            yield f"data: {data}\n\n"
            if state in ["SUCCESS", "FAILURE"]:
                break
            time.sleep(2)

            
    return Response(stream_with_context(event_stream()), mimetype="text/event-stream")


@api.route("/task-status/<task_id>", methods=["GET"])
@handle_errors
def get_task_status(task_id):
    """Get the status of an asynchronous task by its ID."""
    celery = get_celery_instance()
    async_result = celery.AsyncResult(task_id)
    response = {"task_id": task_id, "state": async_result.state}

    if async_result.state == "PROGRESS":
        response.update(async_result.info)
    elif async_result.state == "SUCCESS":
        response["result"] = async_result.result.get("result")
    elif async_result.state == "FAILURE":
        response["error"] = str(async_result.result)

    return jsonify(response)


@api.route("/export", methods=["POST"])
@handle_errors
def export_protocol():
    """Export protocol primers to a TSV file."""
    data = request.get_json()
    primers = data.get("primers", [])
    if not primers:
        return jsonify({"error": "No primers to export"}), 400

    filename = f"primers_{utils.generate_unique_id()}.tsv"
    filepath = f"static/exports/{filename}"
    
    try:
        with open(filepath, "w") as f:
            f.write("Primer Name\tSequence\tAmplicon\n")
            for primer in primers:
                f.write(f"{primer[0]}\t{primer[1]}\t{primer[2]}\n")
        return jsonify({"download_url": f"/static/exports/{filename}"})
    except IOError as e:
        logger.error(f"Failed to write export file: {e}", exc_info=True)
        return jsonify({"error": "Failed to create export file"}), 500


@api.route("/species", methods=["GET"])
@handle_errors
def get_species():
    """Get available species for codon optimization."""
    return jsonify({"species": utils.get_available_species()})


@api.route("/dummy", methods=["GET"])
@handle_errors
def get_dummy_data():
    """Get active application configuration."""
    print(f"current_app.config: {current_app.config}")
    return jsonify(current_app.config.get("PREFILL_DATA", {}))