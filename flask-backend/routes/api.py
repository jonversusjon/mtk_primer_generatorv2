# services/api.py
from flask import Blueprint, request, jsonify, current_app, Response
from models import ProtocolRequest
from services import ProtocolMaker, GoldenGateUtils
from flask_cors import CORS
import json
import time
import threading
from pydantic import BaseModel
from log_utils import logger
from functools import partial

# Cleanup timeout in seconds (e.g. one hour)
JOB_CLEANUP_TIMEOUT = 3600

api = Blueprint("api", __name__, url_prefix="/api")
CORS(api)  # Enable CORS for all routes

utils = GoldenGateUtils()

# Global dictionaries for job status and job creation timestamps.
job_status = {}
job_timestamps = {}
job_status_lock = threading.Lock()


def send_progress_update(job_id, seq_name, **updates):
    """
    Helper to merge a default progress dictionary with update values and then update the job status.
    
    Args:
        job_id (str): ID of the job.
        seq_name (str): Identifier for the sequence.
        updates (dict): Arbitrary keyword arguments (e.g., percentage, message, step, status, data).
    """
    base_progress = {
        "status": "running",
        "percentage": 0,
        "message": "",
        "step": ""
    }
    # Merge the default progress with the provided updates
    progress_data = base_progress | updates
    
    # Add a timestamp to force detection of the update
    progress_data["timestamp"] = time.time()
    
    update_payload = {"sequences": {seq_name: progress_data}}
    
    # Also update the job's overall percentage based on sequence progress
    if "percentage" in updates and updates["percentage"] > 0:
        update_payload["sequence_progress"] = updates["percentage"]
    
    update_job_status(job_id, update_payload)
    
    logger.debug(f"[Job {job_id}] [{seq_name}] {progress_data['percentage']}% - "
                f"{progress_data['message']} ({progress_data['step']})")


def convert_to_serializable(obj):
    """
    Recursively convert Pydantic models and other non-serializable objects to dictionaries.
    """
    if isinstance(obj, BaseModel):
        return obj.model_dump() if hasattr(obj, 'model_dump') else obj.dict()
    elif isinstance(obj, dict):
        return {k: convert_to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_serializable(item) for item in obj]
    elif isinstance(obj, tuple):
        return tuple(convert_to_serializable(item) for item in obj)
    else:
        return obj


def process_sequence(job_id, seq, kozak, max_mut_per_site, verbose_mode, template_sequence, max_results, codon_usage_dict):
    """
    Process a single sequence, updating its own status as well as overall job counts.
    """
    seq_name = getattr(seq, "primerName", None) or f"seq_{id(seq)}"
    
    send_progress_update(job_id,
                         seq_name,
                         status="running",
                         percentage=0,
                         message=f"Started processing {seq_name}",
                         step="initialization")
    
    progress_callback = partial(send_progress_update, job_id, seq_name)
    
    try:
        protocol_maker = ProtocolMaker(
            request_idx=0,
            sequence_to_domesticate=seq,
            codon_usage_dict=codon_usage_dict,
            max_mutations=max_mut_per_site,
            template_seq=template_sequence,
            kozak=kozak,
            max_results=max_results,
            verbose=verbose_mode,
            debug=False,
            job_id=job_id
        )
        result = protocol_maker.create_gg_protocol(progress_callback=progress_callback)
        serializable_result = utils.convert_non_serializable(result)
        
        # Merge primer name if provided.
        if hasattr(seq, 'primerName') and seq.primerName:
            if isinstance(serializable_result, list) and serializable_result:
                serializable_result[0]["primerName"] = seq.primerName
            elif isinstance(serializable_result, dict):
                serializable_result["primerName"] = seq.primerName
        
        progress_callback(percentage=100,
                          status="completed",
                          message=f"Completed processing {seq_name}",
                          step="finished",
                          data={"result": serializable_result})        
    except Exception as e:
        logger.error(f"Error processing sequence [{seq_name}] in job [{job_id}]: {str(e)}", exc_info=True)
        progress_callback(percentage=-1,
                          status="error",
                          message=f"Error: {str(e)}",
                          step="error",
                          data={"error": str(e)})
    finally:
        # Update overall job progress.
        with job_status_lock:
            job_status[job_id]["completed_count"] = job_status[job_id].get("completed_count", 0) + 1
            total = job_status[job_id].get("total_count", 1)
            if job_status[job_id]["completed_count"] >= total:
                job_status[job_id]["status"] = "completed"
                job_status[job_id]["end_time"] = time.time()



# 1. Modify stream_status in api.py to detect sequence-level changes

@api.route("/status/<job_id>", methods=["GET"])
def stream_status(job_id):
    """
    Stream job status updates as Server-Sent Events (SSE).
    """
    def event_stream():
        last_status_json = None  # Store the JSON string instead of the object
        retry_count = 0
        max_retries = 5
        
        while True:
            current_status = None
            with job_status_lock:
                if job_id in job_status:
                    current_status = convert_to_serializable(job_status[job_id])
                    retry_count = 0
                else:
                    retry_count += 1
            
            if retry_count >= max_retries:
                yield f"data: {json.dumps({'error': 'Job not found'})}\n\n"
                break
            
            if current_status:
                # Convert to JSON string for comparison
                current_status_json = json.dumps(current_status)
                
                # Send update if anything changed (including sequence details)
                if current_status_json != last_status_json:
                    logger.debug(f"Sending SSE update for job {job_id}: {current_status_json[:100]}...")
                    yield f"data: {current_status_json}\n\n"
                    last_status_json = current_status_json
                    
                    # Check if job is complete
                    if current_status.get("status") in ["completed", "error"]:
                        # Send a final message to ensure client gets closure notification
                        final_message = {
                            "status": current_status.get("status"),
                            "message": "Job processing complete",
                            "percentage": 100 if current_status.get("status") == "completed" else -1
                        }
                        yield f"data: {json.dumps(final_message)}\n\n"
                        break
            
            # Reduce sleep time to capture updates more frequently
            time.sleep(0.25)  # Check 4 times per second
            
    return Response(event_stream(), mimetype="text/event-stream")

# 3. Update the update_job_status function to better track sequence progress

def update_job_status(job_id, update_data):
    """
    Thread-safe update to the global job_status dictionary.
    
    Uses a deep-update to merge nested dictionaries.
    """
    with job_status_lock:
        if job_id not in job_status:
            job_status[job_id] = {
                "sequences": {},
                "sequence_updates_count": 0  # Track how many sequence updates we've received
            }
        
        # Deep update helper
        def deep_update(d, u):
            for k, v in u.items():
                if isinstance(v, dict) and k in d and isinstance(d[k], dict):
                    deep_update(d[k], v)
                else:
                    d[k] = v
        
        # Track sequence updates
        if "sequences" in update_data:
            job_status[job_id]["sequence_updates_count"] = job_status[job_id].get("sequence_updates_count", 0) + 1
        
        deep_update(job_status[job_id], update_data)
        
        # Update overall job percentage based on completed sequences
        if job_status[job_id].get("completed_count") is not None and job_status[job_id].get("total_count") is not None:
            completed = job_status[job_id]["completed_count"]
            total = job_status[job_id]["total_count"]
            
            # If there's sequence progress info, use a weighted calculation
            if "sequence_progress" in update_data:
                # Factor in the current sequence's progress
                job_status[job_id]["percentage"] = int((completed * 100 + update_data["sequence_progress"]) / total)
            else:
                # Just use completed/total ratio
                job_status[job_id]["percentage"] = int(100 * completed / total)
        
        # Ensure completed jobs always show 100%
        if job_status[job_id].get("status") == "completed":
            job_status[job_id]["percentage"] = 100
            
def cleanup_old_jobs():
    """
    Remove completed or errored jobs older than JOB_CLEANUP_TIMEOUT from memory.
    """
    current_time = time.time()
    with job_status_lock:
        jobs_to_remove = []
        for job_id, timestamp in job_timestamps.items():
            if current_time - timestamp > JOB_CLEANUP_TIMEOUT:
                if job_status.get(job_id, {}).get("status") in ["completed", "error"]:
                    jobs_to_remove.append(job_id)
        for job_id in jobs_to_remove:
            del job_status[job_id]
            del job_timestamps[job_id]
            logger.debug(f"Cleaned up job {job_id} (older than {JOB_CLEANUP_TIMEOUT}s)")


def start_job(job_id: str, req: ProtocolRequest):
    """
    Initialize a new job using data from a ProtocolRequest model, store initial status,
    and spawn a thread for each sequence.
    """
    try:
        # Record job creation time for cleanup.
        job_timestamp = time.time()
        with job_status_lock:
            job_timestamps[job_id] = job_timestamp
            job_status[job_id] = {
                "status": "running",
                "start_time": job_timestamp,
                "total_count": len(req.sequences_to_domesticate),
                "completed_count": 0,
                "sequences": {}
            }
        
        codon_usage_dict = utils.get_codon_usage_dict(req.species)
        # Spawn a thread for each sequence.
        for seq in req.sequences_to_domesticate:
            threading.Thread(
                target=process_sequence,
                args=(
                    job_id,
                    seq,
                    req.kozak,
                    req.max_mut_per_site,
                    req.verbose_mode,
                    req.template_sequence,
                    req.max_results,
                    codon_usage_dict
                ),
                daemon=True
            ).start()
        
        return {
            "jobId": job_id,
            "status": "started",
            "message": f"Job started with {len(req.sequences_to_domesticate)} sequences",
            "total_sequences": len(req.sequences_to_domesticate)
        }
    except Exception as e:
        logger.error(f"Error starting job {job_id}: {str(e)}", exc_info=True)
        update_job_status(job_id, {
            "status": "error",
            "error": str(e),
            "message": f"Failed to start job: {str(e)}"
        })
        return {
            "jobId": job_id,
            "status": "error",
            "message": f"Error starting job: {str(e)}"
        }


@api.route("/generate_protocol", methods=["POST"])
def generate_protocol():
    """
    Start a protocol generation job in the background using a Pydantic model to
    parse incoming request data. Immediately return a jobId and the submitted data.
    """
    try:
        # Use Pydantic to parse and validate the incoming JSON, converting camelCase to snake_case.
        req_data = ProtocolRequest.parse_obj(request.json)
        # Use provided jobId or generate one based on the current time.
        job_id = req_data.job_id or str(int(time.time() * 1000))
        req_data.job_id = job_id  # set the job_id back in the model
        
        result = start_job(job_id, req_data)
        
        # Launch a cleanup thread.
        threading.Thread(target=cleanup_old_jobs, daemon=True).start()
        
        return jsonify(result), 202
    except Exception as e:
        logger.error(f"Error in generate_protocol endpoint: {str(e)}", exc_info=True)
        return jsonify({"error": str(e)}), 500


@api.route("/export", methods=["POST"])
def export_protocol():
    """
    Export protocol results as a TSV file.
    """
    try:
        data = request.get_json()
        primers = data.get("primers", [])

        if not primers:
            return jsonify({"error": "No primers to export"}), 400

        filename = f"primers_{utils.generate_unique_id()}.tsv"
        filepath = f"static/exports/{filename}"

        with open(filepath, "w") as f:
            f.write("Primer Name\tSequence\tAmplicon\n")
            for primer in primers:
                f.write(f"{primer[0]}\t{primer[1]}\t{primer[2]}\n")

        download_url = f"/static/exports/{filename}"
        return jsonify({"download_url": download_url})
    except Exception as e:
        logger.error(f"Error exporting protocol: {str(e)}", exc_info=True)
        return jsonify({"error": str(e)}), 500


@api.route("/species", methods=["GET"])
def get_species():
    """
    Return available species as JSON.
    """
    try:
        species = utils.get_available_species()
        return jsonify({"species": species})
    except Exception as e:
        logger.error(f"Error fetching species: {str(e)}", exc_info=True)
        return jsonify({'error': 'Failed to fetch species'}), 500


@api.route('/config', methods=['GET'])
def get_config():
    """
    Return the currently loaded configuration from app.py.
    """
    config = current_app.config.get("ACTIVE_CONFIG", {})
    return jsonify(config)
