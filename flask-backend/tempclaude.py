# services/api.py
from flask import Blueprint, request, jsonify, current_app, Response
from services.utils import GoldenGateUtils
from services.protocol_maker import ProtocolMaker
from flask_cors import CORS
import json
import time
import threading
from concurrent.futures import ThreadPoolExecutor
from models import SequenceToDomesticate
from pydantic import BaseModel
from log_utils import logger

# Maximum number of worker threads in the pool
MAX_WORKERS = 10

# Timeout for cleaning up completed jobs (in seconds)
JOB_CLEANUP_TIMEOUT = 3600

api = Blueprint("api", __name__, url_prefix="/api")
CORS(api)  # Enable CORS for all routes

utils = GoldenGateUtils()

# Thread pool for running background jobs
thread_pool = ThreadPoolExecutor(max_workers=MAX_WORKERS)

# Global dictionary to store job statuses with thread-safe access
job_status = {}
job_status_lock = threading.Lock()

# Dictionary to track job creation times for cleanup
job_timestamps = {}


def convert_to_serializable(obj):
    """Convert Pydantic models and other non-serializable objects to dictionaries."""
    if isinstance(obj, BaseModel):
        # Use model_dump() for Pydantic v2 or dict() for v1
        return obj.model_dump() if hasattr(obj, 'model_dump') else obj.dict()
    elif isinstance(obj, dict):
        return {k: convert_to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_serializable(item) for item in obj]
    elif isinstance(obj, tuple):
        return tuple(convert_to_serializable(item) for item in obj)
    else:
        return obj


def parse_request_data(data):
    """
    Parse and validate incoming request data for the protocol job.
    
    Args:
        data: Dict containing request data
        
    Returns:
        Tuple of parsed parameters
    """
    sequences_json = data.get("sequencesToDomesticate", [])
    sequences_to_domesticate = [SequenceToDomesticate.model_validate(seq) for seq in sequences_json]

    return (
        sequences_to_domesticate,
        data.get("species", ""),
        data.get("kozak", "MTK"),
        data.get("max_mut_per_site", 3),
        data.get("verbose_mode", True),
        data.get("templateSequence", ""),
        data.get("maxResults", "err"),
    )


def update_job_status(job_id, update_data):
    """
    Thread-safe update to the job status dictionary.
    
    Args:
        job_id: Unique job identifier
        update_data: Dictionary with data to update or add
    """
    with job_status_lock:
        if job_id not in job_status:
            job_status[job_id] = {}
        
        # Deep update that merges nested dictionaries
        def deep_update(d, u):
            for k, v in u.items():
                if isinstance(v, dict) and k in d and isinstance(d[k], dict):
                    deep_update(d[k], v)
                else:
                    d[k] = v
                    
        deep_update(job_status[job_id], update_data)


def process_sequence(job_id, seq, species, kozak, max_mut_per_site, verbose_mode, template_sequence, max_results):
    """
    Process a single sequence and update its status in the job_status dictionary.
    
    Args:
        job_id: Unique job identifier
        seq: SequenceToDomesticate object
        species: Species for codon optimization
        kozak: Kozak sequence
        max_mut_per_site: Maximum mutations per site
        verbose_mode: Enable verbose logging
        template_sequence: Template sequence
        max_results: Maximum results to return
    """
    seq_name = seq.primerName if hasattr(seq, 'primerName') and seq.primerName else f"seq_{id(seq)}"
    
    # Create a progress callback for this sequence
    def progress_callback(percentage, message, step):
        update_job_status(job_id, {
            "sequences": {
                seq_name: {
                    "percentage": percentage,
                    "message": message,
                    "step": step
                }
            }
        })
        logger.debug(f"[Job {job_id}] [{seq_name}] {percentage}% - {message} ({step})")
    
    try:
        # Update initial status
        update_job_status(job_id, {
            "sequences": {
                seq_name: {
                    "status": "running",
                    "percentage": 0,
                    "message": f"Started processing {seq_name}",
                    "step": "initialization"
                }
            }
        })
        
        # Create and run protocol maker
        protocol_maker = ProtocolMaker(
            sequences_to_domesticate=[seq],
            codon_usage_dict=utils.get_codon_usage_dict(species),
            max_mutations=max_mut_per_site,
            template_seq=template_sequence,
            kozak=kozak,
            max_results=max_results,
            verbose=verbose_mode,
            debug=False,
        )
        
        result = protocol_maker.create_gg_protocol(progress_callback=progress_callback)
        serializable_result = utils.convert_non_serializable(result)
        
        # Add primer name if available
        if hasattr(seq, 'primerName') and seq.primerName:
            if isinstance(serializable_result, list) and serializable_result:
                serializable_result[0]["primerName"] = seq.primerName
            elif isinstance(serializable_result, dict):
                serializable_result["primerName"] = seq.primerName
        
        # Update final status
        update_job_status(job_id, {
            "sequences": {
                seq_name: {
                    "status": "completed",
                    "percentage": 100,
                    "message": f"Completed processing {seq_name}",
                    "step": "finished",
                    "result": serializable_result
                }
            }
        })
        
    except Exception as e:
        logger.error(f"Error processing sequence [{seq_name}] in job [{job_id}]: {str(e)}", exc_info=True)
        update_job_status(job_id, {
            "sequences": {
                seq_name: {
                    "status": "error",
                    "percentage": -1,
                    "message": f"Error: {str(e)}",
                    "step": "error",
                    "error": str(e)
                }
            }
        })
    finally:
        # Update completion counter
        with job_status_lock:
            job_status[job_id]["completed_count"] = job_status[job_id].get("completed_count", 0) + 1
            # Check if all sequences are done
            if job_status[job_id]["completed_count"] >= job_status[job_id]["total_count"]:
                job_status[job_id]["status"] = "completed"
                job_status[job_id]["end_time"] = time.time()


def start_job(job_id, data):
    """
    Initialize a new job and start processing sequences in parallel.
    
    Args:
        job_id: Unique job identifier
        data: Request data dictionary
    
    Returns:
        Dict with job status information
    """
    try:
        # Parse request data
        sequences_to_domesticate, species, kozak, max_mut_per_site, verbose_mode, template_sequence, max_results = parse_request_data(data)
        
        # Initialize job status
        job_timestamp = time.time()
        job_timestamps[job_id] = job_timestamp
        
        # Create initial job status
        update_job_status(job_id, {
            "status": "running",
            "start_time": job_timestamp,
            "total_count": len(sequences_to_domesticate),
            "completed_count": 0,
            "sequences": {}
        })
        
        # Submit each sequence for processing
        for seq in sequences_to_domesticate:
            thread_pool.submit(
                process_sequence,
                job_id,
                seq,
                species,
                kozak,
                max_mut_per_site,
                verbose_mode,
                template_sequence,
                max_results
            )
        
        return {
            "jobId": job_id,
            "status": "started",
            "message": f"Job started with {len(sequences_to_domesticate)} sequences",
            "total_sequences": len(sequences_to_domesticate)
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


def cleanup_old_jobs():
    """Remove completed jobs older than JOB_CLEANUP_TIMEOUT from memory."""
    current_time = time.time()
    with job_status_lock:
        jobs_to_remove = []
        for job_id, timestamp in job_timestamps.items():
            if current_time - timestamp > JOB_CLEANUP_TIMEOUT:
                # Only clean up completed or errored jobs
                if job_id in job_status and job_status[job_id].get("status") in ["completed", "error"]:
                    jobs_to_remove.append(job_id)
        
        for job_id in jobs_to_remove:
            del job_status[job_id]
            del job_timestamps[job_id]
            logger.debug(f"Cleaned up job {job_id} (older than {JOB_CLEANUP_TIMEOUT}s)")


@api.route("/generate_protocol", methods=["POST"])
def generate_protocol():
    """
    Start a protocol generation job in the background and immediately
    return a jobId and the submitted data for the frontend to use.
    """
    try:
        data = request.json
        
        # Use jobId from data or generate one
        job_id = data.get("jobId", str(int(time.time() * 1000)))
        data["jobId"] = job_id
        
        # Start the job processing
        result = start_job(job_id, data)
        
        # Clean up old jobs in the background
        threading.Thread(target=cleanup_old_jobs, daemon=True).start()
        
        return jsonify(result), 202
        
    except Exception as e:
        logger.error(f"Error in generate_protocol endpoint: {str(e)}", exc_info=True)
        return jsonify({"error": str(e)}), 500


@api.route("/status/<job_id>", methods=["GET"])
def stream_status(job_id):
    """
    Stream job status updates as SSE (Server-Sent Events) so the frontend
    can receive real-time updates.
    """
    def event_stream():
        last_status = None
        retry_count = 0
        max_retries = 5  # Stop streaming if job not found after retries
        
        while True:
            current_status = None
            with job_status_lock:
                if job_id in job_status:
                    current_status = convert_to_serializable(job_status[job_id])
                    retry_count = 0
                else:
                    retry_count += 1
            
            # If job not found after retries, stop streaming
            if retry_count >= max_retries:
                yield f"data: {json.dumps({'error': 'Job not found'})}\n\n"
                break
                
            # If status changed, send update
            if current_status and current_status != last_status:
                yield f"data: {json.dumps(current_status)}\n\n"
                last_status = current_status
                
                # Break if job completed or errored
                if current_status.get("status") in ["completed", "error"]:
                    break
                    
            time.sleep(1)
            
    return Response(event_stream(), mimetype="text/event-stream")


@api.route("/export", methods=["POST"])
def export_protocol():
    """Export protocol results as a TSV file."""
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
    """Return available species as JSON."""
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
    return jsonify(config)