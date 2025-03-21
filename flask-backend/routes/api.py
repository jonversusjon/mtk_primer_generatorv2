from flask import Blueprint, request, jsonify
from celery_worker import celery
from models import ProtocolRequest
from services import GoldenGateUtils
from log_utils import logger

api = Blueprint("api", __name__, url_prefix="/api")
utils = GoldenGateUtils()

@celery.task(bind=True)
def generate_protocol_task(self, req_dict):
    """Celery task that processes all sequences serially."""
    req = ProtocolRequest.model_validate(req_dict)
    job_id = self.request.id
    total = len(req.sequences_to_domesticate)

    for idx, seq in enumerate(req.sequences_to_domesticate):
        progress = int((idx / total) * 100)
        self.update_state(state="PROGRESS", meta={"step": idx+1, "total": total, "progress": progress})
        
        result = GoldenGateUtils().process_sequence(
            job_id=job_id,
            seq=seq,
            kozak=req.kozak,
            max_mut_per_site=req.max_mut_per_site,
            verbose_mode=req.verbose_mode,
            template_sequence=req.template_sequence,
            max_results=req.max_results,
            codon_usage_dict=utils.get_codon_usage_dict(req.species)
        )
    
    return {"status": "completed", "result": result}

@api.route("/generate_protocol", methods=["POST"])
def generate_protocol():
    data = request.get_json()
    req = ProtocolRequest.model_validate(data)
    task = generate_protocol_task.apply_async(args=[req.model_dump()])
    return jsonify({"task_id": task.id}), 202

@api.route("/task-status/<task_id>", methods=["GET"])
def get_task_status(task_id):
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
def export_protocol():
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
    return jsonify({"download_url": f"/static/exports/{filename}"})

@api.route("/species", methods=["GET"])
def get_species():
    try:
        return jsonify({"species": utils.get_available_species()})
    except Exception as e:
        logger.error(e, exc_info=True)
        return jsonify({'error': 'Failed to fetch species'}), 500

@api.route("/config", methods=["GET"])
def get_config():
    return jsonify(request.app.config.get("ACTIVE_CONFIG", {}))
