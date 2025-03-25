# celery_tasks.py
from celery import shared_task, group
from flask_sse import sse
from flask_backend.services import GoldenGateUtils, ProtocolMaker
from flask_backend.models import ProtocolRequest, DomesticationResult
from flask_backend.logging import logger

def publish_progress(step: str, message: str, progress: float, job_id: str, sequence_idx: int, **kwargs):
    """Publish progress updates via Flask-SSE to a Redis channel.

    Accept arbitrary keyword arguments to be included directly in the payload.
    """
    channel = f"job_{job_id}_{sequence_idx}"
    logger.log_step("Progress Update", f"Step: {step}, Message: {message}, Progress: {progress}")
    logger.debug(f"Publishing progress to channel {channel}: {message}")

    payload = {
        "step": step,
        "message": message,
        "progress": progress,
        "sequenceIdx": sequence_idx,
        **kwargs
    }

    sse.publish(
        payload,
        type="progress",
        channel=channel,
    )



@shared_task(ignore_result=False)
def process_protocol_sequence(req_dict: dict, index: int):
    # No need to call create_app() here because tasks run in an app context.
    req = ProtocolRequest.model_validate(req_dict)
    seq = req.sequences_to_domesticate[index]

    def progress_callback(step: str, message: str, progress: float, sequenceIdx: int = None, **kwargs):
        publish_progress(
            step=step,
            message=message,
            progress=progress,
            job_id=req.job_id,
            sequence_idx=sequenceIdx if sequenceIdx is not None else index,
            **kwargs,
        )


    protocol_maker = ProtocolMaker(
        request_idx=index,
        sequence_to_domesticate=seq,
        codon_usage_dict=GoldenGateUtils().get_codon_usage_dict(req.species),
        max_mutations=req.max_mut_per_site,
        template_seq=req.template_sequence,
        kozak=req.kozak,
        max_results=req.max_results,
        verbose=req.verbose_mode,
        job_id=f"{req.job_id}_{index}",
    )

    # Execute protocol and explicitly serialize result.
    result: DomesticationResult = protocol_maker.create_gg_protocol(progress_callback)
    
    return {"sequenceIdx": index, "result": result.model_dump()}

@shared_task(ignore_result=False)
def generate_protocol_task(req_dict: dict):
    req = ProtocolRequest.model_validate(req_dict)
    total_sequences = len(req.sequences_to_domesticate)

    tasks = group(
        process_protocol_sequence.s(req_dict, idx) for idx in range(total_sequences)
    )
    group_result = tasks.apply_async()

    return {
        "status": "started",
        "group_task_id": group_result.id,
        "total": total_sequences
    }
