# celery_tasks.py
import json
from typing import Optional
from celery import shared_task, group
from flask_sse import sse
from flask_backend.services import GoldenGateUtils, ProtocolMaker
from flask_backend.models import ProtocolRequest, DomesticationResult
from flask_backend.logging import logger

def publish_sse(step: str, message: str, job_id: str, sequence_idx: int, progress: Optional[float] = None, event_type: Optional[str] = None, **kwargs):
    """
    Publish updates via Flask-SSE to a Redis channel.
    
    Parameters:
        step (str): Current step name.
        message (str): Human-readable message.
        job_id (str): Unique identifier for the job.
        sequence_idx (int): Index of the current sequence.
        progress (Optional[float]): Progress percentage, if applicable.
        event_type (Optional[str]): Explicit event type ('progress' or 'data'). If None, inferred from 'progress'.
        **kwargs: Arbitrary additional data to include in the payload.
    """
    channel = f"job_{job_id}_{sequence_idx}"
    payload = {
        "step": step,
        "message": message,
        "sequenceIdx": sequence_idx,
        **kwargs
    }

    if progress is not None:
        payload["progress"] = progress
        inferred_type = "progress"
        logger.log_step("Progress Update", f"Step: {step}, Message: {message}, Progress: {progress}")
    else:
        inferred_type = "data"
        logger.log_step("Data Update", f"Step: {step}, Message: {message}")

    logger.log_step(" ******* SSE Publish ****** ", f"Publishing to channel {channel}: {json.dumps(payload, default=str)}")
    
    event_type = event_type or inferred_type

    sse.publish(
        payload,
        type=event_type,
        channel=channel,
    )

    
@shared_task(ignore_result=False)
def process_protocol_sequence(req_dict: dict, index: int):
    req = ProtocolRequest.model_validate(req_dict)
    seq = req.sequences_to_domesticate[index]

    def progress_callback(step: str, message: str, progress: Optional[float] = None, sequence_idx: Optional[int] = None, event_type: Optional[str] = None, **kwargs):
        logger.log_step("Progress Callback", f"Step: {step}, Message: {message}, Progress: {progress}")
        publish_sse(
            step=step,
            message=message,
            progress=progress,
            job_id=req.job_id,
            sequence_idx=sequence_idx if sequence_idx is not None else index,
            event_type=event_type,
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
