# celery_tasks.py
import json
from functools import partial
from typing import Optional
from celery import shared_task, group
from flask_sse import sse
from flask_backend.services import GoldenGateUtils, ProtocolMaker
from flask_backend.models import ProtocolRequest, DomesticationResult
from flask_backend.logging import logger

def publish_sse(
    job_id: str,
    sequence_idx: int,
    step: str,
    message: str,
    step_progress: Optional[float] = None,
    **kwargs):
    """
    Publish updates via Flask-SSE to a Redis channel.
    
    Parameters:
        step (str): Current step name.
        message (str): Human-readable message.
        job_id (str): Unique identifier for the job.
        sequence_idx (int): Index of the current sequence.
        step_progress (Optional[float]): Progress percentage for the current step (0-100).
                                       If None, no progress bar is shown for this step.
        **kwargs: Arbitrary additional data to include in the payload.
    """
    channel = f"job_{job_id}_{sequence_idx}"
    payload = {
        "jobId": job_id,
        "sequenceIdx": sequence_idx,
        "step": step,
        "message": message,
        **kwargs
    }
    if step_progress is not None:
        payload["stepProgress"] = step_progress

    logger.log_step("SSE Publish", f"Publishing to channel {channel}: {json.dumps(payload, default=str)}")
    
    sse.publish(
        payload,
        channel=channel,
    )

    
@shared_task(ignore_result=False)
def process_protocol_sequence(req_dict: dict, index: int):
    req = ProtocolRequest.model_validate(req_dict)
    seq = req.sequences_to_domesticate[index]

    progress_callback = partial(publish_sse, req.job_id, index)

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
