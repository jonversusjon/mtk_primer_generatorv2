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
    prog: Optional[float] = None,
    **kwargs):
    """
    Publish updates via Flask-SSE to a Redis channel.
    
    Parameters:
        step (str): Current step name.
        message (str): Human-readable message.
        job_id (str): Unique identifier for the job.
        sequence_idx (int): Index of the current sequence.
        prog (Optional[float]): Progress percentage for the current step (0-100).
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
    if prog is not None:
        payload["stepProgress"] = prog

    logger.log_step("SSE Publish", f"Publishing to channel {channel}: {json.dumps(payload, default=str)}")
    try:
        # Attempt to serialize keyword arguments
        json.dumps(kwargs)
    except TypeError as e:
        # Log the error and the problematic data
        print(f"Serialization Error: {e}")
        print(f"Problematic data: {kwargs}")
        # Decide how to handle: maybe send only serializable parts,
        # or use a custom serializer like the one in utils.py
        # For now, maybe just send message and progress
        # Or raise the error to halt execution
        raise e
    
    sse.publish(
        payload,
        channel=channel,
    )

    
@shared_task(ignore_result=False)
def process_protocol_sequence(req_dict: dict, index: int):
    req = ProtocolRequest.model_validate(req_dict)
    seq = req.sequences_to_domesticate[index]

    progress_callback = partial(publish_sse, job_id=req.job_id, sequence_idx=index)

    # # Log all variables sent to ProtocolMaker
    # logger.log_step("ProtocolMaker Input", f"Request Index: {index}")
    # logger.log_step("ProtocolMaker Input", f"Sequence to Domesticate: {seq}")
    # logger.log_step("ProtocolMaker Input", f"Codon Usage Dict: {GoldenGateUtils().get_codon_usage_dict(req.species)}")
    # logger.log_step("ProtocolMaker Input", f"Max Mutations: {req.max_mut_per_site}")
    # logger.log_step("ProtocolMaker Input", f"Template Sequence: {req.template_sequence}")
    # logger.log_step("ProtocolMaker Input", f"Kozak: {req.kozak}")
    # logger.log_step("ProtocolMaker Input", f"Max Results: {req.max_results}")
    # logger.log_step("ProtocolMaker Input", f"Verbose Mode: {req.verbose_mode}")
    # logger.log_step("ProtocolMaker Input", f"Job ID: {req.job_id}_{index}")

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
    result: DomesticationResult = protocol_maker.create_gg_protocol(send_update=progress_callback)
    
    return {"sequenceIdx": index, "result": result.model_dump(by_alias=True)}

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
