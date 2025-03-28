# celery_tasks.py
import json
from functools import partial
from typing import Optional, Any
from celery import shared_task, group
from flask_sse import sse
from flask_backend.services import GoldenGateUtils, ProtocolMaker
from flask_backend.models import ProtocolRequest, DomesticationResult, FrontendFriendly
from flask_backend.logging import logger

from pydantic import BaseModel

# Import your existing to_camel function
from flask_backend.models.base_models import to_camel

def process_payload_values(value: Any) -> Any:
    """
    Recursively process values to:
    1. Convert Pydantic models to dictionaries using model_dump(by_alias=True)
    2. Convert snake_case keys to camelCase in regular dictionaries
    """
    # If it's a Pydantic model, use model_dump with by_alias=True
    if isinstance(value, (BaseModel, FrontendFriendly)):
        return value.model_dump(by_alias=True)
    
    # If it's a list, process each item
    elif isinstance(value, list):
        return [process_payload_values(item) for item in value]
    
    # If it's a dictionary, process keys and values
    elif isinstance(value, dict):
        result = {}
        for k, v in value.items():
            # Convert keys to camelCase if they're snake_case
            camel_key = to_camel(k) if isinstance(k, str) and '_' in k else k
            # Process values recursively
            processed_value = process_payload_values(v)
            result[camel_key] = processed_value
        return result
    
    # Return other types as is
    else:
        return value

def publish_sse(
    job_id: str,
    sequence_idx: int,
    step: str,
    message: str,
    prog: Optional[float] = None,
    **kwargs):
    """
    Publish updates via Flask-SSE to a Redis channel.
    
    Handles a mix of Pydantic models and primitive types, ensuring all keys are camelCase.
    """
    channel = f"job_{job_id}_{sequence_idx}"
    
    # Process kwargs to handle Pydantic models and convert snake_case to camelCase
    processed_kwargs = process_payload_values(kwargs)
    
    # Create the base payload
    payload = {
        "jobId": job_id,
        "sequenceIdx": sequence_idx,
        "step": step,
        "message": message,
        **processed_kwargs  # Use the processed kwargs
    }
    
    # Add progress if provided
    if prog is not None:
        payload["stepProgress"] = prog

    logger.log_step("SSE Publish", f"Publishing to channel {channel}: {json.dumps(payload, default=str)}")
    
    try:
        # Verify the payload is serializable
        json.dumps(payload, default=str)
    except TypeError as e:
        # Log the error and the problematic data
        print(f"Serialization Error: {e}")
        print(f"Problematic payload: {payload}")
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
    if(index == 1):
        logger.log_step("SSE Publish", f"Setting up task for channel 1 job_{req.job_id}_{index}: {json.dumps({'jobId': req.job_id, 'sequenceIdx': index, 'step': 'start', 'message': 'Starting protocol generation'})}")
    # # Log all variables sent to ProtocolMaker
    # logger.log_step("ProtocolMaker Input", f"Request Index: {index}")
    # logger.log_step("ProtocolMaker Input", f"Sequence to Domesticate: {seq}")
    # logger.log_step("ProtocolMaker Input", f"Codon Usage Dict: {GoldenGateUtils().get_codon_usage_dict(req.species)}")
    logger.log_step("ProtocolMaker Input", f"Max Mutations: {req.max_mut_per_site}")
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
