# celery_tasks.py
from celery import shared_task, group
from flask_sse import sse  # This will use the same Redis connection configured in your Flask app
from flask_backend.services import GoldenGateUtils, ProtocolMaker
from flask_backend.models import ProtocolRequest # noqa

@shared_task(ignore_result=False)
def process_protocol_sequence(req_dict, index: int):
    # Validate and parse the request.
    req = ProtocolRequest.model_validate(req_dict)
    
    # Define a progress callback that publishes events via Flask‑SSE.
    def progress_callback(step: str, message: str, progress: float, sequenceIdx: int = None):
        if sequenceIdx is None:
            sequenceIdx = index
        # Publish the progress update to a Redis channel unique to this sequence.
        channel = f"job_{req.job_id}_{index}"
        sse.publish(
            {
                "step": step,
                "message": message,
                "progress": progress,
                "sequenceIdx": sequenceIdx,
            },
            type="progress",
            channel=channel
        )

    # Get the sequence to process.
    seq = req.sequences_to_domesticate[index]
    
    # Create a unique job id per sequence by appending the index.
    protocolMaker = ProtocolMaker(
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
    
    # Run the protocol – progress updates are published via the callback.
    result = protocolMaker.create_gg_protocol(progress_callback)
    return {"sequenceIdx": index, "result": result}

@shared_task(ignore_result=False)
def generate_protocol_task(req_dict):
    # Validate the incoming protocol request.
    req = ProtocolRequest.model_validate(req_dict)
    total = len(req.sequences_to_domesticate)
    
    # Launch a group of tasks—one per sequence.
    job = group(process_protocol_sequence.s(req_dict, idx) for idx in range(total))
    group_result = job.apply_async()
    
    # Return immediately with a status and total sequence count.
    return {"status": "started", "group_task_id": group_result.id, "total": total}
