from celery import shared_task

from flask_backend.services import GoldenGateUtils, ProtocolMaker
from flask_backend.models import ProtocolRequest

@shared_task(bind=True, name="generate_protocol_task")
def generate_protocol_task(self, req_dict):
    req = ProtocolRequest.model_validate(req_dict)
    total = len(req.sequences_to_domesticate)

    def progress_callback(step: str, message: str, progress: float, sequenceIdx: int = None):
        self.update_state(
            state="PROGRESS",
            meta={
                "step": step,
                "message": message,
                "progress": progress,
                "total": total,
                "sequenceIdx": sequenceIdx,
            },
        )


    result = None
    for idx, seq in enumerate(req.sequences_to_domesticate):
        protocolMaker = ProtocolMaker(
            request_idx=idx,
            sequence_to_domesticate=seq,
            codon_usage_dict=GoldenGateUtils().get_codon_usage_dict(req.species),
            max_mutations=req.max_mut_per_site,
            template_seq=req.template_sequence,
            kozak=req.kozak,
            max_results=req.max_results,
            verbose=req.verbose_mode,
            job_id=req.job_id,
        )
        protocolMaker.create_gg_protocol(progress_callback)

    return {"status": "completed", "result": result}

