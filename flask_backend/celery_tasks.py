from celery import shared_task

from flask_backend.services import GoldenGateUtils
from flask_backend.models import ProtocolRequest

@shared_task(bind=True, name="generate_protocol_task")
def generate_protocol_task(self, req_dict):
    """Celery task that processes protocol generation."""
    req = ProtocolRequest.model_validate(req_dict)
    total = len(req.sequences_to_domesticate)
    for idx, seq in enumerate(req.sequences_to_domesticate):
        progress = int((idx / total) * 100)
        self.update_state(state="PROGRESS", meta={"step": idx + 1, "total": total, "progress": progress})
        result = GoldenGateUtils().process_sequence(
            job_id=self.request.id,
            seq=seq,
            kozak=req.kozak,
            max_mut_per_site=req.max_mut_per_site,
            verbose_mode=req.verbose_mode,
            template_sequence=req.template_sequence,
            max_results=req.max_results,
            codon_usage_dict=GoldenGateUtils().get_codon_usage_dict(req.species)
        )
    return {"status": "completed", "result": result}
