from pydantic import BaseModel
from models import Primer


class PCRReaction(BaseModel):
    name: str
    forward_primer: Primer
    reverse_primer: Primer
    amplicon_size: int
