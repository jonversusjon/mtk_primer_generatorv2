from typing import List, Optional
from pydantic import BaseModel
from models import Mutation


class Primer(BaseModel):
    name: str = ""
    sequence: str = ""
    binding_region: Optional[str] = None
    tm: Optional[float] = None
    gc_content: Optional[float] = None
    length: Optional[int] = None


class MutationPrimerPair(BaseModel):
    # forward / reverse for a single restriction site
    site: str
    position: int
    forward: Primer
    reverse: Primer
    mutation: Mutation


class MutationPrimerSet(BaseModel):
    # forward / reverse for all restriction sites
    mut_primer_pairs: List[MutationPrimerPair]
    

class EdgePrimerPair(BaseModel):
    forward: Primer
    reverse: Primer