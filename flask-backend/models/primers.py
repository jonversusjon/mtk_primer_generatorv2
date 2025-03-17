from typing import Dict, List, Optional
from pydantic import BaseModel, RootModel
from models import Mutation

# Primer model
class ReactionDetail(BaseModel):
    forward: str
    reverse: str
    amplicon_size: Optional[int] = None


class Primer(BaseModel):
    name: str = ""
    sequence: str = ""
    binding_region: Optional[str] = None
    tm: Optional[float] = None
    gc_content: Optional[float] = None
    length: Optional[int] = None


class MutationPrimerPair(BaseModel):
    site: str
    position: int
    forward: Primer
    mutation_info: Mutation


class EdgePrimerPair(BaseModel):
    forward: Primer
    reverse: Primer

class PrimerDesignResult(RootModel[Dict[int, List[MutationPrimerPair]]]):
    pass
