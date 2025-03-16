from typing import Dict, List, Optional
from pydantic import BaseModel, RootModel
from .mutations import MutationOption

# --- Primer and Reaction models ---


class ReactionDetail(BaseModel):
    forward: str
    reverse: str
    amplicon_size: Optional[int] = None


class Primer(BaseModel):
    name: str
    sequence: str
    binding_region: Optional[str] = None
    tm: Optional[float] = None
    gc_content: Optional[float] = None
    length: Optional[int] = None


# --- Models for mutation_primers ---


class MutationPrimerPair(BaseModel):
    site: str
    position: int
    forward: Primer
    mutation_info: MutationOption

# --- EdgePrimers model ---


class EdgePrimers(BaseModel):
    forward_primer: Primer
    product_size: int
    reverse_primer: Primer


class PrimerDesignResult(RootModel[Dict[int, List[MutationPrimerPair]]]):
    pass
