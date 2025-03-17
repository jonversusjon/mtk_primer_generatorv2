from pydantic import BaseModel
from typing import List, Optional
from models import Codon, Mutation

# Restriction site model
class RestrictionSite(BaseModel):
    position: int
    frame: int
    codons: List[Codon]
    strand: str
    context_seq: str
    context_rs_indices: List[int]
    context_first_base: int
    context_last_base: int
    recognition_seq: str
    enzyme: str
    mutations: Optional[Mutation] = None