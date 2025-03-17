from pydantic import BaseModel
from typing import List

# Codon model
class Codon(BaseModel):
    amino_acid: str
    context_position: int
    codon_sequence: str
    rs_overlap: List[int]
    usage: float