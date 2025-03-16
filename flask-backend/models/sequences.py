from pydantic import BaseModel
from typing import List, Optional, Tuple


def to_camel(string: str) -> str:
    parts = string.split('_')
    return parts[0] + ''.join(word.capitalize() for word in parts[1:])

class Codon(BaseModel):
    amino_acid: str
    context_position: int
    codon_sequence: str
    rs_overlap: Tuple[int, int, int]


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
        
    
class SequenceToDomesticate(BaseModel):
    primer_name: Optional[str] = None
    sequence: str
    mtk_part_left: str
    mtk_part_right: str
    restriction_sites: Optional[List[RestrictionSite]] = None

    class Config:
        alias_generator = to_camel
        populate_by_name = True


class OverhangOption(BaseModel):
    bottom_overhang: str
    top_overhang: str
    overhang_start_index: int


class Overhangs(BaseModel):
    overhang_options: List[OverhangOption]
