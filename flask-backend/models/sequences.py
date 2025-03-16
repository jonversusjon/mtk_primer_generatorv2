from pydantic import BaseModel
from typing import List, Optional, Any


def to_camel(string: str) -> str:
    parts = string.split('_')
    return parts[0] + ''.join(word.capitalize() for word in parts[1:])


class Codon(BaseModel):
    amino_acid: str
    context_position: int
    codon_sequence: Optional[str] = None

    class Config:
        extra = "allow"


class RestrictionSite(BaseModel):
    frame: int
    context_sequence: str
    codons: List[Codon]
    context_recognition_site_indices: List[int]

    class Config:
        extra = "allow"


class SequenceToDomesticate(BaseModel):
    primer_name: Optional[str] = None
    sequence: str
    mtk_part_left: str
    mtk_part_right: str
    restriction_sites: Optional[List[RestrictionSite]] = None

    class Config:
        alias_generator = to_camel
        allow_population_by_field_name = True


class OverhangOption(BaseModel):
    bottom_overhang: str
    top_overhang: str
    overhang_start_index: int


class Overhangs(BaseModel):
    overhang_options: List[OverhangOption]
