from pydantic import BaseModel, field_validator, Field
from typing import List, Dict, Optional, Any

# --- Bottom-level models ---


class BottomOverhang(BaseModel):
    seq: str
    overhang_start_index: int


class OverhangOption(BaseModel):
    bottom_overhang: BottomOverhang
    overhang_start_index: int


class Overhangs(BaseModel):
    overhang_options: List[OverhangOption]

# --- Mutation option used in both mutations and mutation_primers ---


class MutationOption(BaseModel):
    alternative_codon_sequence: str
    mutated_base_index: int
    original_codon_sequence: Optional[str] = None
    overhangs: Overhangs
    mutation_positions_in_context: List[int]

    class Config:
        extra = "allow"  # allow additional keys as per schema

# --- Mutations container ---


class Mutations(BaseModel):
    all_mutation_options: List[List[MutationOption]]
    # We'll store the nested 4^N matrix as a generic list of unknown nesting
    compatibility: List

    @field_validator("compatibility")
    def validate_compatibility(cls, value):
        def check_nested(item):
            # If item is an int, it's a valid leaf
            if isinstance(item, int):
                return True
            # If it's a list, recursively check all items
            if isinstance(item, list):
                return all(check_nested(i) for i in item)
            return False

        if not isinstance(value, list) or not all(check_nested(x) for x in value):
            raise ValueError(
                "compatibility must be a nested list of integers (matrix shape 4^N)"
            )
        return value


class MutationEntry(BaseModel):
    site: str
    position: int
    original_codon_sequence: str
    alternative_codon_sequence: str
    mutated_base_index: Optional[int]
    overhangs: Overhangs
    mutated_context: str
    mutation_positions_in_context: List[int]


class MutationSet(BaseModel):
    __root__: List[MutationEntry]

# --- Codon and RestrictionSite models ---


class Codon(BaseModel):
    amino_acid: str
    context_position: int
    # Optionally include any additional field (for example, a codon sequence)
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

# --- Primer and Reaction models ---


class ReactionDetail(BaseModel):
    forward: str
    reverse: str


class Primer(BaseModel):
    name: str
    sequence: str
    binding_region: Optional[str] = None
    tm: Optional[float] = None
    gc_content: Optional[float] = None
    length: Optional[int] = None

    class Config:
        extra = "allow"

# --- Models for mutation_primers ---


class MutationPrimer(BaseModel):
    site: str
    position: int
    forward: Primer
    mutation_info: MutationOption

# --- EdgePrimers model ---


class EdgePrimers(BaseModel):
    forward_primer: Primer
    product_size: int
    # Including a reverse_primer if needed (optional)
    reverse_primer: Optional[Primer] = None


class PrimerDesignResult(BaseModel):
    # The final output is {set_index: [MutationPrimer, MutationPrimer, ...], ...}
    __root__: Dict[int, List[MutationPrimer]]

# --- Main ResultData item ---


class ResultDataItem(BaseModel):
    sequence_index: int
    # Based on the snippet, we assume there are three parts for "mtk":
    mtk_part_left: str
    mtk_part_center: str
    mtk_part_right: str
    restriction_sites: List[RestrictionSite]
    mutations: Mutations
    # The PCR_reactions field is nested: first key (int) -> second key (int) -> reaction name (str) -> ReactionDetail
    PCR_reactions: Dict[int, Dict[int, Dict[str, ReactionDetail]]]
    messages: List[str]
    errors: Optional[Any] = None
    mutation_primers: Dict[int, List[MutationPrimer]]
    edge_primers: EdgePrimers

# --- Top-level model renamed to MTKDomesticationProtocol ---


class MTKDomesticationProtocol(BaseModel):
    result_data: Dict[int, ResultDataItem]
