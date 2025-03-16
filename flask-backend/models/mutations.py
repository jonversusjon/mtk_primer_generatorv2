from pydantic import BaseModel, RootModel, field_validator
from typing import List, Optional, Tuple, Any
from .sequences import Overhangs

# --- Mutation option used in both mutations and mutation_primers ---


class AlternativeCodon(BaseModel):
    seq: str  # The alternative codon sequence.
    mutations: Optional[Tuple[int, int, int]] = None
    mutated_context: Optional[str] = None
    sticky_ends: Optional[Overhangs] = None


class CodonMutation(BaseModel):
    original_codon_sequence: str
    context_position: Optional[int]
    amino_acid: Optional[str]
    alternative_codons: List[AlternativeCodon]


class MutationSite(BaseModel):
    position: int
    sequence: str
    frame: int
    strand: str
    enzyme: str
    codons: List[CodonMutation]


class MutationOption(BaseModel):
    alternative_codon_sequence: str
    mutated_base_index: int
    codon_sequence: Optional[str] = None
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
    seq: str
    alternative_codon_sequence: str
    mutated_base_index: Optional[int]
    overhangs: Overhangs
    mutated_context: str
    mutation_positions_in_context: List[int]


class MutationSet(RootModel[List[MutationEntry]]):
    pass
