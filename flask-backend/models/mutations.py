from pydantic import BaseModel, field_validator
from typing import List, Optional
from models import Codon


# Mutation models
class OverhangOption(BaseModel):
    bottom_overhang: str
    top_overhang: str
    overhang_start_index: int
    
    
class Mutation(BaseModel):
    alt_codons: List[Codon] # List of alternative codons for this particular mutation option.
    mut_indices_rs: Optional[List[int]] = None
    mut_indices_codon: Optional[List[int]] = None
    mut_context:str
    first_mut_idx: int
    last_mut_idx: int
    overhang_options: List[OverhangOption]
    
    
class MutationSet(BaseModel):
    mutations: List[Mutation]
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

        if not isinstance(value, list) or not all(check_nested(x)
                                                  for x in value):
            raise ValueError(
                "compatibility must be a nested list of integers (matrix "
                "shape 4^N)"
            )
        return value