import math
import numpy as np
from pydantic import BaseModel, BeforeValidator, PlainSerializer, WithJsonSchema
from typing import Annotated, List, Optional, Any
from models import Codon

def validate_numpy_array(v: Any) -> np.ndarray:
    """
    Convert lists to a NumPy array (with dtype=int) and validate that the total number of elements 
    equals 4^N for some positive integer N.
    """
    if isinstance(v, np.ndarray):
        arr = v
    elif isinstance(v, list):
        arr = np.array(v, dtype=int)
    else:
        raise TypeError("compatibility must be a numpy ndarray or nested list of integers.")
    
    if arr.size == 0:
        raise ValueError("Array must not be empty.")
    
    num_elements = arr.size
    exponent = math.log(num_elements, 4)
    if not exponent.is_integer():
        raise ValueError("The total number of elements must be 4^N for some positive integer N.")
    return arr

def serialize_numpy_array(x: np.ndarray) -> dict:
    """
    Serializes the array into metadata:
      - snippet: first few entries from the flattened array,
      - shape: the shape of the array,
      - ones_percentage: percentage of elements equal to 1.
    """
    snippet_length = 5  # adjust the snippet length as desired
    flattened = x.flatten()
    snippet = flattened[:snippet_length].tolist()
    shape = list(x.shape)
    ones_count = int((x == 1).sum())
    total_count = x.size
    ones_percentage = (ones_count / total_count) * 100
    return {"snippet": snippet, "shape": shape, "ones_percentage": ones_percentage}

# Create a custom Annotated type that uses our validator and serializer.
NumpyArray = Annotated[
    np.ndarray,
    BeforeValidator(validate_numpy_array),
    PlainSerializer(serialize_numpy_array, return_type=dict),
    WithJsonSchema({
        "type": "object",
        "properties": {
            "snippet": {"type": "array", "items": {"type": "integer"}},
            "shape": {"type": "array", "items": {"type": "integer"}},
            "ones_percentage": {"type": "number"}
        },
        "required": ["snippet", "shape", "ones_percentage"]
    })
]

# Define a custom BaseModel that allows arbitrary types
class ConfiguredBaseModel(BaseModel):
    model_config = {"arbitrary_types_allowed": True}

# Mutation models using the custom base model
class OverhangOption(ConfiguredBaseModel):
    bottom_overhang: str
    top_overhang: str
    overhang_start_index: int

class MutationCodon(ConfiguredBaseModel):
    codon: Codon
    nth_codon_in_rs: int

class Mutation(ConfiguredBaseModel):
    mut_codons: List[MutationCodon]
    mut_indices_rs: Optional[List[int]] = None
    mut_indices_codon: Optional[List[int]] = None
    mut_context: str
    first_mut_idx: int
    last_mut_idx: int
    overhang_options: List[OverhangOption]

class MutationSet(ConfiguredBaseModel):
    mutations: List[Mutation]
    compatibility: NumpyArray

class MutationSetCollection(ConfiguredBaseModel):
    sites_to_mutate: List[str]
    sets: List[MutationSet]
    