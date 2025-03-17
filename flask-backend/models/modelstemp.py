# from pydantic import BaseModel, RootModel, ConfigDict, field_validator
# from typing import List, Optional, Dict, Any

# # Codon model
# class Codon(BaseModel):
#     amino_acid: str
#     context_position: int
#     codon_sequence: str
#     rs_overlap: List[int]
#     usage: float
    
    
# # Mutation models
# class OverhangOption(BaseModel):
#     bottom_overhang: str
#     top_overhang: str
#     overhang_start_index: int
    

# class Mutation(BaseModel):
#     alt_codons: List[Codon] # List of alternative codons for this particular mutation option.
#     mut_indices_rs: Optional[List[int]] = None
#     mut_indices_codon: Optional[List[int]] = None
#     mut_context:str
#     first_mut_idx: int
#     last_mut_idx: int
#     overhang_options: List[OverhangOption]
    
# class MutationSet(RootModel[List[Mutation]]):
#     mutations: List[Mutation]
#     compatibility: List

#     @field_validator("compatibility")
#     def validate_compatibility(cls, value):
#         def check_nested(item):
#             # If item is an int, it's a valid leaf
#             if isinstance(item, int):
#                 return True
#             # If it's a list, recursively check all items
#             if isinstance(item, list):
#                 return all(check_nested(i) for i in item)
#             return False

#         if not isinstance(value, list) or not all(check_nested(x)
#                                                   for x in value):
#             raise ValueError(
#                 "compatibility must be a nested list of integers (matrix "
#                 "shape 4^N)"
#             )
#         return value


# # Restriction site model
# class RestrictionSite(BaseModel):
#     position: int
#     frame: int
#     codons: List[Codon]
#     strand: str
#     context_seq: str
#     context_rs_indices: List[int]
#     context_first_base: int
#     context_last_base: int
#     recognition_seq: str
#     enzyme: str
#     mutations: Optional[Mutation] = None
    
# # Sequence models
# def to_camel(string: str) -> str:
#     parts = string.split('_')
#     return parts[0] + ''.join(word.capitalize() for word in parts[1:])
        
    
# class SequenceToDomesticate(BaseModel):
#     primer_name: Optional[str] = None
#     sequence: str
#     mtk_part_left: str
#     mtk_part_right: str
#     restriction_sites: Optional[List[RestrictionSite]] = None

#     model_config = ConfigDict(
#         alias_generator=to_camel,
#         populate_by_alias=True
#     )

# # Primer model
# class ReactionDetail(BaseModel):
#     forward: str
#     reverse: str
#     amplicon_size: Optional[int] = None


# class Primer(BaseModel):
#     name: str = ""
#     sequence: str = ""
#     binding_region: Optional[str] = None
#     tm: Optional[float] = None
#     gc_content: Optional[float] = None
#     length: Optional[int] = None


# class MutationPrimerPair(BaseModel):
#     site: str
#     position: int
#     forward: Primer
#     mutation_info: Mutation


# class EdgePrimerPair(BaseModel):
#     forward: Primer
#     reverse: Primer

# class PrimerDesignResult(RootModel[Dict[int, List[MutationPrimerPair]]]):
#     pass


# # PCR Reaction model
# class PCRReaction(BaseModel):
#     name: str
#     forward_primer: Primer
#     reverse_primer: Primer
#     amplicon_size: int
#     product_size: int
#     melting_temp: float
#     gc_content: float


# # Protocol model
# class DomesticationResult(BaseModel):
#     sequence_index: int = 0
#     processed_sequence: str = ""
#     mtk_part_left: str = ""
#     mtk_part_right: str = ""
#     restriction_sites: List[RestrictionSite] = []
#     mutation_options: List[Mutation] = []
#     edge_primers: List[Primer] = []
#     mutation_primers: List[Primer] = []
#     PCR_reactions: List[PCRReaction] = []
#     messages: List[str] = []
#     errors: Optional[Any] = None


# class MTKDomesticationProtocol(BaseModel):
#     result_data: Dict[int, DomesticationResult]



    


