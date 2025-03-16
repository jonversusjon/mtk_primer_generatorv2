from pydantic import BaseModel
from typing import List, Dict, Optional, Any
from .mutations import Mutations
from .primers import EdgePrimers, MutationPrimerPair
from .sequences import RestrictionSite
from .pcr_reactions import PCRReaction

# --- Main ResultData item ---


class ResultDataItem(BaseModel):
    sequence_index: int
    mtk_part_left: str
    mtk_part_right: str
    restriction_sites: List[RestrictionSite]
    mutations: Mutations
    PCR_reactions: List[PCRReaction]
    messages: List[str]
    errors: Optional[Any] = None
    mutation_primers: Dict[int, List[MutationPrimerPair]]
    edge_primers: EdgePrimers

# --- Top-level model renamed to MTKDomesticationProtocol ---


class MTKDomesticationProtocol(BaseModel):
    result_data: Dict[int, ResultDataItem]
