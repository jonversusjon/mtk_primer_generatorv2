from pydantic import BaseModel
from typing import List, Dict, Optional, Any
from models import Mutation, Primer, RestrictionSite, PCRReaction, EdgePrimerPair, MutationPrimerSet


# Protocol model
class DomesticationResult(BaseModel):
    sequence_index: int = -1
    max_results: int = 1
    processed_sequence: str = ""
    mtk_part_left: str = ""
    mtk_part_right: str = ""
    restriction_sites: List[RestrictionSite] = []
    mutation_options: List[Mutation] = []
    edge_primers: EdgePrimerPair = EdgePrimerPair(forward=Primer(), reverse=Primer())
    mut_primers: List[MutationPrimerSet] = []
    PCR_reactions: List[PCRReaction] = []
    messages: List[str] = []
    errors: Optional[Any] = None


class MTKDomesticationProtocol(BaseModel):
    result_data: Dict[int, DomesticationResult]
