from pydantic import BaseModel
from typing import List, Dict, Optional, Any
from models import Mutations, Primer, RestrictionSite, PCRReaction


class DomesticationResult(BaseModel):
    sequence_index: int = 0
    processed_sequence: str = ""
    mtk_part_left: str = ""
    mtk_part_right: str = ""
    restriction_sites: List[RestrictionSite] = []
    mutation_options: List[Mutations] = []
    edge_primers: List[Primer] = []
    mutation_primers: List[Primer] = []
    PCR_reactions: List[PCRReaction] = []
    messages: List[str] = []
    errors: Optional[Any] = None


class MTKDomesticationProtocol(BaseModel):
    result_data: Dict[int, DomesticationResult]
