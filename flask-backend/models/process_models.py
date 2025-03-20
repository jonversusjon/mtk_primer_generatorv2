from typing import List, Optional
from pydantic import BaseModel, ConfigDict
from models import Primer, Mutation, RestrictionSite, to_camel


class SequenceToDomesticate(BaseModel):
    primer_name: Optional[str] = None
    sequence: str
    mtk_part_left: str
    mtk_part_right: str
    restriction_sites: Optional[List[RestrictionSite]] = None

    model_config = ConfigDict(
        alias_generator=to_camel,
        populate_by_alias=True
    )
    
    
class MutationPrimerPair(BaseModel):
    # forward / reverse for a single restriction site
    site: str
    position: int
    forward: Primer
    reverse: Primer
    mutation: Mutation


class MutationPrimerSet(BaseModel):
    # forward / reverse for all restriction sites
    mut_primer_pairs: List[MutationPrimerPair]