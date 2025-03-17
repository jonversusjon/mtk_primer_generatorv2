from pydantic import BaseModel, ConfigDict
from typing import List, Optional
from models import RestrictionSite

# Sequence models
def to_camel(string: str) -> str:
    parts = string.split('_')
    return parts[0] + ''.join(word.capitalize() for word in parts[1:])
        
    
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
