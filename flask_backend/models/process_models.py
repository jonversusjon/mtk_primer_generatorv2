from typing import List, Optional
from flask_backend.models import Primer, Mutation, RestrictionSite, FrontendFriendly


class SequenceToDomesticate(FrontendFriendly):
    primer_name: Optional[str] = None
    sequence: str
    mtk_part_left: str
    mtk_part_right: str
    restriction_sites: Optional[List[RestrictionSite]] = None


class ProtocolRequest(FrontendFriendly):
    sequences_to_domesticate: List[SequenceToDomesticate]
    species: str = ""
    kozak: str = "MTK"
    max_mut_per_site: int = 3
    verbose_mode: bool = True
    template_sequence: str = ""
    max_results: str = "err"
    job_id: Optional[str] = None

        
class MutationPrimerPair(FrontendFriendly):
    # forward / reverse for a single restriction site
    site: str
    position: int
    forward: Primer
    reverse: Primer
    mutation: Mutation


class MutationPrimerSet(FrontendFriendly):
    # forward / reverse for all restriction sites
    mut_primer_pairs: List[MutationPrimerPair]