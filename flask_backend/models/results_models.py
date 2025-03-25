from typing import List, Dict, Optional, Any

from flask_backend.models import Mutation, Primer, RestrictionSite, MutationPrimerSet, NumpyArray, FrontendFriendly, FrontendNumpyFriendly


### Results ###
class MutationSet(FrontendNumpyFriendly):
    mutations: List[Mutation]
    compatibility: NumpyArray
    mut_primer_sets: List[MutationPrimerSet] = []


class MutationSetCollection(FrontendFriendly):
    sites_to_mutate: List[str]
    sets: List[MutationSet]
    

class PCRReaction(FrontendFriendly):
    name: str
    forward_primer: Primer
    reverse_primer: Primer
    amplicon_size: int

class EdgePrimerPair(FrontendFriendly):
    forward: Primer
    reverse: Primer


# Protocol model
class DomesticationResult(FrontendFriendly):
    sequence_index: int = -1
    max_results: str = "one"
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


class MTKDomesticationProtocol(FrontendFriendly):
    result_data: Dict[int, DomesticationResult]

    


    
