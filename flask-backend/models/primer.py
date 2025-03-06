# models/primer.py
from dataclasses import dataclass, field
from typing import Optional, List

__all__ = ['Primer', 'MutationPrimer', 'EdgePrimer', 'PrimerSet']


@dataclass
class Primer:
    name: str
    sequence: str
    binding_region: Optional[str] = None
    tm: Optional[float] = None
    gc_content: Optional[float] = None
    length: Optional[int] = None


@dataclass
class MutationPrimer:
    site: str
    position: int  # Position is used for ordering
    forward: Primer
    reverse: Primer
    mutation_info: dict


@dataclass
class EdgePrimer:
    forward: Primer
    reverse: Primer
    product_size: Optional[int] = None


@dataclass
class PrimerSet:
    edge: EdgePrimer
    mutations: List[MutationPrimer] = field(default_factory=list)
