"""
Golden Gate Assembly Services Package

This package provides a collection of services for Golden Gate Assembly,
including sequence preparation, restriction site detection, mutation analysis,
optimization, primer design, and protocol generation.
"""

# Import version information
__version__ = "0.1.0"

# Sequence preparation
from .sequence_preparator import SequencePreparator

# Restriction site detection
from .rs_detector import RestrictionSiteDetector

# Mutation analysis and optimization
from .mut_analyzer import MutationAnalyzer
from .mut_optimizer import MutationOptimizer

# Primer design and selection
from .primer_designer import PrimerDesigner

# Reaction organization
from .reaction_organizer import ReactionOrganizer

# Protocol generation
from .protocol_maker import ProtocolMaker

# Utilities
from .utils import GoldenGateUtils

# Define what gets imported with 'from services import *'
__all__ = [
    'SequencePreparator',
    'RestrictionSiteDetector',
    'MutationAnalyzer',
    'MutationOptimizer',
    'PrimerDesigner',
    'PrimerSelector',
    'ReactionOrganizer',
    'ProtocolMaker',
    'GoldenGateUtils',
]