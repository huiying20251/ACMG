#!/usr/bin/env python3
"""
Literature Retrieval Module

Provides literature search and evidence extraction for ACMG variant classification:
- NCBI PubMed literature retrieval
- PubTator API variant annotation retrieval
- Evidence routing based on inheritance pattern

Usage:
    from literature_retrieval import LiteratureRetriever, EvidenceRouter

    # Search literature for a variant
    retriever = LiteratureRetriever()
    articles = retriever.search_variant_literature(variant_info)

    # Route evidence based on inheritance pattern
    router = EvidenceRouter()
    evidence = router.route_evidence(articles, inheritance_pattern="AD")
"""

from .literature_utils import (
    Article,
    VariantLiterature,
    CaseReport,
    EvidenceResult,
    InheritancePattern,
    IndividualVariantObservation,
)
from .ncbi_retriever import NCBILiteratureRetriever
from .pubtator_api import PubTatorClient
from .evidence_router import EvidenceRouter
from .literature_annot import annotate_literature, get_annotate_literature
from .pm3_evaluator import (
    PM3Strength,
    PM3Assessment,
    PM3LLMEvaluator,
    assess_pm3_from_literature,
)
from .ps2_evaluator import (
    PS2Strength,
    PS2Assessment,
    PS2LLMEvaluator,
    DeNovoStatus,
    assess_ps2_from_literature,
)
from .segregation_evaluator import (
    SegregationStrength,
    SegregationStatus,
    SegregationAssessment,
    FamilyMember,
    PedigreeData,
    SegregationEvaluator,
    assess_pp1_from_pedigree,
    assess_bs4_from_pedigree,
)
from .bp2_evaluator import (
    BP2Strength,
    BP2Assessment,
    TransPartnerVariant,
    BP2Evaluator,
    assess_bp2_from_partners,
)
from .ps4_evaluator import (
    PS4Strength,
    DiseaseRarityCategory,
    PS4Assessment,
    PS4LLMEvaluator,
    assess_ps4_from_cases,
)
from .evidence_cache import (
    EvidenceCacheDB,
    EvidenceCacheConfig,
    cache_evidence,
    get_cached_evidence,
    has_cached_evidence,
    get_cache,
)

__all__ = [
    "Article",
    "VariantLiterature",
    "CaseReport",
    "EvidenceResult",
    "InheritancePattern",
    "IndividualVariantObservation",
    "NCBILiteratureRetriever",
    "PubTatorClient",
    "EvidenceRouter",
    "annotate_literature",
    "get_annotate_literature",
    "PM3Strength",
    "PM3Assessment",
    "PM3LLMEvaluator",
    "assess_pm3_from_literature",
    "PS2Strength",
    "PS2Assessment",
    "PS2LLMEvaluator",
    "DeNovoStatus",
    "assess_ps2_from_literature",
    # Segregation (PP1/BS4)
    "SegregationStrength",
    "SegregationStatus",
    "SegregationAssessment",
    "FamilyMember",
    "PedigreeData",
    "SegregationEvaluator",
    "assess_pp1_from_pedigree",
    "assess_bs4_from_pedigree",
    # BP2
    "BP2Strength",
    "BP2Assessment",
    "TransPartnerVariant",
    "BP2Evaluator",
    "assess_bp2_from_partners",
    # PS4
    "PS4Strength",
    "DiseaseRarityCategory",
    "PS4Assessment",
    "PS4LLMEvaluator",
    "assess_ps4_from_cases",
    # Evidence Cache
    "EvidenceCacheDB",
    "EvidenceCacheConfig",
    "cache_evidence",
    "get_cached_evidence",
    "has_cached_evidence",
    "get_cache",
]
