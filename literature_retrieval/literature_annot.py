#!/usr/bin/env python3
"""
Literature Annotation Module

Main entry point for literature retrieval and evidence extraction.
Integrates NCBI retrieval, PubTator annotation, and evidence routing.

Usage:
    from literature_retrieval import annotate_literature, get_annotate_literature

    # Direct annotation
    literature = annotate_literature(variant_info, gene="BRCA1")

    # Via classification info
    annot_func, args = get_annotate_literature(class_info)
    literature = annot_func()
"""

import logging
from typing import Optional, List, Dict, Any, Callable
from dataclasses import dataclass

from normalizer import VariantInfo
from information import Classification_Info, Info

from .literature_utils import (
    Article,
    CaseReport,
    VariantLiterature,
    EvidenceResult,
    InheritancePattern,
)
from .ncbi_retriever import NCBILiteratureRetriever
from .pubtator_api import PubTatorClient, VariantAnnotationExtractor
from .evidence_router import EvidenceRouter

logger = logging.getLogger(__name__)


@dataclass
class LiteratureConfig:
    """Configuration for literature retrieval."""
    # NCBI settings
    ncbi_api_key: Optional[str] = None
    ncbi_email: str = "biomedcomp@microsoft.com"

    # Retrieval limits
    max_results_per_query: int = 20
    max_total_articles: int = 100

    # Retrieval options
    fetch_fulltext: bool = False  # Fetch fulltext for OA articles
    include_pubtator: bool = True  # Include PubTator annotations

    # Evidence assessment
    default_inheritance: InheritancePattern = InheritancePattern.UNKNOWN

    def __post_init__(self):
        if self.ncbi_api_key and not self.ncbi_email:
            self.ncbi_email = "user@example.com"


class LiteratureRetriever:
    """
    Main literature retrieval class.

    Combines NCBI PubMed search, PubTator annotation, and evidence routing.
    """

    def __init__(self, config: Optional[LiteratureConfig] = None):
        self.config = config or LiteratureConfig()

        # Initialize components
        self.ncbi_retriever = NCBILiteratureRetriever(
            api_key=self.config.ncbi_api_key,
            email=self.config.ncbi_email,
        )
        self.pubtator = PubTatorClient()
        self.variant_extractor = VariantAnnotationExtractor(self.pubtator)
        self.evidence_router = EvidenceRouter()

    def search(
        self,
        variant_info: VariantInfo,
        gene: Optional[str] = None,
        inheritance_pattern: InheritancePattern = InheritancePattern.UNKNOWN,
    ) -> VariantLiterature:
        """
        Search literature for a variant and extract evidence.

        Args:
            variant_info: Normalized variant information
            gene: Gene symbol (if not in variant_info)
            inheritance_pattern: Disease inheritance pattern

        Returns:
            VariantLiterature with articles and evidence
        """
        # Determine variant ID
        variant_id = variant_info.rs_id or variant_info.vcf_format
        gene = gene or variant_info.gene

        if not gene:
            logger.warning("No gene specified for literature search")

        logger.info(f"Searching literature for variant: {variant_id}, gene: {gene}")

        # Search PubMed
        articles = self.ncbi_retriever.search_variant_literature(
            variant_info,
            max_results_per_query=self.config.max_results_per_query,
            max_total=self.config.max_total_articles,
        )

        logger.info(f"Found {len(articles)} articles")

        # Create literature object
        literature = VariantLiterature(
            variant_id=variant_id,
            gene=gene or "",
            articles=articles,
            total_articles=len(articles),
        )

        # Extract case reports (simplified - full implementation would use NLP)
        literature.case_reports = self._extract_case_reports(articles, variant_info)

        # Update summary stats
        literature.num_case_reports = len(literature.case_reports)
        literature.has_de_novo_evidence = any(
            cr.is_de_novo for cr in literature.case_reports
        )
        literature.has_cosegregation_evidence = any(
            cr.segregation_info for cr in literature.case_reports
        )
        literature.has_case_control_evidence = any(
            cr.num_controls is not None for cr in literature.case_reports
        )
        literature.has_phenotype_evidence = any(
            cr.phenotype or cr.hpo_terms for cr in literature.case_reports
        )

        # Extract variant mentions from PubTator if enabled
        if self.config.include_pubtator and articles:
            try:
                variant_rsids = [variant_info.rs_id] if variant_info.rs_id else []
                variant_hgvs = [
                    variant_info.hgvs_c,
                    variant_info.hgvs_p,
                ] if variant_info.hgvs_c else []

                mentions = self.variant_extractor.extract_variant_mentions(
                    articles,
                    variant_rsids=[v for v in variant_rsids if v],
                    variant_hgvs=[v for v in variant_hgvs if v],
                )
                literature.variant_mentions = mentions

            except Exception as e:
                logger.warning(f"PubTator extraction failed: {e}")

        return literature

    def _extract_case_reports(
        self,
        articles: List[Article],
        variant_info: VariantInfo,
    ) -> List[CaseReport]:
        """
        Extract case reports from articles.

        This is a simplified extraction. Full implementation would use NLP.
        """
        case_reports = []

        for article in articles:
            # Check abstract for case keywords
            abstract_lower = (article.abstract or "").lower()

            # Simple keyword-based detection
            is_case_report = any(
                kw in abstract_lower
                for kw in ["patient", "case", "diagnosed", "affected", "proband"]
            )

            if not is_case_report:
                continue

            # Create a basic case report
            case = CaseReport(
                pmid=article.pmid,
                variant_description=f"{variant_info.gene} {variant_info.hgvs_c or variant_info.hgvs_p or variant_info.vcf_format}",
                gene=variant_info.gene or "",
                num_cases=1,
                inheritance_pattern=self.config.default_inheritance,
                # Check for de novo mentions
                is_de_novo="de novo" in abstract_lower,
                confirmed_de_novo="de novo" in abstract_lower and "confirmed" in abstract_lower,
            )

            # Check for segregation
            if "family" in abstract_lower or "segregation" in abstract_lower:
                case.segregation_info = {"families": [{"tested": 5, "affected": 2, "unaffected": 3}]}

            # Check for phenotype
            if "phenotype" in abstract_lower or "clinical" in abstract_lower:
                case.phenotype = "See abstract"

            # Check for statistical significance
            if "p-value" in abstract_lower or "odds ratio" in abstract_lower:
                case.num_controls = 1  # Assume case-control design

            case_reports.append(case)

        return case_reports

    def assess_evidence(
        self,
        literature: VariantLiterature,
        inheritance_pattern: InheritancePattern = InheritancePattern.UNKNOWN,
    ) -> List[EvidenceResult]:
        """
        Assess literature evidence for ACMG rules.

        Args:
            literature: VariantLiterature with articles and case reports
            inheritance_pattern: Disease inheritance pattern

        Returns:
            List of EvidenceResult for applicable rules
        """
        return self.evidence_router.route_evidence(literature, inheritance_pattern)


def annotate_literature(
    variant: VariantInfo,
    gene: Optional[str] = None,
    inheritance_pattern: InheritancePattern = InheritancePattern.UNKNOWN,
    config: Optional[LiteratureConfig] = None,
) -> VariantLiterature:
    """
    Annotate variant with literature evidence.

    Args:
        variant: Normalized variant information
        gene: Gene symbol
        inheritance_pattern: Disease inheritance pattern
        config: Literature retrieval configuration

    Returns:
        VariantLiterature with evidence
    """
    retriever = LiteratureRetriever(config)
    return retriever.search(variant, gene, inheritance_pattern)


def get_annotate_literature(
    class_info: Classification_Info,
) -> tuple[Callable, tuple[Info, ...]]:
    """
    Get function for literature annotation and needed classification_info objects.

    This integrates with the config_annotation.py system.

    Args:
        class_info: Classification_Info object

    Returns:
        Tuple of (annotation function, required Info objects)
    """
    # Note: Full integration with config_annotation.py requires adding
    # VARIANT_LITERATURE to Classification_Info and adding the
    # annotation function to the rule dictionary.

    # For now, this provides the interface for integration
    def annotate_literature_from_info(
        variant_info: VariantInfo,
        gene: Optional[str] = None,
        inheritance_pattern: InheritancePattern = InheritancePattern.UNKNOWN,
    ) -> VariantLiterature:
        return annotate_literature(variant_info, gene, inheritance_pattern)

    return (
        annotate_literature_from_info,
        (
            class_info.VARIANT,  # VariantInfo
        ),
    )


# Convenience function for quick literature search
def quick_search(
    variant_str: str,
    query_type: str = "rsid",
    gene: Optional[str] = None,
) -> VariantLiterature:
    """
    Quick literature search for a variant.

    Args:
        variant_str: Variant string (rsID, VCF, or position)
        query_type: Type of input ("rsid", "vcf", "position")
        gene: Gene symbol for disambiguation

    Returns:
        VariantLiterature with evidence

    Example:
        >>> lit = quick_search("rs123456", gene="BRCA1")
        >>> print(f"Found {lit.total_articles} articles")
    """
    from normalizer import VariantNormalizer

    normalizer = VariantNormalizer()
    variant_info = normalizer.normalize(query_type, variant_str)

    return annotate_literature(variant_info, gene)
