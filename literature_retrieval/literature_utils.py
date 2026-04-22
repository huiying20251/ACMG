#!/usr/bin/env python3
"""
Literature Retrieval Utilities

Data classes for literature retrieval and evidence extraction.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, List, Dict, Any


class InheritancePattern(Enum):
    """Inheritance pattern for evidence routing."""
    AUTOSOMAL_DOMINANT = "AD"
    AUTOSOMAL_RECESSIVE = "AR"
    X_LINKED = "X"
    UNKNOWN = "UNKNOWN"


@dataclass
class Article:
    """Represents a PubMed article."""
    pmid: str
    title: str
    abstract: Optional[str] = None
    journal: Optional[str] = None
    first_author: Optional[str] = None
    pub_year: Optional[int] = None
    doi: Optional[str] = None
    pmcid: Optional[str] = None
    citation: Optional[str] = None
    link: Optional[str] = None
    # Open access info
    can_access: bool = False
    is_oa: bool = False
    license: Optional[str] = None
    # Full text (if available)
    fulltext_xml: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pmid": self.pmid,
            "title": self.title,
            "abstract": self.abstract,
            "journal": self.journal,
            "first_author": self.first_author,
            "pub_year": self.pub_year,
            "doi": self.doi,
            "pmcid": self.pmcid,
            "citation": self.citation,
            "link": self.link,
            "can_access": self.can_access,
            "is_oa": self.is_oa,
            "license": self.license,
        }


@dataclass
class CaseReport:
    """
    Represents a case report extracted from literature.

    Used for evidence rules:
    - PS4: Case control studies (multiple affected individuals)
    - PS2: de novo variants (confirmed in proband but not parents)
    - PP1: Cosegregation with disease in family
    - PP4: Phenotype specific to variant
    """
    pmid: str
    variant_description: str
    gene: str

    # Case counts
    num_cases: int = 1
    num_controls: Optional[int] = None  # For PS4

    # Inheritance
    inheritance: Optional[str] = None  # AD, AR, etc.
    inheritance_pattern: InheritancePattern = InheritancePattern.UNKNOWN

    # Segregation data
    segregation_info: Optional[Dict[str, Any]] = None
    # {"families": [{"affected": 5, "unaffected": 3, "tested": 8}]}

    # Phenotype
    phenotype: Optional[str] = None
    hpo_terms: List[str] = field(default_factory=list)

    # de novo status
    is_de_novo: bool = False
    confirmed_de_novo: bool = False  # Confirmed by parent testing

    # Statistical significance (if reported)
    p_value: Optional[float] = None
    odds_ratio: Optional[float] = None

    # Confidence level
    confidence: str = "medium"  # high, medium, low

    @property
    def applicable_rules(self) -> List[str]:
        """Get list of potentially applicable ACMG rules based on case data."""
        rules = []

        if self.is_de_novo:
            rules.append("PS2")

        if self.num_controls is not None and self.num_cases > 0:
            rules.append("PS4")

        if self.segregation_info:
            rules.append("PP1")

        if self.phenotype or self.hpo_terms:
            rules.append("PP4")

        return rules


@dataclass
class IndividualVariantObservation:
    """
    Individual-level variant observation (per healthfutures-evagg design).

    Represents a single individual's observation of a specific variant.
    This is the atomic unit for ACMG evidence assessment.

    For PS2/PM3/PP1 rules, each observation is assessed individually,
    then aggregated for evidence strength determination.
    """
    # Identification
    pmid: str
    individual_id: str           # Patient identifier: "patient 1", "IV-1", "proband", etc.
    variant_description: str      # Original variant text from paper
    gene: str

    # Variant inheritance and zygosity
    variant_inheritance: str = "unknown"  # "inherited", "de novo", "unknown"
    parental_testing: bool = False        # Whether parents were tested
    zygosity: str = "unknown"            # "homozygous", "heterozygous", "compound_heterozygous", "unknown"

    # Variant classification from ClinVar (string format for output)
    clinvar_status: str = ""    # e.g., "Pathogenic", "Likely pathogenic", "VUS", "Likely benign", "Benign"
    clinvar_significance: str = ""  # e.g., "single hit", "double hit" for AR genes

    # Phenotype information
    phenotype: str = ""
    hpo_terms: List[str] = field(default_factory=list)

    # Segregation data
    segregation_data: Optional[Dict[str, Any]] = None
    # {"families": [{"id": "family 1", "affected_carriers": 5, "unaffected_carriers": 3, "total_tested": 8, "lod_score": 2.5}]}

    # Confidence level
    confidence: str = "medium"  # high, medium, low

    def to_string(self) -> str:
        """
        Convert observation to a human-readable string for LLM analysis.

        Returns:
            Formatted string representation of this observation.
        """
        parts = [
            f"PMID:{self.pmid}",
            f"Individual:{self.individual_id}",
            f"Variant:{self.variant_description}",
            f"Gene:{self.gene}",
            f"Inheritance:{self.variant_inheritance}",
            f"ParentalTesting:{'Yes' if self.parental_testing else 'No'}",
            f"Zygosity:{self.zygosity}",
        ]

        if self.clinvar_status:
            parts.append(f"ClinVar:{self.clinvar_status}")
        if self.clinvar_significance:
            parts.append(f"ClinVarSignificance:{self.clinvar_significance}")

        if self.phenotype:
            parts.append(f"Phenotype:{self.phenotype}")
        if self.hpo_terms:
            parts.append(f"HPO:{'|'.join(self.hpo_terms)}")

        if self.segregation_data:
            seg_str = self._format_segregation()
            parts.append(f"Segregation:{seg_str}")

        parts.append(f"Confidence:{self.confidence}")

        return "; ".join(parts)

    def _format_segregation(self) -> str:
        """Format segregation data as string."""
        if not self.segregation_data:
            return ""

        families = self.segregation_data.get("families", [])
        parts = []
        for fam in families:
            fam_parts = []
            if "id" in fam:
                fam_parts.append(f"Family:{fam['id']}")
            if "affected_carriers" in fam:
                fam_parts.append(f"Affected:{fam['affected_carriers']}")
            if "unaffected_carriers" in fam:
                fam_parts.append(f"Unaffected:{fam['unaffected_carriers']}")
            if "total_tested" in fam:
                fam_parts.append(f"Tested:{fam['total_tested']}")
            if "lod_score" in fam:
                fam_parts.append(f"LOD:{fam['lod_score']}")
            parts.append(";".join(fam_parts))
        return "|".join(parts)

    @property
    def is_de_novo(self) -> bool:
        """Check if variant is reported as de novo."""
        return self.variant_inheritance.lower() == "de novo"

    @property
    def is_confirmed_de_novo(self) -> bool:
        """Check if de novo is confirmed by parental testing."""
        return self.is_de_novo and self.parental_testing

    @property
    def applicable_rules(self) -> List[str]:
        """Get list of potentially applicable ACMG rules."""
        rules = []
        if self.is_de_novo:
            rules.append("PS2")
        if self.zygosity in ("homozygous", "compound_heterozygous"):
            rules.append("PM3")
        if self.segregation_data:
            rules.append("PP1")
        if self.phenotype or self.hpo_terms:
            rules.append("PP4")
        return rules


@dataclass
class VariantLiterature:
    """
    Aggregated literature data for a variant.

    Contains all retrieved articles and extracted case reports.
    """
    variant_id: str  # rsID or VCF format
    gene: str

    # Retrieved articles
    articles: List[Article] = field(default_factory=list)
    total_articles: int = 0

    # Extracted case reports (article-level, legacy)
    case_reports: List[CaseReport] = field(default_factory=list)

    # Extracted individual-level observations (per patient)
    individual_observations: List[IndividualVariantObservation] = field(default_factory=list)

    # Summary statistics
    num_case_reports: int = 0
    num_functional_studies: int = 0
    num_review_articles: int = 0

    # Evidence summary
    has_de_novo_evidence: bool = False
    has_cosegregation_evidence: bool = False
    has_case_control_evidence: bool = False
    has_phenotype_evidence: bool = False

    def to_dict(self) -> Dict[str, Any]:
        return {
            "variant_id": self.variant_id,
            "gene": self.gene,
            "total_articles": self.total_articles,
            "case_reports": [
                {
                    "pmid": cr.pmid,
                    "variant_description": cr.variant_description,
                    "num_cases": cr.num_cases,
                    "inheritance": cr.inheritance,
                    "applicable_rules": cr.applicable_rules,
                }
                for cr in self.case_reports
            ],
            "individual_observations": [
                obs.to_string() for obs in self.individual_observations
            ],
            "evidence_flags": {
                "de_novo": self.has_de_novo_evidence,
                "cosegregation": self.has_cosegregation_evidence,
                "case_control": self.has_case_control_evidence,
                "phenotype": self.has_phenotype_evidence,
            }
        }


@dataclass
class EvidenceResult:
    """
    Evidence result from literature analysis.

    Used by ACMG rules to score variants.
    """
    rule: str  # e.g., "PS4", "PP1", "PS2", "PP4"
    applicable: bool
    strength: Optional[str] = None  # STRONG, MODERATE, SUPPORTING

    # Evidence details
    num_cases: Optional[int] = None
    num_controls: Optional[int] = None
    segregation_data: Optional[Dict[str, Any]] = None
    phenotype: Optional[str] = None

    # Source
    pmids: List[str] = field(default_factory=list)
    confidence: str = "medium"

    # Comment for classification report
    comment: Optional[str] = None

    def to_rule_result(self) -> Dict[str, Any]:
        """Convert to rule result format."""
        return {
            "rule": self.rule,
            "applicable": self.applicable,
            "strength": self.strength,
            "num_cases": self.num_cases,
            "pmids": self.pmids,
            "confidence": self.confidence,
            "comment": self.comment,
        }
