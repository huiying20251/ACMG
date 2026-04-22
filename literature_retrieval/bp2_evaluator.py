#!/usr/bin/env python3
"""
BP2 Evaluator - Benign variant in trans with a known pathogenic variant

BP2 Definition:
- For autosomal recessive diseases
- A variant is found in trans (compound heterozygous) with a KNOWN BENIGN or LIKELY BENIGN variant
- This suggests the variant in question is also benign (not causing disease)

Evidence Strength: Supporting (-1 point)

This differs from PM3 which deals with pathogenic variants in trans.
BP2 applies when the OTHER variant is known benign, suggesting the gene can tolerate
variation at this position without disease.
"""

import logging
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
from enum import Enum

from llm_prompts import get_prompt

logger = logging.getLogger(__name__)


class BP2Strength(Enum):
    """BP2 evidence strength levels."""
    SUPPORTING = "P"  # Standard BP2
    STRONG = "S"     # Multiple lines of evidence
    NOT_APPLICABLE = "NA"


@dataclass
class TransPartnerVariant:
    """Represents a variant in trans configuration."""
    variant_description: str
    clinical_significance: str  # Benign, Likely Benign
    evidence_source: str  # ClinVar, literature, etc.
    pmid: Optional[str] = None


def _get_bp2_prompt() -> str:
    return get_prompt("bp2_assessment")


@dataclass
class BP2Assessment:
    """Result of BP2 assessment."""
    applicable: bool

    # Configuration details
    is_trans_config: bool = False  # Compound heterozygous
    is_in_cis: bool = False       # Same chromosome (not relevant for BP2)
    is_homozygous: bool = False   # Homozygous (argues against BP2)

    # Partner variant info
    partner_variant: Optional[TransPartnerVariant] = None

    # Assessment details
    inheritance_pattern: str = "autosomal recessive"
    confirmed_by_testing: bool = False  # Phase confirmed by parental testing

    # Evidence
    strength: BP2Strength = BP2Strength.NOT_APPLICABLE

    # Case info
    num_cases: int = 0
    case_details: List[Dict[str, Any]] = field(default_factory=list)

    # Source
    pmids: List[str] = field(default_factory=list)
    confidence: str = "medium"

    # Reasoning
    reasoning: Optional[str] = None
    comment: Optional[str] = None


# BP2 LLM Prompt Template - now loaded from config/llm_prompts.yaml via get_prompt("bp2_assessment")


class BP2Evaluator:
    """
    BP2 evidence assessor.

    Evaluates whether a variant is found in trans with a known benign variant,
    supporting the classification of the variant as benign.
    """

    def __init__(self, llm_api_key: Optional[str] = None):
        self.llm_api_key = llm_api_key

    def assess(
        self,
        gene: str,
        variant_description: str,
        partner_variants: List[TransPartnerVariant],
        case_reports: Optional[List[Dict[str, Any]]] = None,
        inheritance_pattern: str = "autosomal recessive",
    ) -> BP2Assessment:
        """
        Assess BP2 evidence.

        Args:
            gene: Gene symbol
            variant_description: HGVS or VCF description
            partner_variants: List of known variants in trans configuration
            case_reports: List of case reports with phase information
            inheritance_pattern: Expected inheritance pattern

        Returns:
            BP2Assessment with determination
        """
        # Check inheritance pattern
        if inheritance_pattern.lower() not in ["autosomal recessive", "ar", "recessive"]:
            return BP2Assessment(
                applicable=False,
                confidence="high",
                reasoning="BP2 only applies to autosomal recessive inheritance"
            )

        # Filter for benign/likely benign partner variants
        benign_partners = [
            pv for pv in partner_variants
            if pv.clinical_significance.lower() in ["benign", "likely benign"]
        ]

        if not benign_partners:
            return BP2Assessment(
                applicable=False,
                confidence="high",
                reasoning="No benign or likely benign partner variants found in trans"
            )

        # Select the most significant benign partner
        best_partner = benign_partners[0]

        # Build case reports text
        case_text = self._format_case_reports(case_reports or [])

        # Build partner variants text
        partner_text = self._format_partner_variants(benign_partners)

        # Construct prompt
        prompt = _get_bp2_prompt().format(
            gene=gene,
            variant=variant_description,
            inheritance=inheritance_pattern,
            partner_variants=partner_text,
            case_reports=case_text,
        )

        # Call LLM or use heuristic fallback
        assessment = self._call_llm(prompt, best_partner, case_reports or [])

        return assessment

    def _format_partner_variants(self, partners: List[TransPartnerVariant]) -> str:
        """Format partner variants for prompt."""
        if not partners:
            return "No partner variants provided."

        lines = []
        for i, pv in enumerate(partners, 1):
            lines.append(f"""
Partner Variant {i}:
- Variant: {pv.variant_description}
- Clinical Significance: {pv.clinical_significance}
- Source: {pv.evidence_source}
- PMID: {pv.pmid or 'Not provided'}
""")
        return "\n".join(lines)

    def _format_case_reports(self, case_reports: List[Dict[str, Any]]) -> str:
        """Format case reports for prompt."""
        if not case_reports:
            return "No case reports available."

        lines = []
        for i, cr in enumerate(case_reports, 1):
            lines.append(f"""
Case {i}:
- PMID: {cr.get('pmid', 'Not provided')}
- Description: {cr.get('description', 'Not provided')}
- Phase Confirmed: {cr.get('phase_confirmed', False)}
- Phenotype: {cr.get('phenotype', 'Not provided')}
""")
        return "\n".join(lines)

    def _call_llm(
        self,
        prompt: str,
        partner: TransPartnerVariant,
        case_reports: List[Dict[str, Any]],
    ) -> BP2Assessment:
        """Call LLM API with prompt and parse response."""
        logger.warning("LLM API not configured - using heuristic fallback")
        return self._heuristic_assessment(prompt, partner, case_reports)

    def _heuristic_assessment(
        self,
        prompt: str,
        partner: TransPartnerVariant,
        case_reports: List[Dict[str, Any]],
    ) -> BP2Assessment:
        """Fallback heuristic assessment when LLM is not available."""
        import re

        prompt_lower = prompt.lower()

        # Check for phase confirmation
        phase_confirmed = any(kw in prompt_lower for kw in [
            "phase confirmed", "parental testing", "trans configuration confirmed",
            "compound heterozygous confirmed"
        ])

        # Check for molecular confirmation
        confirmed_by_testing = any(kw in prompt_lower for kw in [
            "molecular testing", "genetic testing", "sequencing confirmed"
        ]) or phase_confirmed

        # Determine applicability based on benign partner
        applicable = partner.clinical_significance.lower() in ["benign", "likely benign"]

        if not applicable:
            return BP2Assessment(
                applicable=False,
                confidence="medium",
                reasoning="Partner variant is not classified as benign or likely benign"
            )

        # Count cases
        num_cases = len(case_reports) if case_reports else 0

        # Look for additional cases in literature
        case_matches = re.findall(r'\b(\d+)\s+(?:cases?|patients?)\b', prompt_lower)
        if case_matches:
            num_cases = max(num_cases, max([int(m) for m in case_matches]))

        # BP2 is typically Supporting strength
        strength = BP2Strength.SUPPORTING

        # Check for enhanced evidence
        if phase_confirmed and num_cases >= 2:
            # Multiple cases with confirmed phase = Strong
            strength = BP2Strength.STRONG

        return BP2Assessment(
            applicable=applicable,
            is_trans_config=True,  # Assuming trans if benign partner exists
            partner_variant=partner,
            confirmed_by_testing=confirmed_by_testing,
            num_cases=num_cases,
            case_details=case_reports,
            strength=strength,
            confidence="low",
            reasoning=f"Heuristic: benign partner={partner.variant_description}, phase_confirmed={phase_confirmed}",
            comment=(
                f"BP2 applicable: variant in trans with known {partner.clinical_significance} variant "
                f"{partner.variant_description} ({partner.evidence_source})"
            )
        )


def assess_bp2_from_partners(
    gene: str,
    variant_description: str,
    partner_variants: List[TransPartnerVariant],
    case_reports: Optional[List[Dict[str, Any]]] = None,
    inheritance_pattern: str = "autosomal recessive",
    llm_api_key: Optional[str] = None,
) -> BP2Assessment:
    """
    Convenience function to assess BP2 from partner variant data.

    Args:
        gene: Gene symbol
        variant_description: Variant description
        partner_variants: List of TransPartnerVariant objects
        case_reports: List of case reports
        inheritance_pattern: Inheritance pattern
        llm_api_key: LLM API key

    Returns:
        BP2Assessment
    """
    evaluator = BP2Evaluator(llm_api_key=llm_api_key)
    return evaluator.assess(
        gene=gene,
        variant_description=variant_description,
        partner_variants=partner_variants,
        case_reports=case_reports,
        inheritance_pattern=inheritance_pattern,
    )