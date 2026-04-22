#!/usr/bin/env python3
"""
PM3 LLM Evaluator

Uses LLM to assess PM3 evidence based on:
1. PM3 ClinGen SVI specification rules
2. Literature case reports for the variant

PM3 is applicable for autosomal recessive inheritance when:
- A variant is found in trans (compound heterozygous) with a known pathogenic variant
- Or when a variant is homozygous/hemizygous consistent with disease inheritance

Evidence Strength (based on PM3 SVI proposal):
- Very Strong (VS): 4+ probands with trans/homozygous configuration
- Strong (S): 3 probands OR trans + known pathogenic + functional study
- Moderate (M): 2 probands OR trans + likely pathogenic
- Supporting (P): 1 proband with trans/homozygous/hemizygous

The base strength is determined by:
1. Type of configuration (trans vs homozygous/hemizygous)
2. Significance of the second variant (Pathogenic > Likely Pathogenic)

Then the strength can be UPGRADED based on the number of independent probands.
"""

import logging
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
from enum import Enum

from .literature_utils import Article, CaseReport, InheritancePattern
from llm_prompts import get_prompt

logger = logging.getLogger(__name__)


class PM3Strength(Enum):
    """PM3 evidence strength levels."""
    VERY_STRONG = "VS"
    STRONG = "S"
    MODERATE = "M"
    SUPPORTING = "P"
    NOT_APPLICABLE = "NA"

    @classmethod
    def from_value(cls, value: str) -> "PM3Strength":
        """Create from string value."""
        for strength in cls:
            if strength.value == value.upper():
                return strength
        return cls.NOT_APPLICABLE


# PM3 Evidence strength thresholds based on number of probands
PM3_PROBAND_THRESHOLDS = {
    # proband_count: base_strength_provided
    # If base < threshold, upgrade to threshold
    1: PM3Strength.SUPPORTING,      # 1 proband → Supporting
    2: PM3Strength.MODERATE,         # 2 probands → Moderate
    3: PM3Strength.STRONG,           # 3 probands → Strong
    4: PM3Strength.VERY_STRONG,     # 4+ probands → Very Strong
}


def _get_pm3_prompt() -> str:
    return get_prompt("pm3_assessment")


@dataclass
class PM3Assessment:
    """Result of PM3 assessment."""
    applicable: bool
    strength: PM3Strength = PM3Strength.NOT_APPLICABLE

    # Evidence details
    is_trans_config: bool = False  # Compound heterozygous
    is_homozygous: bool = False   # Homozygous
    is_hemizygous: bool = False  # Hemizygous (X-linked)

    # Known variant info
    known_pathogenic_variant: Optional[str] = None
    known_variant_significance: Optional[str] = None

    # Functional study support
    has_lof_functional_study: bool = False

    # Confirmed by testing
    confirmed_by_testing: bool = False  # Confirmed trans configuration

    # Case count info
    num_independent_probands: int = 0
    proband_details: List[Dict[str, Any]] = field(default_factory=list)

    # Source
    pmids: List[str] = field(default_factory=list)
    confidence: str = "medium"  # high, medium, low

    # LLM reasoning
    reasoning: Optional[str] = None
    comment: Optional[str] = None


# PM3 LLM Prompt Template - now loaded from config/llm_prompts.yaml via get_prompt("pm3_assessment")


class PM3LLMEvaluator:
    """
    LLM-based PM3 evidence assessor.

    Uses LLM to analyze literature case reports and determine
    if PM3 criteria are met and at what strength level.

    The strength is determined by:
    1. Base criteria (trans configuration, pathogenicity of second variant)
    2. Number of independent probands (can upgrade strength)
    """

    def __init__(self, llm_api_key: Optional[str] = None):
        self.llm_api_key = llm_api_key

    def assess(
        self,
        gene: str,
        variant_description: str,
        case_reports: List[CaseReport],
        inheritance_pattern: InheritancePattern = InheritancePattern.AUTOSOMAL_RECESSIVE,
        known_pathogenic_variants: Optional[List[Dict[str, Any]]] = None,
    ) -> PM3Assessment:
        """
        Assess PM3 evidence using LLM.

        Args:
            gene: Gene symbol
            variant_description: HGVS or VCF description
            case_reports: List of case reports from literature
            inheritance_pattern: Expected inheritance pattern
            known_pathogenic_variants: List of known pathogenic variants in this gene

        Returns:
            PM3Assessment with strength determination
        """
        # Check if inheritance pattern supports PM3
        if inheritance_pattern != InheritancePattern.AUTOSOMAL_RECESSIVE:
            return PM3Assessment(
                applicable=False,
                confidence="high",
                reasoning="PM3 only applies to autosomal recessive inheritance"
            )

        # Build case reports text
        case_text = self._format_case_reports(case_reports)

        # Build known variants text
        known_var_text = self._format_known_variants(known_pathogenic_variants)

        # Construct prompt
        prompt = _get_pm3_prompt().format(
            gene=gene,
            variant=variant_description,
            inheritance=inheritance_pattern.value,
            case_reports=case_text,
            known_variants=known_var_text,
        )

        # Call LLM (placeholder - would integrate with actual LLM API)
        assessment = self._call_llm(prompt)

        # Apply strength adjustment based on proband count
        if assessment.applicable:
            assessment = self._adjust_strength_by_proband_count(assessment)

        return assessment

    def _adjust_strength_by_proband_count(
        self,
        assessment: PM3Assessment,
    ) -> PM3Assessment:
        """
        Adjust PM3 strength based on number of independent probands.

        Per PM3 SVI proposal:
        - 1 proband: Supporting (P)
        - 2 probands: Moderate (M)
        - 3 probands: Strong (S)
        - 4+ probands: Very Strong (VS)

        The final strength is the MAX of base strength and proband-based strength.
        """
        num_probands = assessment.num_independent_probands

        if num_probands <= 0:
            return assessment

        # Determine strength based on proband count
        if num_probands >= 4:
            proband_strength = PM3Strength.VERY_STRONG
        elif num_probands == 3:
            proband_strength = PM3Strength.STRONG
        elif num_probands == 2:
            proband_strength = PM3Strength.MODERATE
        else:
            proband_strength = PM3Strength.SUPPORTING

        # Take the stronger of base strength and proband-based strength
        base_value = self._strength_to_value(assessment.strength)
        proband_value = self._strength_to_value(proband_strength)

        if proband_value > base_value:
            assessment.strength = proband_strength
            assessment.comment = (
                f"Strength upgraded from {assessment.strength.value} to {proband_strength.value} "
                f"based on {num_probands} independent probands"
            )
        else:
            assessment.comment = (
                f"Base strength {assessment.strength.value} maintained "
                f"(strength from {num_probands} proband(s) does not upgrade)"
            )

        return assessment

    def _strength_to_value(self, strength: PM3Strength) -> int:
        """Convert PM3Strength to numeric value for comparison."""
        strength_values = {
            PM3Strength.NOT_APPLICABLE: 0,
            PM3Strength.SUPPORTING: 1,
            PM3Strength.MODERATE: 2,
            PM3Strength.STRONG: 3,
            PM3Strength.VERY_STRONG: 4,
        }
        return strength_values.get(strength, 0)

    def _format_case_reports(self, case_reports: List[CaseReport]) -> str:
        """Format case reports for prompt."""
        if not case_reports:
            return "No case reports available."

        lines = []
        for i, cr in enumerate(case_reports, 1):
            lines.append(f"""
Case {i}:
- PMID: {cr.pmid}
- Variant: {cr.variant_description}
- Gene: {cr.gene}
- Number of cases: {cr.num_cases}
- Inheritance: {cr.inheritance or 'Not specified'}
- Inheritance pattern: {cr.inheritance_pattern.value}
- Phenotype: {cr.phenotype or 'Not specified'}
- de novo: {cr.is_de_novo} (confirmed: {cr.confirmed_de_novo})
- Segregation data: {cr.segregation_info or 'Not available'}
""")

        return "\n".join(lines)

    def _format_known_variants(
        self,
        known_variants: Optional[List[Dict[str, Any]]]
    ) -> str:
        """Format known pathogenic variants for prompt."""
        if not known_variants:
            return "No known pathogenic variants provided."

        lines = []
        for var in known_variants:
            lines.append(f"""
- Variant: {var.get('description', 'Unknown')}
- Significance: {var.get('significance', 'Unknown')}
- PMID: {var.get('pmid', 'Unknown')}
""")

        return "\n".join(lines)

    def _call_llm(self, prompt: str) -> PM3Assessment:
        """
        Call LLM API with prompt and parse response.

        This is a placeholder - integrate with actual LLM API
        (DeepSeek, OpenAI, Claude, etc.)
        """
        # Placeholder implementation
        # In production, this would call the actual LLM API

        logger.warning("LLM API not configured - using heuristic fallback")

        # Fallback to heuristic assessment
        return self._heuristic_assessment(prompt)

    def _heuristic_assessment(self, prompt: str) -> PM3Assessment:
        """
        Fallback heuristic assessment when LLM is not available.

        Analyzes the prompt text for keywords indicating PM3 applicability.
        """
        prompt_lower = prompt.lower()

        # Check for trans configuration keywords
        has_trans = any(kw in prompt_lower for kw in [
            "compound heterozyg", "trans", "opposite chromosomes",
            "biallelic", "two different variants"
        ])

        # Check for homozygous keywords
        has_homozygous = any(kw in prompt_lower for kw in [
            "homozygous", "homozygosity", "homozygote"
        ])

        # Check for hemizygous keywords
        has_hemizygous = any(kw in prompt_lower for kw in [
            "hemizygous", "hemizygote", "x-linked"
        ])

        # Check for known pathogenic variant
        has_known_pathogenic = any(kw in prompt_lower for kw in [
            "pathogenic", "likely pathogenic", "P/LP"
        ])

        # Check for functional study
        has_functional = any(kw in prompt_lower for kw in [
            "functional study", "lof", "loss of function",
            "nonsense", "frameshift", "splice"
        ])

        # Check for confirmation by testing
        has_testing_confirmation = any(kw in prompt_lower for kw in [
            "confirmed by testing", "parental testing", "molecular confirmation"
        ])

        # Count independent probands (heuristic: count "proband" mentions and "family" mentions)
        proband_count = 0
        # Count proband mentions
        import re
        proband_matches = re.findall(r'\bproband[s]?\b', prompt_lower)
        proband_count = len(proband_matches)

        # Also try to count independent families (look for family-related case groupings)
        family_matches = re.findall(r'\bfamily [a-z0-9]+\b', prompt_lower)
        proband_count = max(proband_count, len(family_matches))

        # If no explicit proband count, look for case numbers
        if proband_count == 0:
            # Look for "n cases" or "n patients"
            case_matches = re.findall(r'\b(\d+)\s+(?:cases?|patients?|families?)\b', prompt_lower)
            if case_matches:
                proband_count = max([int(m) for m in case_matches])

        # Determine applicability and strength
        applicable = has_trans or has_homozygous or has_hemizygous

        if not applicable:
            return PM3Assessment(
                applicable=False,
                confidence="medium",
                reasoning="No evidence of trans configuration or homozygous state found"
            )

        # Determine base strength
        base_strength = PM3Strength.SUPPORTING

        if has_trans and has_known_pathogenic:
            if has_functional and has_testing_confirmation:
                base_strength = PM3Strength.VERY_STRONG
            elif has_testing_confirmation:
                base_strength = PM3Strength.STRONG
            else:
                base_strength = PM3Strength.MODERATE
        elif has_homozygous or has_hemizygous:
            base_strength = PM3Strength.SUPPORTING

        # Determine proband-based strength
        if proband_count >= 4:
            proband_strength = PM3Strength.VERY_STRONG
        elif proband_count == 3:
            proband_strength = PM3Strength.STRONG
        elif proband_count == 2:
            proband_strength = PM3Strength.MODERATE
        elif proband_count == 1:
            proband_strength = PM3Strength.SUPPORTING
        else:
            proband_strength = base_strength  # No upgrade if count unknown

        # Take the stronger strength
        base_value = self._strength_to_value(base_strength)
        proband_value = self._strength_to_value(proband_strength)
        final_strength = base_strength if base_value >= proband_value else proband_strength

        # Build proband details
        proband_details = []
        if proband_count > 0:
            proband_details.append({
                "count": proband_count,
                "source": "heuristic_extraction"
            })

        return PM3Assessment(
            applicable=True,
            strength=final_strength,
            is_trans_config=has_trans,
            is_homozygous=has_homozygous,
            is_hemizygous=has_hemizygous,
            known_pathogenic_variant="detected" if has_known_pathogenic else None,
            has_lof_functional_study=has_functional,
            confirmed_by_testing=has_testing_confirmation,
            num_independent_probands=proband_count,
            proband_details=proband_details,
            confidence="low",
            reasoning=f"Heuristic: base={base_strength.value}, probands={proband_count}→{proband_strength.value}"
        )


def assess_pm3_from_literature(
    literature,
    gene: str,
    variant_description: str,
    inheritance_pattern: InheritancePattern = InheritancePattern.AUTOSOMAL_RECESSIVE,
    known_variants: Optional[List[Dict[str, Any]]] = None,
    llm_api_key: Optional[str] = None,
) -> PM3Assessment:
    """
    Convenience function to assess PM3 from literature data.

    Args:
        literature: VariantLiterature object
        gene: Gene symbol
        variant_description: Variant description
        inheritance_pattern: Inheritance pattern
        known_variants: Known pathogenic variants
        llm_api_key: LLM API key

    Returns:
        PM3Assessment
    """
    evaluator = PM3LLMEvaluator(llm_api_key=llm_api_key)

    return evaluator.assess(
        gene=gene,
        variant_description=variant_description,
        case_reports=literature.case_reports if hasattr(literature, 'case_reports') else [],
        inheritance_pattern=inheritance_pattern,
        known_pathogenic_variants=known_variants,
    )
