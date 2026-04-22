#!/usr/bin/env python3
"""
PS2 LLM Evaluator

Uses LLM to assess PS2 evidence based on:
1. PS2 ClinGen SVI specification rules (de novo criteria)
2. Literature case reports for the variant

PS2 is applicable for autosomal dominant inheritance when:
- A variant is confirmed to be de novo (not inherited from either parent)

Evidence Strength:
- Very Strong (VS): 4+ probands with confirmed de novo + high phenotype consistency
- Strong (S): 3 probands OR confirmed de novo + functional study
- Moderate (M): 2 probands with confirmed de novo
- Supporting (P): 1 proband with confirmed de novo

PM6 (Assumed de novo, not confirmed):
- 0.5 points if de novo is inferred but not confirmed by parental testing
"""

import logging
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
from enum import Enum

from .literature_utils import CaseReport, InheritancePattern
from llm_prompts import get_prompt

logger = logging.getLogger(__name__)


class PS2Strength(Enum):
    """PS2 evidence strength levels."""
    VERY_STRONG = "VS"
    STRONG = "S"
    MODERATE = "M"
    SUPPORTING = "P"
    NOT_APPLICABLE = "NA"


class DeNovoStatus(Enum):
    """De novo confirmation status."""
    CONFIRMED = "confirmed"       # Confirmed by parental testing
    ASSUMED = "assumed"          # Inferred but not confirmed
    NOT_DEC_NOVO = "not_de_novo"  # Inherited variant


def _get_ps2_prompt() -> str:
    return get_prompt("ps2_assessment")


@dataclass
class PS2Assessment:
    """Result of PS2 assessment."""
    applicable: bool
    rule: str = "NA"  # PS2, PM6, or NA
    strength: PS2Strength = PS2Strength.NOT_APPLICABLE

    # De novo status
    de_novo_status: DeNovoStatus = DeNovoStatus.NOT_DEC_NOVO
    has_parental_testing: bool = False

    # Case count info
    num_independent_probands: int = 0
    num_families: int = 0
    probands: List[Dict[str, Any]] = field(default_factory=list)

    # Phenotype info
    phenotype_consistency: str = "unknown"  # high, medium, low
    phenotype_details: str = ""
    phenotype_analysis: Optional[str] = None  # explanation of heterogeneity assessment

    # Evidence details
    is_x_linked: bool = False
    is_autosomal_dominant: bool = False

    # Source
    pmids: List[str] = field(default_factory=list)
    confidence: str = "medium"  # high, medium, low

    # LLM reasoning
    reasoning: Optional[str] = None
    comment: Optional[str] = None


class PS2LLMEvaluator:
    """
    LLM-based PS2 evidence assessor.

    Uses LLM to analyze literature case reports and determine
    if PS2 criteria are met and at what strength level.

    Strength determination:
    1. Base strength from de novo status (confirmed vs assumed)
    2. Upgrade based on number of independent probands
    3. Adjust for phenotypic consistency
    """

    # Strength thresholds based on proband count
    PROBAND_THRESHOLDS = {
        1: PS2Strength.SUPPORTING,
        2: PS2Strength.MODERATE,
        3: PS2Strength.STRONG,
        4: PS2Strength.VERY_STRONG,
    }

    def __init__(self, llm_api_key: Optional[str] = None):
        self.llm_api_key = llm_api_key

    def assess(
        self,
        gene: str,
        variant_description: str,
        case_reports: List[CaseReport],
        inheritance_pattern: InheritancePattern = InheritancePattern.AUTOSOMAL_DOMINANT,
    ) -> PS2Assessment:
        """
        Assess PS2 evidence using LLM.

        Args:
            gene: Gene symbol
            variant_description: HGVS or VCF description
            case_reports: List of case reports from literature
            inheritance_pattern: Expected inheritance pattern

        Returns:
            PS2Assessment with strength determination
        """
        # Check if inheritance pattern supports PS2
        if inheritance_pattern not in [
            InheritancePattern.AUTOSOMAL_DOMINANT,
            InheritancePattern.X_LINKED,
        ]:
            return PS2Assessment(
                applicable=False,
                confidence="high",
                reasoning="PS2 only applies to autosomal dominant or X-linked inheritance"
            )

        # Build case reports text
        case_text = self._format_case_reports(case_reports)

        # Construct prompt
        prompt = _get_ps2_prompt().format(
            gene=gene,
            variant=variant_description,
            inheritance=inheritance_pattern.value,
            case_reports=case_text,
        )

        # Call LLM (placeholder - would integrate with actual LLM API)
        assessment = self._call_llm(prompt)

        # Apply strength adjustment based on proband count
        if assessment.applicable:
            assessment = self._adjust_strength_by_proband_count(assessment)
            assessment = self._adjust_for_phenotype(assessment)

        return assessment

    def _adjust_strength_by_proband_count(
        self,
        assessment: PS2Assessment,
    ) -> PS2Assessment:
        """
        Adjust PS2 strength based on number of independent probands.

        Per PS2 SVI proposal:
        - 1 proband: Supporting
        - 2 probands: Moderate
        - 3 probands: Strong
        - 4+ probands: Very Strong
        """
        num_probands = assessment.num_independent_probands

        if num_probands <= 0:
            return assessment

        # Determine strength based on proband count
        if num_probands >= 4:
            proband_strength = PS2Strength.VERY_STRONG
        elif num_probands == 3:
            proband_strength = PS2Strength.STRONG
        elif num_probands == 2:
            proband_strength = PS2Strength.MODERATE
        else:
            proband_strength = PS2Strength.SUPPORTING

        # Take the stronger of base strength and proband-based strength
        base_value = self._strength_to_value(assessment.strength)
        proband_value = self._strength_to_value(proband_strength)

        if proband_value > base_value:
            assessment.strength = proband_strength
            assessment.comment = (
                f"Strength upgraded to {proband_strength.value} "
                f"based on {num_probands} independent probands"
            )
        else:
            assessment.comment = (
                f"Strength {assessment.strength.value} maintained "
                f"(proband count {num_probands} does not upgrade further)"
            )

        return assessment

    def _adjust_for_phenotype(
        self,
        assessment: PS2Assessment,
    ) -> PS2Assessment:
        """
        Adjust PS2 strength based on phenotypic consistency.

        High phenotypic consistency (≥3 probands with same phenotype):
        - Strength may be maintained/upgraded

        Low phenotypic consistency (different phenotypes):
        - May reduce strength or not apply PS2
        """
        if assessment.phenotype_consistency == "low":
            # High heterogeneity - reduce strength or don't apply
            if assessment.num_independent_probands < 2:
                assessment.applicable = False
                assessment.reasoning = (
                    "PS2 not applicable due to high phenotypic heterogeneity "
                    "and insufficient case numbers"
                )
                assessment.comment = (
                    f"PS2 not applied: phenotypes are from different disease spectra "
                    f"({assessment.phenotype_details}). Low heterogeneity requires ≥2 independent probands."
                )
            else:
                # Reduce strength by one level
                current_value = self._strength_to_value(assessment.strength)
                if current_value > 1:
                    new_strength = self._value_to_strength(current_value - 1)
                    assessment.strength = new_strength
                    assessment.comment = (
                        f"Strength reduced to {new_strength.value} due to "
                        f"high phenotypic heterogeneity ({assessment.phenotype_details})"
                    )
                else:
                    assessment.comment = (
                        f"Phenotypic heterogeneity noted ({assessment.phenotype_details}). "
                        "Strength cannot be reduced further."
                    )

        elif assessment.phenotype_consistency == "high":
            # High consistency - support upgrade
            assessment.comment = (
                f"High phenotypic consistency supports {assessment.strength.value} strength"
            ) if assessment.comment is None else assessment.comment

        return assessment

    def _strength_to_value(self, strength: PS2Strength) -> int:
        """Convert PS2Strength to numeric value for comparison."""
        strength_values = {
            PS2Strength.NOT_APPLICABLE: 0,
            PS2Strength.SUPPORTING: 1,
            PS2Strength.MODERATE: 2,
            PS2Strength.STRONG: 3,
            PS2Strength.VERY_STRONG: 4,
        }
        return strength_values.get(strength, 0)

    def _value_to_strength(self, value: int) -> PS2Strength:
        """Convert numeric value to PS2Strength."""
        value_to_strength = {
            0: PS2Strength.NOT_APPLICABLE,
            1: PS2Strength.SUPPORTING,
            2: PS2Strength.MODERATE,
            3: PS2Strength.STRONG,
            4: PS2Strength.VERY_STRONG,
        }
        return value_to_strength.get(value, PS2Strength.NOT_APPLICABLE)

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
- de novo: {cr.is_de_novo} (confirmed: {cr.confirmed_de_novo})
- Phenotype: {cr.phenotype or 'Not specified'}
- HPO terms: {', '.join(cr.hpo_terms[:5]) if cr.hpo_terms else 'None'}
- Segregation data: {cr.segregation_info or 'Not available'}
""")

        return "\n".join(lines)

    def _call_llm(self, prompt: str) -> PS2Assessment:
        """
        Call LLM API with prompt and parse response.

        This is a placeholder - integrate with actual LLM API.
        """
        logger.warning("LLM API not configured - using heuristic fallback")

        # Fallback to heuristic assessment
        return self._heuristic_assessment(prompt)

    def _heuristic_assessment(self, prompt: str) -> PS2Assessment:
        """
        Fallback heuristic assessment when LLM is not available.

        Analyzes the prompt text for keywords indicating PS2 applicability.
        """
        import re
        prompt_lower = prompt.lower()

        # Check for de novo keywords
        confirmed_de_novo = any(kw in prompt_lower for kw in [
            "confirmed de novo", "both parents negative",
            "de novo confirmed", "parental testing confirmed"
        ])

        assumed_de_novo = any(kw in prompt_lower for kw in [
            "assumed de novo", "presumed de novo", "likely de novo",
            "de novo (not tested)", "inherited from affected"
        ]) and not confirmed_de_novo

        # Check for parental testing
        has_parental_testing = any(kw in prompt_lower for kw in [
            "parental testing", "both parents tested", "maternal", "paternal"
        ])

        # Count de novo mentions
        de_novo_matches = re.findall(r'\bde novo\b', prompt_lower)
        num_de_novo = len(de_novo_matches)

        # Count probands (heuristic)
        proband_matches = re.findall(r'\bproband[s]?\b', prompt_lower)
        num_probands = max(len(proband_matches), num_de_novo)

        # If no explicit count, look for case numbers
        if num_probands == 0:
            case_matches = re.findall(r'\b(\d+)\s+(?:cases?|patients?|families?)\b', prompt_lower)
            if case_matches:
                num_probands = max([int(m) for m in case_matches])

        # Determine de novo status
        if confirmed_de_novo:
            de_novo_status = DeNovoStatus.CONFIRMED
            rule = "PS2"
        elif assumed_de_novo:
            de_novo_status = DeNovoStatus.ASSUMED
            rule = "PM6"
        else:
            de_novo_status = DeNovoStatus.NOT_DEC_NOVO
            rule = "NA"

        # Determine applicability
        applicable = confirmed_de_novo or assumed_de_novo

        if not applicable:
            return PS2Assessment(
                applicable=False,
                confidence="medium",
                reasoning="No evidence of de novo variant found"
            )

        # Determine base strength based on de novo status
        if de_novo_status == DeNovoStatus.CONFIRMED:
            base_strength = PS2Strength.MODERATE  # 1 point
        else:
            base_strength = PS2Strength.SUPPORTING  # 0.5 points

        # Determine phenotype consistency (heuristic)
        phenotype_consistency = "unknown"
        phenotype_details = ""

        # Look for phenotype keywords
        phenotype_mentions = re.findall(r'(?:phenotype|disease|manifestation):\s*([^\n]+)', prompt_lower)
        unique_phenotypes = set(phenotype_mentions)

        if len(unique_phenotypes) >= 3:
            phenotype_consistency = "low"
            phenotype_details = f"High heterogeneity: {len(unique_phenotypes)} different phenotypes"
            phenotype_analysis = (
                f"Found {len(unique_phenotypes)} distinct phenotypes, indicating low consistency. "
                "Cases show features from different disease systems or unrelated conditions."
            )
        elif len(unique_phenotypes) == 1 and len(unique_phenotypes) > 0:
            phenotype_consistency = "high"
            phenotype_details = f"High consistency: {list(unique_phenotypes)[0]}"
            phenotype_analysis = (
                "All cases share the same primary phenotype, indicating high consistency. "
                "Cases are consistent with the same disease spectrum."
            )
        elif len(unique_phenotypes) > 1:
            phenotype_consistency = "medium"
            phenotype_details = f"Moderate heterogeneity: {len(unique_phenotypes)} phenotypes"
            phenotype_analysis = (
                f"Found {len(unique_phenotypes)} different phenotypes. "
                "Some phenotypic variation but may still be part of the same disease spectrum with variable expressivity."
            )

        return PS2Assessment(
            applicable=applicable,
            rule=rule,
            strength=base_strength,
            de_novo_status=de_novo_status,
            has_parental_testing=has_parental_testing,
            num_independent_probands=num_probands,
            num_families=num_probands,
            phenotype_consistency=phenotype_consistency,
            phenotype_details=phenotype_details,
            phenotype_analysis=phenotype_analysis,
            confidence="low",
            reasoning=f"Heuristic: de_novo_status={de_novo_status.value}, probands={num_probands}",
        )


def assess_ps2_from_literature(
    literature,
    gene: str,
    variant_description: str,
    inheritance_pattern: InheritancePattern = InheritancePattern.AUTOSOMAL_DOMINANT,
    llm_api_key: Optional[str] = None,
) -> PS2Assessment:
    """
    Convenience function to assess PS2 from literature data.

    Args:
        literature: VariantLiterature object
        gene: Gene symbol
        variant_description: Variant description
        inheritance_pattern: Inheritance pattern
        llm_api_key: LLM API key

    Returns:
        PS2Assessment
    """
    evaluator = PS2LLMEvaluator(llm_api_key=llm_api_key)

    return evaluator.assess(
        gene=gene,
        variant_description=variant_description,
        case_reports=literature.case_reports if hasattr(literature, 'case_reports') else [],
        inheritance_pattern=inheritance_pattern,
    )
