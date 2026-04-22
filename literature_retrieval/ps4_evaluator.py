#!/usr/bin/env python3
"""
PS4 Rule - Prevalence in cases

PS4: Variant enriched in cases compared to controls, or multiple unrelated cases
with the same variant and similar phenotype.

Uses literature evidence from variant.variant_literature["evidence"]["ps4"].

Evidence Strength (based on number of probands and disease rarity):

1. Standard dominant genetic diseases (e.g., familial hypercholesterolemia, Marfan, hereditary deafness):
   - Very Strong (VS): 15+ independent probands
   - Strong (S): 10+ independent probands
   - Moderate (M): 6+ independent probands
   - Supporting (P): 2+ independent probands

2. Very rare developmental delay disorders (dominant):
   - Strong (S): 4+ independent probands
   - Moderate (M): 2+ independent probands
   - Supporting (P): 1+ independent proband

3. Intermediate rare diseases (rarer than FH/Marfan/deafness, more common than Dravet):
   - Strong (S): 7+ independent probands
   - Moderate (M): 4+ independent probands
   - Supporting (P): 2+ independent probands

Note: All proband counts refer to UNRELATED individuals, not from the same family.
"""

import logging
from enum import Enum
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional

from .literature_utils import CaseReport, InheritancePattern
from llm_prompts import get_prompt

logger = logging.getLogger(__name__)


class PS4Strength(Enum):
    """PS4 evidence strength levels."""
    VERY_STRONG = "VS"
    STRONG = "S"
    MODERATE = "M"
    SUPPORTING = "P"
    NOT_APPLICABLE = "NA"


class DiseaseRarityCategory(Enum):
    """Disease rarity categories affecting PS4 thresholds."""
    STANDARD = "standard"  # Standard dominant diseases
    VERY_RARE_DD = "very_rare_developmental_delay"  # Very rare developmental delay disorders
    INTERMEDIATE_RARE = "intermediate_rare"  # Rarer than FH/Marfan/deafness, more common than Dravet


# PS4 thresholds based on disease rarity category and evidence strength
# Format: {category: {strength: min_cases_needed}}
PS4_THRESHOLDS = {
    DiseaseRarityCategory.STANDARD: {
        PS4Strength.VERY_STRONG: 10,
        PS4Strength.STRONG: 6,
        PS4Strength.MODERATE: 2,
        PS4Strength.SUPPORTING: 1,
    },
    DiseaseRarityCategory.VERY_RARE_DD: {
        PS4Strength.VERY_STRONG: 4,
        PS4Strength.STRONG: 4,
        PS4Strength.MODERATE: 2,
        PS4Strength.SUPPORTING: 1,
    },
    DiseaseRarityCategory.INTERMEDIATE_RARE: {
        PS4Strength.VERY_STRONG: 7,
        PS4Strength.STRONG: 4,
        PS4Strength.MODERATE: 2,
        PS4Strength.SUPPORTING: 2,
    },
}


def _get_ps4_prompt() -> str:
    return get_prompt("ps4_assessment")


@dataclass
class PS4Assessment:
    """Result of PS4 assessment."""
    applicable: bool
    rule: str = "NA"  # PS4 or NA
    strength: PS4Strength = PS4Strength.NOT_APPLICABLE

    # Case count info
    num_independent_probands: int = 0
    probands: List[Dict[str, Any]] = field(default_factory=list)

    # Disease category
    disease_category: DiseaseRarityCategory = DiseaseRarityCategory.STANDARD

    # Phenotype info
    phenotype_consistency: str = "unknown"  # high, medium, low
    phenotype_details: str = ""

    # Inheritance (must be dominant for PS4)
    inheritance: InheritancePattern = InheritancePattern.UNKNOWN
    is_de_novo: bool = False  # If True, should use PS2 instead

    # Evidence details
    case_control_enrichment: Optional[float] = None  # e.g., odds ratio
    p_value: Optional[float] = None

    # Source
    pmids: List[str] = field(default_factory=list)
    confidence: str = "medium"  # high, medium, low

    # LLM reasoning
    reasoning: Optional[str] = None
    comment: Optional[str] = None


class PS4LLMEvaluator:
    """
    LLM-based PS4 evidence assessor.

    Uses LLM to analyze literature case reports and determine
    if PS4 criteria are met and at what strength level.

    PS4 applies when:
    - NOT de novo (if de novo, use PS2 instead)
    - IS dominant inheritance
    - Multiple unrelated cases with the same variant and similar phenotype

    Strength determination:
    1. Count independent probands with the variant
    2. Determine disease rarity category
    3. Apply threshold based on category

    Note: All proband counts refer to UNRELATED individuals, not from the same family.
    """

    # Standard dominant genetic diseases (e.g., familial hypercholesterolemia, Marfan, hereditary deafness)
    # Thresholds: ≥2 → Supporting, ≥6 → Moderate, ≥10 → Strong, ≥15 → Very Strong
    STANDARD_PROBAND_THRESHOLDS = {
        2: PS4Strength.SUPPORTING,
        6: PS4Strength.MODERATE,
        10: PS4Strength.STRONG,
        15: PS4Strength.VERY_STRONG,
    }

    # Very rare developmental delay disorders (dominant) - more lenient thresholds
    # Thresholds: ≥1 → Supporting, ≥2 → Moderate, ≥4 → Strong
    DD_PROBAND_THRESHOLDS = {
        1: PS4Strength.SUPPORTING,
        2: PS4Strength.MODERATE,
        4: PS4Strength.STRONG,
    }

    # Intermediate rare diseases (rarer than FH/Marfan/deafness, more common than Dravet)
    # More lenient than standard but stricter than very rare DD
    # Thresholds: ≥2 → Supporting, ≥4 → Moderate, ≥7 → Strong
    INTERMEDIATE_PROBAND_THRESHOLDS = {
        2: PS4Strength.SUPPORTING,
        4: PS4Strength.MODERATE,
        7: PS4Strength.STRONG,
    }

    def __init__(self, llm_api_key: Optional[str] = None):
        self.llm_api_key = llm_api_key

    def _strength_to_value(self, strength: PS4Strength) -> int:
        """Convert PS4Strength to numeric value for comparison."""
        strength_values = {
            PS4Strength.NOT_APPLICABLE: 0,
            PS4Strength.SUPPORTING: 1,
            PS4Strength.MODERATE: 2,
            PS4Strength.STRONG: 3,
            PS4Strength.VERY_STRONG: 4,
        }
        return strength_values.get(strength, 0)

    def _value_to_strength(self, value: int) -> PS4Strength:
        """Convert numeric value to PS4Strength."""
        value_to_strength = {
            0: PS4Strength.NOT_APPLICABLE,
            1: PS4Strength.SUPPORTING,
            2: PS4Strength.MODERATE,
            3: PS4Strength.STRONG,
            4: PS4Strength.VERY_STRONG,
        }
        return value_to_strength.get(value, PS4Strength.NOT_APPLICABLE)

    def evaluate_case_report(self, case: CaseReport) -> PS4Assessment:
        """
        Evaluate a single case report for PS4 evidence.

        Returns PS4Assessment with strength based on:
        - Number of independent probands
        - Disease rarity category
        - Phenotype consistency
        """
        assessment = PS4Assessment(
            applicable=False,
            reasoning="PS4 requires dominant inheritance and multiple cases without de novo status",
        )

        # PS4 does not apply to de novo variants (use PS2 instead)
        if case.is_de_novo or case.de_novo_status == "confirmed":
            assessment.comment = "Variant is de novo; PS2 should be used instead of PS4"
            return assessment

        # PS4 requires dominant inheritance
        if case.inheritance != InheritancePattern.AUTOSOMAL_DOMINANT:
            assessment.comment = f"PS4 requires autosomal dominant inheritance, got {case.inheritance}"
            return assessment

        # Check if case has phenotype data
        if not case.case_description:
            assessment.comment = "No case description available for PS4 assessment"
            return assessment

        assessment.applicable = True
        assessment.num_independent_probands = max(1, case.num_cases or 1)
        assessment.inheritance = case.inheritance
        assessment.is_de_novo = case.is_de_novo
        assessment.phenotype_consistency = getattr(case, 'phenotype_consistency', 'medium')
        assessment.pmids = [case.pmid] if case.pmid else []

        # Determine disease category (would typically come from OMIM or similar)
        disease_category = getattr(case, 'disease_category', DiseaseRarityCategory.STANDARD)
        assessment.disease_category = disease_category

        # Determine strength based on proband count and disease category
        strength = self._determine_strength(
            assessment.num_independent_probands,
            disease_category
        )
        assessment.strength = strength

        return assessment

    def _determine_strength(
        self,
        num_probands: int,
        disease_category: DiseaseRarityCategory = DiseaseRarityCategory.STANDARD
    ) -> PS4Strength:
        """
        Determine PS4 strength based on number of probands and disease category.

        Disease rarity categories:
        - STANDARD: Standard dominant diseases (thresholds: 2, 6, 10, 15)
        - VERY_RARE_DD: Very rare developmental delay (thresholds: 1, 2, 4)
        - INTERMEDIATE_RARE: Rarer than FH/Marfan but more common than Dravet (thresholds: 2, 4, 7)

        Note: All proband counts refer to UNRELATED individuals.
        """
        if disease_category == DiseaseRarityCategory.VERY_RARE_DD:
            thresholds = self.DD_PROBAND_THRESHOLDS
        elif disease_category == DiseaseRarityCategory.INTERMEDIATE_RARE:
            thresholds = self.INTERMEDIATE_PROBAND_THRESHOLDS
        else:
            thresholds = self.STANDARD_PROBAND_THRESHOLDS

        # Find the highest strength level met
        applicable_strengths = []
        for count, strength in thresholds.items():
            if num_probands >= count:
                applicable_strengths.append((count, strength))

        if not applicable_strengths:
            return PS4Strength.NOT_APPLICABLE

        # Return the highest strength (highest count threshold met)
        applicable_strengths.sort(key=lambda x: x[0], reverse=True)
        return applicable_strengths[0][1]

    def aggregate_assessments(
        self,
        assessments: List[PS4Assessment]
    ) -> PS4Assessment:
        """
        Aggregate multiple PS4 assessments into a single result.

        Takes the highest strength from all assessments and
        sums the proband counts.
        """
        if not assessments:
            return PS4Assessment(
                applicable=False,
                comment="No PS4 assessments available"
            )

        # Filter to applicable assessments
        applicable = [a for a in assessments if a.applicable]
        if not applicable:
            return PS4Assessment(
                applicable=False,
                comment="No applicable PS4 evidence found"
            )

        # Sum proband counts
        total_probands = sum(a.num_independent_probands for a in applicable)

        # Take the strongest strength
        strength_values = [self._strength_to_value(a.strength) for a in applicable]
        max_strength_value = max(strength_values)
        strongest_strength = self._value_to_strength(max_strength_value)

        # Determine disease category (prefer more rare if mixed)
        categories = [a.disease_category for a in applicable]
        if DiseaseRarityCategory.VERY_RARE_DD in categories:
            disease_category = DiseaseRarityCategory.VERY_RARE_DD
        elif DiseaseRarityCategory.INTERMEDIATE_RARE in categories:
            disease_category = DiseaseRarityCategory.INTERMEDIATE_RARE
        else:
            disease_category = DiseaseRarityCategory.STANDARD

        # Collect all PMIDs
        all_pmids = []
        for a in applicable:
            all_pmids.extend(a.pmids)

        # Phenotype consistency: require high consistency for higher strengths
        phenotype_consistencies = [a.phenotype_consistency for a in applicable]
        high_consistency_count = phenotype_consistencies.count("high")
        if high_consistency_count >= len(applicable) * 0.5:
            phenotype_consistency = "high"
        elif high_consistency_count > 0:
            phenotype_consistency = "medium"
        else:
            phenotype_consistency = "low"

        return PS4Assessment(
            applicable=True,
            rule="PS4",
            strength=strongest_strength,
            num_independent_probands=total_probands,
            disease_category=disease_category,
            phenotype_consistency=phenotype_consistency,
            inheritance=InheritancePattern.AUTOSOMAL_DOMINANT,
            pmids=all_pmids,
            confidence="high" if len(applicable) >= 3 else "medium",
            comment=f"PS4 based on {total_probands} independent probands with {phenotype_consistency} phenotype consistency"
        )

    def assess_from_literature(
        self,
        literature: Dict[str, Any]
    ) -> PS4Assessment:
        """
        Assess PS4 from literature evidence dictionary.

        Expected structure:
        {
            "evidence": {
                "ps4": {
                    "num_cases": int,
                    "phenotype_consistency": str,
                    "disease_category": str,
                    "pmids": list,
                    ...
                }
            }
        }
        """
        ps4_evidence = literature.get("evidence", {}).get("ps4")

        if not ps4_evidence:
            return PS4Assessment(
                applicable=False,
                comment="No PS4 evidence in literature"
            )

        # Check if de novo (should use PS2 instead)
        is_de_novo = ps4_evidence.get("is_de_novo", False)
        if is_de_novo:
            return PS4Assessment(
                applicable=False,
                comment="Variant is de novo; PS2 should be used instead of PS4"
            )

        # Check inheritance
        inheritance = ps4_evidence.get("inheritance", "")
        if inheritance.upper() not in ("AD", "AUTOSOMAL_DOMINANT", "DOMINANT"):
            return PS4Assessment(
                applicable=False,
                comment=f"PS4 requires autosomal dominant inheritance, got {inheritance}"
            )

        num_cases = ps4_evidence.get("num_cases", 0)
        if num_cases < 1:
            return PS4Assessment(
                applicable=False,
                comment="No cases reported for PS4"
            )

        # Determine disease category
        category_str = ps4_evidence.get("disease_category", "standard")
        disease_category_map = {
            "very_rare_dd": DiseaseRarityCategory.VERY_RARE_DD,
            "intermediate_rare": DiseaseRarityCategory.INTERMEDIATE_RARE,
            "standard": DiseaseRarityCategory.STANDARD,
        }
        disease_category = disease_category_map.get(
            category_str.lower(),
            DiseaseRarityCategory.STANDARD
        )

        # Determine strength
        strength = self._determine_strength(num_cases, disease_category)

        return PS4Assessment(
            applicable=True,
            rule="PS4",
            strength=strength,
            num_independent_probands=num_cases,
            disease_category=disease_category,
            phenotype_consistency=ps4_evidence.get("phenotype_consistency", "medium"),
            inheritance=InheritancePattern.AUTOSOMAL_DOMINANT,
            pmids=ps4_evidence.get("pmids", []),
            confidence=ps4_evidence.get("confidence", "medium"),
            reasoning=ps4_evidence.get("reasoning"),
            comment=ps4_evidence.get("comment")
        )


def assess_ps4_from_cases(
    cases: List[CaseReport],
    disease_category: DiseaseRarityCategory = DiseaseRarityCategory.STANDARD
) -> PS4Assessment:
    """
    Convenience function to assess PS4 from a list of case reports.

    Excludes de novo cases (for PS2) and requires dominant inheritance.
    """
    evaluator = PS4LLMEvaluator()

    # Filter to non-de novo dominant cases
    eligible_cases = [
        c for c in cases
        if not c.is_de_novo
        and c.inheritance == InheritancePattern.AUTOSOMAL_DOMINANT
    ]

    if not eligible_cases:
        return PS4Assessment(
            applicable=False,
            comment="No eligible cases for PS4 (requires non-de novo dominant cases)"
        )

    assessments = []
    for case in eligible_cases:
        assessment = evaluator.evaluate_case_report(case)
        if assessment.applicable:
            assessments.append(assessment)

    if not assessments:
        return PS4Assessment(
            applicable=False,
            comment="No applicable PS4 evidence from cases"
        )

    # Update disease category if specified
    for a in assessments:
        a.disease_category = disease_category

    return evaluator.aggregate_assessments(assessments)