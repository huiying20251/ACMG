#!/usr/bin/env python3
"""
PP1/BS4 LLM Evaluator

Uses LLM to assess co-segregation (PP1) or lack of co-segregation (BS4) evidence based on:
1. Family pedigree data
2. Number of affected/unaffected carriers vs non-carriers
3. Calculated LOD score

PP1: Co-segregation with disease in multiple affected family members
- Segregation count determines evidence strength (NOT LOD score)
- Higher segregation count = stronger evidence

BS4: Lack of segregation in affected members of a family
- Used when variant found in healthy carriers

Evidence Strength (PP1 - based on segregation count):

**For Dominant or X-linked inheritance:**
| Segregations | Strength | Points |
|--------------|----------|--------|
| 2-3 | Supporting (PP1) | 1分 |
| 4 | Moderate (PP1_Moderate) | 2分 |
| 5+ | Strong (PP1_Strong) | 4分 |

**For Recessive inheritance:**
| Evidence | Strength |
|----------|----------|
| 1 family with 2 patients both carrying variant = 1 segregation | Supporting (PP1) |
| 1 family with carrier parents, 3 affected carriers + 1 unaffected non-carrier | Moderate (PP1_Moderate) |
| 1 family with carrier parents, 2 affected carriers + 2 unaffected non-carriers | Moderate (PP1_Moderate) |
| 1 family with 2 affected carriers + 3 unaffected non-carriers | Strong (PP1_Strong) |

Note: "Segregation" = one meiosis in which the variant segregates with disease.
For dominant: affected carrier relative = 1 segregation
For recessive: affected carrier = 1 segregation, confirmed by unaffected non-carrier = additional evidence
"""

import logging
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
from enum import Enum
import math

from llm_prompts import get_prompt

logger = logging.getLogger(__name__)


class SegregationStrength(Enum):
    """Segregation evidence strength levels."""
    VERY_STRONG = "VS"
    STRONG = "S"
    MODERATE = "M"
    SUPPORTING = "P"
    NOT_APPLICABLE = "NA"


class SegregationStatus(Enum):
    """Co-segregation status."""
    SEGREGATION = "segregation"         # Variant segregates with disease
    NO_SEGREGATION = "no_segregation"   # Variant does not segregate
    INCONCLUSIVE = "inconclusive"        # Cannot determine


# LOD Score Thresholds for PP1 (legacy, based on ClinGen SVI)
# DEPRECATED: Now using segregation count instead of LOD score
PP1_LOD_THRESHOLDS = {
    2.0: SegregationStrength.VERY_STRONG,
    1.5: SegregationStrength.STRONG,
    1.0: SegregationStrength.MODERATE,
    0.5: SegregationStrength.SUPPORTING,
}

# PP1 Segregation Count Thresholds (new logic)
# For dominant or X-linked inheritance
PP1_SEGREGATION_THRESHOLDS_DOMINANT = {
    5: SegregationStrength.STRONG,       # 5+ segregations
    4: SegregationStrength.MODERATE,     # 4 segregations
    2: SegregationStrength.SUPPORTING,   # 2-3 segregations
}

# PP1 Segregation Count Thresholds (new logic)
# For recessive inheritance
# Format: (affected_carriers, unaffected_non_carriers) minimum required
PP1_SEGREGATION_THRESHOLDS_RECESSIVE = {
    (2, 3): SegregationStrength.STRONG,       # 2 affected + 3 unaffected non-carriers
    (3, 1): SegregationStrength.MODERATE,     # 3 affected + 1 unaffected non-carrier
    (2, 2): SegregationStrength.MODERATE,     # 2 affected + 2 unaffected non-carriers
    (2, 0): SegregationStrength.SUPPORTING,   # 2 affected carriers (minimum)
}


# Negative LOD thresholds for BS4
BS4_LOD_THRESHOLDS = {
    -2.0: SegregationStrength.VERY_STRONG,  # Strong lack of segregation
    -1.5: SegregationStrength.STRONG,
    -1.0: SegregationStrength.MODERATE,
    -0.5: SegregationStrength.SUPPORTING,
}


@dataclass
class FamilyMember:
    """Represents a single family member in the pedigree."""
    individual_id: str
    affected: bool
    genotype: str  # "carrier", "non-carrier", "unknown"
    relationship: str  # "proband", "parent", "sibling", "offspring", "grandparent", "other"
    sex: str = "unknown"  # "male", "female", "unknown"


@dataclass
class PedigreeData:
    """Family pedigree data for segregation analysis."""
    family_id: str
    proband_id: str
    members: List[FamilyMember]

    # Pre-calculated LOD score (if provided)
    lod_score: Optional[float] = None

    # Parsed from literature
    num_affected_carriers: int = 0
    num_affected_non_carriers: int = 0
    num_unaffected_carriers: int = 0
    num_unaffected_non_carriers: int = 0

    @classmethod
    def from_family_members(cls, family_id: str, proband_id: str, members: List[FamilyMember]) -> "PedigreeData":
        """Create PedigreeData from list of FamilyMember objects."""
        pd = cls(family_id=family_id, proband_id=proband_id, members=members)
        pd._calculate_counts()
        return pd

    def _calculate_counts(self):
        """Calculate carrier/non-carrier counts."""
        for member in self.members:
            if member.genotype == "carrier":
                if member.affected:
                    self.num_affected_carriers += 1
                else:
                    self.num_unaffected_carriers += 1
            elif member.genotype == "non-carrier":
                if member.affected:
                    self.num_affected_non_carriers += 1
                else:
                    self.num_unaffected_non_carriers += 1


@dataclass
class SegregationAssessment:
    """Result of segregation assessment (PP1 or BS4)."""
    applicable: bool
    rule: str  # PP1 or BS4

    # Segregation status
    segregation_status: SegregationStatus = SegregationStatus.INCONCLUSIVE
    lod_score: float = 0.0

    # Evidence strength
    strength: SegregationStrength = SegregationStrength.NOT_APPLICABLE

    # Case details
    num_affected_carriers: int = 0
    num_affected_non_carriers: int = 0
    num_unaffected_carriers: int = 0
    num_unaffected_non_carriers: int = 0

    # Inheritance pattern (for PP1 strength determination)
    inheritance_pattern: str = "unknown"  # "autosomal dominant", "x-linked", "autosomal recessive"

    # Number of segregations (for PP1 strength calculation)
    # For dominant: count of affected carrier relatives (excluding proband)
    # For recessive: count of affected carriers with unaffected non-carrier confirmation
    num_segregations: int = 0

    # Family info
    families: List[str] = field(default_factory=list)
    family_details: List[Dict[str, Any]] = field(default_factory=list)

    # Source
    pmids: List[str] = field(default_factory=list)
    confidence: str = "medium"  # high, medium, low

    # LLM reasoning
    reasoning: Optional[str] = None
    comment: Optional[str] = None


def _get_segregation_prompt() -> str:
    return get_prompt("segregation_assessment")


class SegregationEvaluator:
    """
    LLM-based segregation evidence assessor.

    Evaluates PP1 (co-segregation) or BS4 (lack of co-segregation)
    based on family pedigree data.
    """

    def __init__(self, llm_api_key: Optional[str] = None):
        self.llm_api_key = llm_api_key

    def assess_pp1(
        self,
        gene: str,
        variant_description: str,
        pedigree_data: List[PedigreeData],
        inheritance_pattern: str = "autosomal dominant",
        published_lod: Optional[float] = None,
    ) -> SegregationAssessment:
        """
        Assess PP1 evidence (co-segregation).

        Args:
            gene: Gene symbol
            variant_description: HGVS or VCF description
            pedigree_data: List of family pedigrees
            inheritance_pattern: Expected inheritance pattern
            published_lod: Pre-calculated LOD score from literature

        Returns:
            SegregationAssessment with PP1 determination
        """
        return self._assess(
            gene=gene,
            variant_description=variant_description,
            pedigree_data=pedigree_data,
            inheritance_pattern=inheritance_pattern,
            published_lod=published_lod,
            rule="PP1",
        )

    def assess_bs4(
        self,
        gene: str,
        variant_description: str,
        pedigree_data: List[PedigreeData],
        inheritance_pattern: str = "autosomal dominant",
        published_lod: Optional[float] = None,
    ) -> SegregationAssessment:
        """
        Assess BS4 evidence (lack of co-segregation).

        Args:
            gene: Gene symbol
            variant_description: HGVS or VCF description
            pedigree_data: List of family pedigrees
            inheritance_pattern: Expected inheritance pattern
            published_lod: Pre-calculated LOD score from literature

        Returns:
            SegregationAssessment with BS4 determination
        """
        return self._assess(
            gene=gene,
            variant_description=variant_description,
            pedigree_data=pedigree_data,
            inheritance_pattern=inheritance_pattern,
            published_lod=published_lod,
            rule="BS4",
        )

    def _assess(
        self,
        gene: str,
        variant_description: str,
        pedigree_data: List[PedigreeData],
        inheritance_pattern: str,
        published_lod: Optional[float],
        rule: str,
    ) -> SegregationAssessment:
        """Internal assessment method."""
        # Build pedigree text
        pedigree_text = self._format_pedigree_data(pedigree_data)

        # Determine rule name for prompt
        if rule == "PP1":
            rule_description = "PP1 (Co-segregation with disease)"
        else:
            rule_description = "BS4 (Lack of co-segregation)"

        # Construct prompt
        prompt = _get_segregation_prompt().format(
            gene=gene,
            variant=variant_description,
            inheritance=inheritance_pattern,
            pedigree_data=pedigree_text,
            published_lod=f"{published_lod:.2f}" if published_lod is not None else "Not provided",
        )

        # Call LLM or use heuristic fallback
        assessment = self._call_llm(prompt, rule)

        # Store inheritance pattern for strength calculation
        assessment.inheritance_pattern = inheritance_pattern

        # Calculate segregation count
        assessment.num_segregations = self._calculate_segregation_count(
            pedigree_data, inheritance_pattern
        )

        # Post-process: use published LOD if available and more reliable
        if published_lod is not None and assessment.confidence == "low":
            assessment.lod_score = published_lod

        # For PP1: use segregation count for strength determination
        if rule == "PP1":
            assessment = self._calculate_pp1_strength_from_segregation(assessment)
        else:
            # For BS4: use LOD-based strength (legacy)
            assessment = self._calculate_strength_from_lod(assessment, rule)

        return assessment

    def _format_pedigree_data(self, pedigree_data: List[PedigreeData]) -> str:
        """Format pedigree data for prompt."""
        if not pedigree_data:
            return "No pedigree data available."

        lines = []
        for i, family in enumerate(pedigree_data, 1):
            lines.append(f"\nFamily {i} (ID: {family.family_id}):")
            lines.append(f"  Proband: {family.proband_id}")
            lines.append(f"  Members:")

            for member in family.members:
                status = "AFFECTED" if member.affected else "UNAFFECTED"
                genotype = member.genotype.upper()
                rel = member.relationship
                lines.append(f"    - {member.individual_id} ({rel}): {status}, {genotype}")

            # Summary counts
            lines.append(f"  Summary:")
            lines.append(f"    Affected carriers: {family.num_affected_carriers}")
            lines.append(f"    Affected non-carriers: {family.num_affected_non_carriers}")
            lines.append(f"    Unaffected carriers: {family.num_unaffected_carriers}")
            lines.append(f"    Unaffected non-carriers: {family.num_unaffected_non_carriers}")

            if family.lod_score is not None:
                lines.append(f"  LOD Score: {family.lod_score:.2f}")

        return "\n".join(lines)

    def _call_llm(self, prompt: str, rule: str) -> SegregationAssessment:
        """Call LLM API with prompt and parse response."""
        logger.warning("LLM API not configured - using heuristic fallback")
        return self._heuristic_assessment(prompt, rule)

    def _heuristic_assessment(self, prompt: str, rule: str) -> SegregationAssessment:
        """Fallback heuristic assessment when LLM is not available."""
        import re

        prompt_lower = prompt.lower()

        # Extract inheritance pattern
        inheritance_pattern = "unknown"
        if "autosomal dominant" in prompt_lower or "ad" in prompt_lower:
            inheritance_pattern = "autosomal dominant"
        elif "x-linked" in prompt_lower or "xld" in prompt_lower:
            inheritance_pattern = "x-linked"
        elif "autosomal recessive" in prompt_lower or "ar" in prompt_lower:
            inheritance_pattern = "autosomal recessive"

        # Extract LOD score if mentioned
        lod_match = re.search(r'lod\s*[:=]?\s*(-?\d+\.?\d*)', prompt_lower)
        lod_score = float(lod_match.group(1)) if lod_match else None

        # Extract counts
        affected_carriers = len(re.findall(r'affected\s+carrier', prompt_lower))
        unaffected_carriers = len(re.findall(r'unaffected\s+carrier', prompt_lower))
        affected_non_carriers = len(re.findall(r'affected\s+non-carrier', prompt_lower))
        unaffected_non_carriers = len(re.findall(r'unaffected\s+non-carrier', prompt_lower))

        # Calculate LOD if we have counts
        if lod_score is None:
            lod_score = self._calculate_lod_from_counts(
                affected_carriers,
                unaffected_carriers,
                affected_non_carriers,
                rule
            )

        # Determine segregation status
        if lod_score is not None:
            if lod_score > 0.5:
                segregation_status = SegregationStatus.SEGREGATION
            elif lod_score < -0.5:
                segregation_status = SegregationStatus.NO_SEGREGATION
            else:
                segregation_status = SegregationStatus.INCONCLUSIVE
        else:
            segregation_status = SegregationStatus.INCONCLUSIVE

        # Determine if applicable
        if rule == "PP1":
            applicable = lod_score is not None and lod_score >= 0.5
        else:  # BS4
            applicable = lod_score is not None and lod_score <= -0.5

        # Calculate strength
        assessment = SegregationAssessment(
            applicable=applicable,
            rule=rule,
            segregation_status=segregation_status,
            lod_score=lod_score or 0.0,
            num_affected_carriers=affected_carriers,
            num_unaffected_carriers=unaffected_carriers,
            num_affected_non_carriers=affected_non_carriers,
            num_unaffected_non_carriers=unaffected_non_carriers,
            inheritance_pattern=inheritance_pattern,
            confidence="low",
            reasoning=f"Heuristic: LOD={lod_score}, segregation={segregation_status.value}",
        )

        # For PP1: use segregation count for strength determination
        if rule == "PP1":
            assessment = self._calculate_pp1_strength_from_segregation(assessment)
        else:
            # For BS4: use LOD-based strength (legacy)
            assessment = self._calculate_strength_from_lod(assessment, rule)

        return assessment

    def _calculate_lod_from_counts(
        self,
        affected_carriers: int,
        unaffected_carriers: int,
        affected_non_carriers: int,
        rule: str,
    ) -> float:
        """
        Calculate approximate LOD score from carrier counts.

        Simplified calculation for autosomal dominant:
        LOD = (affected_carriers × log10(2)) - log10(2) × affected_non_carriers
        """
        if affected_carriers == 0:
            return -1.0 if rule == "BS4" else 0.0

        # Base LOD from affected carriers (assuming θ = 0)
        # For each affected carrier: +log10(2) ≈ +0.301
        lod = affected_carriers * math.log10(2)

        # Penalty for affected non-carriers (unexpected)
        # For each affected non-carrier: -log10(2) ≈ -0.301
        lod -= affected_non_carriers * math.log10(2)

        # Penalty for unaffected carriers (possible phenocopy)
        # For each unaffected carrier: -log10(2) ≈ -0.301
        lod -= unaffected_carriers * math.log10(2)

        return lod

    def _calculate_strength_from_lod(
        self,
        assessment: SegregationAssessment,
        rule: str,
    ) -> SegregationAssessment:
        """Determine strength from LOD score (legacy method)."""
        lod = assessment.lod_score

        if rule == "PP1":
            if lod >= 2.0:
                assessment.strength = SegregationStrength.VERY_STRONG
            elif lod >= 1.5:
                assessment.strength = SegregationStrength.STRONG
            elif lod >= 1.0:
                assessment.strength = SegregationStrength.MODERATE
            elif lod >= 0.5:
                assessment.strength = SegregationStrength.SUPPORTING
            else:
                assessment.strength = SegregationStrength.NOT_APPLICABLE
        else:  # BS4
            if lod <= -2.0:
                assessment.strength = SegregationStrength.VERY_STRONG
            elif lod <= -1.5:
                assessment.strength = SegregationStrength.STRONG
            elif lod <= -1.0:
                assessment.strength = SegregationStrength.MODERATE
            elif lod <= -0.5:
                assessment.strength = SegregationStrength.SUPPORTING
            else:
                assessment.strength = SegregationStrength.NOT_APPLICABLE

        return assessment

    def _calculate_pp1_strength_from_segregation(
        self,
        assessment: SegregationAssessment,
    ) -> SegregationAssessment:
        """
        Determine PP1 strength based on segregation count.

        For Dominant or X-linked inheritance:
        | Segregations | Strength |
        |--------------|----------|
        | 5+ | Strong |
        | 4 | Moderate |
        | 2-3 | Supporting |
        | <2 | Not applicable |

        For Recessive inheritance:
        | Affected Carriers | Unaffected Non-carriers | Strength |
        |-------------------|-------------------------|----------|
        | 2+ | 3+ | Strong |
        | 3+ | 1+ | Moderate |
        | 2+ | 0+ | Supporting |
        """
        inheritance = assessment.inheritance_pattern.lower() if assessment.inheritance_pattern else "unknown"

        if inheritance in ("autosomal_dominant", "ad", "x_linked", "xld", "dominant"):
            # Dominant/X-linked: count affected carrier relatives (excluding proband)
            num_segregations = assessment.num_segregations

            if num_segregations >= 5:
                assessment.strength = SegregationStrength.STRONG
            elif num_segregations == 4:
                assessment.strength = SegregationStrength.MODERATE
            elif num_segregations >= 2:
                assessment.strength = SegregationStrength.SUPPORTING
            else:
                assessment.strength = SegregationStrength.NOT_APPLICABLE
                assessment.applicable = False

        elif inheritance in ("autosomal_recessive", "ar", "recessive"):
            # Recessive: affected carriers + unaffected non-carriers confirmation
            num_affected = assessment.num_affected_carriers
            num_unaffected_non_carriers = assessment.num_unaffected_non_carriers

            if num_affected >= 2 and num_unaffected_non_carriers >= 3:
                assessment.strength = SegregationStrength.STRONG
            elif num_affected >= 3 and num_unaffected_non_carriers >= 1:
                assessment.strength = SegregationStrength.MODERATE
            elif num_affected >= 2 and num_unaffected_non_carriers >= 0:
                assessment.strength = SegregationStrength.SUPPORTING
            else:
                assessment.strength = SegregationStrength.NOT_APPLICABLE
                assessment.applicable = False
        else:
            # Unknown inheritance: use LOD-based fallback
            assessment.strength = SegregationStrength.NOT_APPLICABLE
            assessment.applicable = False

        return assessment

    def _calculate_segregation_count(
        self,
        pedigree_data: List[PedigreeData],
        inheritance_pattern: str,
    ) -> int:
        """
        Calculate number of segregations from pedigree data.

        For dominant/X-linked: count affected carrier relatives (excluding proband)
        For recessive: count affected carriers (each affected carrier = 1 segregation)

        Returns:
            Number of segregations
        """
        inheritance = inheritance_pattern.lower() if inheritance_pattern else "unknown"

        if inheritance in ("autosomal_dominant", "ad", "x_linked", "xld", "dominant"):
            # Count affected carriers excluding the proband
            total = 0
            for family in pedigree_data:
                for member in family.members:
                    if member.genotype == "carrier" and member.affected:
                        # Exclude proband from count
                        if member.relationship.lower() != "proband":
                            total += 1
            return total

        elif inheritance in ("autosomal_recessive", "ar", "recessive"):
            # For recessive: each affected carrier counts as 1 segregation
            total = 0
            for family in pedigree_data:
                total += family.num_affected_carriers
            return total

        return 0


def assess_pp1_from_pedigree(
    gene: str,
    variant_description: str,
    pedigree_data: List[PedigreeData],
    inheritance_pattern: str = "autosomal dominant",
    published_lod: Optional[float] = None,
    llm_api_key: Optional[str] = None,
) -> SegregationAssessment:
    """
    Convenience function to assess PP1 from pedigree data.

    Args:
        gene: Gene symbol
        variant_description: Variant description
        pedigree_data: List of PedigreeData objects
        inheritance_pattern: Inheritance pattern
        published_lod: Pre-calculated LOD score
        llm_api_key: LLM API key

    Returns:
        SegregationAssessment
    """
    evaluator = SegregationEvaluator(llm_api_key=llm_api_key)
    return evaluator.assess_pp1(
        gene=gene,
        variant_description=variant_description,
        pedigree_data=pedigree_data,
        inheritance_pattern=inheritance_pattern,
        published_lod=published_lod,
    )


def assess_bs4_from_pedigree(
    gene: str,
    variant_description: str,
    pedigree_data: List[PedigreeData],
    inheritance_pattern: str = "autosomal dominant",
    published_lod: Optional[float] = None,
    llm_api_key: Optional[str] = None,
) -> SegregationAssessment:
    """
    Convenience function to assess BS4 from pedigree data.

    Args:
        gene: Gene symbol
        variant_description: Variant description
        pedigree_data: List of PedigreeData objects
        inheritance_pattern: Inheritance pattern
        published_lod: Pre-calculated LOD score
        llm_api_key: LLM API key

    Returns:
        SegregationAssessment
    """
    evaluator = SegregationEvaluator(llm_api_key=llm_api_key)
    return evaluator.assess_bs4(
        gene=gene,
        variant_description=variant_description,
        pedigree_data=pedigree_data,
        inheritance_pattern=inheritance_pattern,
        published_lod=published_lod,
    )