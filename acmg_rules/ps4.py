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

from typing import Callable, Optional

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    rule_type,
    evidence_type,
)
from information import Classification_Info, Info


class Ps4(abstract_rule):
    """
    PS4: Prevalence in cases evidence (case-control enrichment).

    Evidence Strength (based on number of probands and disease rarity):

    Standard dominant genetic diseases (e.g., familial hypercholesterolemia, Marfan, hereditary deafness):
    - Very Strong (VS): 15+ independent probands
    - Strong (S): 10+ independent probands
    - Moderate (M): 6+ independent probands
    - Supporting (P): 2+ independent probands

    Note: PS4 does not apply to de novo variants (use PS2 instead).
    Note: All proband counts refer to UNRELATED individuals.
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.VARIANT_LITERATURE,),
        )

    @classmethod
    def assess_rule(
        cls,
        variant_literature: Optional[dict] = None,
    ) -> RuleResult:
        """
        Assess PS4 from literature evidence.

        Args:
            variant_literature: dict with "evidence" key containing literature assessments
        """
        # Try literature evidence
        if variant_literature and "evidence" in variant_literature:
            ps4_literature = variant_literature["evidence"].get("ps4")
            if ps4_literature and hasattr(ps4_literature, 'applicable') and ps4_literature.applicable:
                return cls._assess_from_literature(ps4_literature)

        # No PS4 evidence available
        return RuleResult(
            "PS4",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            False,
            evidence_strength.SUPPORTING,
            "No case-control prevalence evidence available from literature.",
        )

    @classmethod
    def _assess_from_literature(cls, ps4_evidence) -> RuleResult:
        """Assess PS4 from literature EvidenceResult."""
        from literature_retrieval.ps4_evaluator import PS4Strength

        applicable = getattr(ps4_evidence, 'applicable', False)
        strength_obj = getattr(ps4_evidence, 'strength', PS4Strength.SUPPORTING)
        comment = getattr(ps4_evidence, 'comment', '') or ''
        pmids = getattr(ps4_evidence, 'pmids', []) or []
        num_probands = getattr(ps4_evidence, 'num_independent_probands', 0)

        # Convert PS4Strength to evidence_strength
        # PS4Strength: VS=Very Strong, S=Strong, M=Moderate, P=Supporting
        ps4_to_evidence = {
            PS4Strength.VERY_STRONG: evidence_strength.VERY_STRONG,
            PS4Strength.STRONG: evidence_strength.STRONG,
            PS4Strength.MODERATE: evidence_strength.MODERATE,
            PS4Strength.SUPPORTING: evidence_strength.SUPPORTING,
            PS4Strength.NOT_APPLICABLE: evidence_strength.SUPPORTING,
        }

        if isinstance(strength_obj, PS4Strength):
            strength = ps4_to_evidence.get(strength_obj, evidence_strength.SUPPORTING)
        else:
            # Fallback for string or other types
            try:
                strength = evidence_strength(strength_obj.lower())
            except (ValueError, AttributeError):
                strength = evidence_strength.SUPPORTING

        # Add PMIDs to comment
        if pmids:
            pmid_str = ", ".join([f"PMID:{p}" for p in pmids])
            if comment:
                comment = f"{comment} ({pmid_str})"
            else:
                comment = f"Evidence from {pmid_str}"

        # Add proband count to comment
        if num_probands > 0:
            disease_category = getattr(ps4_evidence, 'disease_category', 'standard')
            if comment:
                comment = f"{comment} [{num_probands} probands, category: {disease_category.value}]"
            else:
                comment = f"{num_probands} independent probands with consistent phenotype"

        return RuleResult(
            "PS4",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            applicable,
            strength,
            comment,
        )