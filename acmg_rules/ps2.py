#!/usr/bin/env python3
"""
PS2 Rule - de novo variants

PS2: Confirmed de novo origin of the variant (both parents negative).

Uses literature evidence from variant.variant_literature["evidence"]["ps2"].
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


class Ps2(abstract_rule):
    """
    PS2: Confirmed de novo variant (both parents negative).

    Evidence Strength (based on number of probands):
    - Very Strong (VS): 4+ independent probands with confirmed de novo
    - Strong (S): 3 independent probands OR confirmed de novo + functional study
    - Moderate (M): 2 independent probands with confirmed de novo
    - Supporting (P): 1 proband with confirmed de novo

    PM6: Assumed de novo (not confirmed by parental testing)
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
        Assess PS2 from literature evidence.

        Args:
            variant_literature: dict with "evidence" key containing literature assessments
        """
        # Try literature evidence
        if variant_literature and "evidence" in variant_literature:
            ps2_literature = variant_literature["evidence"].get("ps2")
            if ps2_literature and hasattr(ps2_literature, 'applicable') and ps2_literature.applicable:
                return cls._assess_from_literature(ps2_literature)

        # No PS2 evidence available
        return RuleResult(
            "PS2",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            False,
            evidence_strength.SUPPORTING,
            "No de novo evidence available from literature.",
        )

    @classmethod
    def _assess_from_literature(cls, ps2_evidence) -> RuleResult:
        """Assess PS2 from literature EvidenceResult."""
        from literature_retrieval.ps2_evaluator import PS2Strength

        applicable = getattr(ps2_evidence, 'applicable', False)
        strength_obj = getattr(ps2_evidence, 'strength', PS2Strength.SUPPORTING)
        comment = getattr(ps2_evidence, 'comment', '') or ''
        pmids = getattr(ps2_evidence, 'pmids', []) or []

        # Convert PS2Strength to evidence_strength
        # PS2Strength: VS=Very Strong, S=Strong, M=Moderate, P=Supporting
        ps2_to_evidence = {
            PS2Strength.VERY_STRONG: evidence_strength.VERY_STRONG,
            PS2Strength.STRONG: evidence_strength.STRONG,
            PS2Strength.MODERATE: evidence_strength.MODERATE,
            PS2Strength.SUPPORTING: evidence_strength.SUPPORTING,
            PS2Strength.NOT_APPLICABLE: evidence_strength.SUPPORTING,
        }

        if isinstance(strength_obj, PS2Strength):
            strength = ps2_to_evidence.get(strength_obj, evidence_strength.SUPPORTING)
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

        # Check for de novo status
        de_novo_status = getattr(ps2_evidence, 'de_novo_status', None)
        if de_novo_status:
            # PS2Assessment has de_novo_status attribute
            num_probands = getattr(ps2_evidence, 'num_independent_probands', 1)
            if num_probands > 1:
                comment = f"{comment} ({num_probands} independent probands)"

        return RuleResult(
            "PS2",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            applicable,
            strength,
            comment,
        )
