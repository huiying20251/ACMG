#!/usr/bin/env python3
"""
BP2 Rule - Variant in trans with a pathogenic variant (benign evidence)

BP2: For autosomal recessive inheritance - variant found in trans with a known pathogenic variant
in a gene definitively known to cause autosomal recessive disease.

This is a benign evidence criterion - observing a candidate benign variant in trans with a
known pathogenic variant in a recessive disease gene supports the benign classification.

Uses literature evidence from variant.variant_literature["evidence"]["bp2"].
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


class Bp2(abstract_rule):
    """
    BP2: Benign variant in trans with pathogenic variant for autosomal recessive disease.

    Evidence Strength:
    - Supporting: 1 proband with variant in trans
    - Moderate: 2+ probands with variant in trans
    - Strong: Variant observed in multiple unrelated families in trans configuration
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
        Assess BP2 from literature evidence.

        Args:
            variant_literature: dict with "evidence" key containing literature assessments
        """
        # Try literature evidence
        if variant_literature and "evidence" in variant_literature:
            bp2_literature = variant_literature["evidence"].get("bp2")
            if bp2_literature and hasattr(bp2_literature, 'applicable') and bp2_literature.applicable:
                return cls._assess_from_literature(bp2_literature)

        # No BP2 evidence available
        return RuleResult(
            "BP2",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            False,
            evidence_strength.SUPPORTING,
            "No trans configuration evidence available from literature.",
        )

    @classmethod
    def _assess_from_literature(cls, bp2_evidence) -> RuleResult:
        """Assess BP2 from literature EvidenceResult."""
        from literature_retrieval.bp2_evaluator import BP2Strength

        applicable = getattr(bp2_evidence, 'applicable', False)
        strength_obj = getattr(bp2_evidence, 'strength', BP2Strength.SUPPORTING)
        comment = getattr(bp2_evidence, 'comment', '') or ''
        pmids = getattr(bp2_evidence, 'pmids', []) or []

        # Convert BP2Strength to evidence_strength
        # BP2Strength: S=Strong, P=Supporting
        bp2_to_evidence = {
            BP2Strength.STRONG: evidence_strength.STRONG,
            BP2Strength.SUPPORTING: evidence_strength.SUPPORTING,
            BP2Strength.NOT_APPLICABLE: evidence_strength.SUPPORTING,
        }

        if isinstance(strength_obj, BP2Strength):
            strength = bp2_to_evidence.get(strength_obj, evidence_strength.SUPPORTING)
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

        # Get partner variant info if available
        partner_variant = getattr(bp2_evidence, 'partner_variant', None)
        if partner_variant:
            comment = f"{comment} [Partner: {partner_variant}]"

        return RuleResult(
            "BP2",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            applicable,
            strength,
            comment,
        )
