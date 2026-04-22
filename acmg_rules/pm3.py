#!/usr/bin/env python3
"""
PM3 Rule - Biallelic variants (compound heterozygous/homozygous/hemizygous)

PM3: For autosomal recessive inheritance - variant found in trans with a known pathogenic variant,
or homozygous/hemizygous consistent with disease inheritance.

Uses literature evidence from variant.variant_literature["evidence"]["pm3"].
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


class Pm3(abstract_rule):
    """
    PM3: Biallelic inheritance evidence (compound heterozygous/homozygous/hemizygous).

    Evidence Strength (based on number of probands and configuration):
    - Very Strong (VS): 4+ independent probands with trans configuration
    - Strong (S): 3 independent probands OR trans + known pathogenic + functional study
    - Moderate (M): 2 independent probands OR trans + likely pathogenic
    - Supporting (P): 1 proband with trans/homozygous/hemizygous configuration
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
        Assess PM3 from literature evidence.

        Args:
            variant_literature: dict with "evidence" key containing literature assessments
        """
        # Try literature evidence
        if variant_literature and "evidence" in variant_literature:
            pm3_literature = variant_literature["evidence"].get("pm3")
            if pm3_literature and hasattr(pm3_literature, 'applicable') and pm3_literature.applicable:
                return cls._assess_from_literature(pm3_literature)

        # No PM3 evidence available
        return RuleResult(
            "PM3",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            False,
            evidence_strength.SUPPORTING,
            "No biallelic inheritance evidence available from literature.",
        )

    @classmethod
    def _assess_from_literature(cls, pm3_evidence) -> RuleResult:
        """Assess PM3 from literature EvidenceResult."""
        from literature_retrieval.pm3_evaluator import PM3Strength

        applicable = getattr(pm3_evidence, 'applicable', False)
        strength_obj = getattr(pm3_evidence, 'strength', PM3Strength.SUPPORTING)
        comment = getattr(pm3_evidence, 'comment', '') or ''
        pmids = getattr(pm3_evidence, 'pmids', []) or []

        # Convert PM3Strength to evidence_strength
        # PM3Strength: VS=Very Strong, S=Strong, M=Moderate, P=Supporting
        pm3_to_evidence = {
            PM3Strength.VERY_STRONG: evidence_strength.VERY_STRONG,
            PM3Strength.STRONG: evidence_strength.STRONG,
            PM3Strength.MODERATE: evidence_strength.MODERATE,
            PM3Strength.SUPPORTING: evidence_strength.SUPPORTING,
            PM3Strength.NOT_APPLICABLE: evidence_strength.SUPPORTING,
        }

        if isinstance(strength_obj, PM3Strength):
            strength = pm3_to_evidence.get(strength_obj, evidence_strength.SUPPORTING)
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

        # Check for configuration type
        is_trans = getattr(pm3_evidence, 'is_trans_config', False)
        is_homozygous = getattr(pm3_evidence, 'is_homozygous', False)
        is_hemizygous = getattr(pm3_evidence, 'is_hemizygous', False)

        config_parts = []
        if is_trans:
            config_parts.append("compound heterozygous")
        if is_homozygous:
            config_parts.append("homozygous")
        if is_hemizygous:
            config_parts.append("hemizygous")

        if config_parts:
            comment = f"{comment} [Configuration: {', '.join(config_parts)}]"

        return RuleResult(
            "PM3",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            applicable,
            strength,
            comment,
        )
