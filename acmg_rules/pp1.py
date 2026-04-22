#!/usr/bin/env python3

from typing import Callable, Optional, Any

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    rule_type,
    evidence_type,
)
from information import Classification_Info, Info
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds
from variant import MultifactorialLikelihood


class Pp1(abstract_rule):
    """
    PP1: Co-segregation with disease in multiple affected family members in a gene definitively known to cause the disease.
    Here used for Co-segregation

    Uses literature evidence from variant.variant_literature["evidence"]["pp1"] when available.
    Falls back to ClinVar multifactorial evidence if literature evidence is not available.
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_MULTIFACTORIAL_LIKELIHOOD,
                class_info.THRESHOLD_LIKELIHOOD_PATHOGENIC,
                class_info.VARIANT_LITERATURE,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        multifactorial_likelihood: MultifactorialLikelihood,
        thresholds: Threshold,
        variant_literature: Optional[dict] = None,
    ) -> RuleResult:
        # Try literature evidence first
        if variant_literature and "evidence" in variant_literature:
            pp1_literature = variant_literature["evidence"].get("pp1")
            if pp1_literature and hasattr(pp1_literature, 'applicable') and pp1_literature.applicable:
                # Use literature evidence
                return cls._assess_from_literature(pp1_literature)

        # Fall back to ClinVar multifactorial evidence
        return cls._assess_from_clinvar(multifactorial_likelihood, thresholds)

    @classmethod
    def _assess_from_literature(cls, pp1_evidence) -> RuleResult:
        """Assess PP1 from literature evidence."""
        from literature_retrieval.segregation_evaluator import SegregationStrength

        strength_obj = getattr(pp1_evidence, 'strength', None)
        if strength_obj is None:
            strength = evidence_strength.SUPPORTING
        elif isinstance(strength_obj, SegregationStrength):
            # Convert SegregationStrength to evidence_strength
            segregation_to_evidence = {
                SegregationStrength.VERY_STRONG: evidence_strength.VERY_STRONG,
                SegregationStrength.STRONG: evidence_strength.STRONG,
                SegregationStrength.MODERATE: evidence_strength.MODERATE,
                SegregationStrength.SUPPORTING: evidence_strength.SUPPORTING,
                SegregationStrength.NOT_APPLICABLE: evidence_strength.SUPPORTING,
            }
            strength = segregation_to_evidence.get(strength_obj, evidence_strength.SUPPORTING)
        else:
            # Fallback for string or other types
            try:
                strength = evidence_strength(strength_obj.lower())
            except (ValueError, AttributeError):
                strength = evidence_strength.SUPPORTING

        comment = getattr(pp1_evidence, 'comment', '') or ''
        pmids = getattr(pp1_evidence, 'pmids', []) or []

        if pmids and comment:
            pmid_str = ", ".join([f"PMID:{p}" for p in pmids])
            comment = f"{comment} ({pmid_str})"

        applicable = getattr(pp1_evidence, 'applicable', False)

        return RuleResult(
            "PP1",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            applicable,
            strength,
            comment,
        )

    @classmethod
    def _assess_from_clinvar(
        cls,
        multifactorial_likelihood: MultifactorialLikelihood,
        thresholds: Threshold,
    ) -> RuleResult:
        """Assess PP1 from ClinVar multifactorial evidence."""
        num_thresholds_met = assess_thresholds(
            thresholds, multifactorial_likelihood.co_segregation
        )
        if num_thresholds_met is None:
            result = False
            comment = "No co-segregation data available for variant."
            strength = evidence_strength.SUPPORTING

        elif num_thresholds_met == 0:
            result = False
            strength = evidence_strength.SUPPORTING
            comment = f"Co-segregation of {multifactorial_likelihood.co_segregation} given for the variant meets no threshold for pathogenic evidence."
        else:
            result = True
            strength = thresholds.strengths[num_thresholds_met - 1]
            comment = f"Co-segregation of {multifactorial_likelihood.co_segregation} given for variant meets threshold for {strength.value} pathogenic evidence ({thresholds.thresholds[num_thresholds_met -1]})."
        return RuleResult(
            "PP1",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )
