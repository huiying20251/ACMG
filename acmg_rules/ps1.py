#!/usr/bin/env python3

from typing import Callable

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from information import Classification_Info, Info
from clinvar_utils import ClinVar_Type, ClinVar


class Ps1_protein(abstract_rule):
    """
    PS1: Position is classified as pathogenic
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.VARIANT_CLINVAR,),
        )

    @classmethod
    def assess_rule(cls, clinvar_result: dict[ClinVar_Type, ClinVar]) -> RuleResult:
        clinvar_same_aa = clinvar_result[ClinVar_Type.SAME_AA_CHANGE]
        if clinvar_same_aa.pathogenic:
            comment = f"The following ClinVar entries show the same amino acid change as pathogenic: {clinvar_same_aa.ids}."
            result = True
            if clinvar_same_aa.associated_ids:
                comment = (
                    comment
                    + f" The following ClinVar entries show the same amino acid change as likely pathogenic: {clinvar_same_aa.associated_ids}."
                )
        else:
            comment = "No ClinVar entries found that show the same amino acid change as pathogenic."
            result = False
        return RuleResult(
            "PS1",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.STRONG,
            comment,
        )


class Ps1_splicing(abstract_rule):
    """
    PS1 for splicing: Splice variant in same position has been show to be pathogenic

    Evidence strength depends on:
    - SAME_NUCLEOTIDE (same position): PS1 or PS1_Moderate
    - SAME_SPLICE_SITE (same motif, different position): PS1_Moderate or PS1_Supporting
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.VARIANT_CLINVAR,),
        )

    @classmethod
    def assess_rule(cls, clinvar_result: dict[ClinVar_Type, ClinVar]) -> RuleResult:
        clinvar_same_nucleotide = clinvar_result[ClinVar_Type.SAME_NUCLEOTIDE]
        clinvar_same_splice_site = clinvar_result[ClinVar_Type.SAME_SPLICE_SITE]

        # Check same nucleotide (same position, different alt)
        has_same_nuc = clinvar_same_nucleotide.pathogenic
        has_same_nuc_lp = bool(clinvar_same_nucleotide.associated_ids)

        # Check same splice site (different position, same motif)
        has_same_motif = clinvar_same_splice_site.pathogenic
        has_same_motif_lp = bool(clinvar_same_splice_site.associated_ids)

        # Determine result and strength
        result = has_same_nuc or has_same_nuc_lp or has_same_motif or has_same_motif_lp

        if not result:
            return RuleResult(
                "PS1",
                rule_type.SPLICING,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.STRONG,
                "No ClinVar entries found that show splice variants at the same nucleotide position or same splice motif as pathogenic.",
            )

        # Build comment
        comment_parts = []
        if has_same_nuc:
            comment_parts.append(f"same nucleotide position P: {clinvar_same_nucleotide.ids}")
        if has_same_nuc_lp:
            comment_parts.append(f"same nucleotide position LP: {clinvar_same_nucleotide.associated_ids}")
        if has_same_motif:
            comment_parts.append(f"same splice motif P: {clinvar_same_splice_site.ids}")
        if has_same_motif_lp:
            comment_parts.append(f"same splice motif LP: {clinvar_same_splice_site.associated_ids}")

        comment = "ClinVar entries show: " + "; ".join(comment_parts)

        # Determine evidence strength
        # 相同位点: P → STRONG, LP → MODERATE
        # 同基序不同位点: P/LP → SUPPORTING
        if has_same_nuc:
            strength = evidence_strength.STRONG if has_same_nuc else evidence_strength.MODERATE
        elif has_same_nuc_lp:
            strength = evidence_strength.MODERATE
        elif has_same_motif:
            strength = evidence_strength.MODERATE
        else:
            strength = evidence_strength.SUPPORTING

        return RuleResult(
            "PS1",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            True,
            strength,
            comment,
        )
