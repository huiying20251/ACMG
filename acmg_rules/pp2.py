#!/usr/bin/env python3

from typing import Callable, Optional
import pathlib

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    rule_type,
    evidence_type,
)
from information import Classification_Info, Info
from var_type import VARTYPE_GROUPS
from variant import TranscriptInfo, VariantInfo, PopulationDatabases_gnomAD


class Pp2(abstract_rule):
    """
    PP2: Missense variant in a gene that has a low rate of benign missense variation
    and where missense variants are a common mechanism of disease.

    Uses gnomAD missense z-score to determine if gene has low benign missense rate.
    z-score > threshold (default 2.09) indicates gene is intolerant to missense variants.
    """

    DEFAULT_ZSCORE_THRESHOLD = 2.09

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.VARIANT, class_info.TRANSCRIPT, class_info.VARIANT_GNOMAD_POPMAX),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        transcripts: list[TranscriptInfo],
        gnomad: Optional[PopulationDatabases_gnomAD],
        zscore_threshold: float = DEFAULT_ZSCORE_THRESHOLD,
    ) -> RuleResult:
        """
        Assess PP2 using gnomAD missense z-score.

        Args:
            variant: VariantInfo
            transcripts: List of transcript info
            gnomad: gnomAD population database with missense_zscore
            zscore_threshold: Minimum z-score threshold (default 2.09)
        """
        # In case one disease variant transcript is defined, use type of variant in that transcript
        # Otherwise use all variant types defined for variant
        if len(transcripts) == 1:
            variant_types = transcripts[0].var_type
        else:
            variant_types = variant.var_type

        # Check if variant is missense
        if not any(var_type in VARTYPE_GROUPS.MISSENSE.value for var_type in variant_types):
            return RuleResult(
                "PP2",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                f"PP2 does not apply to this variant, as PP2 does not apply to variant types {', '.join([var_type.value for var_type in variant.var_type])}.",
            )

        # Check gnomAD missense z-score
        if gnomad is None:
            return RuleResult(
                "PP2",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                f"PP2 cannot be assessed for {variant.gene_name} - no gnomAD data available.",
            )

        missense_zscore = getattr(gnomad, 'missense_zscore', None)

        if missense_zscore is None:
            return RuleResult(
                "PP2",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                f"PP2 cannot be assessed for {variant.gene_name} - missense z-score not available.",
            )

        if missense_zscore > zscore_threshold:
            return RuleResult(
                "PP2",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                True,
                evidence_strength.SUPPORTING,
                f"Missense variants in {variant.gene_name} are a known mechanism of disease (missense z-score: {missense_zscore:.2f} > {zscore_threshold}).",
            )
        else:
            return RuleResult(
                "PP2",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                f"PP2 does not apply - {variant.gene_name} missense z-score ({missense_zscore:.2f}) is below threshold ({zscore_threshold}).",
            )
