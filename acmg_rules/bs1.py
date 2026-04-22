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
from information import Info, Classification_Info
from variant import PopulationDatabases_gnomAD, VariantInfo


def _get_default_clingen_db_path() -> pathlib.Path:
    """Get the default ClinGen rules database path."""
    return pathlib.Path(__file__).parent.parent.parent / "clingen_rules.db"


class Bs1(abstract_rule):
    """
    BS1: Frequency of variant higher in population than expected based on disease frequency.

    Uses ClinGen gene-specific frequency cutoffs when available, otherwise falls back
    to config thresholds.
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT,
                class_info.VARIANT_GNOMAD_POPMAX,
                class_info.THRESHOLD_BS1,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        gnomad: PopulationDatabases_gnomAD,
        threshold_bs1: float,
    ) -> RuleResult:
        # Try ClinGen gene-specific cutoff first
        gene_name = variant.gene_name
        clingen_db_path = _get_default_clingen_db_path()

        if clingen_db_path and clingen_db_path.exists():
            result = cls._assess_with_clingen(gene_name, gnomad, clingen_db_path)
            if result is not None:
                return result

        # Fall back to config threshold
        if gnomad.subpopulation_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif gnomad.subpopulation_frequency > threshold_bs1:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
        if gnomad.subpopulation == "None":
            comment = f"Variant does not occur in gnomAD, allele frequency in gnomAD is assumed to be 0."
        if gnomad.subpopulation == "ALL":
            comment = f"Variant has no subpopulation. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        return RuleResult(
            "BS1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.STRONG,
            comment,
        )

    @classmethod
    def _assess_with_clingen(
        cls,
        gene_name: str,
        gnomad: PopulationDatabases_gnomAD,
        clingen_db_path: pathlib.Path,
    ) -> Optional[RuleResult]:
        """Assess BS1 using ClinGen gene-specific frequency cutoffs."""
        if gnomad.subpopulation_frequency is None:
            return RuleResult(
                "BS1",
                rule_type.GENERAL,
                evidence_type.BENIGN,
                False,
                evidence_strength.STRONG,
                "gnomAD frequency is None.",
            )

        try:
            from check_clingen_rules import ClinGenRulesDB

            with ClinGenRulesDB(clingen_db_path) as db:
                applies, applied_strength, reason = db.check_frequency_applies(
                    gene_name, "BS1", gnomad.subpopulation_frequency
                )

                if applies:
                    if applied_strength == "Stand Alone":
                        strength = evidence_strength.STAND_ALONE
                    elif applied_strength == "Strong":
                        strength = evidence_strength.STRONG
                    elif applied_strength == "Supporting":
                        strength = evidence_strength.SUPPORTING
                    else:
                        strength = evidence_strength.STRONG

                    return RuleResult(
                        "BS1",
                        rule_type.GENERAL,
                        evidence_type.BENIGN,
                        True,
                        strength,
                        f"ClinGen BS1: {gnomad.subpopulation_frequency} meets {applied_strength} threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )
                else:
                    return RuleResult(
                        "BS1",
                        rule_type.GENERAL,
                        evidence_type.BENIGN,
                        False,
                        evidence_strength.STRONG,
                        f"ClinGen BS1: {gnomad.subpopulation_frequency} does not meet threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )

        except Exception:
            return None  # Fall back to config


class Bs1_faf(abstract_rule):
    """
    BS1: Frequency of variant higher in population than expected based on disease frequency
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            Bs1.assess_rule,
            (
                class_info.VARIANT,
                class_info.VARIANT_GNOMAD_FAF,
                class_info.THRESHOLD_BS1,
            ),
        )


class Bs1_with_absolute(abstract_rule):
    """
    BS1: Frequency of variant higher in population than expected based on disease frequency
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT,
                class_info.VARIANT_GNOMAD_POPMAX,
                class_info.THRESHOLD_BS1,
                class_info.THRESHOLD_BS1_ABSOLUTE,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        gnomad: PopulationDatabases_gnomAD,
        threshold_bs1: float,
        threshold_bs1_absolute: int,
    ) -> RuleResult:
        # Try ClinGen gene-specific cutoff first
        gene_name = variant.gene_name
        clingen_db_path = _get_default_clingen_db_path()

        if clingen_db_path and clingen_db_path.exists():
            result = cls._assess_with_clingen(gene_name, gnomad, clingen_db_path)
            if result is not None:
                return result

        # Fall back to config threshold
        if gnomad.subpopulation_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif (
            gnomad.subpopulation_frequency > threshold_bs1
            and gnomad.subpopulation_allele_count >= threshold_bs1_absolute
        ):
            comment = f"Variant occures with a frequeny of {gnomad.subpopulation_frequency} and a total of {gnomad.subpopulation_allele_count} times in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} and a total of {gnomad.subpopulation_allele_count} times in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
        if gnomad.subpopulation == "None":
            comment = f"Variant does not occur in gnomAD, allele frequency in gnomAD is assumed to be 0."
        if gnomad.subpopulation == "ALL":
            comment = f"Variant has no entry for subopulation. Variant occurs with {gnomad.subpopulation_frequency} and a total of {gnomad.subpopulation_allele_count} times in gnomAD."
        return RuleResult(
            "BS1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.STRONG,
            comment,
        )

    @classmethod
    def _assess_with_clingen(
        cls,
        gene_name: str,
        gnomad: PopulationDatabases_gnomAD,
        clingen_db_path: pathlib.Path,
    ) -> Optional[RuleResult]:
        """Assess BS1 using ClinGen gene-specific frequency cutoffs."""
        if gnomad.subpopulation_frequency is None:
            return RuleResult(
                "BS1",
                rule_type.GENERAL,
                evidence_type.BENIGN,
                False,
                evidence_strength.STRONG,
                "gnomAD frequency is None.",
            )

        try:
            from check_clingen_rules import ClinGenRulesDB

            with ClinGenRulesDB(clingen_db_path) as db:
                applies, applied_strength, reason = db.check_frequency_applies(
                    gene_name, "BS1", gnomad.subpopulation_frequency
                )

                if applies:
                    if applied_strength == "Stand Alone":
                        strength = evidence_strength.STAND_ALONE
                    elif applied_strength == "Strong":
                        strength = evidence_strength.STRONG
                    elif applied_strength == "Supporting":
                        strength = evidence_strength.SUPPORTING
                    else:
                        strength = evidence_strength.STRONG

                    return RuleResult(
                        "BS1",
                        rule_type.GENERAL,
                        evidence_type.BENIGN,
                        True,
                        strength,
                        f"ClinGen BS1: {gnomad.subpopulation_frequency} meets {applied_strength} threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )
                else:
                    return RuleResult(
                        "BS1",
                        rule_type.GENERAL,
                        evidence_type.BENIGN,
                        False,
                        evidence_strength.STRONG,
                        f"ClinGen BS1: {gnomad.subpopulation_frequency} does not meet threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )

        except Exception:
            return None  # Fall back to config


class Bs1_with_supporting(abstract_rule):
    """
    BS1: Frequency of variant higher in population than expected based on disease frequency
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT,
                class_info.VARIANT_GNOMAD_POPMAX,
                class_info.THRESHOLD_BS1,
                class_info.THRESHOLD_BS1_SUPPORTING,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        gnomad: PopulationDatabases_gnomAD,
        threshold_bs1: float,
        threshold_bs1_supporting: float,
    ) -> RuleResult:
        # Try ClinGen gene-specific cutoff first
        gene_name = variant.gene_name
        clingen_db_path = _get_default_clingen_db_path()

        if clingen_db_path and clingen_db_path.exists():
            result = cls._assess_with_clingen(gene_name, gnomad, clingen_db_path)
            if result is not None:
                return result

        # Fall back to config threshold
        if gnomad.subpopulation_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif gnomad.subpopulation_frequency > threshold_bs1:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
            strength = evidence_strength.STRONG
        elif gnomad.subpopulation_frequency > threshold_bs1_supporting:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
            strength = evidence_strength.SUPPORTING
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
            strength = evidence_strength.STRONG
        if gnomad.subpopulation == "None":
            comment = f"Variant does not occur in gnomAD, allele frequency in gnomAD is assumed to be 0."
        if gnomad.subpopulation == "ALL":
            comment = f"Variant has no entry for subpopulation. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        return RuleResult(
            "BS1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            strength,
            comment,
        )

    @classmethod
    def _assess_with_clingen(
        cls,
        gene_name: str,
        gnomad: PopulationDatabases_gnomAD,
        clingen_db_path: pathlib.Path,
    ) -> Optional[RuleResult]:
        """Assess BS1 using ClinGen gene-specific frequency cutoffs."""
        if gnomad.subpopulation_frequency is None:
            return RuleResult(
                "BS1",
                rule_type.GENERAL,
                evidence_type.BENIGN,
                False,
                evidence_strength.STRONG,
                "gnomAD frequency is None.",
            )

        try:
            from check_clingen_rules import ClinGenRulesDB

            with ClinGenRulesDB(clingen_db_path) as db:
                applies, applied_strength, reason = db.check_frequency_applies(
                    gene_name, "BS1", gnomad.subpopulation_frequency
                )

                if applies:
                    if applied_strength == "Stand Alone":
                        strength = evidence_strength.STAND_ALONE
                    elif applied_strength == "Strong":
                        strength = evidence_strength.STRONG
                    elif applied_strength == "Supporting":
                        strength = evidence_strength.SUPPORTING
                    else:
                        strength = evidence_strength.STRONG

                    return RuleResult(
                        "BS1",
                        rule_type.GENERAL,
                        evidence_type.BENIGN,
                        True,
                        strength,
                        f"ClinGen BS1: {gnomad.subpopulation_frequency} meets {applied_strength} threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )
                else:
                    return RuleResult(
                        "BS1",
                        rule_type.GENERAL,
                        evidence_type.BENIGN,
                        False,
                        evidence_strength.STRONG,
                        f"ClinGen BS1: {gnomad.subpopulation_frequency} does not meet threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )

        except Exception:
            return None  # Fall back to config


class Bs1_with_supporting_faf(abstract_rule):
    """
    BS1: Frequency of variant higher in population than expected based on disease frequency
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            Bs1_with_supporting.assess_rule,
            (
                class_info.VARIANT,
                class_info.VARIANT_GNOMAD_FAF,
                class_info.THRESHOLD_BS1,
                class_info.THRESHOLD_BS1_SUPPORTING,
            ),
        )


class Bs1_clingen(abstract_rule):
    """
    BS1 based on ClinGen gene-specific frequency cutoffs.

    Uses ClinGen_rules_frequency_cutoffs database to check if a variant's
    population frequency exceeds the gene-specific threshold for BS1.

    Evidence strength follows ClinGen specification:
    - Strong: frequency >= gene-specific BS1 threshold
    - Supporting: if defined with lower threshold

    Reference: ClinGen frequency cutoffs specifications
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT,
                class_info.VARIANT_GNOMAD_POPMAX,
                class_info.THRESHOLD_BS1,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        gnomad: PopulationDatabases_gnomAD,
        threshold_bs1: float,
    ) -> RuleResult:
        """BS1 assess_rule with ClinGen gene-specific cutoffs."""
        gene_name = variant.gene_name
        clingen_db_path = _get_default_clingen_db_path()

        if clingen_db_path and clingen_db_path.exists():
            result = cls._assess_with_clingen(gene_name, gnomad, clingen_db_path)
            if result is not None:
                return result

        # Fall back to config threshold
        if gnomad.subpopulation_frequency is None:
            raise ValueError(
                f"The gnomAD allele frequency is None. Please check variant import."
            )
        elif gnomad.subpopulation_frequency > threshold_bs1:
            comment = f"Variant occurs with {gnomad.subpopulation_frequency} in gnomAD (>{threshold_bs1})."
            result = True
        else:
            comment = f"Variant frequency {gnomad.subpopulation_frequency} does not exceed threshold {threshold_bs1}."
            result = False
        if gnomad.subpopulation == "None":
            comment = f"Variant does not occur in gnomAD."
        return RuleResult(
            "BS1",
            rule_type.GENERAL,
            evidence_type.BENIGN,
            result,
            evidence_strength.STRONG,
            comment,
        )

    @classmethod
    def _assess_with_clingen(
        cls,
        gene_name: str,
        gnomad: PopulationDatabases_gnomAD,
        clingen_db_path: pathlib.Path,
    ) -> Optional[RuleResult]:
        """Assess BS1 using ClinGen gene-specific frequency cutoffs."""
        if gnomad.subpopulation_frequency is None:
            return RuleResult(
                "BS1",
                rule_type.GENERAL,
                evidence_type.BENIGN,
                False,
                evidence_strength.STRONG,
                "gnomAD frequency is None.",
            )

        try:
            from check_clingen_rules import ClinGenRulesDB

            with ClinGenRulesDB(clingen_db_path) as db:
                applies, applied_strength, reason = db.check_frequency_applies(
                    gene_name, "BS1", gnomad.subpopulation_frequency
                )

                if applies:
                    if applied_strength == "Stand Alone":
                        strength = evidence_strength.STAND_ALONE
                    elif applied_strength == "Strong":
                        strength = evidence_strength.STRONG
                    elif applied_strength == "Supporting":
                        strength = evidence_strength.SUPPORTING
                    else:
                        strength = evidence_strength.STRONG

                    return RuleResult(
                        "BS1",
                        rule_type.GENERAL,
                        evidence_type.BENIGN,
                        True,
                        strength,
                        f"ClinGen BS1: {gnomad.subpopulation_frequency} meets {applied_strength} threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )
                else:
                    return RuleResult(
                        "BS1",
                        rule_type.GENERAL,
                        evidence_type.BENIGN,
                        False,
                        evidence_strength.STRONG,
                        f"ClinGen BS1: {gnomad.subpopulation_frequency} does not meet threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )

        except Exception:
            return None  # Fall back to config
