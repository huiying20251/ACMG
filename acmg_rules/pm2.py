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
from variant import PopulationDatabases_gnomAD, VariantInfo


def _get_default_clingen_db_path() -> pathlib.Path:
    """Get the default ClinGen rules database path."""
    return pathlib.Path(__file__).parent.parent.parent / "clingen_rules.db"


class Pm2(abstract_rule):
    """
    PM2: Variant is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate.

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
                class_info.THRESHOLD_PM2,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        gnomad: PopulationDatabases_gnomAD,
        threshold_pm2: float,
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
        elif gnomad.subpopulation_frequency <= threshold_pm2:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
        if gnomad.subpopulation == "None":
            comment = f"Variant does not occur in gnomAD, allele frequency in gnomAD is assumed to be 0."
        if gnomad.subpopulation == "ALL":
            comment = f"Variant has no entry for subpopulation. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        return RuleResult(
            "PM2",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )

    @classmethod
    def _assess_with_clingen(
        cls,
        gene_name: str,
        gnomad: PopulationDatabases_gnomAD,
        clingen_db_path: pathlib.Path,
    ) -> Optional[RuleResult]:
        """Assess PM2 using ClinGen gene-specific frequency cutoffs."""
        if gnomad.subpopulation_frequency is None:
            return RuleResult(
                "PM2",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.MODERATE,
                "gnomAD frequency is None.",
            )

        try:
            from check_clingen_rules import ClinGenRulesDB

            with ClinGenRulesDB(clingen_db_path) as db:
                applies, applied_strength, reason = db.check_frequency_applies(
                    gene_name, "PM2", gnomad.subpopulation_frequency
                )

                if applies:
                    # Map ClinGen strength to evidence_strength enum
                    if applied_strength == "Very Strong":
                        strength = evidence_strength.VERY_STRONG
                    elif applied_strength == "Strong":
                        strength = evidence_strength.STRONG
                    elif applied_strength == "Moderate":
                        strength = evidence_strength.MODERATE
                    elif applied_strength == "Supporting":
                        strength = evidence_strength.SUPPORTING
                    else:
                        strength = evidence_strength.MODERATE

                    return RuleResult(
                        "PM2",
                        rule_type.GENERAL,
                        evidence_type.PATHOGENIC,
                        True,
                        strength,
                        f"ClinGen PM2: {gnomad.subpopulation_frequency} meets {applied_strength} threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )
                else:
                    return RuleResult(
                        "PM2",
                        rule_type.GENERAL,
                        evidence_type.PATHOGENIC,
                        False,
                        evidence_strength.MODERATE,
                        f"ClinGen PM2: {gnomad.subpopulation_frequency} does not meet threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )

        except Exception:
            return None  # Fall back to config


class Pm2_supporting(abstract_rule):
    """
    PM2: Variant is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate
    Default strength of PM2 is set to supporting following SVI recommendations
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
                class_info.THRESHOLD_PM2,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        gnomad: PopulationDatabases_gnomAD,
        threshold_pm2: float,
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
        elif gnomad.subpopulation_frequency <= threshold_pm2:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
        if gnomad.subpopulation == "None":
            comment = f"Variant does not occur in gnomAD, allele frequency in gnomAD is assumed to be 0."
        if gnomad.subpopulation == "ALL":
            comment = f"Variant has no entry for subpopulation. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        return RuleResult(
            "PM2",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )

    @classmethod
    def _assess_with_clingen(
        cls,
        gene_name: str,
        gnomad: PopulationDatabases_gnomAD,
        clingen_db_path: pathlib.Path,
    ) -> Optional[RuleResult]:
        """Assess PM2 using ClinGen gene-specific frequency cutoffs."""
        if gnomad.subpopulation_frequency is None:
            return RuleResult(
                "PM2",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                "gnomAD frequency is None.",
            )

        try:
            from check_clingen_rules import ClinGenRulesDB

            with ClinGenRulesDB(clingen_db_path) as db:
                applies, applied_strength, reason = db.check_frequency_applies(
                    gene_name, "PM2", gnomad.subpopulation_frequency
                )

                if applies:
                    if applied_strength == "Very Strong":
                        strength = evidence_strength.VERY_STRONG
                    elif applied_strength == "Strong":
                        strength = evidence_strength.STRONG
                    elif applied_strength == "Moderate":
                        strength = evidence_strength.MODERATE
                    elif applied_strength == "Supporting":
                        strength = evidence_strength.SUPPORTING
                    else:
                        strength = evidence_strength.SUPPORTING

                    return RuleResult(
                        "PM2",
                        rule_type.GENERAL,
                        evidence_type.PATHOGENIC,
                        True,
                        strength,
                        f"ClinGen PM2: {gnomad.subpopulation_frequency} meets {applied_strength} threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )
                else:
                    return RuleResult(
                        "PM2",
                        rule_type.GENERAL,
                        evidence_type.PATHOGENIC,
                        False,
                        evidence_strength.SUPPORTING,
                        f"ClinGen PM2: {gnomad.subpopulation_frequency} does not meet threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )

        except Exception:
            return None  # Fall back to config


class Pm2_supporting_faf(abstract_rule):
    """
    PM2: Variant is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate
    Default strength of PM2 is set to supporting following SVI recommendations
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            Pm2_supporting.assess_rule,
            (
                class_info.VARIANT,
                class_info.VARIANT_GNOMAD_FAF,
                class_info.THRESHOLD_PM2,
            ),
        )


class Pm2_supporting_less(abstract_rule):
    """
    PM2: Varinat is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate
    Default strength of PM2 is set to supporting following SVI recommendations
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
                class_info.THRESHOLD_PM2,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        gnomad: PopulationDatabases_gnomAD,
        threshold_pm2: float,
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
        elif gnomad.subpopulation_frequency < threshold_pm2:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
        if gnomad.subpopulation == "None":
            comment = f"Variant does not occur in gnomAD, allele frequency in gnomAD is assumed to be 0."
        if gnomad.subpopulation == "ALL":
            comment = f"Variant has no entry for subpopulation. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        return RuleResult(
            "PM2",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )

    @classmethod
    def _assess_with_clingen(
        cls,
        gene_name: str,
        gnomad: PopulationDatabases_gnomAD,
        clingen_db_path: pathlib.Path,
    ) -> Optional[RuleResult]:
        """Assess PM2 using ClinGen gene-specific frequency cutoffs."""
        if gnomad.subpopulation_frequency is None:
            return RuleResult(
                "PM2",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                "gnomAD frequency is None.",
            )

        try:
            from check_clingen_rules import ClinGenRulesDB

            with ClinGenRulesDB(clingen_db_path) as db:
                applies, applied_strength, reason = db.check_frequency_applies(
                    gene_name, "PM2", gnomad.subpopulation_frequency
                )

                if applies:
                    if applied_strength == "Very Strong":
                        strength = evidence_strength.VERY_STRONG
                    elif applied_strength == "Strong":
                        strength = evidence_strength.STRONG
                    elif applied_strength == "Moderate":
                        strength = evidence_strength.MODERATE
                    elif applied_strength == "Supporting":
                        strength = evidence_strength.SUPPORTING
                    else:
                        strength = evidence_strength.SUPPORTING

                    return RuleResult(
                        "PM2",
                        rule_type.GENERAL,
                        evidence_type.PATHOGENIC,
                        True,
                        strength,
                        f"ClinGen PM2: {gnomad.subpopulation_frequency} meets {applied_strength} threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )
                else:
                    return RuleResult(
                        "PM2",
                        rule_type.GENERAL,
                        evidence_type.PATHOGENIC,
                        False,
                        evidence_strength.SUPPORTING,
                        f"ClinGen PM2: {gnomad.subpopulation_frequency} does not meet threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )

        except Exception:
            return None  # Fall back to config


class Pm2_supporting_less_faf(abstract_rule):
    """
    PM2: Varinat is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate
    Default strength of PM2 is set to supporting following SVI recommendations
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            Pm2_supporting_less.assess_rule,
            (
                class_info.VARIANT,
                class_info.VARIANT_GNOMAD_FAF,
                class_info.THRESHOLD_PM2,
            ),
        )


class Pm2_supporting_no_ins_del_indel(abstract_rule):
    """
    PM2: Varinat is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate
    Default strength of PM2 is set to supporting following SVI recommendations
    PM2 does not apply to insertions, deletions or indels
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
                class_info.THRESHOLD_PM2,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant: VariantInfo,
        gnomad: PopulationDatabases_gnomAD,
        threshold_pm2: float,
    ) -> RuleResult:
        # Check if variant is indel first
        if len(variant.var_ref) != 1 or len(variant.var_obs) != 1:
            comment = f"PM2 does not apply to insertions, deletions or delins."
            result = False
            return RuleResult(
                "PM2",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                result,
                evidence_strength.SUPPORTING,
                comment,
            )

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
        elif gnomad.subpopulation_frequency <= threshold_pm2:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = True
        else:
            comment = f"Variant occures with {gnomad.subpopulation_frequency} in gnomAD subpopulation {gnomad.subpopulation}."
            result = False
        if gnomad.subpopulation == "None":
            comment = f"Variant does not occur in gnomAD, allele frequency in gnomAD is assumed to be 0."
        if gnomad.subpopulation == "ALL":
            comment = f"Variant has no entry for subpopulation. Variant occurs with {gnomad.subpopulation_frequency} in gnomAD."
        return RuleResult(
            "PM2",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )

    @classmethod
    def _assess_with_clingen(
        cls,
        gene_name: str,
        gnomad: PopulationDatabases_gnomAD,
        clingen_db_path: pathlib.Path,
    ) -> Optional[RuleResult]:
        """Assess PM2 using ClinGen gene-specific frequency cutoffs."""
        if gnomad.subpopulation_frequency is None:
            return RuleResult(
                "PM2",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                "gnomAD frequency is None.",
            )

        try:
            from check_clingen_rules import ClinGenRulesDB

            with ClinGenRulesDB(clingen_db_path) as db:
                applies, applied_strength, reason = db.check_frequency_applies(
                    gene_name, "PM2", gnomad.subpopulation_frequency
                )

                if applies:
                    if applied_strength == "Very Strong":
                        strength = evidence_strength.VERY_STRONG
                    elif applied_strength == "Strong":
                        strength = evidence_strength.STRONG
                    elif applied_strength == "Moderate":
                        strength = evidence_strength.MODERATE
                    elif applied_strength == "Supporting":
                        strength = evidence_strength.SUPPORTING
                    else:
                        strength = evidence_strength.SUPPORTING

                    return RuleResult(
                        "PM2",
                        rule_type.GENERAL,
                        evidence_type.PATHOGENIC,
                        True,
                        strength,
                        f"ClinGen PM2: {gnomad.subpopulation_frequency} meets {applied_strength} threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )
                else:
                    return RuleResult(
                        "PM2",
                        rule_type.GENERAL,
                        evidence_type.PATHOGENIC,
                        False,
                        evidence_strength.SUPPORTING,
                        f"ClinGen PM2: {gnomad.subpopulation_frequency} does not meet threshold for {gene_name}. {reason[:200] if reason else ''}",
                    )

        except Exception:
            return None  # Fall back to config


class Pm2_supporting_no_ins_del_indel_faf(abstract_rule):
    """
    PM2: Varinat is absent from control population
    In case of recessive disorders: Variant occurres less than expected carrier rate
    Default strength of PM2 is set to supporting following SVI recommendations
    PM2 does not apply to insertions, deletions or indels
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            Pm2_supporting_no_ins_del_indel.assess_rule,
            (
                class_info.VARIANT,
                class_info.VARIANT_GNOMAD_FAF,
                class_info.THRESHOLD_PM2,
            ),
        )
