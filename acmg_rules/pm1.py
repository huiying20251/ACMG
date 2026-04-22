#!/usr/bin/env python3

from typing import Callable, Optional, Tuple, List
import pathlib

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    evidence_type,
    rule_type,
)
from information import Classification_Info, Info
from variant import PopulationDatabases
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds


def check_bed_file(
    chrom: str,
    pos: int,
    bed_path: pathlib.Path,
) -> bool:
    """
    Check if a genomic position falls within any interval in a BED file.

    Args:
        chrom: Chromosome (e.g., "17" or "chr17")
        pos: Genomic position
        bed_path: Path to BED file

    Returns:
        True if position is in any interval, False otherwise
    """
    import csv

    chrom_clean = chrom.replace("chr", "")
    with open(bed_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if row[0].replace("chr", "") == chrom_clean:
                start = int(row[1])
                end = int(row[2])
                if start <= pos <= end:
                    return True
    return False


def get_uniprot_domains(
    gene_name: str,
    protein_position: int,
    uniprot_domain_path: pathlib.Path,
) -> Optional[dict]:
    """
    Get UniProt protein domain for a given protein position.

    Args:
        gene_name: Gene name (e.g., "BRCA1")
        protein_position: Protein position (amino acid number)
        uniprot_domain_path: Path to UniProt domain BED file

    Returns:
        Dict with domain info or None
    """
    import csv

    with open(uniprot_domain_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("gene_name") == gene_name or row.get("gene") == gene_name:
                domain_start = int(row["protein_start"])
                domain_end = int(row["protein_end"])
                if domain_start <= protein_position <= domain_end:
                    return {
                        "domain_name": row.get("domain_name", row.get("name", "unknown")),
                        "start": domain_start,
                        "end": domain_end,
                        "source": row.get("source", "UniProt"),
                    }
    return None


def check_clinvar_density_in_domain(
    clinvar_client,
    gene_name: str,
    domain_start: int,
    domain_end: int,
    min_pathogenic: int = 3,
) -> Tuple[bool, int]:
    """
    Check ClinVar pathogenic variant density within a protein domain.

    Args:
        clinvar_client: ClinVarEutilsClient instance
        gene_name: Gene name
        domain_start: Domain start position (protein)
        domain_end: Domain end position (protein)
        min_pathogenic: Minimum pathogenic variants to apply PM1

    Returns:
        (applies_pm1, count)
    """
    count = clinvar_client.count_pathogenic_in_gene(
        gene_name=gene_name,
        domain_start=domain_start,
        domain_end=domain_end,
    )
    return count >= min_pathogenic, count


def check_clinvar_density_in_window(
    clinvar_client,
    chrom: str,
    pos: int,
    window_bp: int = 25,
    min_pathogenic: int = 4,
) -> Tuple[bool, int]:
    """
    Check ClinVar pathogenic variant density within a genomic window.

    Args:
        clinvar_client: ClinVarEutilsClient instance
        chrom: Chromosome
        pos: Genomic position
        window_bp: Window size in bp (default 25)
        min_pathogenic: Minimum pathogenic variants to apply PM1

    Returns:
        (applies_pm1, count)
    """
    variants = clinvar_client.get_pathogenic_variants_in_region(
        gene_name="",
        chr=chrom,
        start=pos,
        end=pos,
        window_bp=window_bp,
    )
    count = len([v for v in variants if is_pathogenic_or_likely_pathogenic(v)])
    return count >= min_pathogenic, count


def is_pathogenic_or_likely_pathogenic(variant: dict) -> bool:
    """
    Check if a ClinVar variant is pathogenic or likely pathogenic.

    Args:
        variant: ClinVar variant dict

    Returns:
        True if pathogenic/likely pathogenic
    """
    significance = variant.get("clinical_significance", "")
    if isinstance(significance, list) and len(significance) > 0:
        significance = significance[0]
    if isinstance(significance, dict):
        significance = significance.get("description", "")
    significance_lower = significance.lower() if significance else ""
    return "pathogenic" in significance_lower or "likely pathogenic" in significance_lower


class Pm1_enhanced(abstract_rule):
    """
    PM1 Enhanced: Variant located in mutational hot spot or critical protein region.

    Enhanced logic:
    1. Check predefined hotspot BED file
    2. If not found, query UniProt for protein domains
    3. Check ClinVar density (≥3 pathogenic in domain OR ±25bp ≥4 pathogenic)
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_HOTSPOT_ANNOTATION,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant_in_hotspot: bool,
    ) -> RuleResult:
        if variant_in_hotspot:
            comment = "Variant in predefined mutational hotspot."
            result = True
        else:
            comment = "Variant not in predefined hotspot BED file."
            result = False
        return RuleResult(
            "PM1",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )

    @classmethod
    def assess_rule_with_clinvar(
        cls,
        chrom: str,
        pos: int,
        gene_name: str,
        protein_position: int,
        hotspot_bed_path: Optional[pathlib.Path] = None,
        uniprot_domain_path: Optional[pathlib.Path] = None,
        clinvar_client=None,
        domain_min_pathogenic: int = 3,
        window_min_pathogenic: int = 4,
        window_bp: int = 25,
    ) -> RuleResult:
        """
        Enhanced PM1 assessment with ClinVar density check.

        Args:
            chrom: Chromosome
            pos: Genomic position
            gene_name: Gene name
            protein_position: Protein position (amino acid number)
            hotspot_bed_path: Path to predefined hotspot BED file
            uniprot_domain_path: Path to UniProt domain BED file
            clinvar_client: ClinVarEutilsClient instance
            domain_min_pathogenic: Min pathogenic in domain to apply PM1 (default 3)
            window_min_pathogenic: Min pathogenic in window to apply PM1 (default 4)
            window_bp: Genomic window size in bp (default 25)

        Returns:
            RuleResult for PM1
        """
        # Step 1: Check predefined hotspot BED
        if hotspot_bed_path and hotspot_bed_path.exists():
            if check_bed_file(chrom, pos, hotspot_bed_path):
                return RuleResult(
                    "PM1",
                    rule_type.GENERAL,
                    evidence_type.PATHOGENIC,
                    True,
                    evidence_strength.MODERATE,
                    "Variant in predefined mutational hotspot BED file.",
                )

        # Step 2: Check UniProt domain
        if uniprot_domain_path and uniprot_domain_path.exists():
            domain_info = get_uniprot_domains(gene_name, protein_position, uniprot_domain_path)
            if domain_info:
                # Step 2a: Check domain ClinVar density
                if clinvar_client:
                    applies, count = check_clinvar_density_in_domain(
                        clinvar_client,
                        gene_name,
                        domain_info["start"],
                        domain_info["end"],
                        min_pathogenic=domain_min_pathogenic,
                    )
                    if applies:
                        return RuleResult(
                            "PM1",
                            rule_type.GENERAL,
                            evidence_type.PATHOGENIC,
                            True,
                            evidence_strength.MODERATE,
                            f"Domain {domain_info['domain_name']} has {count} pathogenic variants (≥{domain_min_pathogenic}).",
                        )

                # Step 2b: Check ±window_bp genomic range density
                if clinvar_client:
                    applies, count = check_clinvar_density_in_window(
                        clinvar_client,
                        chrom,
                        pos,
                        window_bp=window_bp,
                        min_pathogenic=window_min_pathogenic,
                    )
                    if applies:
                        return RuleResult(
                            "PM1",
                            rule_type.GENERAL,
                            evidence_type.PATHOGENIC,
                            True,
                            evidence_strength.MODERATE,
                            f"{count} pathogenic variants within ±{window_bp}bp (≥{window_min_pathogenic}).",
                        )

        # No PM1 applies
        return RuleResult(
            "PM1",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            False,
            evidence_strength.MODERATE,
            "No significant domain or ClinVar density found for PM1.",
        )


class Pm1(abstract_rule):
    """
    PM1: Variant located in mutational hot spot or citical protein region
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.VARIANT_HOTSPOT_ANNOTATION,),
        )

    @classmethod
    def assess_rule(cls, variant_in_hotspot: bool) -> RuleResult:
        if variant_in_hotspot:
            comment = f"Variant in mutational hotspot."
            result = True
        else:
            comment = "Variant not in mutational hotspot"
            result = False
        return RuleResult(
            "PM1",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )


class Pm1_supporting(abstract_rule):
    """
    PM1: Variant located in mutational hot spot or citical protein region
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.VARIANT_HOTSPOT_ANNOTATION,),
        )

    @classmethod
    def assess_rule(cls, variant_in_hotspot: bool) -> RuleResult:
        if variant_in_hotspot:
            comment = f"Variant in mutational hotspot."
            result = True
        else:
            comment = "Variant not in mutational hotspot"
            result = False
        return RuleResult(
            "PM1",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.SUPPORTING,
            comment,
        )


class Pm1_tp53(abstract_rule):
    """
    PM1: Variant located in mutational hot spot or citical protein region
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_HOTSPOT_ANNOTATION,
                class_info.VARIANT_CANCERHOTSPOTS,
                class_info.THRESHOLD_CANCERHOTSPOTS_AC,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        variant_in_hotspot: bool,
        cancerhotspots: PopulationDatabases,
        threshold_cancerhotspots_ac: float,
        prediction_dict: dict[str, float],
        threshold: Threshold,
    ) -> RuleResult:
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)
        if variant_in_hotspot:
            comment = f"Variant in defined mutational hotspot."
            result = True
        elif cancerhotspots.count is None:
            comment = f"Variant not located in defined mutational hotspot and no Cancer Hotspots entry database for the variant."
            result = False
        elif cancerhotspots.count >= threshold_cancerhotspots_ac:
            comment = f"Variant occurs in Cancer Hotspots {cancerhotspots.count} times."
            result = True
        else:
            comment = f"Variant not in mutational hotspot and occurence of variant in Cancer Hotspots ({cancerhotspots.count}) does not meet threshold ({threshold_cancerhotspots_ac})."
            result = False
        if num_thresholds_met is None:
            comment = (
                "ATTENTION: No splicing prediction available for this variant. "
                + comment
            )
        elif num_thresholds_met == 0:
            comment = f"Variant is not predicted to not affect splicing, therefore PM1 does not apply."
            result = False
        else:
            comment = (
                comment
                + f" {threshold.name} predicts no effect on splicing for the variant (threshold: {threshold.thresholds[num_thresholds_met -1]}, value: {prediction_value})."
            )
        return RuleResult(
            "PM1",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )


class Pm1_clingen(abstract_rule):
    """
    PM1 based on ClinGen gene-specific PM1 domain specifications.

    Uses ClinGen_rules_pm1_details database to check if a variant's
    protein position falls within a ClinGen-defined critical domain.

    Evidence strength follows ClinGen specification:
    - Moderate: variant in defined critical domain
    - Supporting: variant in defined domain (if specified as Supporting)

    Reference: ClinGen PM1 domain specifications
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.VARIANT_HOTSPOT_ANNOTATION,),
        )

    @classmethod
    def assess_rule(cls, variant_in_hotspot: bool) -> RuleResult:
        """
        Simplified assess_rule using predefined hotspot annotation.
        For full ClinGen-based assessment, use assess_rule_with_clingen().
        """
        if variant_in_hotspot:
            comment = "Variant in ClinGen-defined mutational hotspot."
            result = True
        else:
            comment = "Variant not in ClinGen-defined hotspot"
            result = False
        return RuleResult(
            "PM1",
            rule_type.GENERAL,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.MODERATE,
            comment,
        )

    @classmethod
    def assess_rule_with_clingen(
        cls,
        gene_name: str,
        protein_position: int,
        clingen_db_path: Optional[pathlib.Path] = None,
    ) -> RuleResult:
        """
        Assess PM1 using ClinGen gene-specific domain specifications.

        Args:
            gene_name: Gene symbol (e.g., "BRCA1", "MYH7")
            protein_position: Protein amino acid position (e.g., 1850)
            clingen_db_path: Path to clingen_rules.db

        Returns:
            RuleResult for PM1 with ClinGen-based strength
        """
        if clingen_db_path is None or not clingen_db_path.exists():
            return RuleResult(
                "PM1",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.MODERATE,
                "ClinGen PM1 database not configured.",
            )

        try:
            from check_clingen_rules import ClinGenRulesDB

            with ClinGenRulesDB(clingen_db_path) as db:
                pm1_detail = db.get_pm1_for_position(gene_name, protein_position)

                if pm1_detail is None:
                    return RuleResult(
                        "PM1",
                        rule_type.GENERAL,
                        evidence_type.PATHOGENIC,
                        False,
                        evidence_strength.MODERATE,
                        f"Protein position {protein_position} not in ClinGen-defined PM1 domain for {gene_name}.",
                    )

                # Determine evidence strength from ClinGen specification
                if pm1_detail.strength == "Moderate":
                    strength = evidence_strength.MODERATE
                elif pm1_detail.strength == "Supporting":
                    strength = evidence_strength.SUPPORTING
                else:
                    strength = evidence_strength.MODERATE

                domain_info = f"{pm1_detail.domain} ({pm1_detail.protein_positions})" if pm1_detail.domain else pm1_detail.protein_positions

                return RuleResult(
                    "PM1",
                    rule_type.GENERAL,
                    evidence_type.PATHOGENIC,
                    True,
                    strength,
                    f"Variant in ClinGen-defined domain: {domain_info} (strength: {pm1_detail.strength}).",
                )

        except Exception as e:
            return RuleResult(
                "PM1",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.MODERATE,
                f"ClinGen PM1 assessment error: {str(e)}",
            )
