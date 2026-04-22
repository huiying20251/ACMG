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
from acmg_rules.functional_splicing_assay_utils import (
    summarise_func_data,
)
from information import Info, Classification_Info
from variant import FunctionalData
from check_splicevardb import SpliceVarDB


class Ps3(abstract_rule):
    """
    PS3: Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product.
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.FUNCTIONAL_ASSAY,),
        )

    @classmethod
    def assess_rule(
        cls,
        func_data: list[FunctionalData],
    ) -> RuleResult:
        pathogenic_count, benign_count, uncertain_count = summarise_func_data(func_data)
        if not func_data:
            result = False
            comment = f"No functional assay performed."
        elif pathogenic_count == 0:
            result = False
            comment = f"None of the {len(func_data)} functional assay(s) performed indicate pathogenicity of the variant."
        elif benign_count > 0 and pathogenic_count > 0:
            result = False
            comment = f"{pathogenic_count} of the {len(func_data)} perforemd assays indicate that the variant is pathogenic and {benign_count} of the {len(func_data)} indicate that the variant is benign. Due to this conflicting evidenced PS3 can not be applied."
        elif pathogenic_count > 0:
            result = True
            comment = f"{pathogenic_count} of the {len(func_data)} performed assays indicate that the variant is pathogenic."
            if uncertain_count != 0:
                comment = (
                    comment
                    + f" ATTENTION: {uncertain_count} of the {len(func_data)} performed assays show no clear result."
                )
        else:
            result = False
            comment = (
                "Functional assay performed but results do not indicate pathogenicity."
            )
        return RuleResult(
            "PS3",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            evidence_strength.STRONG,
            comment,
        )


def check_functional_region_bed(
    chrom: str,
    pos: int,
    functional_bed_path: pathlib.Path,
) -> bool:
    """
    Check if a genomic position falls within a functional assay region.

    Args:
        chrom: Chromosome (e.g., "17" or "chr17")
        pos: Genomic position
        functional_bed_path: Path to functional assay regions BED file

    Returns:
        True if position is in a functional region, False otherwise
    """
    import csv

    chrom_clean = chrom.replace("chr", "")
    with open(functional_bed_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if row[0].replace("chr", "") == chrom_clean:
                start = int(row[1])
                end = int(row[2])
                if start <= pos <= end:
                    return True
    return False


class Ps3_functional_region(abstract_rule):
    """
    PS3/BS3: Based on functional assay regions BED file.

    - Variant in functional assay region → PS3_Supporting (supports pathogenic)
    - Variant NOT in functional assay region → BS3_Supporting (supports benign)

    This is used when functional assay data needs to be interpreted based on
    whether the variant falls within a known functional domain/region.
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.FUNCTIONAL_ASSAY,),
        )

    @classmethod
    def assess_rule(
        cls,
        func_data: list[FunctionalData],
    ) -> RuleResult:
        """
        Assess PS3/BS3 based on functional assay regions.

        Note: This is a simplified version that doesn't use BED file checking.
        Use assess_rule_with_bed() for full functionality with BED file.
        """
        # If no functional data, rule doesn't apply
        if not func_data:
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                "No functional assay data available.",
            )

        # Summarize functional data
        pathogenic_count, benign_count, uncertain_count = summarise_func_data(func_data)

        # If any functional assay shows pathogenic effect
        if pathogenic_count > 0:
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                True,
                evidence_strength.SUPPORTING,
                f"{pathogenic_count} of {len(func_data)} assays indicate pathogenic effect.",
            )

        # If any functional assay shows benign effect
        if benign_count > 0:
            return RuleResult(
                "BS3",
                rule_type.PROTEIN,
                evidence_type.BENIGN,
                True,
                evidence_strength.SUPPORTING,
                f"{benign_count} of {len(func_data)} assays indicate benign effect.",
            )

        # Uncertain results
        return RuleResult(
            "PS3",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            False,
            evidence_strength.SUPPORTING,
            f"No conclusive functional assay results ({uncertain_count} uncertain).",
        )

    @classmethod
    def assess_rule_with_bed(
        cls,
        chrom: str,
        pos: int,
        func_data: list[FunctionalData],
        functional_bed_path: Optional[pathlib.Path] = None,
    ) -> RuleResult:
        """
        Assess PS3/BS3 with functional region BED file checking.

        Args:
            chrom: Chromosome
            pos: Genomic position
            func_data: Functional assay data
            functional_bed_path: Path to functional assay regions BED file

        Returns:
            RuleResult for PS3 or BS3 depending on region
        """
        # First check functional data
        if func_data:
            pathogenic_count, benign_count, uncertain_count = summarise_func_data(func_data)

            # Conflicting results - can't apply
            if pathogenic_count > 0 and benign_count > 0:
                return RuleResult(
                    "PS3",
                    rule_type.PROTEIN,
                    evidence_type.PATHOGENIC,
                    False,
                    evidence_strength.SUPPORTING,
                    "Conflicting functional assay results.",
                )

            # Pathogenic functional data
            if pathogenic_count > 0:
                # Check if in functional region
                if functional_bed_path and functional_bed_path.exists():
                    if check_functional_region_bed(chrom, pos, functional_bed_path):
                        return RuleResult(
                            "PS3",
                            rule_type.PROTEIN,
                            evidence_type.PATHOGENIC,
                            True,
                            evidence_strength.SUPPORTING,
                            f"Variant in functional assay region, {pathogenic_count} assays support pathogenicity.",
                        )
                    else:
                        return RuleResult(
                            "BS3",
                            rule_type.PROTEIN,
                            evidence_type.BENIGN,
                            True,
                            evidence_strength.SUPPORTING,
                            f"Variant NOT in functional assay region, {pathogenic_count} assays support pathogenicity.",
                        )
                else:
                    return RuleResult(
                        "PS3",
                        rule_type.PROTEIN,
                        evidence_type.PATHOGENIC,
                        True,
                        evidence_strength.SUPPORTING,
                        f"{pathogenic_count} assays support pathogenicity.",
                    )

            # Benign functional data
            if benign_count > 0:
                return RuleResult(
                    "BS3",
                    rule_type.PROTEIN,
                    evidence_type.BENIGN,
                    True,
                    evidence_strength.SUPPORTING,
                    f"{benign_count} assays support benign effect.",
                )

        # No functional data - check BED file only
        if functional_bed_path and functional_bed_path.exists():
            in_functional_region = check_functional_region_bed(chrom, pos, functional_bed_path)
            if in_functional_region:
                return RuleResult(
                    "PS3",
                    rule_type.PROTEIN,
                    evidence_type.PATHOGENIC,
                    True,
                    evidence_strength.SUPPORTING,
                    "Variant in functional assay region (no functional data available).",
                )
            else:
                return RuleResult(
                    "BS3",
                    rule_type.PROTEIN,
                    evidence_type.BENIGN,
                    True,
                    evidence_strength.SUPPORTING,
                    "Variant NOT in functional assay region (no functional data available).",
                )

        # No data available
        return RuleResult(
            "PS3",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            False,
            evidence_strength.SUPPORTING,
            "No functional assay data or functional region information available.",
        )


class Ps3_splicevardb(abstract_rule):
    """
    PS3/BS3 based on SpliceVarDB functional evidence database.

    For non-canonical splice variants:
    - Splice-altering in SpliceVarDB → PS3_Supporting
    - Normal in SpliceVarDB → BS3_Supporting
    - Low-frequency → no evidence applied
    - Not found → no evidence applied
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.FUNCTIONAL_ASSAY,),
        )

    @classmethod
    def assess_rule(
        cls,
        func_data: list[FunctionalData],
    ) -> RuleResult:
        """
        Simplified assess_rule - requires manual SpliceVarDB lookup.

        Use assess_rule_with_splicevardb() for SpliceVarDB-based assessment.
        """
        if not func_data:
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                "No functional assay data available.",
            )

        pathogenic_count, benign_count, uncertain_count = summarise_func_data(func_data)

        if pathogenic_count > 0:
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                True,
                evidence_strength.SUPPORTING,
                f"{pathogenic_count} assays indicate pathogenic effect.",
            )

        if benign_count > 0:
            return RuleResult(
                "BS3",
                rule_type.PROTEIN,
                evidence_type.BENIGN,
                True,
                evidence_strength.SUPPORTING,
                f"{benign_count} assays indicate benign effect.",
            )

        return RuleResult(
            "PS3",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            False,
            evidence_strength.SUPPORTING,
            f"No conclusive functional assay results ({uncertain_count} uncertain).",
        )

    @classmethod
    def assess_rule_with_splicevardb(
        cls,
        gene: str,
        hgvs: Optional[str] = None,
        cDNA_pos: Optional[int] = None,
        chrom: Optional[str] = None,
        pos: Optional[int] = None,
        transcript: Optional[str] = None,
        splicevardb_path: Optional[pathlib.Path] = None,
        is_classic_splice: bool = False,
    ) -> RuleResult:
        """
        Assess PS3/BS3 using SpliceVarDB database.

        Args:
            gene: Gene symbol
            hgvs: HGVS cDNA notation
            cDNA_pos: cDNA position
            chrom: Chromosome
            pos: Genomic position
            transcript: Transcript ID
            splicevardb_path: Path to splicevardb_with_pmid.csv
            is_classic_splice: If True, variant is classic splice site (±1, ±2)

        Returns:
            RuleResult for PS3 (Splice-altering) or BS3 (Normal) or no evidence
            For classic splice variants, returns PVS1 instead of PS3 when
            SpliceVarDB shows splice-altering effect.
        """
        if splicevardb_path is None or not splicevardb_path.exists():
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                "SpliceVarDB not configured.",
            )

        try:
            db = SpliceVarDB(splicevardb_path)
            has_pathogenic, has_benign, reason = db.get_functional_evidence(
                gene=gene,
                hgvs=hgvs,
                cDNA_pos=cDNA_pos,
                chrom=chrom,
                pos=pos,
                transcript=transcript,
            )

            if has_pathogenic:
                # For classic splice variants, return PVS1 not PS3
                if is_classic_splice:
                    return RuleResult(
                        "PVS1",
                        rule_type.SPLICING,
                        evidence_type.PATHOGENIC,
                        True,
                        evidence_strength.VERY_STRONG,
                        f"SpliceVarDB: {reason}",
                    )
                return RuleResult(
                    "PS3",
                    rule_type.PROTEIN,
                    evidence_type.PATHOGENIC,
                    True,
                    evidence_strength.SUPPORTING,
                    f"SpliceVarDB: {reason}",
                )

            if has_benign:
                return RuleResult(
                    "BS3",
                    rule_type.PROTEIN,
                    evidence_type.BENIGN,
                    True,
                    evidence_strength.SUPPORTING,
                    f"SpliceVarDB: {reason}",
                )

            # Low-frequency or not found - no evidence
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                f"SpliceVarDB: {reason}",
            )

        except Exception as e:
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                f"SpliceVarDB error: {str(e)}",
            )


class Ps3_functional_db(abstract_rule):
    """
    PS3/BS3 based on functional_evidence.db (MAVE functional assay database).

    For variants with functional evidence in the database:
    - Loss-of-function or Gain-of-function → PS3
    - Functionally normal → BS3
    - Conflicting evidence (classification_review = "冲突") → no evidence applied
    - Not found → no evidence applied

    Evidence strength rules:
    - >=6 publications with consistent results → Moderate (2 points)
    - <6 publications with consistent results → Supporting (1 point)

    The database contains pre-aggregated functional evidence from multiple publications,
    deduplicated by (transcript_id, cDNA position).
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.FUNCTIONAL_ASSAY,),
        )

    @classmethod
    def assess_rule(
        cls,
        func_data: list[FunctionalData],
    ) -> RuleResult:
        """
        Simplified assess_rule - requires manual database lookup.

        Use assess_rule_with_db() for database-based assessment.
        """
        if not func_data:
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                "No functional assay data available.",
            )

        pathogenic_count, benign_count, uncertain_count = summarise_func_data(func_data)

        if pathogenic_count > 0:
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                True,
                evidence_strength.SUPPORTING,
                f"{pathogenic_count} assays indicate pathogenic effect.",
            )

        if benign_count > 0:
            return RuleResult(
                "BS3",
                rule_type.PROTEIN,
                evidence_type.BENIGN,
                True,
                evidence_strength.SUPPORTING,
                f"{benign_count} assays indicate benign effect.",
            )

        return RuleResult(
            "PS3",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            False,
            evidence_strength.SUPPORTING,
            f"No conclusive functional assay results ({uncertain_count} uncertain).",
        )

    @classmethod
    def assess_rule_with_db(
        cls,
        transcript_id: str,
        cdna_pos: str,
        gene_name: str,
        position: int,
        chrom: str,
        functional_db_path: Optional[pathlib.Path] = None,
        var_type: Optional[list] = None,
        ref_alt: Optional[str] = None,
    ) -> RuleResult:
        """
        Assess PS3/BS3 using functional_evidence.db database.

        Evidence strength rules:
        - >=6 publications with consistent results → Moderate (2 points)
        - <6 publications with consistent results → Supporting (1 point)
        - Conflicting evidence → no evidence applied

        Lookup strategy:
        - For splice/intronic variants: exact match with (transcript, cdna_pos, ref_alt)
        - For missense/other variants: match with (transcript, cdna_pos)

        Args:
            transcript_id: Transcript ID (e.g., "NM_005157.6")
            cdna_pos: cDNA position (e.g., "c.80-4")
            gene_name: Gene symbol
            position: Genomic position
            chrom: Chromosome
            functional_db_path: Path to functional_evidence.db
            var_type: List of VARTYPE (optional, for determining lookup strategy)
            ref_alt: Ref/Alt string like "C/T" (optional, for splice variant exact match)

        Returns:
            RuleResult for PS3 (Loss-of-function/Gain-of-function) or
            BS3 (Functionally normal) or no evidence
        """
        if functional_db_path is None or not functional_db_path.exists():
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                "Functional evidence DB not configured.",
            )

        # Check if variant is splice/intronic type (needs exact ref/alt match)
        is_splice_variant = False
        if var_type:
            from var_type import VARTYPE, VARTYPE_GROUPS
            splice_types = {
                VARTYPE.SPLICE_ACCEPTOR,
                VARTYPE.SPLICE_DONOR,
                VARTYPE.SPLICE_ACCEPTOR_VARIANT,
                VARTYPE.SPLICE_DONOR_VARIANT,
            }
            for vt in var_type:
                if vt in splice_types or vt in VARTYPE_GROUPS.INTRONIC.value:
                    is_splice_variant = True
                    break

        try:
            from functional_evidence_db import FunctionalEvidenceDB

            with FunctionalEvidenceDB(functional_db_path) as db:
                result = None

                if is_splice_variant and ref_alt:
                    # Splice variants: exact match with ref/alt
                    result = db.lookup_by_transcript_cdna_refalt(transcript_id, cdna_pos, ref_alt)
                    if result is None:
                        result = db.lookup_by_gene_cdna_refalt(gene_name, cdna_pos, ref_alt)
                else:
                    # Missense/other variants: match without ref/alt
                    result = db.lookup_by_transcript_cdna(transcript_id, cdna_pos)
                    if result is None:
                        result = db.lookup_by_gene_cdna(gene_name, cdna_pos)
                    if result is None:
                        result = db.lookup_by_gene_position(gene_name, position, chrom)

                if result is None:
                    return RuleResult(
                        "PS3",
                        rule_type.PROTEIN,
                        evidence_type.PATHOGENIC,
                        False,
                        evidence_strength.SUPPORTING,
                        "Variant not found in functional evidence database.",
                    )

                # Check for conflicting evidence
                if result.is_conflicting:
                    return RuleResult(
                        "PS3",
                        rule_type.PROTEIN,
                        evidence_type.PATHOGENIC,
                        False,
                        evidence_strength.SUPPORTING,
                        f"Conflicting functional evidence (multiple publications with different conclusions: {result.unique_publications}).",
                    )

                # Determine evidence strength based on publication count
                pub_count = len(result.unique_publications.split('|')) if result.unique_publications else 1
                if pub_count >= 6:
                    strength = evidence_strength.MODERATE
                else:
                    strength = evidence_strength.SUPPORTING

                evidence_type_val = evidence_type.PATHOGENIC
                rule_name = "PS3"

                # Pathogenic evidence
                if result.is_pathogenic:
                    return RuleResult(
                        rule_name,
                        rule_type.PROTEIN,
                        evidence_type_val,
                        True,
                        strength,
                        f"Functional evidence: {result.classification} (method: {result.mave_technique}, {pub_count} publication(s), PMID: {result.unique_publications}).",
                    )

                # Benign evidence
                if result.is_benign:
                    return RuleResult(
                        "BS3",
                        rule_type.PROTEIN,
                        evidence_type.BENIGN,
                        True,
                        strength,
                        f"Functional evidence: {result.classification} (method: {result.mave_technique}, {pub_count} publication(s), PMID: {result.unique_publications}).",
                    )

                # Should not reach here, but handle gracefully
                return RuleResult(
                    "PS3",
                    rule_type.PROTEIN,
                    evidence_type.PATHOGENIC,
                    False,
                    evidence_strength.SUPPORTING,
                    f"Unknown classification: {result.classification}",
                )

        except Exception as e:
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                f"Functional evidence DB error: {str(e)}",
            )


# Strength to score mapping for combining PS3 from multiple sources
PS3_STRENGTH_SCORES = {
    "supporting": 1,
    "moderate": 2,
    "strong": 4,
    "very_strong": 8,
}


def get_ps3_score_from_strength(strength: evidence_strength) -> int:
    """
    Convert evidence_strength to Bayes score.

    Args:
        strength: evidence_strength enum value

    Returns:
        Score: supporting=1, moderate=2, strong=4, very_strong=8
    """
    if strength == evidence_strength.VERY_STRONG:
        return 8
    elif strength == evidence_strength.STRONG:
        return 4
    elif strength == evidence_strength.MODERATE:
        return 2
    else:
        return 1  # SUPPORTING or default


class Ps3_combined(abstract_rule):
    """
    PS3 combined from multiple sources:

    1. SpliceVarDB functional evidence
    2. functional_evidence.db (MAVE functional assay data)
    3. Literature retrieval (functional studies)

    Takes the MAXIMUM score from all applicable sources.

    Evidence strength determination:
    - SpliceVarDB: Always SUPPORTING (1 point) unless modified
    - functional_evidence.db: >=6 publications → MODERATE (2 points), <6 → SUPPORTING (1 point)
    - Literature: high confidence → MODERATE (2 points), medium/low → SUPPORTING (1 point)

    Score mapping:
    - VERY_STRONG: 8 points
    - STRONG: 4 points
    - MODERATE: 2 points
    - SUPPORTING: 1 point
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (class_info.FUNCTIONAL_ASSAY, class_info.VARIANT_LITERATURE),
        )

    @classmethod
    def assess_rule(
        cls,
        func_data: list[FunctionalData],
        variant_literature: Optional[dict] = None,
    ) -> RuleResult:
        """
        Simplified assess_rule - requires manual combined assessment.

        Use assess_rule_with_all_sources() for full combined assessment.
        """
        if not func_data and not variant_literature:
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                "No functional evidence available.",
            )

        # Summarize functional data
        pathogenic_count, benign_count, uncertain_count = summarise_func_data(func_data)

        # If we have functional data, return basic result
        if pathogenic_count > 0:
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                True,
                evidence_strength.SUPPORTING,
                f"{pathogenic_count} assays indicate pathogenic effect.",
            )

        if benign_count > 0:
            return RuleResult(
                "BS3",
                rule_type.PROTEIN,
                evidence_type.BENIGN,
                True,
                evidence_strength.SUPPORTING,
                f"{benign_count} assays indicate benign effect.",
            )

        return RuleResult(
            "PS3",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            False,
            evidence_strength.SUPPORTING,
            f"No conclusive functional evidence ({uncertain_count} uncertain).",
        )

    @classmethod
    def assess_rule_with_all_sources(
        cls,
        transcript_id: str,
        cdna_pos: str,
        gene_name: str,
        position: int,
        chrom: str,
        functional_db_path: Optional[pathlib.Path] = None,
        splicevardb_path: Optional[pathlib.Path] = None,
        variant_literature: Optional[dict] = None,
        var_type: Optional[list] = None,
        ref_alt: Optional[str] = None,
    ) -> RuleResult:
        """
        Assess PS3/BS3 combining all available sources.

        Takes the MAXIMUM score from:
        1. SpliceVarDB
        2. functional_evidence.db
        3. Literature retrieval

        Args:
            transcript_id: Transcript ID
            cdna_pos: cDNA position (e.g., "c.80-4")
            gene_name: Gene symbol
            position: Genomic position
            chrom: Chromosome
            functional_db_path: Path to functional_evidence.db
            splicevardb_path: Path to splicevardb_with_pmid.csv
            variant_literature: Literature evidence dict with "ps3_literature" key
            var_type: List of VARTYPE (optional)
            ref_alt: Ref/Alt string (optional)

        Returns:
            RuleResult with combined PS3/BS3 from best source
        """
        # Collect all applicable PS3 results with scores
        ps3_results: list[tuple[RuleResult, int]] = []
        bs3_results: list[tuple[RuleResult, int]] = []

        # Check if variant is classic splice site (±1, ±2)
        is_classic = False
        if hgvs or cdna_pos:
            from null_variant_utils import is_classic_splice
            check_str = hgvs if hgvs else cdna_pos
            is_classic = is_classic_splice(check_str)

        # Source 1: SpliceVarDB
        if splicevardb_path and splicevardb_path.exists():
            splicevardb_result = cls._assess_splicevardb(
                gene=gene_name,
                hgvs=None,
                cDNA_pos=cdna_pos,
                chrom=chrom,
                pos=position,
                transcript=transcript_id,
                splicevardb_path=splicevardb_path,
                is_classic_splice=is_classic,
            )
            if splicevardb_result.status:
                score = get_ps3_score_from_strength(splicevardb_result.strength)
                if splicevardb_result.name == "BS3":
                    bs3_results.append((splicevardb_result, score))
                else:
                    ps3_results.append((splicevardb_result, score))

        # Source 2: functional_evidence.db
        if functional_db_path and functional_db_path.exists():
            func_db_result = cls._assess_functional_db(
                transcript_id=transcript_id,
                cdna_pos=cdna_pos,
                gene_name=gene_name,
                position=position,
                chrom=chrom,
                functional_db_path=functional_db_path,
                var_type=var_type,
                ref_alt=ref_alt,
            )
            if func_db_result.status:
                score = get_ps3_score_from_strength(func_db_result.strength)
                if func_db_result.name == "BS3":
                    bs3_results.append((func_db_result, score))
                else:
                    ps3_results.append((func_db_result, score))

        # Source 3: Literature
        if variant_literature and "ps3_literature" in variant_literature:
            lit_result = variant_literature["ps3_literature"]
            if lit_result and lit_result.get("applicable"):
                score = lit_result.get("score", 1)
                lit_rule_result = RuleResult(
                    "PS3",
                    rule_type.PROTEIN,
                    evidence_type.PATHOGENIC,
                    True,
                    cls._score_to_strength(score),
                    lit_result.get("reason", "Literature functional study support"),
                )
                ps3_results.append((lit_rule_result, score))

        # Determine best PS3 result (highest score)
        best_ps3 = None
        best_ps3_score = 0
        best_ps3_comment = ""

        for result, score in ps3_results:
            if score > best_ps3_score:
                best_ps3 = result
                best_ps3_score = score
                best_ps3_comment = result.comment

        # Determine best BS3 result (highest score)
        best_bs3 = None
        best_bs3_score = 0
        best_bs3_comment = ""

        for result, score in bs3_results:
            if score > best_bs3_score:
                best_bs3 = result
                best_bs3_score = score
                best_bs3_comment = result.comment

        # Return the best result (PS3 or BS3 with higher score)
        # PS3 and BS3 are mutually exclusive in ACMG
        if best_ps3 and best_bs3:
            # Both exist - conflicting evidence, don't apply
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                f"Conflicting functional evidence: PS3 ({best_ps3_comment}) vs BS3 ({best_bs3_comment}).",
            )

        if best_ps3:
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                True,
                best_ps3.strength,
                f"Combined PS3 (score={best_ps3_score}): {best_ps3_comment}",
            )

        if best_bs3:
            return RuleResult(
                "BS3",
                rule_type.PROTEIN,
                evidence_type.BENIGN,
                True,
                best_bs3.strength,
                f"Combined BS3 (score={best_bs3_score}): {best_bs3_comment}",
            )

        # No evidence found
        sources_tried = []
        if functional_db_path and functional_db_path.exists():
            sources_tried.append("functional_evidence.db")
        if splicevardb_path and splicevardb_path.exists():
            sources_tried.append("SpliceVarDB")
        if variant_literature:
            sources_tried.append("literature")

        return RuleResult(
            "PS3",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            False,
            evidence_strength.SUPPORTING,
            f"No functional evidence found in: {', '.join(sources_tried) if sources_tried else 'any source'}.",
        )

    @staticmethod
    def _score_to_strength(score: int) -> evidence_strength:
        """Convert numeric score to evidence_strength."""
        if score >= 8:
            return evidence_strength.VERY_STRONG
        elif score >= 4:
            return evidence_strength.STRONG
        elif score >= 2:
            return evidence_strength.MODERATE
        else:
            return evidence_strength.SUPPORTING

    @classmethod
    def _assess_splicevardb(
        cls,
        gene: str,
        hgvs: Optional[str],
        cDNA_pos: Optional[str],
        chrom: Optional[str],
        pos: Optional[int],
        transcript: Optional[str],
        splicevardb_path: pathlib.Path,
        is_classic_splice: bool = False,
    ) -> RuleResult:
        """Assess SpliceVarDB and return RuleResult.

        Args:
            gene: Gene symbol
            hgvs: HGVS cDNA notation
            cDNA_pos: cDNA position
            chrom: Chromosome
            pos: Genomic position
            transcript: Transcript ID
            splicevardb_path: Path to SpliceVarDB
            is_classic_splice: If True, variant is classic splice site (±1, ±2)

        For classic splice variants, PVS1 should be used instead of PS3
        when SpliceVarDB shows splice-altering effect.
        """
        try:
            db = SpliceVarDB(splicevardb_path)
            has_pathogenic, has_benign, reason = db.get_functional_evidence(
                gene=gene,
                hgvs=hgvs,
                cDNA_pos=cDNA_pos,
                chrom=chrom,
                pos=pos,
                transcript=transcript,
            )

            if has_pathogenic:
                # For classic splice variants, return PVS1 not PS3
                if is_classic_splice:
                    return RuleResult(
                        "PVS1",
                        rule_type.SPLICING,
                        evidence_type.PATHOGENIC,
                        True,
                        evidence_strength.VERY_STRONG,
                        f"SpliceVarDB: {reason}",
                    )
                return RuleResult(
                    "PS3",
                    rule_type.PROTEIN,
                    evidence_type.PATHOGENIC,
                    True,
                    evidence_strength.SUPPORTING,
                    f"SpliceVarDB: {reason}",
                )

            if has_benign:
                return RuleResult(
                    "BS3",
                    rule_type.PROTEIN,
                    evidence_type.BENIGN,
                    True,
                    evidence_strength.SUPPORTING,
                    f"SpliceVarDB: {reason}",
                )

            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                f"SpliceVarDB: {reason}",
            )

        except Exception as e:
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                f"SpliceVarDB error: {str(e)}",
            )

    @classmethod
    def _assess_functional_db(
        cls,
        transcript_id: str,
        cdna_pos: str,
        gene_name: str,
        position: int,
        chrom: str,
        functional_db_path: pathlib.Path,
        var_type: Optional[list] = None,
        ref_alt: Optional[str] = None,
    ) -> RuleResult:
        """Assess functional_evidence.db and return RuleResult."""
        try:
            from functional_evidence_db import FunctionalEvidenceDB

            # Check if variant is splice/intronic type
            is_splice_variant = False
            if var_type:
                from var_type import VARTYPE, VARTYPE_GROUPS
                splice_types = {
                    VARTYPE.SPLICE_ACCEPTOR,
                    VARTYPE.SPLICE_DONOR,
                    VARTYPE.SPLICE_ACCEPTOR_VARIANT,
                    VARTYPE.SPLICE_DONOR_VARIANT,
                }
                for vt in var_type:
                    if vt in splice_types or vt in VARTYPE_GROUPS.INTRONIC.value:
                        is_splice_variant = True
                        break

            with FunctionalEvidenceDB(functional_db_path) as db:
                result = None

                if is_splice_variant and ref_alt:
                    result = db.lookup_by_transcript_cdna_refalt(transcript_id, cdna_pos, ref_alt)
                    if result is None:
                        result = db.lookup_by_gene_cdna_refalt(gene_name, cdna_pos, ref_alt)
                else:
                    result = db.lookup_by_transcript_cdna(transcript_id, cdna_pos)
                    if result is None:
                        result = db.lookup_by_gene_cdna(gene_name, cdna_pos)
                    if result is None:
                        result = db.lookup_by_gene_position(gene_name, position, chrom)

                if result is None:
                    return RuleResult(
                        "PS3",
                        rule_type.PROTEIN,
                        evidence_type.PATHOGENIC,
                        False,
                        evidence_strength.SUPPORTING,
                        "Variant not found in functional evidence database.",
                    )

                if result.is_conflicting:
                    return RuleResult(
                        "PS3",
                        rule_type.PROTEIN,
                        evidence_type.PATHOGENIC,
                        False,
                        evidence_strength.SUPPORTING,
                        f"Conflicting functional evidence: {result.unique_publications}.",
                    )

                # Determine strength based on publication count
                pub_count = len(result.unique_publications.split('|')) if result.unique_publications else 1
                strength = evidence_strength.MODERATE if pub_count >= 6 else evidence_strength.SUPPORTING

                if result.is_pathogenic:
                    return RuleResult(
                        "PS3",
                        rule_type.PROTEIN,
                        evidence_type.PATHOGENIC,
                        True,
                        strength,
                        f"Functional evidence: {result.classification} ({pub_count} pub(s): {result.unique_publications}).",
                    )

                if result.is_benign:
                    return RuleResult(
                        "BS3",
                        rule_type.PROTEIN,
                        evidence_type.BENIGN,
                        True,
                        strength,
                        f"Functional evidence: {result.classification} ({pub_count} pub(s): {result.unique_publications}).",
                    )

                return RuleResult(
                    "PS3",
                    rule_type.PROTEIN,
                    evidence_type.PATHOGENIC,
                    False,
                    evidence_strength.SUPPORTING,
                    f"Unknown classification: {result.classification}",
                )

        except Exception as e:
            return RuleResult(
                "PS3",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.SUPPORTING,
                f"Functional evidence DB error: {str(e)}",
            )
