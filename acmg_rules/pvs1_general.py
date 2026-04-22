#!/usr/bin/env python3
"""
General PVS1 rule for LOF variant assessment.

This PVS1 implementation is gene-agnostic and can be applied to any gene.
It determines PVS1 applicability based on:
1. Is the variant a NULL/LOF variant (frameshift, splice_site, start_loss, stop_gained)?
2. Does it cause NMD (nonsense-mediated decay)?
3. Is it in a critical protein region?

Evidence strength follows ACMG guidelines:
- NULL variant + NMD → PVS1 Very Strong
- NULL variant + critical region → PVS1 Strong
- NULL variant + non-critical region → PVS1 Moderate
"""

import pathlib
from typing import Callable, Optional

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    evidence_type,
    rule_type,
    abstract_rule,
    summarise_results_per_transcript,
)
from information import Classification_Info, Info
from variant import RNAData, TranscriptInfo, VariantInfo
from transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
    TranscriptInfo_exonic_inframe,
)
from var_type import VARTYPE_GROUPS
from null_variant_utils import (
    is_null_variant_by_hgvs,
    get_null_variant_description as get_null_variant_description_by_hgvs,
)
from acmg_rules.computation_evidence_utils import Threshold, assess_thresholds
from acmg_rules.functional_splicing_assay_utils import (
    adjust_strength_according_to_rna_data_pvs1,
)


def is_null_variant(variant: VariantInfo) -> bool:
    """
    Check if variant is a NULL/LOF variant.

    NULL variants include:
    - stop_gained (nonsense): p.Trp123*, p.Trp123X, p.Trp123Ter
    - frameshift: p.Gly123fs*12
    - classic splice site: c.321+1G>C, c.123-2G>A (±1, ±2)
    - start_lost: p.Met1?, p.Met1X
    - stop_lost: p.*123Gly

    判断优先级:
    1. 先尝试 HGVS 正则判断 (更精确)
    2. 如果 HGVS 不可用，回退到 VEP var_type 判断

    Args:
        variant: VariantInfo object

    Returns:
        True if variant is NULL/LOF type
    """
    # 方法1: 通过 HGVS 正则判断 (优先)
    if variant.hgvs_protein:
        if is_null_variant_by_hgvs(variant.hgvs_protein):
            return True

    # 方法2: 通过 VEP var_type 判断 (回退)
    if variant.var_type:
        null_types = set(VARTYPE_GROUPS.NULL_LOF.value)
        for var_type in variant.var_type:
            if var_type in null_types:
                return True
    return False


def get_null_variant_description(variant: VariantInfo) -> str:
    """
    Get human-readable description of NULL variant type.

    Args:
        variant: VariantInfo object

    Returns:
        Description string like "frameshift", "stop_gained", "splice site", etc.
    """
    # 方法1: 通过 HGVS 正则获取描述 (优先)
    if variant.hgvs_protein:
        hgvs_desc = get_null_variant_description_by_hgvs(variant.hgvs_protein)
        if hgvs_desc != "non-NULL":
            return hgvs_desc

    # 方法2: 通过 VEP var_type 获取描述 (回退)
    if variant.var_type:
        null_types = VARTYPE_GROUPS.NULL_LOF.value
        for var_type in variant.var_type:
            if var_type in null_types:
                return var_type.value.replace("_variant", "").replace("_", " ")
    return "unknown"


class Pvs1_general(abstract_rule):
    """
    General PVS1 rule applicable to any gene.

    This implements ACMG PVS1 without gene-specific configurations:
    - NULL variant + NMD → PVS1 Very Strong
    - NULL variant + critical region (disease relevant) → PVS1 Strong
    - NULL variant + not in critical region → PVS1 Moderate
    - Non-NULL variant → PVS1 does not apply
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.ANNOTATED_TRANSCRIPT_LIST,
                class_info.VARIANT,
                class_info.THRESHOLD_DIFF_LEN_PROT_PERCENT,
                class_info.SPLICING_ASSAY,
                class_info.VARIANT_PREDICTION,
                class_info.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC,
                class_info.MANE_TRANSCRIPT_LIST_PATH,
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        annotated_transcripts: list[TranscriptInfo],
        variant: VariantInfo,
        threshold_diff_len_prot_percent: float,
        splice_assay: Optional[list[RNAData]],
        prediction_dict: dict[str, float],
        threshold: Threshold,
        mane_path: pathlib.Path,
    ) -> RuleResult:
        """
        Assess PVS1 for general gene application.

        Args:
            annotated_transcripts: List of annotated transcripts
            variant: Variant information
            threshold_diff_len_prot_percent: Threshold for protein length change
            splice_assay: RNA splicing assay data
            prediction_dict: Splice prediction scores
            threshold: Threshold for splice prediction
            mane_path: Path to MANE transcript list

        Returns:
            RuleResult for PVS1
        """
        # First check if variant is a NULL/LOF variant
        if not is_null_variant(variant):
            return RuleResult(
                "PVS1",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.VERY_STRONG,
                f"PVS1 does not apply - variant type '{get_null_variant_description(variant)}' is not a NULL/LOF variant.",
            )

        # Process each transcript
        results = {}
        for transcript in annotated_transcripts:
            if isinstance(transcript, TranscriptInfo_exonic):
                result = cls.assess_pvs1_frameshift_general(
                    transcript, threshold_diff_len_prot_percent
                )
                results[transcript.transcript_id] = result
            elif isinstance(transcript, TranscriptInfo_intronic):
                result = cls.assess_pvs1_splice_general(
                    transcript,
                    prediction_dict,
                    threshold,
                    threshold_diff_len_prot_percent,
                )
                if splice_assay:
                    result = adjust_strength_according_to_rna_data_pvs1(
                        splice_assay, result
                    )
                results[transcript.transcript_id] = result
            elif isinstance(transcript, TranscriptInfo_start_loss):
                result = cls.assess_pvs1_start_loss_general(transcript)
                results[transcript.transcript_id] = result

        if len(results) == 0:
            return RuleResult(
                "PVS1",
                rule_type.GENERAL,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.VERY_STRONG,
                f"PVS1 does not apply to this variant, as PVS1 does not apply to variant types {', '.join([v.value for v in variant.var_type])}.",
            )

        return summarise_results_per_transcript(results, "PVS1", mane_path)

    @classmethod
    def assess_pvs1_frameshift_general(
        cls,
        transcript: TranscriptInfo_exonic,
        threshold_diff_len_prot_percent: float,
    ) -> RuleResult:
        """
        Assess PVS1 for frameshift/null variants.

        Logic:
        - NMD + critical region → Very Strong
        - NMD only → Very Strong
        - Not NMD + critical region → Strong
        - Not NMD + not critical region → Moderate
        """
        null_type = get_null_variant_description_from_transcript(transcript)

        if transcript.is_NMD:
            # NMD applies - Very Strong evidence
            if transcript.is_truncated_region_disease_relevant:
                comment = (
                    f"Transcript {transcript.transcript_id}: {null_type} variant causes NMD. "
                    f"Truncated region is disease relevant. "
                    f"{transcript.comment_truncated_region}"
                )
            else:
                comment = (
                    f"Transcript {transcript.transcript_id}: {null_type} variant causes NMD. "
                    f"Truncated region is not disease relevant."
                )
            return RuleResult(
                "PVS1",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                True,
                evidence_strength.VERY_STRONG,
                comment,
            )
        else:
            # Not NMD - check critical region
            if transcript.is_truncated_region_disease_relevant:
                comment = (
                    f"Transcript {transcript.transcript_id}: {null_type} variant does not cause NMD. "
                    f"Truncated region is disease relevant. "
                    f"{transcript.comment_truncated_region}"
                )
                strength = evidence_strength.STRONG
            else:
                # Check protein length change
                if transcript.diff_len_protein_percent > threshold_diff_len_prot_percent:
                    comment = (
                        f"Transcript {transcript.transcript_id}: {null_type} variant does not cause NMD. "
                        f"Protein length change: {transcript.diff_len_protein_percent}% (>{threshold_diff_len_prot_percent}% threshold)."
                    )
                    strength = evidence_strength.MODERATE
                else:
                    comment = (
                        f"Transcript {transcript.transcript_id}: {null_type} variant does not cause NMD. "
                        f"Protein length change: {transcript.diff_len_protein_percent}% (≤{threshold_diff_len_prot_percent}% threshold)."
                    )
                    strength = evidence_strength.MODERATE

            return RuleResult(
                "PVS1",
                rule_type.PROTEIN,
                evidence_type.PATHOGENIC,
                True,
                strength,
                comment,
            )

    @classmethod
    def assess_pvs1_splice_general(
        cls,
        transcript: TranscriptInfo_intronic,
        prediction_dict: dict[str, float],
        threshold: Threshold,
        threshold_diff_len_prot_percent: float,
    ) -> RuleResult:
        """
        Assess PVS1 for splice variants.

        Logic:
        - Predicted splicing alteration + NMD → Very Strong
        - Predicted splicing alteration + critical region → Strong
        - Predicted splicing alteration + not critical region → Moderate
        - No splicing alteration predicted → Does not apply
        """
        prediction_value = prediction_dict.get(threshold.name, None)
        num_thresholds_met = assess_thresholds(threshold, prediction_value)

        if not transcript.are_exons_skipped or not num_thresholds_met:
            return RuleResult(
                "PVS1",
                rule_type.SPLICING,
                evidence_type.PATHOGENIC,
                False,
                evidence_strength.VERY_STRONG,
                f"Transcript {transcript.transcript_id}: No splicing alteration predicted.",
            )

        if transcript.is_NMD:
            comment = (
                f"Transcript {transcript.transcript_id}: Splicing alteration causes NMD. "
                f"All exons are disease relevant."
            )
            return RuleResult(
                "PVS1",
                rule_type.SPLICING,
                evidence_type.PATHOGENIC,
                True,
                evidence_strength.VERY_STRONG,
                comment,
            )

        if transcript.is_truncated_region_disease_relevant:
            comment = (
                f"Transcript {transcript.transcript_id}: Splicing alteration does not cause NMD. "
                f"Skipped exon is disease relevant. "
                f"{transcript.comment_truncated_region}"
            )
            strength = evidence_strength.STRONG
        else:
            if not transcript.is_reading_frame_preserved:
                if transcript.diff_len_protein_percent > threshold_diff_len_prot_percent:
                    comment = (
                        f"Transcript {transcript.transcript_id}: Splicing alteration does not cause NMD. "
                        f"Reading frame not preserved. Protein length change: {transcript.diff_len_protein_percent}%."
                    )
                    strength = evidence_strength.MODERATE
                else:
                    comment = (
                        f"Transcript {transcript.transcript_id}: Splicing alteration does not cause NMD. "
                        f"Reading frame not preserved but protein length change minimal."
                    )
                    strength = evidence_strength.MODERATE
            else:
                comment = (
                    f"Transcript {transcript.transcript_id}: Splicing alteration does not cause NMD. "
                    f"Reading frame preserved."
                )
                strength = evidence_strength.MODERATE

        return RuleResult(
            "PVS1",
            rule_type.SPLICING,
            evidence_type.PATHOGENIC,
            True,
            strength,
            comment,
        )

    @classmethod
    def assess_pvs1_start_loss_general(
        cls,
        transcript: TranscriptInfo_start_loss,
    ) -> RuleResult:
        """
        Assess PVS1 for start loss variants.

        Original ACMG logic:
        - No alternative start codon → Moderate
        - Alternative start codon + critical region → Moderate
        - Alternative start codon + not critical region → Supporting
        """
        if not transcript.exists_alternative_start_codon:
            comment = (
                f"Transcript {transcript.transcript_id}: No alternative start codons were detected in transcript."
            )
            result = True
            strength = evidence_strength.MODERATE
        else:
            comment = f"Transcript {transcript.transcript_id}: Alternative start codon observed."
            if transcript.is_truncated_region_disease_relevant:
                comment = (
                    comment
                    + f" Alternative start codon leads to the exclusion of a disease relevant protein region. "
                    + transcript.comment_truncated_region
                )
                result = True
                strength = evidence_strength.MODERATE
            else:
                comment = (
                    comment
                    + f" No pathogenic variant detected between start codon and alternative start codon."
                )
                result = True
                strength = evidence_strength.SUPPORTING

        return RuleResult(
            "PVS1",
            rule_type.PROTEIN,
            evidence_type.PATHOGENIC,
            result,
            strength,
            comment,
        )


def get_null_variant_description_from_transcript(transcript: TranscriptInfo) -> str:
    """
    Get description of NULL variant type from transcript.

    判断优先级:
    1. 先尝试 HGVS 正则判断 (更精确)
    2. 如果 HGVS 不可用，回退到 VEP var_type 判断
    """
    # 方法1: 通过 HGVS 正则获取描述 (优先)
    if transcript.var_protein:
        hgvs_desc = get_null_variant_description_by_hgvs(transcript.var_protein)
        if hgvs_desc != "non-NULL":
            return hgvs_desc

    # 方法2: 通过 VEP var_type 获取描述 (回退)
    if transcript.var_type:
        null_types = VARTYPE_GROUPS.NULL_LOF.value
        for var_type in transcript.var_type:
            if var_type in null_types:
                return var_type.value.replace("_variant", "").replace("_", " ")
    return "unknown"
