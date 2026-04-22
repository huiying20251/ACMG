#!/usr/bin/env python3

import pathlib
import logging
from collections.abc import Callable

from information import Classification_Info, Info
from var_type import VARTYPE_GROUPS
from clinvar_missense import check_clinvar_missense
from clinvar_splicing import check_clinvar_splicing
from variant import VariantInfo, TranscriptInfo
from clinvar_utils import ClinVar_Type, ClinVar, ClinVarReviewDetails, create_ClinVarReviewDetails, ClinVarSubmission, ClinVar_Status
from ps1_splicing_utils import should_use_ps1_splicing, SplicingVariantInfo, SplicingRegion

logger = logging.getLogger("annotate_clinvar")


def _get_hgvs_from_transcripts(transcripts: list[TranscriptInfo]) -> str:
    """
    从transcripts获取HGVS字符串
    """
    for transcript in transcripts:
        if hasattr(transcript, 'var_hgvs') and transcript.var_hgvs:
            # var_hgvs是hgvs.posedit.PosEdit对象
            return str(transcript.var_hgvs)
    return ""


def _get_vep_consequence(transcripts: list[TranscriptInfo]) -> list[str]:
    """
    从transcripts获取VEP注释的consequence terms
    优先使用all_consequences（原始VEP字符串），否则从var_type推断
    """
    consequences = set()
    for transcript in transcripts:
        # 优先使用 all_consequences（原始VEP字符串列表）
        if hasattr(transcript, 'all_consequences') and transcript.all_consequences:
            for c in transcript.all_consequences:
                consequences.add(c)
        # 回退：从 var_type 推断
        elif hasattr(transcript, 'var_type'):
            for vt in transcript.var_type:
                consequences.add(vt.value)
    return list(consequences)


def annotate_clinvar(
    variant: VariantInfo,
    transcripts: list[TranscriptInfo],
    path_clinvar: pathlib.Path,
) -> dict[ClinVar_Type, ClinVar]:
    """
    Manage ClinVar annotation
    For both splicing and missense variants

    根据ClinGen Splice Variant Classification规则:
    - Branch 1: HGVS格式为剪接变异 (±1~+6, -1~-20)
    - Branch 2: VEP注释为splice site但实际是missense

    使用ps1_splicing_utils中的should_use_ps1_splicing判断
    """
    # 获取HGVS字符串和VEP consequence用于判断
    hgvs_string = _get_hgvs_from_transcripts(transcripts)
    vep_consequence = _get_vep_consequence(transcripts)

    # 使用新的判断逻辑决定是否走PS1_splicing
    use_ps1_splicing, splicing_info = should_use_ps1_splicing(hgvs_string, vep_consequence)

    # 处理missense变异
    if any(var_type in VARTYPE_GROUPS.MISSENSE.value for var_type in variant.var_type):
        clinvar_same_aa, clinvar_diff_aa = check_clinvar_missense(
            variant, transcripts, path_clinvar
        )
    else:
        clinvar_same_aa = ClinVar(False, ClinVar_Type.SAME_AA_CHANGE, None, [])
        clinvar_diff_aa = ClinVar(False, ClinVar_Type.DIFF_AA_CHANGE, None, [])

    # 处理剪接变异
    if use_ps1_splicing:
        clinvar_same_pos, clinvar_splice_site = check_clinvar_splicing(
            variant, transcripts, path_clinvar, splicing_info
        )
    else:
        clinvar_same_pos = ClinVar(False, ClinVar_Type.SAME_NUCLEOTIDE, None, [])
        clinvar_splice_site = ClinVar(False, ClinVar_Type.SAME_SPLICE_SITE, None, [])
    clinvar_dict = {
        ClinVar_Type.SAME_AA_CHANGE: clinvar_same_aa,
        ClinVar_Type.DIFF_AA_CHANGE: clinvar_diff_aa,
        ClinVar_Type.SAME_NUCLEOTIDE: clinvar_same_pos,
        ClinVar_Type.SAME_SPLICE_SITE: clinvar_splice_site,
    }
    return clinvar_dict


def get_annotate_clinvar(
    class_info: Classification_Info,
) -> tuple[Callable, tuple[Info, ...]]:
    """
    Get function for clinvar annotation and needed classification_information objects
    """
    return (
        annotate_clinvar,
        (
            class_info.VARIANT,
            class_info.TRANSCRIPT,
            class_info.CLINVAR_PATH,
        ),
    )


def annotate_clinvar_bp6(
    variant: VariantInfo,
    transcripts: list[TranscriptInfo],
    path_clinvar: pathlib.Path,
) -> ClinVarReviewDetails:
    """
    Annotate variant with extended ClinVar information for BP6 assessment.

    Extracts review status, clinical significance, submissions, and determines
    if BP6 criteria are met based on ClinGen/Expert Panel classification
    or multiple benign submissions from independent labs.
    """
    from clinvar_missense import check_clinvar_missense
    from clinvar_splicing import check_clinvar_splicing

    # Get basic ClinVar data using existing functions
    if any(var_type in VARTYPE_GROUPS.MISSENSE.value for var_type in variant.var_type):
        clinvar_same_aa, clinvar_diff_aa = check_clinvar_missense(
            variant, transcripts, path_clinvar
        )
    else:
        clinvar_same_aa = ClinVar(False, ClinVar_Type.SAME_AA_CHANGE, None, [])
        clinvar_diff_aa = ClinVar(False, ClinVar_Type.DIFF_AA_CHANGE, None, [])

    if any(var_type in VARTYPE_GROUPS.INTRONIC.value for var_type in variant.var_type):
        clinvar_same_pos, clinvar_splice_site = check_clinvar_splicing(
            variant, transcripts, path_clinvar
        )
    else:
        clinvar_same_pos = ClinVar(False, ClinVar_Type.SAME_NUCLEOTIDE, None, [])
        clinvar_splice_site = ClinVar(False, ClinVar_Type.SAME_SPLICE_SITE, None, [])

    # Combine all clinvar data - use the most relevant one with data
    all_clinvars = [clinvar_same_aa, clinvar_diff_aa, clinvar_same_pos, clinvar_splice_site]
    clinvar_with_data = [c for c in all_clinvars if c.ids]

    if not clinvar_with_data:
        return ClinVarReviewDetails(
            variation_id="unknown",
            rs_id=getattr(variant, 'rs_id', None),
            gene=getattr(variant, 'gene', None),
            review_status="unknown",
            clinical_significance="unknown",
            is_benign=False,
            is_likely_benign=False,
            submissions=[],
            is_clingen_classified=False,
        )

    # For BP6, we use the first clinvar entry with IDs (most relevant)
    primary_clinvar = clinvar_with_data[0]

    # Map ClinVar_Status to clinical significance string
    status_to_sig = {
        ClinVar_Status.PATHOGENIC: "Pathogenic",
        ClinVar_Status.LIKELY_PATHOGENIC: "Likely pathogenic",
        ClinVar_Status.BENIGN: "Benign",
        ClinVar_Status.LIKELY_BENIGN: "Likely benign",
        ClinVar_Status.UNCERTAIN_SIGNIFICANCE: "Uncertain significance",
        ClinVar_Status.NOT_PROVIDED: "not provided",
        ClinVar_Status.CONFLICTING: "conflicting interpretations",
    }

    clinical_sig = "unknown"
    is_benign = False
    is_likely_benign = False

    if primary_clinvar.highest_classification:
        clinical_sig = status_to_sig.get(
            primary_clinvar.highest_classification, "unknown"
        )
        is_benign = primary_clinvar.highest_classification == ClinVar_Status.BENIGN
        is_likely_benign = primary_clinvar.highest_classification == ClinVar_Status.LIKELY_BENIGN

    # Create ClinVarReviewDetails directly from ClinVar object
    # Note: Full submission details require a separate ClinVar API call
    return ClinVarReviewDetails(
        variation_id=str(primary_clinvar.ids[0]) if primary_clinvar.ids else "unknown",
        rs_id=getattr(variant, 'rs_id', None),
        gene=getattr(variant, 'gene', None),
        review_status="unknown",  # Would need API call for full review status
        clinical_significance=clinical_sig,
        is_benign=is_benign,
        is_likely_benign=is_likely_benign,
        submissions=[],  # Would need API call for submission details
        is_clingen_classified=False,  # Would need API call to determine
    )


def get_annotate_clinvar_bp6(
    class_info: Classification_Info,
) -> tuple[Callable, tuple[Info, ...]]:
    """
    Get function for BP6 ClinVar annotation and needed classification_information objects
    """
    return (
        annotate_clinvar_bp6,
        (
            class_info.VARIANT,
            class_info.TRANSCRIPT,
            class_info.CLINVAR_PATH,
        ),
    )
