#!/usr/bin/env python3

import logging
import pandas as pd
from typing import Generator, Optional
from dataclasses import dataclass, field
from enum import Enum
from collections.abc import Iterable

import pyensembl

from variant import TranscriptInfo
from var_type import VARTYPE_GROUPS
from ensembl import ensembl
from custom_exceptions import No_transcript_with_var_type_found

logger = logging.getLogger("GenOtoScope_Classify.genotoscope_clinvar")

# path_clinvar = pathlib.Path("/home/katzkean/clinvar/clinvar_20230730.vcf.gz")


class ClinVar_Type(Enum):
    SAME_AA_CHANGE = "same_aa_change"
    DIFF_AA_CHANGE = "diff_aa_change"
    SAME_NUCLEOTIDE = "same_nucleotide"
    SAME_SPLICE_SITE = "same_splice_site"
    REGION = "region"


class ClinVar_Status(Enum):
    # Pathogenic statuses
    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely pathogenic"

    # Benign statuses
    BENIGN = "Benign"
    LIKELY_BENIGN = "Likely benign"

    # Uncertain significance
    UNCERTAIN_SIGNIFICANCE = "Uncertain significance"

    # Other
    NOT_PROVIDED = "not provided"
    CONFLICTING = "conflicting interpretations"


@dataclass
class ClinVar:
    pathogenic: bool
    type: ClinVar_Type
    highest_classification: Optional[ClinVar_Status]
    ids: list[str] = field(default_factory=list)
    associated_ids: list[str] = field(default_factory=list)
    # Extended fields for frontend output
    review_status: Optional[str] = None  # e.g., "criteria provided, single submitter"
    conditions: list[str] = field(default_factory=list)  # e.g., ["Autosomal recessive early-onset Parkinson disease 6"]
    submission_count: int = 0  # Total number of submissions
    pathogenic_submission_count: int = 0  # Number of pathogenic submissions
    is_clingen_reviewed: bool = False  # True if reviewed by ClinGen/Expert Panel


@dataclass
class ClinVarSubmission:
    """Represents a single ClinVar submission (SCV)."""
    scv: str  # Submission accession (e.g., "SCV000123456")
    clinical_significance: str  # "Benign", "Likely benign", etc.
    lab: str  # Submitting laboratory
    review_status: str  # "criteria provided", "no assertion", etc.
    submission_date: Optional[str] = None


@dataclass
class ClinVarReviewDetails:
    """
    Extended ClinVar information for BP6 assessment.

    Contains detailed information about ClinVar classifications,
    submissions, and review status.
    """
    variation_id: str
    rs_id: Optional[str]
    gene: Optional[str]

    # Overall review status
    review_status: str  # e.g., "criteria provided, single submitter", "reviewed by expert panel"
    clinical_significance: str  # Current clinical significance (from RCV/VCV)

    # Is it benign?
    is_benign: bool  # True if B or LB
    is_likely_benign: bool  # True if LB

    # Submission details
    submissions: list[ClinVarSubmission] = field(default_factory=list)

    # ClinGen specific
    is_clingen_classified: bool  # True if reviewed by ClinGen/Expert Panel
    clingen_review_type: Optional[str] = None  # e.g., "Expert Panel Review"

    @property
    def benign_submissions(self) -> list[ClinVarSubmission]:
        """Get submissions that are Benign or Likely Benign."""
        benign_terms = {"benign", "likely benign"}
        return [s for s in self.submissions
                if s.clinical_significance.lower() in benign_terms]

    @property
    def unique_labs_with_benign(self) -> int:
        """Count unique labs that submitted benign/likely benign classifications."""
        return len(set(s.lab for s in self.benign_submissions))

    @property
    def has_clingen_benign_classification(self) -> bool:
        """Check if ClinGen/Expert Panel has classified as B or LB."""
        return self.is_clingen_classified and self.is_benign

    @property
    def meets_bp6_criteria(self) -> bool:
        """
        Check if BP6 evidence criteria are met.

        BP6 applies when:
        1. ClinGen/Expert Panel has classified as B or LB, OR
        2. Multiple labs (>=3) have submitted only B/LB classifications
        """
        # Criteria 1: ClinGen review panel classification
        if self.has_clingen_benign_classification:
            return True

        # Criteria 2: Multiple labs with benign-only submissions
        benign_subs = self.benign_submissions
        if (len(benign_subs) >= 3 and  # At least 3 submissions
            self.unique_labs_with_benign >= 3 and  # From at least 3 labs
                len([s for s in self.submissions
                     if s.clinical_significance.lower() not in {"benign", "likely benign"}]) == 0):
            # All submissions are benign/likely benign
            return True

        return False


def get_affected_transcript(
    transcripts: Iterable[TranscriptInfo], var_types: VARTYPE_GROUPS
) -> tuple[TranscriptInfo, pyensembl.transcript.Transcript]:
    for transcript in transcripts:
        if any(var_type in transcript.var_type for var_type in var_types.value):
            try:
                ref_transcript = ensembl.transcript_by_id(transcript.transcript_id)
                ref_transcript.coding_sequence
            except ValueError or AttributeError:
                continue
            return transcript, ref_transcript
    raise No_transcript_with_var_type_found


def convert_vcf_gen_to_df(vcf_generator: Generator) -> pd.DataFrame:
    """
    Covnerts cyvcf generator into a pd.DataFrame
    """
    names = ["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info"]
    df = pd.DataFrame(columns=names)
    for entry in vcf_generator:
        clinvar_split_str = str(entry).split("\t")
        clinvar_dict = dict(zip(names, clinvar_split_str))
        df = pd.concat([df, pd.DataFrame([clinvar_dict])], axis=0, ignore_index=True)
    df_format_info = format_info(df)
    return df_format_info


def format_info(data: pd.DataFrame) -> pd.DataFrame:
    """
    Format the info column from ClinVar.vcf file as depicted in cyvcf
    """
    info = data["info"]
    info_split = [entry.split("=") for entry in info]
    info_split = [entry.split("\n")[0].split(";") for entry in info]
    processed_info = [[item.split("=") for item in entry] for entry in info_split]
    dict_info = [
        {item[0]: item[1] for item in entry if len(item) > 1}
        for entry in processed_info
    ]
    return pd.concat([data, pd.DataFrame(dict_info)], axis=1)


def create_ClinVar(clinvar: pd.DataFrame, type: ClinVar_Type) -> ClinVar:
    """
    From clinvar entries, get highest classification and IDs of ClinVar entries with that classification

    Checks for both pathogenic and benign classifications.
    """
    is_pathogenic = False
    highest_classification = None
    clinvar_ids = []
    clinvar_associated_ids = []

    if not clinvar.empty:
        # Check for pathogenic classifications
        if any(clinvar.CLNSIG == "Pathogenic"):
            is_pathogenic = True
            highest_classification = ClinVar_Status.PATHOGENIC
            clinvar_ids = list(clinvar[clinvar.CLNSIG == "Pathogenic"].id)
            clinvar_associated_ids = list(
                clinvar[clinvar.CLNSIG.str.contains("Likely_pathogenic")].id
            )
        elif any(clinvar.CLNSIG == "Likely_pathogenic") or any(
            clinvar.CLNSIG == "Pathogenic/Likely_pathogenic"
        ):
            is_pathogenic = True
            highest_classification = ClinVar_Status.LIKELY_PATHOGENIC
            clinvar_ids = list(
                clinvar[clinvar.CLNSIG.str.contains("Likely_pathogenic")].id
            )

    return ClinVar(
        is_pathogenic, type, highest_classification, clinvar_ids, clinvar_associated_ids
    )


def create_ClinVarReviewDetails(
    clinvar: pd.DataFrame,
    variation_id: str,
    rs_id: Optional[str] = None,
    gene: Optional[str] = None,
) -> ClinVarReviewDetails:
    """
    Create extended ClinVar details for BP6 assessment.

    Extracts review status, clinical significance, and submission details.
    """
    if clinvar.empty:
        return ClinVarReviewDetails(
            variation_id=variation_id,
            rs_id=rs_id,
            gene=gene,
            review_status="unknown",
            clinical_significance="unknown",
            is_benign=False,
            is_likely_benign=False,
            submissions=[],
            is_clingen_classified=False,
        )

    # Get primary classification from first row (VCV level)
    clinical_sig = clinvar.iloc[0].get("CLNSIG", "unknown")
    review_status = clinvar.iloc[0].get("CLNSTATUS", "unknown")

    # Determine if benign
    is_benign = clinical_sig == "Benign"
    is_likely_benign = clinical_sig == "Likely benign"

    # Check if ClinGen/Expert Panel classified this
    # (review_status contains "expert panel" or "clingen")
    clingen_keywords = ["expert panel", "clin gen", "gcep", "scg"]
    is_clingen = any(kw in str(review_status).lower() for kw in clingen_keywords)

    # Extract submissions if available
    # Note: Full submission details require separate API call
    # This extracts what's available in the VCF INFO field
    submissions = []

    # Try to extract submission accessions and labs from the data
    if "SUBS" in clinvar.columns:
        # Submissions are typically comma-separated in CLNACC
        accession_str = clinvar.iloc[0].get("CLNACC", "")
        if accession_str and isinstance(accession_str, str):
            # CLNACC format: "RCV000123456,RCV000123457"
            for i, acc in enumerate(accession_str.split(",")):
                if acc:
                    # Extract lab from CLNORIGIN if available, or use "unknown"
                    lab = clinvar.iloc[min(i, len(clinvar) - 1)].get("CLNORIGIN", "unknown")
                    submissions.append(ClinVarSubmission(
                        scv=acc.strip(),
                        clinical_significance=clinical_sig,
                        lab=str(lab),
                        review_status=review_status,
                    ))

    return ClinVarReviewDetails(
        variation_id=variation_id,
        rs_id=rs_id,
        gene=gene,
        review_status=str(review_status),
        clinical_significance=str(clinical_sig),
        is_benign=is_benign,
        is_likely_benign=is_likely_benign,
        submissions=submissions,
        is_clingen_classified=is_clingen,
    )


def filter_gene(clinvar: pd.DataFrame, gene: str) -> pd.DataFrame:
    """
    Filter out ClinVar entries that don't contain gene in GENEINFO
    """
    clinvar_filtered = clinvar[clinvar.GENEINFO.str.contains(gene)]
    return clinvar_filtered


def summarise_ClinVars(clinvars: list[ClinVar], type: ClinVar_Type) -> ClinVar:
    """
    Summarise a list of ClinVars into one ClinVar object
    """
    if len(clinvars) == 0:
        return ClinVar(pathogenic=False, highest_classification=None, ids=[], type=type)
    elif len(clinvars) == 1:
        return clinvars[0]
    else:
        pathogenic = False
        highest_classification = None
        ids = []
        for entry in clinvars:
            if entry.pathogenic:
                pathogenic = True
                ids.append(entry.ids)
                if not highest_classification == "Pathogenic":
                    highest_classification = entry.highest_classification
        return ClinVar(pathogenic, type, highest_classification, ids)
