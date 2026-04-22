#!/usr/bin/env python3
"""
BP6: Benign variant from reputable source

BP6 Definition:
- If ClinGen variant review小组已定级为B(Benign)或LB(Likely Benign)，则使用BP6
- If ClinGen中仅有B或LB的提交记录，且超过2家实验室提交，则使用BP6

Evidence strength: Supporting (-1 point in Bayes scoring)
"""

from typing import Callable, Optional, List, Dict, Any

from acmg_rules.utils import (
    RuleResult,
    evidence_strength,
    abstract_rule,
    rule_type,
    evidence_type,
)
from information import Classification_Info, Info


@dataclass
class ClinVarSubmission:
    """Represents a single ClinVar submission."""
    scv: str  # Submission accessio (e.g., "SCV000123456")
    clinical_significance: str  # "Benign", "Likely benign", etc.
    lab: str  # Submitting laboratory
    review_status: str  # "criteria provided", "no assertion", etc.
    submission_date: Optional[str] = None


@dataclass
class ClinVarReviewDetails:
    """Extended ClinVar information for BP6 assessment."""
    variation_id: str
    rs_id: Optional[str]
    gene: Optional[str]

    # Review status from ClinGen/Expert Panel
    review_status: str  # "criteria provided, single submitter", "reviewed by expert panel", etc.
    clinical_significance: str  # Current clinical significance classification
    is_benign: bool  # True if B or LB
    is_likely_benign: bool  # True if LB

    # Submission details
    submissions: List[ClinVarSubmission]

    # Gene-specific classification (if from ClinGen Expert Panel)
    is_clingen_classified: bool
    clingen_review_status: Optional[str]  # e.g., "Expert Panel Review", "General Consentation"

    @property
    def benign_submissions(self) -> List[ClinVarSubmission]:
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
        """Check if BP6 evidence criteria are met."""
        # Criteria 1: ClinGen review小组已定级为B或LB
        if self.has_clingen_benign_classification:
            return True

        # Criteria 2: 仅有B/LB提交，且超过2家实验室
        if (len(self.benign_submissions) > 0 and
            self.unique_labs_with_benign >= 3 and  # >2 labs
            len([s for s in self.submissions
                 if s.clinical_significance.lower() not in {"benign", "likely benign"}]) == 0):
            return True

        return False


from dataclasses import dataclass


class Bp6_clingen(abstract_rule):
    """
    BP6: Reputable source (benign classification)

    Uses ClinVar data to determine if this variant has been classified as
    Benign (B) or Likely Benign (LB) by:

    1. ClinGen variant review小组 (Expert Panel review)
    2. Multiple independent lab submissions (>2 labs) with only B/LB classifications

    Evidence strength: Supporting
    """

    @classmethod
    def get_assess_rule(
        cls, class_info: Classification_Info
    ) -> tuple[Callable, tuple[Info, ...]]:
        return (
            cls.assess_rule,
            (
                class_info.VARIANT_CLINVAR_BP6,  # Extended ClinVar info for BP6
            ),
        )

    @classmethod
    def assess_rule(
        cls,
        clinvar_details: ClinVarReviewDetails,
    ) -> RuleResult:
        """
        Assess BP6 evidence based on ClinVar data.

        BP6 is applicable when:
        1. ClinGen/Expert Panel has classified variant as B or LB
        2. OR multiple labs (>2) have submitted only B/LB classifications
        """
        if clinvar_details is None:
            result = False
            comment = "No ClinVar data available for BP6 assessment."
            strength = evidence_strength.SUPPORTING

        elif not clinvar_details.meets_bp6_criteria:
            result = False
            comment = cls._build_no_bp6_comment(clinvar_details)
            strength = evidence_strength.SUPPORTING

        else:
            result = True

            # Determine if using ClinGen classification or multi-lab criteria
            if clinvar_details.has_clingen_benign_classification:
                comment = (
                    f"ClinGen variant review classified as {clinvar_details.clinical_significance} "
                    f"(Review status: {clinvar_details.review_status}). "
                    f"BP6 (Supporting) evidence applied."
                )
            else:
                benign_count = len(clinvar_details.benign_submissions)
                lab_count = clinvar_details.unique_labs_with_benign
                comment = (
                    f"{benign_count} B/LB submissions from {lab_count} independent laboratories. "
                    f"BP6 (Supporting) evidence applied."
                )

            strength = evidence_strength.SUPPORTING

        return RuleResult(
            "BP6",
            rule_type.CLINICAL_DATA,
            evidence_type.BENIGN,
            result,
            strength,
            comment,
        )

    @classmethod
    def _build_no_bp6_comment(cls, details: ClinVarReviewDetails) -> str:
        """Build comment explaining why BP6 is not applicable."""
        if details.submissions:
            benign_count = len(details.benign_submissions)
            lab_count = details.unique_labs_with_benign
            total_count = len(details.submissions)

            if benign_count == 0:
                return (
                    f"No benign/likely benign submissions found in ClinVar. "
                    f"BP6 not applicable."
                )
            elif lab_count <= 2:
                return (
                    f"Only {lab_count} lab(s) with B/LB submissions out of {total_count} total submissions. "
                    f"BP6 requires >2 independent labs for benign-only submissions."
                )
            else:
                # Has benign submissions but also has pathogenic ones
                return (
                    f"Found {benign_count} B/LB submissions from {lab_count} labs, but variant also has "
                    f"non-benign submissions in ClinVar. BP6 not applicable."
                )

        return (
            f"ClinVar review status: {details.review_status or 'unknown'}. "
            f"Clinical significance: {details.clinical_significance or 'unknown'}. "
            f"BP6 not applicable without benign classifications."
        )
