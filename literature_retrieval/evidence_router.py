#!/usr/bin/env python3
"""
Evidence Router

Routes literature evidence to appropriate ACMG rules based on inheritance pattern.

Autosomal Dominant → PS2, PP1, PS4, PP4 experts
Autosomal Recessive → PM3, PP1, PP4 experts
"""

import logging
from typing import List, Dict, Any, Optional

from .literature_utils import (
    Article,
    CaseReport,
    VariantLiterature,
    EvidenceResult,
    InheritancePattern,
    IndividualVariantObservation,
)

logger = logging.getLogger(__name__)


# Evidence expert prompts for different rules
EVIDENCE_PROMPTS = {
    "PS2": """Evaluate if the following case represents a confirmed de novo variant:
Variant: {variant}
Case: {case_description}

Consider:
- Was the variant confirmed to be de novo through parental testing?
- Is the variant absent from both parents?
- Is this a rare variant in a disease-relevant gene?

If confirmed de novo, this supports PS2 (strong pathogenic evidence).""",

    "PS4": """Evaluate case-control evidence for the following variant:
Variant: {variant}
Case: {case_description}

Consider:
- Number of affected cases vs controls
- Statistical significance (p-value, odds ratio)
- Functional studies if available
- Allele frequency in cases vs controls

If significant enrichment in cases, this supports PS4 (strong pathogenic evidence).""",

    "PP1": """Evaluate cosegregation evidence for the following variant:
Variant: {variant}
Family data: {segregation_data}

Consider:
- Number of affected family members with the variant
- Number of unaffected family members without the variant
- LOD score if reported
- Inheritance pattern consistency

If variant cosegregates with disease, this supports PP1 (supporting pathogenic evidence).""",

    "PP4": """Evaluate phenotype-specific evidence for the following variant:
Variant: {variant}
Phenotype: {phenotype}

Consider:
- Does the phenotype match the known disease phenotype?
- Are there HPO terms that support the association?
- Is the variant in a gene known to cause this specific phenotype?

If phenotype matches, this supports PP4 (supporting pathogenic evidence).""",

    "PM3": """Evaluate for autosomal recessive inheritance:
Variant: {variant}
Case: {case_description}

Consider:
- Is this a compound heterozygous case?
- Are there two different pathogenic variants on opposite chromosomes?
- Is there a known pathogenic variant in trans?

If confirmed compound heterozygous with known pathogenic variant, this supports PM3.""",
}


class EvidenceRouter:
    """
    Routes literature evidence to appropriate ACMG rules based on inheritance pattern.

    Evidence Routing:
    - Autosomal Dominant → PS2, PP1, PS4, PP4
    - Autosomal Recessive → PM3, PP1, PP4
    - X-linked → Similar to AD/AR with gender considerations
    """

    # Rules applicable per inheritance pattern
    INHERITANCE_RULES = {
        InheritancePattern.AUTOSOMAL_DOMINANT: ["PS2", "PS4", "PP1", "PP4"],
        InheritancePattern.AUTOSOMAL_RECESSIVE: ["PM3", "PP1", "PP4"],
        InheritancePattern.X_LINKED: ["PS2", "PP1", "PP4", "PM3"],
        InheritancePattern.UNKNOWN: ["PP1", "PP4"],  # Only rules that work without knowing inheritance
    }

    def __init__(self, llm_api_key: Optional[str] = None):
        self.llm_api_key = llm_api_key
        self.evidence_prompts = EVIDENCE_PROMPTS

    def route_evidence(
        self,
        literature: VariantLiterature,
        inheritance_pattern: InheritancePattern = InheritancePattern.UNKNOWN,
    ) -> List[EvidenceResult]:
        """
        Route literature evidence to applicable ACMG rules.

        Args:
            literature: Aggregated literature data for the variant
            inheritance_pattern: Inheritance pattern of the disease

        Returns:
            List of EvidenceResult objects for applicable rules
        """
        applicable_rules = self.INHERITANCE_RULES.get(
            inheritance_pattern,
            self.INHERITANCE_RULES[InheritancePattern.UNKNOWN]
        )

        logger.info(f"Routing evidence for {literature.gene} with {inheritance_pattern.value} pattern")
        logger.info(f"Applicable rules: {applicable_rules}")

        results = []

        # Process case reports
        for case in literature.case_reports:
            # Check de novo status (PS2)
            if "PS2" in applicable_rules and case.is_de_novo:
                result = self._assess_ps2(case, literature)
                if result:
                    results.append(result)

            # Check case-control evidence (PS4)
            if "PS4" in applicable_rules and case.num_controls is not None:
                result = self._assess_ps4(case, literature)
                if result:
                    results.append(result)

            # Check cosegregation (PP1)
            if "PP1" in applicable_rules and case.segregation_info:
                result = self._assess_pp1(case, literature)
                if result:
                    results.append(result)

            # Check phenotype match (PP4)
            if "PP4" in applicable_rules and (case.phenotype or case.hpo_terms):
                result = self._assess_pp4(case, literature)
                if result:
                    results.append(result)

            # Check recessive inheritance (PM3)
            if "PM3" in applicable_rules:
                result = self._assess_pm3(case, literature)
                if result:
                    results.append(result)

        return results

    def _assess_ps2(self, case: CaseReport, literature: VariantLiterature) -> Optional[EvidenceResult]:
        """Assess PS2 (de novo) evidence."""
        if not case.confirmed_de_novo:
            return None

        # Determine strength based on evidence
        strength = "STRONG" if case.confirmed_de_novo else "MODERATE"

        return EvidenceResult(
            rule="PS2",
            applicable=True,
            strength=strength,
            num_cases=1,
            pmids=[case.pmid],
            confidence="high" if case.confirmed_de_novo else "medium",
            comment=f"Confirmed de novo variant in {case.gene} (PMID: {case.pmid}). "
                    f"Inheritance: {case.inheritance or 'not specified'}.",
        )

    def _assess_ps4(self, case: CaseReport, literature: VariantLiterature) -> Optional[EvidenceResult]:
        """Assess PS4 (case-control) evidence."""
        if case.num_controls is None or case.num_controls == 0:
            return None

        # Determine strength based on statistical significance and case count
        strength = "STRONG"
        if case.p_value:
            if case.p_value > 0.001:
                strength = "MODERATE"
            if case.p_value > 0.01:
                strength = "SUPPORTING"

        return EvidenceResult(
            rule="PS4",
            applicable=True,
            strength=strength,
            num_cases=case.num_cases,
            num_controls=case.num_controls,
            pmids=[case.pmid],
            confidence="high" if case.p_value and case.p_value < 0.001 else "medium",
            comment=f"Case-control evidence: {case.num_cases} cases vs {case.num_controls} controls. "
                    f"OR: {case.odds_ratio}, p-value: {case.p_value}. (PMID: {case.pmid})",
        )

    def _assess_pp1(self, case: CaseReport, literature: VariantLiterature) -> Optional[EvidenceResult]:
        """Assess PP1 (cosegregation) evidence."""
        if not case.segregation_info:
            return None

        families = case.segregation_info.get("families", [])
        if not families:
            return None

        # Count total segregation
        total_tested = sum(f.get("tested", 0) for f in families)
        total_affected = sum(f.get("affected", 0) for f in families)
        total_unaffected = sum(f.get("unaffected", 0) for f in families)

        # Determine strength based on number of informative meioses
        # PP1 can be supporting, moderate, or strong based on LOD score
        strength = "SUPPORTING"
        if total_tested >= 10:
            strength = "MODERATE"
        if total_tested >= 20:
            strength = "STRONG"

        return EvidenceResult(
            rule="PP1",
            applicable=True,
            strength=strength,
            segregation_data=case.segregation_info,
            pmids=[case.pmid],
            confidence="high" if total_tested >= 10 else "medium",
            comment=f"Cosegregation: {total_affected} affected with variant, "
                    f"{total_unaffected} unaffected without variant. "
                    f"Total tested: {total_tested}. (PMID: {case.pmid})",
        )

    def _assess_pp4(self, case: CaseReport, literature: VariantLiterature) -> Optional[EvidenceResult]:
        """Assess PP4 (phenotype match) evidence."""
        if not case.phenotype and not case.hpo_terms:
            return None

        phenotype_str = case.phenotype or ", ".join(case.hpo_terms[:5])

        return EvidenceResult(
            rule="PP4",
            applicable=True,
            strength="SUPPORTING",
            phenotype=case.phenotype,
            pmids=[case.pmid],
            confidence="medium",
            comment=f"Phenotype: {phenotype_str}. "
                    f"HPO terms: {', '.join(case.hpo_terms[:3])}. (PMID: {case.pmid})",
        )

    def _assess_pm3(self, case: CaseReport, literature: VariantLiterature) -> Optional[EvidenceResult]:
        """Assess PM3 (recessive - compound heterozygous) evidence."""
        # This requires detailed compound het analysis
        # For now, return a basic assessment

        if case.inheritance_pattern != InheritancePattern.AUTOSOMAL_RECESSIVE:
            return None

        # PM3 typically requires knowing the other variant
        # This is a simplified assessment
        return EvidenceResult(
            rule="PM3",
            applicable=True,
            strength="MODERATE",
            pmids=[case.pmid],
            confidence="low",  # Requires more detailed analysis
            comment=f"Autosomal recessive inheritance reported. "
                    f"Detailed compound het analysis needed. (PMID: {case.pmid})",
        )

    def get_evidence_summary(
        self,
        results: List[EvidenceResult],
    ) -> Dict[str, Any]:
        """
        Generate a summary of evidence results.

        Args:
            results: List of EvidenceResult objects

        Returns:
            Summary dict with evidence counts and details
        """
        summary = {
            "total_evidence": len(results),
            "rules_triggered": {},
            "high_confidence": [],
            "medium_confidence": [],
            "low_confidence": [],
        }

        for result in results:
            # Count by rule
            if result.rule not in summary["rules_triggered"]:
                summary["rules_triggered"][result.rule] = []
            summary["rules_triggered"][result.rule].append(result.to_rule_result())

            # Categorize by confidence
            if result.confidence == "high":
                summary["high_confidence"].append(result.rule)
            elif result.confidence == "medium":
                summary["medium_confidence"].append(result.rule)
            else:
                summary["low_confidence"].append(result.rule)

        return summary

    def route_individual_observations(
        self,
        observations: List[IndividualVariantObservation],
        inheritance_pattern: InheritancePattern = InheritancePattern.UNKNOWN,
    ) -> List[EvidenceResult]:
        """
        Route individual-level observations to appropriate ACMG rules.

        Args:
            observations: List of IndividualVariantObservation objects
            inheritance_pattern: Inheritance pattern of the disease

        Returns:
            List of EvidenceResult objects for applicable rules
        """
        applicable_rules = self.INHERITANCE_RULES.get(
            inheritance_pattern,
            self.INHERITANCE_RULES[InheritancePattern.UNKNOWN]
        )

        logger.info(f"Routing {len(observations)} individual observations for {inheritance_pattern.value} pattern")

        results = []

        for obs in observations:
            # PS2: de novo assessment
            if "PS2" in applicable_rules and obs.is_de_novo:
                strength = "STRONG" if obs.is_confirmed_de_novo else "MODERATE"
                results.append(EvidenceResult(
                    rule="PS2",
                    applicable=True,
                    strength=strength,
                    num_cases=1,
                    pmids=[obs.pmid],
                    confidence="high" if obs.is_confirmed_de_novo else "medium",
                    comment=f"De novo variant in {obs.gene} individual {obs.individual_id}. "
                            f"Zygosity: {obs.zygosity}. "
                            f"ClinVar: {obs.clinvar_status or 'Not specified'}. "
                            f"PMID: {obs.pmid}",
                ))

            # PM3: biallelic (homozygous or compound heterozygous)
            if "PM3" in applicable_rules and obs.zygosity in ("homozygous", "compound_heterozygous"):
                results.append(EvidenceResult(
                    rule="PM3",
                    applicable=True,
                    strength="MODERATE",
                    num_cases=1,
                    pmids=[obs.pmid],
                    confidence="medium",
                    comment=f"Biallelic inheritance in {obs.gene} individual {obs.individual_id}. "
                            f"Zygosity: {obs.zygosity}. "
                            f"ClinVar: {obs.clinvar_status or 'Not specified'}. "
                            f"PMID: {obs.pmid}",
                ))

            # PP1: cosegregation
            if "PP1" in applicable_rules and obs.segregation_data:
                families = obs.segregation_data.get("families", [])
                total_tested = sum(f.get("total_tested", 0) for f in families)
                strength = "SUPPORTING"
                if total_tested >= 10:
                    strength = "MODERATE"
                if total_tested >= 20:
                    strength = "STRONG"

                results.append(EvidenceResult(
                    rule="PP1",
                    applicable=True,
                    strength=strength,
                    segregation_data=obs.segregation_data,
                    pmids=[obs.pmid],
                    confidence="high" if total_tested >= 10 else "medium",
                    comment=f"Cosegregation in {obs.gene} individual {obs.individual_id}. "
                            f"Families: {len(families)}. "
                            f"PMID: {obs.pmid}",
                ))

            # PP4: phenotype match
            if "PP4" in applicable_rules and (obs.phenotype or obs.hpo_terms):
                phenotype_str = obs.phenotype or ", ".join(obs.hpo_terms[:5])
                results.append(EvidenceResult(
                    rule="PP4",
                    applicable=True,
                    strength="SUPPORTING",
                    phenotype=obs.phenotype,
                    pmids=[obs.pmid],
                    confidence="medium",
                    comment=f"Phenotype match for {obs.gene} individual {obs.individual_id}: {phenotype_str}. "
                            f"HPO: {', '.join(obs.hpo_terms[:3])}. "
                            f"ClinVar: {obs.clinvar_status or 'Not specified'}. "
                            f"PMID: {obs.pmid}",
                ))

        return results
