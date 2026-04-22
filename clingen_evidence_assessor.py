#!/usr/bin/env python3
"""
ClinGen Evidence Assessor - integrates ClinGen gene-specific rules into ACMG classification.

This module provides:
1. Gene-specific evidence strength lookup from ClinGen_rules_rules
2. Combined PM1/BA1/BS1/PM2 assessment with ClinGen specifications
3. RAG-ready text chunks for LLM decision-making

Usage:
    assessor = ClinGenEvidenceAssessor(Path("clingen_rules.db"))
    rule = assessor.get_evidence_rule("BRCA1", "PP4")
    print(f"PP4 for BRCA1: {rule.strength}")
"""

import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple


@dataclass
class ClinGenEvidenceRule:
    """ClinGen gene-specific evidence rule."""
    gene: str
    rule_code: str
    strength: str
    strength_order: int
    description: str
    application_criteria: str


class ClinGenEvidenceAssessor:
    """
    Assessor for ClinGen gene-specific ACMG evidence rules.

    Provides methods to:
    - Get evidence strength for a gene/rule combination
    - Check if a variant meets ClinGen criteria
    - Generate RAG text chunks for LLM decision-making
    """

    def __init__(self, db_path: Path):
        self.db_path = db_path
        self._conn: Optional[sqlite3.Connection] = None

    def _get_connection(self) -> sqlite3.Connection:
        if self._conn is None:
            self._conn = sqlite3.connect(str(self.db_path))
            self._conn.row_factory = sqlite3.Row
        return self._conn

    def close(self):
        if self._conn:
            self._conn.close()
            self._conn = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def get_evidence_strength(
        self,
        gene: str,
        rule_code: str,
        default_strength: Optional[str] = None,
    ) -> Optional[str]:
        """
        Get ClinGen-specified evidence strength for a gene/rule.

        Args:
            gene: Gene symbol (e.g., "BRCA1")
            rule_code: ACMG rule code (e.g., "PM1", "PP4")
            default_strength: Fallback strength if no ClinGen rule found

        Returns:
            Evidence strength string or None
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        # Get the highest strength (most significant) rule for this gene/rule
        cursor.execute(
            """SELECT strength FROM clingen_rules
               WHERE gene = ? AND rule_code = ?
               ORDER BY strength_order DESC
               LIMIT 1""",
            (gene, rule_code)
        )

        row = cursor.fetchone()
        if row:
            return row["strength"]

        return default_strength

    def get_evidence_rule(
        self,
        gene: str,
        rule_code: str,
    ) -> Optional[ClinGenEvidenceRule]:
        """
        Get full ClinGen rule for a gene/rule combination.

        Args:
            gene: Gene symbol
            rule_code: ACMG rule code

        Returns:
            ClinGenEvidenceRule if found, None otherwise
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        cursor.execute(
            """SELECT * FROM clingen_rules
               WHERE gene = ? AND rule_code = ?
               ORDER BY strength_order DESC
               LIMIT 1""",
            (gene, rule_code)
        )

        row = cursor.fetchone()
        if row is None:
            return None

        return ClinGenEvidenceRule(
            gene=row["gene"],
            rule_code=row["rule_code"],
            strength=row["strength"],
            strength_order=row["strength_order"],
            description=row["description"] or "",
            application_criteria=row["application_criteria"] or "",
        )

    def get_all_rules_for_gene(self, gene: str) -> List[ClinGenEvidenceRule]:
        """
        Get all ClinGen rules for a gene.

        Args:
            gene: Gene symbol

        Returns:
            List of ClinGenEvidenceRule
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        cursor.execute(
            """SELECT * FROM clingen_rules
               WHERE gene = ?
               ORDER BY rule_code, strength_order DESC""",
            (gene,)
        )

        rules = []
        for row in cursor.fetchall():
            rules.append(ClinGenEvidenceRule(
                gene=row["gene"],
                rule_code=row["rule_code"],
                strength=row["strength"],
                strength_order=row["strength_order"],
                description=row["description"] or "",
                application_criteria=row["application_criteria"] or "",
            ))

        return rules

    def get_applicable_rules_for_gene(
        self,
        gene: str,
    ) -> Dict[str, ClinGenEvidenceRule]:
        """
        Get applicable rules for a gene, one per rule code (highest strength).

        Args:
            gene: Gene symbol

        Returns:
            Dict mapping rule_code to ClinGenEvidenceRule
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        cursor.execute(
            """SELECT DISTINCT rule_code FROM clingen_rules WHERE gene = ?""",
            (gene,)
        )
        rule_codes = [row["rule_code"] for row in cursor.fetchall()]

        result = {}
        for rc in rule_codes:
            rule = self.get_evidence_rule(gene, rc)
            if rule:
                result[rc] = rule

        return result

    def generate_rag_context(
        self,
        gene: str,
        rule_codes: Optional[List[str]] = None,
    ) -> str:
        """
        Generate RAG-ready text context for a gene.

        Args:
            gene: Gene symbol
            rule_codes: Optional list of specific rule codes to include

        Returns:
            Formatted text suitable for RAG retrieval
        """
        rules = self.get_all_rules_for_gene(gene)

        if rule_codes:
            rules = [r for r in rules if r.rule_code in rule_codes]

        if not rules:
            return f"No ClinGen rules available for gene {gene}."

        lines = [f"ClinGen ACMG Rules for {gene}", "=" * 50]

        # Group by rule code
        by_code: Dict[str, List[ClinGenEvidenceRule]] = {}
        for rule in rules:
            if rule.rule_code not in by_code:
                by_code[rule.rule_code] = []
            by_code[rule.rule_code].append(rule)

        for rc, code_rules in by_code.items():
            lines.append(f"\n{rc}:")
            for cr in code_rules:
                lines.append(f"  - Strength: {cr.strength}")
                if cr.description:
                    lines.append(f"    Description: {cr.description[:200]}...")
                if cr.application_criteria:
                    lines.append(f"    Criteria: {cr.application_criteria[:200]}...")

        return "\n".join(lines)

    def generate_variant_assessment_context(
        self,
        gene: str,
        protein_position: Optional[int] = None,
        frequency: Optional[float] = None,
        is_splice_variant: bool = False,
    ) -> str:
        """
        Generate context for variant assessment including relevant ClinGen rules.

        Args:
            gene: Gene symbol
            protein_position: Optional protein position for PM1 lookup
            frequency: Optional population frequency for BA1/BS1/PM2
            is_splice_variant: Whether variant is splice type

        Returns:
            Formatted text with relevant ClinGen rule information
        """
        lines = [f"ClinGen Evidence Assessment for {gene}", "=" * 50]

        # Get applicable rules
        applicable = self.get_applicable_rules_for_gene(gene)
        lines.append(f"\nApplicable ClinGen rules ({len(applicable)}):")

        for rc, rule in applicable.items():
            lines.append(f"  {rc}: {rule.strength}")

        # Add PM1 info if protein position provided
        if protein_position:
            lines.append(f"\nPM1 domain check for protein position {protein_position}:")
            cursor = self._get_connection().cursor()
            cursor.execute(
                """SELECT * FROM pm1_details WHERE gene = ?""",
                (gene,)
            )
            pm1_rows = cursor.fetchall()
            if pm1_rows:
                for row in pm1_rows:
                    lines.append(f"  - {row['strength']}: {row['protein_positions']}")
                    if row['domain']:
                        lines.append(f"    Domain: {row['domain']}")
            else:
                lines.append("  No ClinGen PM1 domains defined.")

        # Add frequency cutoff info if frequency provided
        if frequency is not None:
            lines.append(f"\nFrequency cutoff check for frequency={frequency}:")
            for rule_code in ["BA1", "BS1", "PM2"]:
                cutoff = self.get_frequency_cutoff(gene, rule_code)
                if cutoff:
                    lines.append(f"  {rule_code}: {cutoff['cutoff_operator']} {cutoff['cutoff_value']} ({cutoff['strength']})")

        return "\n".join(lines)

    def get_frequency_cutoff(
        self,
        gene: str,
        rule_code: str,
    ) -> Optional[Dict[str, Any]]:
        """Get frequency cutoff for a gene/rule."""
        conn = self._get_connection()
        cursor = conn.cursor()

        cursor.execute(
            """SELECT * FROM frequency_cutoffs
               WHERE gene = ? AND rule_code = ?
               ORDER BY CASE strength
                   WHEN 'Stand Alone' THEN 1
                   WHEN 'Strong' THEN 2
                   WHEN 'Supporting' THEN 3
                   ELSE 4
               END
               LIMIT 1""",
            (gene, rule_code)
        )

        row = cursor.fetchone()
        if row is None:
            return None

        return dict(row)


def get_clingen_evidence_summary(
    gene: str,
    db_path: Path,
) -> Dict[str, str]:
    """
    Convenience function to get evidence strength summary for a gene.

    Args:
        gene: Gene symbol
        db_path: Path to clingen_rules.db

    Returns:
        Dict mapping rule_code to evidence strength
    """
    with ClinGenEvidenceAssessor(db_path) as assessor:
        applicable = assessor.get_applicable_rules_for_gene(gene)
        return {rc: rule.strength for rc, rule in applicable.items()}
