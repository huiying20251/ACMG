#!/usr/bin/env python3
"""
Query module for ClinGen rules database (clingen_rules.db).

Provides lookup functions for:
- PM1 domain details
- Frequency cutoffs for BA1, BS1, PM2
- Gene-specific ACMG rules
"""

import sqlite3
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Tuple


@dataclass
class PM1Detail:
    """PM1 domain detail from ClinGen."""
    gene: str
    strength: str  # "Moderate" or "Supporting"
    domain: str
    protein_positions: str  # e.g., "p.1-p.40; p.214-p.223"
    cdna_positions: str
    transcript: str
    notes: str
    application_criteria: str


@dataclass
class FrequencyCutoff:
    """Frequency cutoff for BA1, BS1, PM2 from ClinGen."""
    gene: str
    rule_code: str  # BA1, BS1, PM2
    strength: str  # "Stand Alone", "Strong", "Supporting"
    cutoff_operator: str  # "<", ">", "<=", ">="
    cutoff_value: float
    frequency_type: str  # "popmax", "total MAF", "allele frequency"
    database: str  # "gnomAD", "unknown"
    population_group: str
    notes: str


@dataclass
class ClinGenRule:
    """Gene-specific ACMG rule from ClinGen."""
    gene: str
    rule_code: str
    strength: str
    strength_order: int
    description: str
    application_criteria: str


class ClinGenRulesDB:
    """
    Query interface for clingen_rules.db.

    Usage:
        db = ClinGenRulesDB(db_path)
        pm1 = db.get_pm1_for_gene("BRCA1", protein_position=1850)
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

    def get_pm1_for_position(
        self,
        gene: str,
        protein_position: int,
    ) -> Optional[PM1Detail]:
        """
        Get PM1 detail for a specific protein position.

        Args:
            gene: Gene symbol (e.g., "BRCA1")
            protein_position: Protein amino acid position (e.g., 1850)

        Returns:
            PM1Detail if position is in a defined domain, None otherwise
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        cursor.execute(
            "SELECT * FROM pm1_details WHERE gene = ?",
            (gene,)
        )

        for row in cursor.fetchall():
            protein_pos_str = row["protein_positions"]
            if not protein_pos_str:
                continue

            # Parse protein positions (e.g., "p.1-p.40; p.214-p.223; p.Arg329")
            if self._position_in_protein_range(protein_position, protein_pos_str):
                return PM1Detail(
                    gene=row["gene"],
                    strength=row["strength"],
                    domain=row["domain"] or "",
                    protein_positions=protein_pos_str,
                    cdna_positions=row["cdna_positions"] or "",
                    transcript=row["transcript"] or "",
                    notes=row["notes"] or "",
                    application_criteria=row["application_criteria"] or "",
                )

        return None

    def _position_in_protein_range(self, pos: int, protein_range_str: str) -> bool:
        """
        Check if a protein position is within a protein range string.

        Range formats supported:
        - "p.1-p.40" (range)
        - "p.Arg329" (single position)
        - "p.1-p.40; p.214-p.223; p.Arg329" (multiple)
        """
        # Remove "p." prefix and split by semicolon
        parts = protein_range_str.replace("p.", "").split(";")

        for part in parts:
            part = part.strip()
            if not part:
                continue

            # Try range format "Arg329" or "1-40"
            match = re.match(r'([A-Za-z]+)?(\d+)[-]?(\d+)?', part)
            if not match:
                continue

            prefix = match.group(1)  # e.g., "Arg" or None
            start = int(match.group(2))
            end = match.group(3)

            if end:
                # Range format
                if start <= pos <= int(end):
                    return True
            else:
                # Single position format
                if pos == start:
                    return True

        return False

    def get_all_pm1_for_gene(self, gene: str) -> List[PM1Detail]:
        """
        Get all PM1 domains for a gene.

        Args:
            gene: Gene symbol

        Returns:
            List of PM1Detail for the gene
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        cursor.execute(
            "SELECT * FROM pm1_details WHERE gene = ?",
            (gene,)
        )

        results = []
        for row in cursor.fetchall():
            results.append(PM1Detail(
                gene=row["gene"],
                strength=row["strength"],
                domain=row["domain"] or "",
                protein_positions=row["protein_positions"] or "",
                cdna_positions=row["cdna_positions"] or "",
                transcript=row["transcript"] or "",
                notes=row["notes"] or "",
                application_criteria=row["application_criteria"] or "",
            ))

        return results

    def get_frequency_cutoff(
        self,
        gene: str,
        rule_code: str,
        strength: Optional[str] = None,
    ) -> Optional[FrequencyCutoff]:
        """
        Get frequency cutoff for a gene/rule/strength combination.

        Args:
            gene: Gene symbol
            rule_code: Rule code (BA1, BS1, PM2)
            strength: Optional strength filter (Stand Alone, Strong, Supporting)

        Returns:
            FrequencyCutoff if found, None otherwise
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        if strength:
            cursor.execute(
                """SELECT * FROM frequency_cutoffs
                   WHERE gene = ? AND rule_code = ? AND strength = ?
                   LIMIT 1""",
                (gene, rule_code, strength)
            )
        else:
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

        return FrequencyCutoff(
            gene=row["gene"],
            rule_code=row["rule_code"],
            strength=row["strength"],
            cutoff_operator=row["cutoff_operator"],
            cutoff_value=row["cutoff_value"],
            frequency_type=row["frequency_type"] or "",
            database=row["database"] or "",
            population_group=row["population_group"] or "",
            notes=row["notes"] or "",
        )

    def get_clingen_rule(
        self,
        gene: str,
        rule_code: str,
        strength: Optional[str] = None,
    ) -> Optional[ClinGenRule]:
        """
        Get ClinGen rule for a gene/rule/strength combination.

        Args:
            gene: Gene symbol
            rule_code: Rule code (e.g., "PM1", "BS3", "PP4")
            strength: Optional strength filter

        Returns:
            ClinGenRule if found, None otherwise
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        if strength:
            cursor.execute(
                """SELECT * FROM clingen_rules
                   WHERE gene = ? AND rule_code = ? AND strength = ?
                   LIMIT 1""",
                (gene, rule_code, strength)
            )
        else:
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

        return ClinGenRule(
            gene=row["gene"],
            rule_code=row["rule_code"],
            strength=row["strength"],
            strength_order=row["strength_order"],
            description=row["description"] or "",
            application_criteria=row["application_criteria"] or "",
        )

    def get_all_rules_for_gene(self, gene: str) -> List[ClinGenRule]:
        """
        Get all ClinGen rules for a gene.

        Args:
            gene: Gene symbol

        Returns:
            List of ClinGenRule for the gene
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        cursor.execute(
            """SELECT * FROM clingen_rules
               WHERE gene = ?
               ORDER BY rule_code, strength_order DESC""",
            (gene,)
        )

        results = []
        for row in cursor.fetchall():
            results.append(ClinGenRule(
                gene=row["gene"],
                rule_code=row["rule_code"],
                strength=row["strength"],
                strength_order=row["strength_order"],
                description=row["description"] or "",
                application_criteria=row["application_criteria"] or "",
            ))

        return results

    def check_frequency_applies(
        self,
        gene: str,
        rule_code: str,
        frequency: float,
        strength: Optional[str] = None,
    ) -> Tuple[bool, Optional[str], str]:
        """
        Check if a frequency value triggers a BA1/BS1/PM2 rule.

        Args:
            gene: Gene symbol
            rule_code: Rule code (BA1, BS1, PM2)
            frequency: Population frequency value
            strength: Optional specific strength to check

        Returns:
            Tuple of (applies, applied_strength, reason)
        """
        cutoff = self.get_frequency_cutoff(gene, rule_code, strength)
        if cutoff is None:
            return False, None, f"No {rule_code} cutoff defined for {gene}"

        # Check if frequency meets the cutoff
        # Handle Unicode operators: ≥ (>=), ≤ (<=), > (>), < (<)
        applies = False
        op = cutoff.cutoff_operator
        if op in (">=", "≥"):
            applies = frequency >= cutoff.cutoff_value
        elif op in ("<=", "≤"):
            applies = frequency <= cutoff.cutoff_value
        elif op == ">":
            applies = frequency > cutoff.cutoff_value
        elif op == "<":
            applies = frequency < cutoff.cutoff_value
        elif cutoff.cutoff_operator == "<=":
            applies = frequency <= cutoff.cutoff_value

        if applies:
            return True, cutoff.strength, cutoff.notes
        else:
            return False, None, f"Frequency {frequency} does not meet {cutoff.cutoff_operator} {cutoff.cutoff_value}"


def get_pm1_for_variant(
    gene: str,
    protein_position: int,
    db_path: Path,
) -> Optional[PM1Detail]:
    """
    Convenience function to get PM1 detail for a variant.

    Args:
        gene: Gene symbol
        protein_position: Protein amino acid position
        db_path: Path to clingen_rules.db

    Returns:
        PM1Detail if position is in a ClinGen defined domain
    """
    if not db_path or not db_path.exists():
        return None

    with ClinGenRulesDB(db_path) as db:
        return db.get_pm1_for_position(gene, protein_position)
