#!/usr/bin/env python3
"""
Functional Evidence Database query module.

Queries functional_evidence.db for MAVE functional assay data.
Provides lookup by transcript + cDNA position or by gene + position.
"""

import sqlite3
import logging
from dataclasses import dataclass
from typing import Optional, Tuple, List
from pathlib import Path

logger = logging.getLogger("GenOtoScope_Classify.functional_evidence_db")


@dataclass
class FunctionalEvidenceResult:
    """Result from functional evidence database lookup."""
    transcript_id: str
    gene_name: str
    chrom: Optional[str]
    position: int
    cdna_pos: str
    classification: str  # Loss-of-function, Functionally normal, Gain-of-function
    classification_review: str  # "" for consistent, "冲突" for conflicting
    cross_assay_hits: int
    mave_technique: str
    unique_publications: str  # Pipe-separated PMIDs
    score: float
    phenotype: str
    functional_description: str

    @property
    def is_pathogenic(self) -> bool:
        """Return True if classification supports pathogenicity."""
        return self.classification in ("Loss-of-function", "Gain-of-function")

    @property
    def is_benign(self) -> bool:
        """Return True if classification supports benign."""
        return self.classification == "Functionally normal"

    @property
    def is_conflicting(self) -> bool:
        """Return True if there are conflicting evidence from different publications."""
        return self.classification_review == "冲突"

    def to_functional_data(self) -> Tuple[bool, bool]:
        """
        Convert to FunctionalData tuple (pathogenic, benign).

        For conflicting evidence, returns (True, True) to indicate conflicting.
        """
        if self.is_conflicting:
            # Conflicting evidence - both flags set to indicate conflict
            return (True, True)
        return (self.is_pathogenic, self.is_benign)


class FunctionalEvidenceDB:
    """
    Query interface for functional_evidence.db.

    Usage:
        db = FunctionalEvidenceDB(db_path)
        result = db.lookup(transcript_id="NM_005157.6", cdna_pos="c.80-4")
        if result:
            is_patho, is_benign = result.to_functional_data()
    """

    def __init__(self, db_path: Path):
        """
        Initialize with path to functional_evidence.db.

        Args:
            db_path: Path to functional_evidence.db
        """
        self.db_path = db_path
        self._conn: Optional[sqlite3.Connection] = None

    def _get_connection(self) -> sqlite3.Connection:
        """Get or create database connection."""
        if self._conn is None:
            self._conn = sqlite3.connect(str(self.db_path))
            self._conn.row_factory = sqlite3.Row
        return self._conn

    def close(self):
        """Close database connection."""
        if self._conn:
            self._conn.close()
            self._conn = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def lookup_by_transcript_cdna(
        self,
        transcript_id: str,
        cdna_pos: str,
    ) -> Optional[FunctionalEvidenceResult]:
        """
        Lookup by transcript ID and cDNA position.

        Args:
            transcript_id: e.g., "NM_005157.6"
            cdna_pos: e.g., "c.80-4" or "c.1092"

        Returns:
            FunctionalEvidenceResult if found, None otherwise
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        # Normalize transcript ID (remove version for matching)
        base_transcript = transcript_id.split(".")[0] if "." in transcript_id else transcript_id

        cursor.execute(
            """
            SELECT * FROM functional_evidence
            WHERE transcript_id LIKE ?
            AND cdna_pos = ?
            LIMIT 1
            """,
            (f"{base_transcript}%", cdna_pos),
        )

        row = cursor.fetchone()
        if row is None:
            return None

        return self._row_to_result(row)

    def lookup_by_gene_position(
        self,
        gene_name: str,
        position: int,
        chrom: Optional[str] = None,
    ) -> Optional[FunctionalEvidenceResult]:
        """
        Lookup by gene name and genomic position.

        Args:
            gene_name: Gene symbol, e.g., "ABL1"
            position: Genomic position
            chrom: Optional chromosome (e.g., "9" or "chr9")

        Returns:
            FunctionalEvidenceResult if found, None otherwise
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        if chrom:
            chrom_clean = chrom.replace("chr", "")
            cursor.execute(
                """
                SELECT * FROM functional_evidence
                WHERE gene_name = ? AND position = ? AND chrom = ?
                LIMIT 1
                """,
                (gene_name, position, chrom_clean),
            )
        else:
            cursor.execute(
                """
                SELECT * FROM functional_evidence
                WHERE gene_name = ? AND position = ?
                LIMIT 1
                """,
                (gene_name, position),
            )

        row = cursor.fetchone()
        if row is None:
            return None

        return self._row_to_result(row)

    def lookup_by_gene_cdna(
        self,
        gene_name: str,
        cdna_pos: str,
    ) -> Optional[FunctionalEvidenceResult]:
        """
        Lookup by gene name and cDNA position.

        Args:
            gene_name: Gene symbol, e.g., "ABL1"
            cdna_pos: e.g., "c.80-4"

        Returns:
            FunctionalEvidenceResult if found, None otherwise
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        cursor.execute(
            """
            SELECT * FROM functional_evidence
            WHERE gene_name = ? AND cdna_pos = ?
            LIMIT 1
            """,
            (gene_name, cdna_pos),
        )

        row = cursor.fetchone()
        if row is None:
            return None

        return self._row_to_result(row)

    def lookup_by_transcript_cdna_refalt(
        self,
        transcript_id: str,
        cdna_pos: str,
        ref_alt: str,
    ) -> Optional[FunctionalEvidenceResult]:
        """
        Lookup by transcript ID, cDNA position, and Ref/Alt (exact match for splice variants).

        Args:
            transcript_id: e.g., "NM_005157.6"
            cdna_pos: e.g., "c.80-4" or "c.1092"
            ref_alt: e.g., "C/T" or "G/A"

        Returns:
            FunctionalEvidenceResult if found, None otherwise
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        # Normalize transcript ID (remove version for matching)
        base_transcript = transcript_id.split(".")[0] if "." in transcript_id else transcript_id

        cursor.execute(
            """
            SELECT * FROM functional_evidence
            WHERE transcript_id LIKE ?
            AND cdna_pos = ?
            AND ref_alt = ?
            LIMIT 1
            """,
            (f"{base_transcript}%", cdna_pos, ref_alt),
        )

        row = cursor.fetchone()
        if row is None:
            return None

        return self._row_to_result(row)

    def lookup_by_gene_cdna_refalt(
        self,
        gene_name: str,
        cdna_pos: str,
        ref_alt: str,
    ) -> Optional[FunctionalEvidenceResult]:
        """
        Lookup by gene name, cDNA position, and Ref/Alt (exact match for splice variants).

        Args:
            gene_name: Gene symbol, e.g., "ABL1"
            cdna_pos: e.g., "c.80-4"
            ref_alt: e.g., "C/T" or "G/A"

        Returns:
            FunctionalEvidenceResult if found, None otherwise
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        cursor.execute(
            """
            SELECT * FROM functional_evidence
            WHERE gene_name = ? AND cdna_pos = ? AND ref_alt = ?
            LIMIT 1
            """,
            (gene_name, cdna_pos, ref_alt),
        )

        row = cursor.fetchone()
        if row is None:
            return None

        return self._row_to_result(row)

    def get_all_for_gene(self, gene_name: str) -> List[FunctionalEvidenceResult]:
        """
        Get all functional evidence entries for a gene.

        Args:
            gene_name: Gene symbol

        Returns:
            List of FunctionalEvidenceResult
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        cursor.execute(
            """
            SELECT * FROM functional_evidence
            WHERE gene_name = ?
            ORDER BY position
            """,
            (gene_name,),
        )

        return [self._row_to_result(row) for row in cursor.fetchall()]

    def _row_to_result(self, row: sqlite3.Row) -> FunctionalEvidenceResult:
        """Convert sqlite3.Row to FunctionalEvidenceResult."""
        return FunctionalEvidenceResult(
            transcript_id=row["transcript_id"],
            gene_name=row["gene_name"],
            chrom=row["chrom"],
            position=row["position"],
            cdna_pos=row["cdna_pos"],
            classification=row["classification"],
            classification_review=row["classification_review"],
            cross_assay_hits=row["cross_assay_hits"],
            mave_technique=row["mave_technique"],
            unique_publications=row["unique_publications"],
            score=row["score"],
            phenotype=row["phenotype"],
            functional_description=row["functional_description"],
        )


def get_functional_evidence(
    transcript_id: str,
    cdna_pos: str,
    gene_name: str,
    position: int,
    chrom: str,
    db_path: Path,
) -> Optional[FunctionalEvidenceResult]:
    """
    Convenience function to get functional evidence for a variant.

    Tries multiple lookup strategies:
    1. Transcript + cDNA position (most specific)
    2. Gene + cDNA position
    3. Gene + genomic position

    Args:
        transcript_id: e.g., "NM_005157.6"
        cdna_pos: e.g., "c.80-4"
        gene_name: Gene symbol
        position: Genomic position
        chrom: Chromosome
        db_path: Path to functional_evidence.db

    Returns:
        FunctionalEvidenceResult if found, None otherwise
    """
    if not db_path or not db_path.exists():
        logger.warning(f"Functional evidence DB not found: {db_path}")
        return None

    with FunctionalEvidenceDB(db_path) as db:
        # Try transcript + cDNA first (most specific)
        result = db.lookup_by_transcript_cdna(transcript_id, cdna_pos)
        if result:
            return result

        # Try gene + cDNA
        result = db.lookup_by_gene_cdna(gene_name, cdna_pos)
        if result:
            return result

        # Try gene + position
        result = db.lookup_by_gene_position(gene_name, position, chrom)
        if result:
            return result

        return None
