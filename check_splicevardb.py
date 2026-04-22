#!/usr/bin/env python3
"""
SpliceVarDB functional evidence annotation.
Based on splicevardb_with_pmid.csv - a database of splice-altering variants with experimental evidence.

Classification rules:
- Splice-altering → PS3 (pathogenic functional support)
- Normal → BS3 (benign functional support)
- Low-frequency → no functional evidence applied
"""

import pathlib
import csv
import re
import logging
from typing import Optional, Tuple, Dict, Any

from dataclasses import dataclass

logger = logging.getLogger("GenOtoScope_Classify.check_splicevardb")


@dataclass
class SpliceVarDBEntry:
    """Entry from SpliceVarDB database."""
    variant_id: str
    hg19: str
    hg38: str
    gene: str
    hgvs: str
    method: str
    classification: str  # Splice-altering, Normal, Low-frequency
    location: str  # Exonic, Intronic
    doi: str
    pmid: str


class SpliceVarDB:
    """
    SpliceVarDB functional evidence database.

    Loads from splicevardb_with_pmid.csv and provides lookup methods.
    """

    def __init__(self, db_path: pathlib.Path):
        """
        Initialize SpliceVarDB.

        Args:
            db_path: Path to splicevardb_with_pmid.csv
        """
        self.db_path = db_path
        self.entries: Dict[str, SpliceVarDBEntry] = {}
        self._load_database()

    def _load_database(self) -> None:
        """Load SpliceVarDB from CSV file."""
        if not self.db_path.exists():
            logger.warning(f"SpliceVarDB file not found: {self.db_path}")
            return

        with open(self.db_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    entry = SpliceVarDBEntry(
                        variant_id=row["variant_id"],
                        hg19=row["hg19"],
                        hg38=row["hg38"],
                        gene=row["gene"],
                        hgvs=row["hgvs"],
                        method=row["method"],
                        classification=row["classification"],
                        location=row["location"],
                        doi=row.get("doi", ""),
                        pmid=row.get("PMID", ""),
                    )
                    # Index by gene+hgvs for quick lookup
                    key = self._make_key(entry.gene, entry.hgvs)
                    self.entries[key] = entry
                except (KeyError, ValueError) as e:
                    logger.warning(f"Skipping invalid SpliceVarDB row: {e}")
                    continue

        logger.info(f"Loaded {len(self.entries)} entries from SpliceVarDB")

    def _make_key(self, gene: str, hgvs: str) -> str:
        """Create lookup key from gene and HGVS."""
        return f"{gene.upper()}:{hgvs}"

    def _parse_hgvs_position(self, hgvs: str) -> Optional[Tuple[str, int]]:
        """
        Parse HGVS cDNA notation to extract transcript and position.

        Args:
            hgvs: HGVS notation like NM_194292.3:c.1092A>G

        Returns:
            Tuple of (transcript_id, cDNA_position) or None
        """
        # Match patterns like NM_194292.3:c.1092A>G or NM_194292.3:c.670-1G>T
        pattern = r'(NM_\d+\.\d+):c\.(\d+[-+]\d+|\d+)'
        match = re.search(pattern, hgvs)
        if match:
            transcript = match.group(1)
            pos_str = match.group(2)
            # Handle compound positions (e.g., 670-1) by just taking the first number
            pos = int(pos_str.split('-')[0].split('+')[0])
            return transcript, pos
        return None

    def lookup_by_gene_hgvs(self, gene: str, hgvs: str) -> Optional[SpliceVarDBEntry]:
        """
        Lookup variant by gene and HGVS notation.

        Args:
            gene: Gene symbol
            hgvs: HGVS cDNA notation

        Returns:
            SpliceVarDBEntry if found, None otherwise
        """
        key = self._make_key(gene, hgvs)
        return self.entries.get(key)

    def lookup_by_gene_position(self, gene: str, cDNA_pos: int, transcript: str = None) -> Optional[SpliceVarDBEntry]:
        """
        Lookup variant by gene and cDNA position.

        Args:
            gene: Gene symbol
            cDNA_pos: cDNA position (nucleotide number)
            transcript: Optional transcript ID for disambiguation

        Returns:
            SpliceVarDBEntry if found, None otherwise
        """
        gene_upper = gene.upper()

        for key, entry in self.entries.items():
            if entry.gene.upper() == gene_upper:
                parsed = self._parse_hgvs_position(entry.hgvs)
                if parsed:
                    entry_transcript, entry_pos = parsed
                    if entry_pos == cDNA_pos:
                        if transcript is None or entry_transcript == transcript:
                            return entry
        return None

    def lookup_by_genomic_position(self, chrom: str, pos: int, gene: str, build: str = "hg38") -> Optional[SpliceVarDBEntry]:
        """
        Lookup variant by genomic position.

        Args:
            chrom: Chromosome (e.g., "1" or "chr1")
            pos: Genomic position
            gene: Gene symbol
            build: Genome build ("hg38" or "hg19")

        Returns:
            SpliceVarDBEntry if found, None otherwise
        """
        chrom_clean = chrom.replace("chr", "")
        gene_upper = gene.upper()

        pos_col = "hg38" if build == "hg38" else "hg19"

        for entry in self.entries.values():
            if entry.gene.upper() != gene_upper:
                continue

            # Parse position from entry (format: "1-100573238-T-C")
            parts = getattr(entry, pos_col).split("-")
            if len(parts) >= 2:
                try:
                    entry_chrom = parts[0]
                    entry_pos = int(parts[1])
                    if entry_chrom == chrom_clean and entry_pos == pos:
                        return entry
                except (ValueError, IndexError):
                    continue

        return None

    def get_functional_evidence(self, gene: str, hgvs: str = None, cDNA_pos: int = None,
                                  chrom: str = None, pos: int = None,
                                  transcript: str = None, build: str = "hg38") -> Tuple[bool, bool, str]:
        """
        Get functional evidence recommendation for a variant.

        Args:
            gene: Gene symbol
            hgvs: HGVS cDNA notation (preferred)
            cDNA_pos: cDNA position (alternative)
            chrom: Chromosome (for genomic lookup)
            pos: Genomic position (for genomic lookup)
            transcript: Transcript ID for disambiguation
            build: Genome build ("hg38" or "hg19")

        Returns:
            Tuple of (has_pathogenic_evidence, has_benign_evidence, reason)
            - (True, False, "Splice-altering") for PS3
            - (False, True, "Normal") for BS3
            - (False, False, "Low-frequency or not found") for no evidence
        """
        entry = None

        # Try HGVS lookup first
        if hgvs:
            entry = self.lookup_by_gene_hgvs(gene, hgvs)

        # Try cDNA position lookup
        if entry is None and cDNA_pos is not None:
            entry = self.lookup_by_gene_position(gene, cDNA_pos, transcript)

        # Try genomic lookup
        if entry is None and chrom is not None and pos is not None:
            entry = self.lookup_by_genomic_position(chrom, pos, gene, build)

        if entry is None:
            return False, False, "Variant not found in SpliceVarDB"

        classification = entry.classification.lower()

        if classification == "splice-altering":
            return True, False, f"Splice-altering (method: {entry.method}, PMID:{entry.pmid})"
        elif classification == "normal":
            return False, True, f"Normal splicing (method: {entry.method}, PMID:{entry.pmid})"
        elif classification == "low-frequency":
            return False, False, f"Low-frequency ({entry.method})"
        else:
            return False, False, f"Unknown classification: {entry.classification}"


def get_splicevardb_evidence(
    gene: str,
    hgvs: str = None,
    cDNA_pos: int = None,
    chrom: str = None,
    pos: int = None,
    transcript: str = None,
    splicevardb_path: pathlib.Path = None,
) -> Tuple[bool, bool, str]:
    """
    Convenience function to get SpliceVarDB functional evidence.

    Args:
        gene: Gene symbol
        hgvs: HGVS cDNA notation
        cDNA_pos: cDNA position
        chrom: Chromosome
        pos: Genomic position
        transcript: Transcript ID
        splicevardb_path: Path to splicevardb_with_pmid.csv

    Returns:
        Tuple of (has_pathogenic_evidence, has_benign_evidence, reason)
    """
    if splicevardb_path is None:
        return False, False, "SpliceVarDB path not configured"

    db = SpliceVarDB(splicevardb_path)
    return db.get_functional_evidence(
        gene=gene,
        hgvs=hgvs,
        cDNA_pos=cDNA_pos,
        chrom=chrom,
        pos=pos,
        transcript=transcript,
    )
