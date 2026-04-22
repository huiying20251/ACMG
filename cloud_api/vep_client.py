#!/usr/bin/env python3
"""
VEP REST API client.
Provides access to Ensembl VEP (Variant Effect Predictor) for variant annotation.
"""

import logging
from typing import Optional, Dict, Any, List

from cloud_api.base_client import BaseCloudClient

logger = logging.getLogger(__name__)


class VEPCloudClient(BaseCloudClient):
    """
    VEP REST API client.
    Accesses rest.ensembl.org/vep for variant annotation.
    """

    BASE_URL = "https://rest.ensembl.org"

    def __init__(self, timeout: int = 30, max_retries: int = 3):
        """
        Initialize VEP REST API client.

        Args:
            timeout: Request timeout in seconds
            max_retries: Maximum number of retry attempts
        """
        super().__init__(
            base_url=self.BASE_URL,
            timeout=timeout,
            max_retries=max_retries,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
        )

    def annotate_region(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        species: str = "homo_sapiens",
        phenotypes: bool = True,
        CADD: bool = True,
        SpliceAI: bool = True,
        db_xref: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """
        Annotate a variant by genomic region.

        Args:
            chrom: Chromosome (e.g., "17" or "chr17")
            pos: Position
            ref: Reference allele
            alt: Alternate allele
            species: Species (default: homo_sapiens)
            phenotypes: Include phenotype information
            CADD: Include CADD scores
            SpliceAI: Include SpliceAI scores
            db_xref: Include database cross-references

        Returns:
            VEP annotation result or None
        """
        # Strip chr prefix if present
        chrom = chrom.replace("chr", "")

        endpoint = f"/vep/{species}/region/{chrom}:{pos}:{ref}/{alt}"

        params = {
            "phenotypes": str(phenotypes).lower(),
            "CADD": str(CADD).lower(),
            "SpliceAI": str(SpliceAI).lower(),
            "db_xref": str(db_xref).lower(),
        }

        response = self.get(endpoint=endpoint, params=params)

        if not response.success:
            logger.error(f"VEP region annotation failed: {response.error}")
            return None

        if isinstance(response.data, list) and len(response.data) > 0:
            return response.data[0]
        return response.data

    def annotate_rsid(
        self,
        rsid: str,
        species: str = "homo_sapiens",
        phenotypes: bool = True,
        CADD: bool = True,
        SpliceAI: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """
        Annotate a variant by rsID.

        Args:
            rsid: rsID (e.g., "rs123456")
            species: Species (default: homo_sapiens)
            phenotypes: Include phenotype information
            CADD: Include CADD scores
            SpliceAI: Include SpliceAI scores

        Returns:
            VEP annotation result or None
        """
        endpoint = f"/vep/{species}/id/{rsid}"

        params = {
            "phenotypes": str(phenotypes).lower(),
            "CADD": str(CADD).lower(),
            "SpliceAI": str(SpliceAI).lower(),
        }

        response = self.get(endpoint=endpoint, params=params)

        if not response.success:
            logger.error(f"VEP rsID annotation failed: {response.error}")
            return None

        if isinstance(response.data, list) and len(response.data) > 0:
            return response.data[0]
        return response.data

    def get_population_frequency(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        species: str = "homo_sapiens",
    ) -> Optional[Dict[str, float]]:
        """
        Get population frequencies for a variant.

        Args:
            chrom: Chromosome
            pos: Position
            ref: Reference allele
            alt: Alternate allele
            species: Species

        Returns:
            Dict with population frequencies (gnomAD, 1000G, etc.)
        """
        annotation = self.annotate_region(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            species=species,
            phenotypes=False,
            CADD=False,
            SpliceAI=False,
            db_xref=False,
        )

        if not annotation:
            return None

        frequencies = {}

        # Extract gnomAD frequencies
        colocal_variants = annotation.get("colocal_variants", [])
        for variant in colocal_variants:
            if variant.get("type") == "gnomAD":
                freq = variant.get("allele_frequency")
                if freq is not None:
                    frequencies["gnomAD"] = freq

        # Check frequencies from annotation
        frequencies_data = annotation.get("freq", {})
        if isinstance(frequencies_data, dict):
            for source, data in frequencies_data.items():
                if isinstance(data, dict) and "AF" in data:
                    frequencies[source] = data["AF"]

        return frequencies if frequencies else None

    def get_cadd_score(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        species: str = "homo_sapiens",
    ) -> Optional[float]:
        """
        Get CADD score for a variant.

        Args:
            chrom: Chromosome
            pos: Position
            ref: Reference allele
            alt: Alternate allele
            species: Species

        Returns:
            CADD score or None
        """
        annotation = self.annotate_region(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            species=species,
            phenotypes=False,
            CADD=True,
            SpliceAI=False,
            db_xref=False,
        )

        if not annotation:
            return None

        return annotation.get("cadd_score")

    def get_spliceai_score(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        species: str = "homo_sapiens",
    ) -> Optional[Dict[str, float]]:
        """
        Get SpliceAI scores for a variant.

        Args:
            chrom: Chromosome
            pos: Position
            ref: Reference allele
            alt: Alternate allele
            species: Species

        Returns:
            Dict with SpliceAI delta scores (spliceai_ds, spliceai_dg, etc.) or None
        """
        annotation = self.annotate_region(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            species=species,
            phenotypes=False,
            CADD=False,
            SpliceAI=True,
            db_xref=False,
        )

        if not annotation:
            return None

        spliceai = annotation.get("spliceai")
        if spliceai:
            return {
                "spliceai_ds": spliceai.get("spliceai_ds"),
                "spliceai_dg": spliceai.get("spliceai_dg"),
                "spliceai_dp": spliceai.get("spliceai_dp"),
                "spliceai_consequence": spliceai.get("spliceai_consequence"),
            }

        return None

    def get_transcript_consequences(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        species: str = "homo_sapiens",
    ) -> List[Dict[str, Any]]:
        """
        Get transcript-specific consequences for a variant.

        Args:
            chrom: Chromosome
            pos: Position
            ref: Reference allele
            alt: Alternate allele
            species: Species

        Returns:
            List of transcript consequence objects
        """
        annotation = self.annotate_region(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            species=species,
            phenotypes=False,
            CADD=False,
            SpliceAI=False,
            db_xref=False,
        )

        if not annotation:
            return []

        return annotation.get("transcript_consequences", [])

    def get_mane_transcript(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        species: str = "homo_sapiens",
    ) -> Optional[Dict[str, Any]]:
        """
        Get MANE transcript for a variant.

        Args:
            chrom: Chromosome
            pos: Position
            ref: Reference allele
            alt: Alternate allele
            species: Species

        Returns:
            MANE transcript info or None
        """
        annotation = self.annotate_region(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            species=species,
            phenotypes=False,
            CADD=False,
            SpliceAI=False,
            db_xref=False,
        )

        if not annotation:
            return None

        # Look for MANE Select or MANE Plus Clinical
        transcript_consequences = annotation.get("transcript_consequences", [])
        for tc in transcript_consequences:
            if tc.get("MANE_SELECT") or tc.get("MANE_PLUS_CLINICAL"):
                return tc

        return None