#!/usr/bin/env python3
"""
Ensembl REST API client.
Provides access to general Ensembl REST endpoints (gene info, transcript info, etc.).
"""

import logging
from typing import Optional, Dict, Any, List

from cloud_api.base_client import BaseCloudClient

logger = logging.getLogger(__name__)


class EnsemblCloudClient(BaseCloudClient):
    """
    Ensembl REST API client.
    Accesses rest.ensembl.org for general queries (genes, transcripts, homologs, etc.).
    """

    BASE_URL = "https://rest.ensembl.org"

    def __init__(self, timeout: int = 30, max_retries: int = 3):
        """
        Initialize Ensembl REST API client.

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

    def get_gene_info(
        self,
        gene_id: str,
        species: str = "homo_sapiens",
    ) -> Optional[Dict[str, Any]]:
        """
        Get gene information by gene ID.

        Args:
            gene_id: Gene ID (e.g., "ENSG00000012048")
            species: Species (default: homo_sapiens)

        Returns:
            Gene info dict or None
        """
        endpoint = f"/lookup/symbol/{species}/{gene_id}"

        response = self.get(endpoint=endpoint)

        if not response.success:
            logger.error(f"Ensembl gene lookup failed: {response.error}")
            return None

        return response.data

    def get_gene_by_symbol(
        self,
        symbol: str,
        species: str = "homo_sapiens",
    ) -> Optional[Dict[str, Any]]:
        """
        Get gene information by symbol.

        Args:
            symbol: Gene symbol (e.g., "BRCA1")
            species: Species (default: homo_sapiens)

        Returns:
            Gene info dict or None
        """
        endpoint = f"/lookup/symbol/{species}/{symbol}"

        response = self.get(endpoint=endpoint)

        if not response.success:
            logger.error(f"Ensembl symbol lookup failed: {response.error}")
            return None

        return response.data

    def get_transcript_info(
        self,
        transcript_id: str,
        species: str = "homo_sapiens",
    ) -> Optional[Dict[str, Any]]:
        """
        Get transcript information by transcript ID.

        Args:
            transcript_id: Transcript ID (e.g., "ENST00000352993")
            species: Species (default: homo_sapiens)

        Returns:
            Transcript info dict or None
        """
        endpoint = f"/lookup/id/{transcript_id}"

        response = self.get(endpoint=endpoint)

        if not response.success:
            logger.error(f"Ensembl transcript lookup failed: {response.error}")
            return None

        return response.data

    def get_sequence(
        self,
        region: str,
        species: str = "homo_sapiens",
        coord_system: str = "chromosome",
    ) -> Optional[str]:
        """
        Get sequence for a region.

        Args:
            region: Region string (e.g., "17:43045678:43045678:1")
            species: Species
            coord_system: Coordinate system

        Returns:
            Sequence string or None
        """
        endpoint = f"/sequence/region/{species}/{region}"

        params = {"coord_system": coord_system}

        response = self.get(endpoint=endpoint, params=params)

        if not response.success:
            logger.error(f"Ensembl sequence fetch failed: {response.error}")
            return None

        return response.data.get("seq") if isinstance(response.data, dict) else None

    def get_variants(
        self,
        region: str,
        species: str = "homo_sapiens",
        coordinator_system: str = "chromosome",
    ) -> List[Dict[str, Any]]:
        """
        Get variants in a region.

        Args:
            region: Region string (e.g., "17:43000000:43100000")
            species: Species
            coordinator_system: Coordinate system

        Returns:
            List of variant records
        """
        endpoint = f"/variation/{species}/region"

        params = {
            "region": region,
            "coord_system": coordinator_system,
        }

        response = self.get(endpoint=endpoint, params=params)

        if not response.success:
            logger.error(f"Ensembl variants fetch failed: {response.error}")
            return []

        variants = response.data.get("variants", []) if isinstance(response.data, dict) else []
        return variants

    def get_homologs(
        self,
        gene_id: str,
        species: str = "homo_sapiens",
    ) -> List[Dict[str, Any]]:
        """
        Get homologs for a gene.

        Args:
            gene_id: Gene ID
            species: Species

        Returns:
            List of homolog records
        """
        endpoint = f"/homology/id/{species}/{gene_id}"

        response = self.get(endpoint=endpoint)

        if not response.success:
            logger.error(f"Ensembl homolog lookup failed: {response.error}")
            return []

        data = response.data
        if isinstance(data, dict):
            homologies = data.get("data", [{}])[0].get("homologies", [])
            return homologies

        return []

    def get_domain_info(
        self,
        transcript_id: str,
        species: str = "homo_sapiens",
    ) -> List[Dict[str, Any]]:
        """
        Get protein domain information for a transcript.

        Args:
            transcript_id: Transcript ID
            species: Species

        Returns:
            List of domain records
        """
        endpoint = f"/cloud/domains/{species}/{transcript_id}"

        response = self.get(endpoint=endpoint)

        if not response.success:
            logger.error(f"Ensembl domain lookup failed: {response.error}")
            return []

        if isinstance(response.data, dict):
            return response.data.get("domains", [])

        return []

    def get_mane_transcript(
        self,
        gene_id: str,
        species: str = "homo_sapiens",
    ) -> Optional[Dict[str, Any]]:
        """
        Get MANE transcript for a gene.

        Args:
            gene_id: Gene ID (e.g., "ENSG00000012048")
            species: Species

        Returns:
            MANE transcript info or None
        """
        endpoint = f"/info/mane/{species}/{gene_id}"

        response = self.get(endpoint=endpoint)

        if not response.success:
            logger.error(f"Ensembl MANE lookup failed: {response.error}")
            return None

        return response.data