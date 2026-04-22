#!/usr/bin/env python3
"""
gnomAD API client for gene constraint data.
Provides access to gnomAD's gene-level constraint metrics including missense z-score.
"""

import logging
from typing import Optional, Dict, Any

from cloud_api.base_client import BaseCloudClient

logger = logging.getLogger(__name__)


class GnomADConstraintClient(BaseCloudClient):
    """
    gnomAD gene constraint API client.
    Accesses gnomAD for gene-level constraint metrics.
    """

    BASE_URL = "https://gnomad.broadinstitute.org/api"

    def __init__(self, timeout: int = 30, max_retries: int = 3):
        """
        Initialize gnomAD constraint client.

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

    def get_gene_constraint(
        self,
        gene_symbol: str,
        release_version: str = "v4",
    ) -> Optional[Dict[str, Any]]:
        """
        Get gene constraint metrics for a gene.

        Args:
            gene_symbol: Gene symbol (e.g., "BRCA1", "TP53")
            release_version: gnomAD release version ("v4" or "v2")

        Returns:
            Dict with constraint metrics including missense_zscore, or None if not found
        """
        # GraphQL query for gene constraint
        query = """
        query GeneConstraint($geneSymbol: String!, $releaseVersion: String!) {
            gene_constraint(geneSymbol: $geneSymbol, releaseVersion: $releaseVersion) {
                gene_id
                gene_symbol
                mis_z
                mis_z_rank
                mis_z_approx
                n_precision
                obs_mis
                exp_mis
                oe_mis
                oe_mis_upper
                oe_mis_lower
                oe_mis_std
            }
        }
        """

        variables = {
            "geneSymbol": gene_symbol,
            "releaseVersion": release_version,
        }

        response = self.post(
            endpoint="/",
            data={"query": query, "variables": variables},
        )

        if not response.success:
            logger.error(f"gnomAD constraint lookup failed: {response.error}")
            return None

        data = response.data
        if "data" not in data or "gene_constraint" not in data["data"]:
            return None

        constraint_data = data["data"]["gene_constraint"]
        if not constraint_data:
            return None

        return {
            "gene_id": constraint_data.get("gene_id"),
            "gene_symbol": constraint_data.get("gene_symbol"),
            "missense_zscore": constraint_data.get("mis_z"),
            "missense_z_rank": constraint_data.get("mis_z_rank"),
            "missense_z_approx": constraint_data.get("mis_z_approx"),
            "observed_missense": constraint_data.get("obs_mis"),
            "expected_missense": constraint_data.get("exp_mis"),
            "oe_mis": constraint_data.get("oe_mis"),
            "oe_mis_upper": constraint_data.get("oe_mis_upper"),
            "oe_mis_lower": constraint_data.get("oe_mis_lower"),
        }

    def get_missense_zscore(
        self,
        gene_symbol: str,
        release_version: str = "v4",
    ) -> Optional[float]:
        """
        Get missense z-score for a gene.

        Args:
            gene_symbol: Gene symbol (e.g., "BRCA1", "TP53")
            release_version: gnomAD release version ("v4" or "v2")

        Returns:
            Missense z-score or None if not available
        """
        constraint = self.get_gene_constraint(gene_symbol, release_version)
        if constraint:
            return constraint.get("missense_zscore")
        return None


def load_gnomad_constraint_from_tsv(
    tsv_path: str,
    gene_symbol: str,
) -> Optional[float]:
    """
    Load missense z-score from a precomputed TSV file.

    Args:
        tsv_path: Path to TSV file with gene constraint data
        gene_symbol: Gene symbol to look up

    Returns:
        Missense z-score or None if not found
    """
    import csv

    try:
        with open(tsv_path, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row.get("gene_symbol", "").upper() == gene_symbol.upper():
                    zscore = row.get("mis_z") or row.get("missense_zscore")
                    if zscore:
                        return float(zscore)
        logger.warning(f"Gene {gene_symbol} not found in constraint TSV")
        return None
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Error loading gnomAD constraint from TSV: {e}")
        return None
