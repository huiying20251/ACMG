#!/usr/bin/env python3
"""
PubTator API Client

Retrieves variant annotations from PubTator 3 API.
PubTator provides automated gene, variant, and disease annotations for PubMed articles.

Reference: https://www.ncbi.nlm.nih.gov/research/pubtator3/
"""

import logging
from typing import List, Optional, Dict, Any

import requests

logger = logging.getLogger(__name__)


class PubTatorClient:
    """
    Client for PubTator 3 API.

    Provides:
    - Gene annotations
    - Variant annotations (using rsID or HGVS)
    - Disease annotations
    """

    BASE_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3/api"

    # PubTator 3 endpoints
    ANNOTATE_PMC = "/publications/pmc/{pmcid}/annotations"
    ANNOTATE_PMID = "/publications/pmid/{pmid}/annotations"
    SEARCH_GENE = "/genes/search/{query}"
    SEARCH_VARIANT = "/variants/search/{query}"

    def __init__(self, timeout: int = 30):
        self.timeout = timeout

    def get_annotations_for_pmc(
        self,
        pmcid: str,
        include_abstract: bool = True,
    ) -> Dict[str, Any]:
        """
        Get all annotations for a PMC article.

        Args:
            pmcid: PMC article ID (e.g., "PMC1234567")
            include_abstract: Include abstract in annotations

        Returns:
            Dict with "entities" and "relations" keys
        """
        pmcid_clean = pmcid.replace("PMC", "")
        url = self.BASE_URL + self.ANNOTATE_PMC.format(pmcid=pmcid_clean)

        params = {}
        if include_abstract:
            params["include_abstract"] = "true"

        try:
            response = requests.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            return response.json()

        except requests.exceptions.RequestException as e:
            logger.warning(f"PubTator API failed for PMC{pmcid_clean}: {e}")
            return {"entities": [], "relations": []}

    def get_annotations_for_pmid(
        self,
        pmid: str,
        include_abstract: bool = True,
    ) -> Dict[str, Any]:
        """
        Get all annotations for a PubMed article.

        Args:
            pmid: PubMed article ID
            include_abstract: Include abstract in annotations

        Returns:
            Dict with "entities" and "relations" keys
        """
        url = self.BASE_URL + self.ANNOTATE_PMID.format(pmid=pmid)

        params = {}
        if include_abstract:
            params["include_abstract"] = "true"

        try:
            response = requests.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            return response.json()

        except requests.exceptions.RequestException as e:
            logger.warning(f"PubTator API failed for PMID{pmid}: {e}")
            return {"entities": [], "relations": []}

    def search_genes(self, query: str, max_results: int = 10) -> List[Dict[str, Any]]:
        """
        Search for genes by symbol or name.

        Args:
            query: Gene search query
            max_results: Maximum number of results

        Returns:
            List of gene entities
        """
        url = self.BASE_URL + self.SEARCH_GENE.format(query=query)
        params = {"size": max_results}

        try:
            response = requests.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            return data.get("genes", [])

        except requests.exceptions.RequestException as e:
            logger.warning(f"PubTator gene search failed for '{query}': {e}")
            return []

    def search_variants(self, query: str, max_results: int = 10) -> List[Dict[str, Any]]:
        """
        Search for variants by rsID or HGVS.

        Args:
            query: Variant query (rsID or HGVS notation)
            max_results: Maximum number of results

        Returns:
            List of variant entities
        """
        url = self.BASE_URL + self.SEARCH_VARIANT.format(query=query)
        params = {"size": max_results}

        try:
            response = requests.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            return data.get("variants", [])

        except requests.exceptions.RequestException as e:
            logger.warning(f"PubTator variant search failed for '{query}': {e}")
            return []


class VariantAnnotationExtractor:
    """
    Extract variant annotations from PubTator results.

    Used to find variant mentions in literature for:
    - PS1: Same variant reported pathogenic
    - PM5: Different variant at same position
    - BP6: Variant classified as benign
    """

    def __init__(self, pubtator_client: Optional[PubTatorClient] = None):
        self.pubtator = pubtator_client or PubTatorClient()

    def extract_variant_mentions(
        self,
        articles: List,  # List of Article objects
        variant_rsids: Optional[List[str]] = None,
        variant_hgvs: Optional[List[str]] = None,
    ) -> Dict[str, List[Dict[str, Any]]]:
        """
        Extract variant mentions from a list of articles.

        Args:
            articles: List of Article objects
            variant_rsids: List of rsIDs to look for (e.g., ["rs123456"])
            variant_hgvs: List of HGVS notations to look for

        Returns:
            Dict mapping PMID to list of variant mentions
        """
        mentions = {}

        rsid_set = set(variant_rsids) if variant_rsids else set()
        hgvs_set = set(variant_hgvs) if variant_hgvs else set()

        for article in articles:
            try:
                if article.pmcid:
                    annotations = self.pubtator.get_annotations_for_pmc(article.pmcid)
                elif article.pmid:
                    annotations = self.pubtator.get_annotations_for_pmid(article.pmid)
                else:
                    continue

                # Filter for variant entities only
                variants = [
                    ent for ent in annotations.get("entities", [])
                    if ent.get("type") == "Variant"
                ]

                if variants:
                    # Further filter if specific variants requested
                    if rsid_set or hgvs_set:
                        filtered_variants = []
                        for var in variants:
                            # Check if variant matches any of our targets
                            var_id = var.get("id", "")
                            var_hgvs = var.get("hgvs", "")

                            if var_id in rsid_set or any(h in var_hgvs for h in hgvs_set):
                                filtered_variants.append(var)

                        if filtered_variants:
                            mentions[article.pmid] = filtered_variants
                    else:
                        mentions[article.pmid] = variants

            except Exception as e:
                logger.debug(f"Failed to extract variants from {article.pmid}: {e}")
                continue

        return mentions

    def find_pathogenic_variant_reports(
        self,
        mentions: Dict[str, List[Dict[str, Any]]],
        min_cases: int = 1,
    ) -> List[Dict[str, Any]]:
        """
        Find variant mentions that report pathogenicity.

        Used for PS1 (same variant) and PM5 (different variant same position).

        Args:
            mentions: Dict from extract_variant_mentions
            min_cases: Minimum number of cases to report

        Returns:
            List of pathogenic variant reports
        """
        reports = []

        for pmid, variants in mentions.items():
            for var in variants:
                # Look for clinical significance in annotations
                significance = var.get("clinical_significance", "").lower()

                if significance in ("pathogenic", "likely pathogenic", "uncertain significance"):
                    # Count cases if available
                    num_cases = var.get("num_cases", min_cases)

                    if num_cases >= min_cases:
                        reports.append({
                            "pmid": pmid,
                            "variant": var,
                            "significance": significance,
                            "num_cases": num_cases,
                        })

        return reports

    def find_benign_variant_reports(
        self,
        mentions: Dict[str, List[Dict[str, Any]]],
    ) -> List[Dict[str, Any]]:
        """
        Find variant mentions classified as benign.

        Used for BP6 (reputable source).

        Args:
            mentions: Dict from extract_variant_mentions

        Returns:
            List of benign variant reports
        """
        reports = []

        for pmid, variants in mentions.items():
            for var in variants:
                significance = var.get("clinical_significance", "").lower()

                if significance in ("benign", "likely benign"):
                    reports.append({
                        "pmid": pmid,
                        "variant": var,
                        "significance": significance,
                    })

        return reports
