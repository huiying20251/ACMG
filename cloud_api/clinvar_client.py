#!/usr/bin/env python3
"""
ClinVar E-utilities API client.
Provides access to ClinVar data via NCBI E-utilities (esearch, efetch).
"""

import logging
from typing import Optional, Dict, Any, List

from cloud_api.base_client import BaseCloudClient

logger = logging.getLogger(__name__)


class ClinVarEutilsClient(BaseCloudClient):
    """
    ClinVar E-utilities API client.
    Uses NCBI E-utilities: esearch.fcgi for search, efetch.fcgi for retrieval.
    """

    ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    def __init__(self, timeout: int = 30, max_retries: int = 3):
        """
        Initialize ClinVar E-utilities client.

        Args:
            timeout: Request timeout in seconds
            max_retries: Maximum number of retry attempts
        """
        super().__init__(
            base_url=self.EFETCH_URL,
            timeout=timeout,
            max_retries=max_retries,
        )
        self.search_base_url = self.ESEARCH_URL

    def search_variants(
        self,
        gene_name: Optional[str] = None,
        protein_change: Optional[str] = None,
        clinical_significance: Optional[str] = None,
        max_results: int = 100,
    ) -> List[Dict[str, Any]]:
        """
        Search ClinVar for variants matching criteria.

        Args:
            gene_name: Gene name (e.g., "BRCA1")
            protein_change: Protein change (e.g., "p.Gly12Val")
            clinical_significance: Clinical significance filter
            max_results: Maximum number of results to return

        Returns:
            List of variant records from ClinVar
        """
        params = {
            "db": "clinvar",
            "retmode": "json",
            "retmax": max_results,
            "sort": "relevance",
        }

        # Build search term
        term_parts = []
        if gene_name:
            term_parts.append(f'{gene_name}[Gene Name]')
        if protein_change:
            term_parts.append(f'{protein_change}[Protein Change]')
        if clinical_significance:
            term_parts.append(f'{clinical_significance}[Clinical Significance]')

        if term_parts:
            params["term"] = " AND ".join(term_parts)
        else:
            params["term"] = "variant"

        search_response = self._make_request(
            method="GET",
            endpoint=self.search_base_url,
            params=params,
        )

        if not search_response.success:
            logger.error(f"ClinVar search failed: {search_response.error}")
            return []

        id_list = search_response.data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return []

        return self._fetch_variants_by_ids(id_list)

    def _fetch_variants_by_ids(
        self,
        variant_ids: List[str],
        retmode: str = "json",
    ) -> List[Dict[str, Any]]:
        """
        Fetch variant details by ClinVar IDs.

        Args:
            variant_ids: List of ClinVar variant IDs
            retmode: Return mode (json or xml)

        Returns:
            List of variant records
        """
        if not variant_ids:
            return []

        params = {
            "db": "clinvar",
            "id": ",".join(variant_ids),
            "retmode": retmode,
        }

        response = self.get(
            endpoint=self.EFETCH_URL,
            params=params,
        )

        if not response.success:
            logger.error(f"ClinVar fetch failed: {response.error}")
            return []

        if retmode == "json":
            return response.data.get("result", {}).values()
        return []

    def get_variant_details(self, variant_id: str) -> Optional[Dict[str, Any]]:
        """
        Get details for a specific ClinVar variant.

        Args:
            variant_id: ClinVar variant ID

        Returns:
            Variant details or None if not found
        """
        variants = self._fetch_variants_by_ids([variant_id])
        return variants[0] if variants else None

    def get_pathogenic_variants_in_region(
        self,
        gene_name: str,
        chr: str,
        start: int,
        end: int,
        window_bp: int = 0,
    ) -> List[Dict[str, Any]]:
        """
        Get pathogenic/likely pathogenic variants in a genomic region.

        Args:
            gene_name: Gene name
            chr: Chromosome (e.g., "17")
            start: Start position
            end: End position
            window_bp: Additional window in bp to expand search region

        Returns:
            List of pathogenic ClinVar variants in the region
        """
        expanded_start = start - window_bp
        expanded_end = end + window_bp

        # Search with location filter and pathogenic significance
        params = {
            "db": "clinvar",
            "retmode": "json",
            "retmax": 500,
            "term": f'{gene_name}[Gene Name] AND ({chr}[Chromosome] AND {expanded_start}[Base Position] : {expanded_end}[Base Position]) AND (pathogenic[Clinical Significance] OR likely_pathogenic[Clinical Significance])',
        }

        search_response = self._make_request(
            method="GET",
            endpoint=self.search_base_url,
            params=params,
        )

        if not search_response.success:
            logger.error(f"ClinVar region search failed: {search_response.error}")
            return []

        id_list = search_response.data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return []

        return self._fetch_variants_by_ids(id_list)

    def count_pathogenic_in_gene(
        self,
        gene_name: str,
        domain_start: Optional[int] = None,
        domain_end: Optional[int] = None,
    ) -> int:
        """
        Count pathogenic/likely pathogenic variants in a gene or domain.

        Args:
            gene_name: Gene name
            domain_start: Optional domain start position
            domain_end: Optional domain end position

        Returns:
            Count of pathogenic variants
        """
        term = f'{gene_name}[Gene Name]'
        if domain_start and domain_end:
            term += f' AND ({domain_start}[Base Position] : {domain_end}[Base Position])'
        term += ' AND (pathogenic[Clinical Significance] OR likely_pathogenic[Clinical Significance])'

        params = {
            "db": "clinvar",
            "retmode": "json",
            "term": term,
            "retmax": 0,
        }

        response = self._make_request(
            method="GET",
            endpoint=self.search_base_url,
            params=params,
        )

        if not response.success:
            return 0

        count = response.data.get("esearchresult", {}).get("count", 0)
        return int(count) if count else 0

    def check_variant_significance(
        self,
        gene_name: str,
        protein_change: str,
    ) -> Optional[str]:
        """
        Check if a variant with given protein change exists in ClinVar.

        Args:
            gene_name: Gene name
            protein_change: Protein change (e.g., "p.Gly12Val")

        Returns:
            Clinical significance string or None
        """
        variants = self.search_variants(
            gene_name=gene_name,
            protein_change=protein_change,
            max_results=10,
        )

        for variant in variants:
            rcv = variant.get("rcvaccession", [])
            if rcv:
                significance = rcv[0].get("clinical_significance", {}).get("description")
                if significance:
                    return significance

        return None

    def search_missense_by_gene_and_position(
        self,
        gene_name: str,
        protein_position: int,
        amino_acid_ref: Optional[str] = None,
        amino_acid_alt: Optional[str] = None,
        clinical_significance: Optional[str] = None,
        max_results: int = 500,
    ) -> List[Dict[str, Any]]:
        """
        Search ClinVar for missense variants at a specific codon position.

        This is used for PS1 (same AA change) and PM5 (different AA change at same position).

        Args:
            gene_name: Gene name (e.g., "BRCA1")
            protein_position: Protein position (amino acid number, e.g., 185)
            amino_acid_ref: Reference amino acid (e.g., "Gly" for p.Gly185)
            amino_acid_alt: Alternate amino acid (e.g., "Val" for p.185Val)
            clinical_significance: Filter by clinical significance (default: pathogenic/likely pathogenic)
            max_results: Maximum number of results

        Returns:
            List of ClinVar variant records at the specified position
        """
        # Build search term
        # ClinVar E-utilities supports:
        # - Gene name: GENE[Gene Name]
        # - Protein change: p.XXXnnn or p.nnnXXX format
        # - Variant type: variant type[Variant Type]

        # First, try protein position search with gene
        term_parts = [f"{gene_name}[Gene Name]"]

        # Add protein position search
        # ClinVar uses format like "p.Gly185" or "185"
        if amino_acid_ref:
            term_parts.append(f"p.{amino_acid_ref}{protein_position}[Protein Change]")
        else:
            # Search for any change at this position
            term_parts.append(f"{protein_position}[Protein Position]")

        # Add variant type filter for missense
        term_parts.append("missense[Variant Type]")

        # Default to pathogenic if not specified
        if not clinical_significance:
            term_parts.append("(pathogenic[Clinical Significance] OR likely_pathogenic[Clinical Significance])")
        else:
            term_parts.append(f"{clinical_significance}[Clinical Significance]")

        term = " AND ".join(term_parts)

        params = {
            "db": "clinvar",
            "retmode": "json",
            "retmax": max_results,
            "sort": "relevance",
            "term": term,
        }

        search_response = self._make_request(
            method="GET",
            endpoint=self.search_base_url,
            params=params,
        )

        if not search_response.success:
            logger.error(f"ClinVar missense search failed: {search_response.error}")
            return []

        id_list = search_response.data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return []

        variants = self._fetch_variants_by_ids(id_list)

        # Post-fetch filtering: ensure only missense variants are returned
        # Filter out frameshift, nonsense, splice, etc.
        missense_variants = self._filter_missense_only(variants)

        return missense_variants

    def _filter_missense_only(self, variants: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Filter variants to only include missense variants.

        ClinVar VCF has a "Variant_type" field, but API responses may vary.
        This method filters based on available variant type information.

        Args:
            variants: List of ClinVar variant records

        Returns:
            List of missense variants only
        """
        missense_variants = []

        # Variant types to exclude (non-missense)
        exclude_types = {
            "frameshift", "frameshift_variant", "frameshift deletion",
            "frameshift insertion", "stop_gained", "nonsense", "stop_lost",
            "start_lost", "splice_site", "splice_region", "splice_acceptor",
            "splice_donor", "synonymous", "intron", "UTR", "upstream",
            "downstream", "intergenic", "inframe", "inframe_deletion",
            "inframe_insertion", "inframe_indel", "deletion", "insertion",
            "duplication", "structural_variant", "copy_number",
        }

        for variant in variants:
            # Check various fields that might contain variant type info
            variant_type = None

            # Try to get variant type from different possible fields
            if "variant_type" in variant:
                variant_type = str(variant["variant_type"]).lower()
            elif "type" in variant:
                variant_type = str(variant["type"]).lower()
            elif "rcvaccession" in variant:
                # RCV accessions may contain variant type info
                rcv = variant.get("rcvaccession", [{}])[0] if variant.get("rcvaccession") else {}
                if "variant_type" in rcv:
                    variant_type = str(rcv["variant_type"]).lower()

            # If no variant type info, include the variant (conservative approach)
            # The API search term already filters for missense
            if variant_type is None:
                missense_variants.append(variant)
            elif variant_type not in exclude_types:
                missense_variants.append(variant)
            # If variant_type contains "missense", include it
            elif "missense" in variant_type:
                missense_variants.append(variant)

        return missense_variants

    def search_pathogenic_at_codon(
        self,
        gene_name: str,
        protein_position: int,
        chromosome: str,
        genomic_start: int,
        genomic_end: int,
        window_bp: int = 50,
    ) -> List[Dict[str, Any]]:
        """
        Search ClinVar for pathogenic/likely pathogenic variants at a codon region.

        Falls back to genomic region search if protein position search yields no results.

        Args:
            gene_name: Gene name
            protein_position: Protein position (amino acid number)
            chromosome: Chromosome (e.g., "17")
            genomic_start: Genomic start position of codon
            genomic_end: Genomic end position of codon
            window_bp: Window size in bp to expand search

        Returns:
            List of pathogenic ClinVar variants in the region
        """
        # First try protein position search
        variants = self.search_missense_by_gene_and_position(
            gene_name=gene_name,
            protein_position=protein_position,
            max_results=500,
        )

        if variants:
            return variants

        # Fallback: search by genomic region
        logger.info(f"No results from protein position search, falling back to genomic region search for {gene_name} position {protein_position}")

        return self.get_pathogenic_variants_in_region(
            gene_name=gene_name,
            chr=chromosome,
            start=genomic_start,
            end=genomic_end,
            window_bp=window_bp,
        )

    def search_splicing_by_genomic_position(
        self,
        gene_name: str,
        chromosome: str,
        genomic_start: int,
        genomic_end: int,
        variant_type: Optional[str] = None,
        clinical_significance: Optional[str] = None,
        window_bp: int = 5,
        max_results: int = 500,
    ) -> List[Dict[str, Any]]:
        """
        Search ClinVar for splicing variants at a specific genomic position.

        Used for PS1 (same nucleotide) and splice site variant searches.

        Args:
            gene_name: Gene name
            chromosome: Chromosome (e.g., "17")
            genomic_start: Genomic start position
            genomic_end: Genomic end position
            variant_type: Type of variant (e.g., "splice_site", "splice_region")
            clinical_significance: Filter by clinical significance
            window_bp: Window size in bp around the position
            max_results: Maximum number of results

        Returns:
            List of ClinVar variants at the splice site
        """
        # Expand search window
        expanded_start = genomic_start - window_bp
        expanded_end = genomic_end + window_bp

        # Build search term
        term_parts = [f"{gene_name}[Gene Name]"]

        # Add location filter
        term_parts.append(f"({chromosome}[Chromosome] AND {expanded_start}[Base Position] : {expanded_end}[Base Position])")

        # Add variant type filter for splicing
        if variant_type:
            term_parts.append(f"{variant_type}[Variant Type]")
        else:
            term_parts.append("(splice_site[Variant Type] OR splice_region[Variant Type] OR splicing[Variant Type])")

        # Add clinical significance filter
        if not clinical_significance:
            term_parts.append("(pathogenic[Clinical Significance] OR likely_pathogenic[Clinical Significance])")
        else:
            term_parts.append(f"{clinical_significance}[Clinical Significance]")

        term = " AND ".join(term_parts)

        params = {
            "db": "clinvar",
            "retmode": "json",
            "retmax": max_results,
            "sort": "relevance",
            "term": term,
        }

        search_response = self._make_request(
            method="GET",
            endpoint=self.search_base_url,
            params=params,
        )

        if not search_response.success:
            logger.error(f"ClinVar splicing search failed: {search_response.error}")
            return []

        id_list = search_response.data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return []

        variants = self._fetch_variants_by_ids(id_list)

        # Filter to only splicing-related variants
        splicing_variants = self._filter_splicing_only(variants)

        return splicing_variants

    def _filter_splicing_only(self, variants: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Filter variants to only include splicing-related variants.

        Args:
            variants: List of ClinVar variant records

        Returns:
            List of splicing variants only
        """
        splicing_variants = []

        # Splicing-related variant types to include
        include_types = {
            "splice_site", "splice_region", "splice_acceptor", "splice_donor",
            "splicing", "splice", "intron", "splice_site_variant",
            "splice_region_variant",
        }

        for variant in variants:
            variant_type = None

            # Try to get variant type from different possible fields
            if "variant_type" in variant:
                variant_type = str(variant["variant_type"]).lower()
            elif "type" in variant:
                variant_type = str(variant["type"]).lower()
            elif "rcvaccession" in variant:
                rcv = variant.get("rcvaccession", [{}])[0] if variant.get("rcvaccession") else {}
                if "variant_type" in rcv:
                    variant_type = str(rcv["variant_type"]).lower()

            # If no variant type info, include the variant (conservative)
            if variant_type is None:
                splicing_variants.append(variant)
            else:
                # Check if any include_type is in the variant_type string
                if any(t in variant_type for t in include_types):
                    splicing_variants.append(variant)

        return splicing_variants