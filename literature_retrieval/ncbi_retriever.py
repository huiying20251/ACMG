#!/usr/bin/env python3
"""
NCBI Literature Retriever

Retrieves literature from NCBI PubMed based on variant information.
Uses the VariantNormalizer.build_search_queries() to generate search keywords.

Reference: healthfutures-evagg/lib/evagg/ref/ncbi.py
"""

import logging
import urllib.parse
from typing import List, Optional, Dict, Any, Set
from dataclasses import dataclass

import requests
from defusedxml import ElementTree

from normalizer import VariantNormalizer, VariantInfo
from .literature_utils import Article

logger = logging.getLogger(__name__)


class NcbiSettings:
    """NCBI API settings."""
    EUTILS_HOST = "https://eutils.ncbi.nlm.nih.gov"
    ESEARCH_URL = "/entrez/eutils/esearch.fcgi"
    EFETCH_URL = "/entrez/eutils/efetch.fcgi"

    def __init__(self, api_key: Optional[str] = None, email: str = "biomedcomp@microsoft.com"):
        self.api_key = api_key
        self.email = email

    def get_key_string(self) -> str:
        """Build NCBI API key parameter string."""
        params = []
        if self.email:
            params.append(f"email={urllib.parse.quote(self.email)}")
        if self.api_key:
            params.append(f"api_key={self.api_key}")
        return "&".join(params) if params else ""


@dataclass
class PubMedArticle:
    """PubMed article data extracted from XML."""
    pmid: str
    title: str
    abstract: str
    journal: str
    first_author: str
    pub_year: str
    doi: str
    pmcid: str


class NCBILiteratureRetriever:
    """
    NCBI PubMed literature retriever for variants.

    Uses VariantNormalizer.build_search_queries() to generate multiple
    search keywords from a normalized VariantInfo object.
    """

    ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    PMC_OA_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi"

    def __init__(
        self,
        normalizer: Optional[VariantNormalizer] = None,
        api_key: Optional[str] = None,
        email: str = "biomedcomp@microsoft.com",
        timeout: int = 30,
    ):
        self.normalizer = normalizer or VariantNormalizer()
        self.settings = NcbiSettings(api_key=api_key, email=email)
        self.timeout = timeout

    def search_variant_literature(
        self,
        variant_info: VariantInfo,
        max_results_per_query: int = 20,
        max_total: int = 100,
    ) -> List[Article]:
        """
        Search PubMed for variant-related literature.

        Uses build_search_queries() to generate multiple query strategies:
        1. rsID (most specific)
        2. Gene + cDNA change
        3. Gene + protein change
        4. Chromosome position (chr:pos:ref:alt)

        Each query returns up to max_results_per_query PMIDs,
        total limited to max_total.

        Args:
            variant_info: Normalized variant information
            max_results_per_query: Max results per query strategy
            max_total: Maximum total results to return

        Returns:
            List of Article objects
        """
        # Generate all search queries from normalized variant
        queries = self.normalizer.build_search_queries(variant_info)

        logger.info(f"Generated {len(queries)} search queries for variant")

        # Execute searches and collect unique PMIDs
        all_pmids: Set[str] = set()
        for query in queries:
            try:
                pmids = self._search_pubmed(query, max_results_per_query)
                all_pmids.update(pmids)
                logger.debug(f"Query '{query}' returned {len(pmids)} PMIDs")
            except Exception as e:
                logger.warning(f"Search failed for query '{query}': {e}")

            if len(all_pmids) >= max_total:
                break

        logger.info(f"Total unique PMIDs collected: {len(all_pmids)}")

        # Limit to max_total
        pmids_to_fetch = list(all_pmids)[:max_total]

        # Fetch article details
        articles = self._fetch_articles(pmids_to_fetch)

        return articles

    def _search_pubmed(self, query: str, max_results: int) -> List[str]:
        """
        Execute PubMed search and return PMIDs.

        Args:
            query: Search query string
            max_results: Maximum number of results to return

        Returns:
            List of PMID strings
        """
        params = {
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "sort": "relevance",
            "tool": "variant_classification",
        }

        key_string = self.settings.get_key_string()
        url = self.EFETCH_URL + "?" + "&".join([f"{k}={v}" for k, v in params.items()])
        if key_string:
            url += "&" + key_string

        # Actually need to use ESEARCH first to get IDs
        search_url = self.ESEARCH_URL + "?" + "&".join([f"{k}={v}" for k, v in params.items()])
        if key_string:
            search_url += "&" + key_string

        try:
            response = requests.get(
                f"{self.settings.EUTILS_HOST}{search_url}",
                timeout=self.timeout,
            )
            response.raise_for_status()

            root = ElementTree.fromstring(response.content)
            pmids = [
                pmid.text
                for pmid in root.findall(".//IdList/Id")
                if pmid.text
            ]
            return pmids

        except Exception as e:
            logger.warning(f"NCBI search failed for '{query}': {e}")
            return []

    def _fetch_articles(self, pmids: List[str]) -> List[Article]:
        """
        Fetch article details for a list of PMIDs.

        Args:
            pmids: List of PubMed IDs

        Returns:
            List of Article objects
        """
        articles = []

        if not pmids:
            return articles

        # NCBI allows batch fetching with comma-separated IDs
        id_list = ",".join(pmids[:200])  # NCBI limit is typically 200

        params = {
            "db": "pubmed",
            "id": id_list,
            "retmode": "xml",
            "rettype": "abstract",
        }

        key_string = self.settings.get_key_string()
        url = self.EFETCH_URL + "?" + "&".join([f"{k}={v}" for k, v in params.items()])
        if key_string:
            url += "&" + key_string

        try:
            response = requests.get(
                f"{self.settings.EUTILS_HOST}{url}",
                timeout=self.timeout,
            )
            response.raise_for_status()

            root = ElementTree.fromstring(response.content)

            for article_elem in root.findall(".//PubmedArticle"):
                try:
                    article = self._parse_article(article_elem)
                    if article:
                        articles.append(article)
                except Exception as e:
                    logger.warning(f"Failed to parse article: {e}")
                    continue

        except Exception as e:
            logger.error(f"Failed to fetch articles: {e}")

        return articles

    def _parse_article(self, article_elem) -> Optional[Article]:
        """Parse a PubMed article XML element into an Article object."""
        try:
            # Extract PMID
            pmid_elem = article_elem.find(".//PMID")
            if pmid_elem is None or not pmid_elem.text:
                return None
            pmid = pmid_elem.text

            # Extract title
            title_elem = article_elem.find(".//ArticleTitle")
            title = self._get_text(title_elem) if title_elem is not None else ""

            # Extract abstract
            abstract_elem = article_elem.find(".//Abstract")
            abstract = self._get_text(abstract_elem) if abstract_elem is not None else ""

            # Extract journal
            journal_elem = article_elem.find(".//Journal/ISOAbbreviation")
            journal = self._get_text(journal_elem) if journal_elem is not None else ""

            # Extract first author
            author_elem = article_elem.find(".//AuthorList/Author[1]/LastName")
            first_author = self._get_text(author_elem) if author_elem is not None else ""

            # Extract publication year
            year_elem = article_elem.find(".//Journal/JournalIssue/PubDate/Year")
            pub_year_str = self._get_text(year_elem) if year_elem is not None else None
            pub_year = int(pub_year_str) if pub_year_str and pub_year_str.isdigit() else None

            # Extract DOI
            doi_elem = article_elem.find(".//ArticleId[@IdType='doi']")
            doi = self._get_text(doi_elem) if doi_elem is not None else None

            # Extract PMCID
            pmcid_elem = article_elem.find(".//ArticleId[@IdType='pmc']")
            pmcid = self._get_text(pmcid_elem) if pmcid_elem is not None else None

            # Check PMC OA status if PMCID exists
            can_access = False
            is_oa = False
            license = None
            if pmcid:
                try:
                    can_access, is_oa, license = self._check_pmc_oa(pmcid)
                except Exception:
                    pass

            # Build citation
            citation = f"{first_author} ({pub_year}) {journal}" if first_author and pub_year and journal else None
            link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None

            return Article(
                pmid=pmid,
                title=title,
                abstract=abstract,
                journal=journal,
                first_author=first_author,
                pub_year=pub_year,
                doi=doi,
                pmcid=pmcid,
                citation=citation,
                link=link,
                can_access=can_access,
                is_oa=is_oa,
                license=license,
            )

        except Exception as e:
            logger.warning(f"Error parsing article: {e}")
            return None

    def _get_text(self, elem) -> str:
        """Get all text content from an XML element."""
        if elem is None:
            return ""
        return " ".join(elem.itertext()).strip()

    def _check_pmc_oa(self, pmcid: str) -> tuple:
        """
        Check if a PMC article is open access.

        Returns:
            (can_access, is_oa, license)
        """
        if not pmcid:
            return False, False, None

        try:
            # Remove PMC prefix if present
            pmcid_clean = pmcid.replace("PMC", "")
            url = f"{self.PMC_OA_URL}?id={pmcid_clean}"

            response = requests.get(url, timeout=10)
            if response.status_code != 200:
                return False, False, None

            root = ElementTree.fromstring(response.content)
            record = root.find(f"records/record[@id='PMC{pmcid_clean}']")

            if record is None:
                return False, False, None

            license = record.attrib.get("license", "unknown")
            is_oa = True
            can_access = True

            # Check for no-derivatives license
            if "-ND" in license:
                can_access = False

            return can_access, is_oa, license

        except Exception as e:
            logger.debug(f"PMC OA check failed for {pmcid}: {e}")
            return False, False, None

    def fetch_fulltext(self, article: Article) -> Optional[str]:
        """
        Fetch fulltext XML for an article from PMC.

        Args:
            article: Article with pmcid

        Returns:
            Fulltext XML string or None
        """
        if not article.pmcid or not article.can_access:
            return None

        try:
            pmcid_clean = article.pmcid.replace("PMC", "")
            url = f"https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/{pmcid_clean}/ascii"

            response = requests.get(url, timeout=self.timeout)
            if response.status_code == 200:
                return response.text

        except Exception as e:
            logger.warning(f"Failed to fetch fulltext for {article.pmid}: {e}")

        return None
