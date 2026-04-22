#!/usr/bin/env python3
"""
Literature Evidence Cache Database

Caches extracted literature evidence to SQLite for reuse:
- Articles (PubMed metadata)
- Individual variant observations (per-patient data for PS2/PM3/PP1/PS4)
- Evidence assessment results

Schema:
- articles: PubMed article metadata
- individual_observations: Per-patient variant observations
- evidence_results: ACMG evidence assessment results
"""

import json
import logging
import sqlite3
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple

from .literature_utils import (
    Article,
    IndividualVariantObservation,
    EvidenceResult,
    CaseReport,
)

logger = logging.getLogger(__name__)


# Default database path
DEFAULT_DB_PATH = Path(__file__).parent.parent / "literature_evidence.db"


@dataclass
class EvidenceCacheConfig:
    """Configuration for evidence cache."""
    db_path: Path = DEFAULT_DB_PATH
    # Cache validity in days (None = never expires)
    max_age_days: Optional[int] = 90
    # Enable/disable cache
    enabled: bool = True


class EvidenceCacheDB:
    """
    SQLite database for caching literature evidence.

    Usage:
        cache = EvidenceCacheDB()
        cache.init_db()  # Create tables if not exist

        # Store evidence
        cache.store_article(article)
        cache.store_observation(observation, variant_id)
        cache.store_evidence_result(result, variant_id, gene)

        # Retrieve evidence
        evidence = cache.get_evidence_result(variant_id, gene, rule)
    """

    def __init__(self, db_path: Optional[Path] = None):
        self.db_path = db_path or DEFAULT_DB_PATH
        self._conn: Optional[sqlite3.Connection] = None

    def connect(self) -> sqlite3.Connection:
        """Get database connection."""
        if self._conn is None:
            self._conn = sqlite3.connect(str(self.db_path))
            self._conn.row_factory = sqlite3.Row
        return self._conn

    def close(self):
        """Close database connection."""
        if self._conn:
            self._conn.close()
            self._conn = None

    def init_db(self):
        """Initialize database schema."""
        conn = self.connect()
        cursor = conn.cursor()

        # Articles table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS articles (
                pmid TEXT PRIMARY KEY,
                title TEXT,
                abstract TEXT,
                journal TEXT,
                first_author TEXT,
                pub_year INTEGER,
                doi TEXT,
                pmcid TEXT,
                citation TEXT,
                link TEXT,
                can_access INTEGER DEFAULT 0,
                is_oa INTEGER DEFAULT 0,
                license TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)

        # Individual observations table (per-patient data)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS individual_observations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pmid TEXT NOT NULL,
                variant_id TEXT NOT NULL,
                gene TEXT NOT NULL,
                individual_id TEXT,
                variant_description TEXT,
                variant_inheritance TEXT,
                parental_testing INTEGER DEFAULT 0,
                zygosity TEXT,
                clinvar_status TEXT,
                clinvar_significance TEXT,
                phenotype TEXT,
                hpo_terms TEXT,
                segregation_data TEXT,
                confidence TEXT DEFAULT 'medium',
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (pmid) REFERENCES articles(pmid)
            )
        """)

        # Evidence results table (PS2/PM3/PP1/PS4 assessments)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS evidence_results (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                variant_id TEXT NOT NULL,
                gene TEXT NOT NULL,
                rule TEXT NOT NULL,
                applicable INTEGER DEFAULT 0,
                strength TEXT,
                num_cases INTEGER,
                num_controls INTEGER,
                num_segregations INTEGER,
                inheritance_pattern TEXT,
                disease_category TEXT,
                phenotype_consistency TEXT,
                case_control_enrichment TEXT,
                pmids TEXT,
                comment TEXT,
                confidence TEXT DEFAULT 'medium',
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                UNIQUE(variant_id, gene, rule)
            )
        """)

        # Variant literature summary table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS variant_literature_summary (
                variant_id TEXT PRIMARY KEY,
                gene TEXT NOT NULL,
                total_articles INTEGER DEFAULT 0,
                num_case_reports INTEGER DEFAULT 0,
                num_functional_studies INTEGER DEFAULT 0,
                num_individual_observations INTEGER DEFAULT 0,
                has_de_novo_evidence INTEGER DEFAULT 0,
                has_cosegregation_evidence INTEGER DEFAULT 0,
                has_case_control_evidence INTEGER DEFAULT 0,
                has_phenotype_evidence INTEGER DEFAULT 0,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)

        # Indexes
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_obs_variant ON individual_observations(variant_id)
        """)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_obs_gene ON individual_observations(gene)
        """)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_obs_pmid ON individual_observations(pmid)
        """)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_evidence_variant ON evidence_results(variant_id)
        """)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_evidence_gene ON evidence_results(gene)
        """)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_evidence_rule ON evidence_results(rule)
        """)

        conn.commit()
        logger.info(f"Initialized evidence cache database: {self.db_path}")

    # ==================== Article Methods ====================

    def store_article(self, article: Article) -> bool:
        """
        Store an article to cache.

        Args:
            article: Article dataclass

        Returns:
            True if stored successfully
        """
        conn = self.connect()
        cursor = conn.cursor()

        try:
            cursor.execute("""
                INSERT OR REPLACE INTO articles
                (pmid, title, abstract, journal, first_author, pub_year, doi, pmcid, citation, link, can_access, is_oa, license)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                article.pmid,
                article.title,
                article.abstract,
                article.journal,
                article.first_author,
                article.pub_year,
                article.doi,
                article.pmcid,
                article.citation,
                article.link,
                int(article.can_access),
                int(article.is_oa),
                article.license,
            ))
            conn.commit()
            return True
        except sqlite3.Error as e:
            logger.error(f"Error storing article {article.pmid}: {e}")
            return False

    def get_article(self, pmid: str) -> Optional[Article]:
        """Get article by PMID."""
        conn = self.connect()
        cursor = conn.cursor()

        cursor.execute("SELECT * FROM articles WHERE pmid = ?", (pmid,))
        row = cursor.fetchone()

        if row:
            return Article(
                pmid=row["pmid"],
                title=row["title"],
                abstract=row["abstract"],
                journal=row["journal"],
                first_author=row["first_author"],
                pub_year=row["pub_year"],
                doi=row["doi"],
                pmcid=row["pmcid"],
                citation=row["citation"],
                link=row["link"],
                can_access=bool(row["can_access"]),
                is_oa=bool(row["is_oa"]),
                license=row["license"],
            )
        return None

    def store_articles_batch(self, articles: List[Article]) -> int:
        """Store multiple articles in a batch."""
        count = 0
        for article in articles:
            if self.store_article(article):
                count += 1
        return count

    # ==================== Observation Methods ====================

    def store_observation(
        self,
        observation: IndividualVariantObservation,
        variant_id: str,
    ) -> bool:
        """
        Store an individual variant observation.

        Args:
            observation: IndividualVariantObservation dataclass
            variant_id: Variant identifier (rsID or VCF format)

        Returns:
            True if stored successfully
        """
        conn = self.connect()
        cursor = conn.cursor()

        try:
            cursor.execute("""
                INSERT INTO individual_observations
                (pmid, variant_id, gene, individual_id, variant_description,
                 variant_inheritance, parental_testing, zygosity,
                 clinvar_status, clinvar_significance, phenotype, hpo_terms,
                 segregation_data, confidence)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                observation.pmid,
                variant_id,
                observation.gene,
                observation.individual_id,
                observation.variant_description,
                observation.variant_inheritance,
                int(observation.parental_testing),
                observation.zygosity,
                observation.clinvar_status,
                observation.clinvar_significance,
                observation.phenotype,
                json.dumps(observation.hpo_terms),
                json.dumps(observation.segregation_data) if observation.segregation_data else None,
                observation.confidence,
            ))
            conn.commit()
            return True
        except sqlite3.Error as e:
            logger.error(f"Error storing observation: {e}")
            return False

    def get_observations(
        self,
        variant_id: Optional[str] = None,
        gene: Optional[str] = None,
        pmid: Optional[str] = None,
    ) -> List[IndividualVariantObservation]:
        """
        Get individual observations filtered by variant/gene/pmid.

        Args:
            variant_id: Filter by variant ID
            gene: Filter by gene
            pmid: Filter by PMID

        Returns:
            List of IndividualVariantObservation
        """
        conn = self.connect()
        cursor = conn.cursor()

        query = "SELECT * FROM individual_observations WHERE 1=1"
        params = []

        if variant_id:
            query += " AND variant_id = ?"
            params.append(variant_id)
        if gene:
            query += " AND gene = ?"
            params.append(gene)
        if pmid:
            query += " AND pmid = ?"
            params.append(pmid)

        cursor.execute(query, params)
        rows = cursor.fetchall()

        observations = []
        for row in rows:
            observations.append(IndividualVariantObservation(
                pmid=row["pmid"],
                individual_id=row["individual_id"],
                variant_description=row["variant_description"],
                gene=row["gene"],
                variant_inheritance=row["variant_inheritance"],
                parental_testing=bool(row["parental_testing"]),
                zygosity=row["zygosity"],
                clinvar_status=row["clinvar_status"],
                clinvar_significance=row["clinvar_significance"],
                phenotype=row["phenotype"],
                hpo_terms=json.loads(row["hpo_terms"]) if row["hpo_terms"] else [],
                segregation_data=json.loads(row["segregation_data"]) if row["segregation_data"] else None,
                confidence=row["confidence"],
            ))

        return observations

    def store_observations_batch(
        self,
        observations: List[IndividualVariantObservation],
        variant_id: str,
    ) -> int:
        """Store multiple observations in a batch."""
        count = 0
        for obs in observations:
            if self.store_observation(obs, variant_id):
                count += 1
        return count

    # ==================== Evidence Result Methods ====================

    def store_evidence_result(
        self,
        result: EvidenceResult,
        variant_id: str,
        gene: str,
    ) -> bool:
        """
        Store an evidence assessment result.

        Args:
            result: EvidenceResult dataclass
            variant_id: Variant identifier
            gene: Gene symbol

        Returns:
            True if stored successfully
        """
        conn = self.connect()
        cursor = conn.cursor()

        try:
            cursor.execute("""
                INSERT OR REPLACE INTO evidence_results
                (variant_id, gene, rule, applicable, strength,
                 num_cases, num_controls, num_segregations, inheritance_pattern,
                 disease_category, phenotype_consistency, case_control_enrichment,
                 pmids, comment, confidence)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                variant_id,
                gene,
                result.rule,
                int(result.applicable),
                result.strength,
                result.num_cases,
                result.num_controls,
                getattr(result, 'num_segregations', None),
                getattr(result, 'inheritance_pattern', None),
                getattr(result, 'disease_category', None),
                getattr(result, 'phenotype_consistency', None),
                json.dumps(getattr(result, 'case_control_enrichment', None)),
                json.dumps(result.pmids),
                result.comment,
                result.confidence,
            ))
            conn.commit()
            return True
        except sqlite3.Error as e:
            logger.error(f"Error storing evidence result: {e}")
            return False

    def get_evidence_result(
        self,
        variant_id: str,
        gene: str,
        rule: str,
    ) -> Optional[EvidenceResult]:
        """
        Get cached evidence result for a variant/gene/rule combination.

        Args:
            variant_id: Variant identifier
            gene: Gene symbol
            rule: Rule name (e.g., "PS2", "PM3", "PP1", "PS4")

        Returns:
            EvidenceResult if found and not expired, None otherwise
        """
        conn = self.connect()
        cursor = conn.cursor()

        cursor.execute("""
            SELECT * FROM evidence_results
            WHERE variant_id = ? AND gene = ? AND rule = ?
        """, (variant_id, gene, rule))

        row = cursor.fetchone()
        if not row:
            return None

        # Check if expired
        created_at = datetime.fromisoformat(row["created_at"])
        if self._is_expired(created_at):
            logger.debug(f"Evidence result expired for {variant_id}/{gene}/{rule}")
            return None

        # Parse case_control_enrichment
        cce = None
        if row["case_control_enrichment"]:
            try:
                cce = json.loads(row["case_control_enrichment"])
            except json.JSONDecodeError:
                pass

        return EvidenceResult(
            rule=row["rule"],
            applicable=bool(row["applicable"]),
            strength=row["strength"],
            num_cases=row["num_cases"],
            num_controls=row["num_controls"],
            segregation_data=getattr(row, 'segregation_data', None),
            phenotype=getattr(row, 'phenotype', None),
            pmids=json.loads(row["pmids"]) if row["pmids"] else [],
            confidence=row["confidence"],
            comment=row["comment"],
        )

    def get_all_evidence_results(
        self,
        variant_id: str,
        gene: str,
    ) -> Dict[str, EvidenceResult]:
        """
        Get all cached evidence results for a variant.

        Returns:
            Dict mapping rule name to EvidenceResult
        """
        conn = self.connect()
        cursor = conn.cursor()

        cursor.execute("""
            SELECT * FROM evidence_results
            WHERE variant_id = ? AND gene = ?
        """, (variant_id, gene))

        rows = cursor.fetchall()
        results = {}

        for row in rows:
            if self._is_expired(datetime.fromisoformat(row["created_at"])):
                continue

            results[row["rule"]] = EvidenceResult(
                rule=row["rule"],
                applicable=bool(row["applicable"]),
                strength=row["strength"],
                num_cases=row["num_cases"],
                num_controls=row["num_controls"],
                pmids=json.loads(row["pmids"]) if row["pmids"] else [],
                confidence=row["confidence"],
                comment=row["comment"],
            )

        return results

    # ==================== Variant Summary Methods ====================

    def store_variant_summary(
        self,
        variant_id: str,
        gene: str,
        summary: Dict[str, Any],
    ) -> bool:
        """Store variant literature summary."""
        conn = self.connect()
        cursor = conn.cursor()

        try:
            cursor.execute("""
                INSERT OR REPLACE INTO variant_literature_summary
                (variant_id, gene, total_articles, num_case_reports,
                 num_functional_studies, num_individual_observations,
                 has_de_novo_evidence, has_cosegregation_evidence,
                 has_case_control_evidence, has_phenotype_evidence, updated_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
            """, (
                variant_id,
                gene,
                summary.get("total_articles", 0),
                summary.get("num_case_reports", 0),
                summary.get("num_functional_studies", 0),
                summary.get("num_individual_observations", 0),
                int(summary.get("has_de_novo_evidence", False)),
                int(summary.get("has_cosegregation_evidence", False)),
                int(summary.get("has_case_control_evidence", False)),
                int(summary.get("has_phenotype_evidence", False)),
            ))
            conn.commit()
            return True
        except sqlite3.Error as e:
            logger.error(f"Error storing variant summary: {e}")
            return False

    def get_variant_summary(
        self,
        variant_id: str,
        gene: str,
    ) -> Optional[Dict[str, Any]]:
        """Get variant literature summary if not expired."""
        conn = self.connect()
        cursor = conn.cursor()

        cursor.execute("""
            SELECT * FROM variant_literature_summary
            WHERE variant_id = ? AND gene = ?
        """, (variant_id, gene))

        row = cursor.fetchone()
        if not row:
            return None

        if self._is_expired(datetime.fromisoformat(row["updated_at"])):
            return None

        return {
            "total_articles": row["total_articles"],
            "num_case_reports": row["num_case_reports"],
            "num_functional_studies": row["num_functional_studies"],
            "num_individual_observations": row["num_individual_observations"],
            "has_de_novo_evidence": bool(row["has_de_novo_evidence"]),
            "has_cosegregation_evidence": bool(row["has_cosegregation_evidence"]),
            "has_case_control_evidence": bool(row["has_case_control_evidence"]),
            "has_phenotype_evidence": bool(row["has_phenotype_evidence"]),
        }

    # ==================== Utility Methods ====================

    def _is_expired(self, created_at: datetime) -> bool:
        """Check if a cache entry is expired."""
        if self.config.max_age_days is None:
            return False
        age = datetime.now() - created_at
        return age.days > self.config.max_age_days

    def clear_expired(self) -> int:
        """Clear expired cache entries. Returns number of deleted rows."""
        conn = self.connect()
        cursor = conn.cursor()

        if self.config.max_age_days is None:
            return 0

        cursor.execute("""
            DELETE FROM evidence_results
            WHERE datetime(created_at) < datetime('now', '-' || ? || ' days')
        """, (self.config.max_age_days,))

        cursor.execute("""
            DELETE FROM variant_literature_summary
            WHERE datetime(updated_at) < datetime('now', '-' || ? || ' days')
        """, (self.config.max_age_days,))

        deleted = cursor.rowcount
        conn.commit()
        return deleted

    def clear_all(self):
        """Clear all cached data."""
        conn = self.connect()
        cursor = conn.cursor()

        cursor.execute("DELETE FROM individual_observations")
        cursor.execute("DELETE FROM evidence_results")
        cursor.execute("DELETE FROM variant_literature_summary")
        cursor.execute("DELETE FROM articles")

        conn.commit()

    def get_stats(self) -> Dict[str, int]:
        """Get cache statistics."""
        conn = self.connect()
        cursor = conn.cursor()

        stats = {}
        tables = ["articles", "individual_observations", "evidence_results", "variant_literature_summary"]

        for table in tables:
            cursor.execute(f"SELECT COUNT(*) FROM {table}")
            stats[table] = cursor.fetchone()[0]

        return stats

    def __enter__(self):
        """Context manager entry."""
        self.init_db()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()


# Global config
config = EvidenceCacheConfig()


def get_cache() -> EvidenceCacheDB:
    """Get evidence cache instance."""
    return EvidenceCacheDB(db_path=config.db_path)


# ==================== Convenience Functions ====================

def cache_evidence(
    variant_id: str,
    gene: str,
    literature_data: Dict[str, Any],
) -> bool:
    """
    Cache complete literature evidence for a variant.

    Args:
        variant_id: Variant identifier (rsID or VCF format)
        gene: Gene symbol
        literature_data: Dict containing articles, observations, and results

    Returns:
        True if all data cached successfully
    """
    with get_cache() as cache:
        success = True

        # Store articles
        articles = literature_data.get("articles", [])
        if articles:
            count = cache.store_articles_batch(articles)
            logger.info(f"Cached {count}/{len(articles)} articles for {variant_id}")

        # Store individual observations
        observations = literature_data.get("individual_observations", [])
        if observations:
            count = cache.store_observations_batch(observations, variant_id)
            logger.info(f"Cached {count}/{len(observations)} observations for {variant_id}")

        # Store evidence results
        results = literature_data.get("evidence_results", {})
        for rule, result in results.items():
            if result and isinstance(result, EvidenceResult):
                if cache.store_evidence_result(result, variant_id, gene):
                    logger.debug(f"Cached {rule} result for {variant_id}/{gene}")
                else:
                    success = False

        # Store summary
        summary = literature_data.get("summary", {})
        if summary:
            cache.store_variant_summary(variant_id, gene, summary)

        return success


def get_cached_evidence(
    variant_id: str,
    gene: str,
    rules: Optional[List[str]] = None,
) -> Dict[str, EvidenceResult]:
    """
    Retrieve cached evidence results for a variant.

    Args:
        variant_id: Variant identifier
        gene: Gene symbol
        rules: Optional list of rules to retrieve (e.g., ["PS2", "PM3"])
               If None, retrieves all available rules

    Returns:
        Dict mapping rule name to EvidenceResult
    """
    with get_cache() as cache:
        if rules:
            results = {}
            for rule in rules:
                result = cache.get_evidence_result(variant_id, gene, rule)
                if result:
                    results[rule] = result
            return results
        else:
            return cache.get_all_evidence_results(variant_id, gene)


def has_cached_evidence(variant_id: str, gene: str) -> bool:
    """Check if variant has any cached evidence."""
    with get_cache() as cache:
        results = cache.get_all_evidence_results(variant_id, gene)
        return len(results) > 0
