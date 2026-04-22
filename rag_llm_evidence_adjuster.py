#!/usr/bin/env python3
"""
RAG + LLM Evidence Adjuster for ClinGen gene-specific ACMG rule adjustments.

This module provides:
1. RAG retrieval of ClinGen gene-specific rules
2. LLM API integration for evidence adjustment (supports OpenAI and DeepSeek)
3. Evidence adjustment before final classification

Usage:
    adjuster = EvidenceAdjuster(clingen_db_path)
    adjustment = adjuster.adjust_evidence(
        gene="BRCA1",
        variant_info={...},
        rule_results=[...],
    )

API Key Configuration (unified via llm_config):
    - OPENAI_API_KEY: OpenAI API key (preferred)
    - DEEPSEEK_API_KEY: DeepSeek API key (fallback)
    - LLM_API_KEY: Generic fallback key
"""

import json
import logging
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple
from dataclasses import dataclass
from enum import Enum

from llm_config import get_llm_config, LLMConfig
from llm_prompts import get_prompt

logger = logging.getLogger("GenOtoScope_Classify.rag_llm_adjuster")


class AdjustmentAction(Enum):
    """Evidence adjustment action types."""
    KEEP = "keep"  # Keep original evidence as is
    UPGRADE = "upgrade"  # Upgrade evidence strength
    DOWNGRADE = "downgrade"  # Downgrade evidence strength
    REMOVE = "remove"  # Remove evidence (rule not applicable)
    ADD = "add"  # Add new evidence


@dataclass
class EvidenceAdjustment:
    """Single evidence adjustment recommendation."""
    rule_code: str
    original_strength: str
    adjusted_strength: str
    action: AdjustmentAction
    reason: str
    clingen_criteria: str


@dataclass
class AdjustmentResult:
    """Result of LLM evidence adjustment."""
    gene: str
    adjustments: List[EvidenceAdjustment]
    summary: str
    raw_llm_response: Optional[str] = None


class LLMAPI:
    """Unified LLM API client supporting OpenAI and DeepSeek."""

    def __init__(self, api_key: Optional[str] = None, provider: str = None, model: str = None):
        """
        Initialize LLM API client.

        Args:
            api_key: Optional API key override
            provider: "openai" or "deepseek" (auto-detected if not provided)
            model: Model name (auto-selected based on provider if not provided)
        """
        if api_key:
            # Explicit API key provided
            llm_config = get_llm_config(provider)
            llm_config.api_key = api_key
        else:
            # Use unified config
            llm_config = get_llm_config(provider)

        if not llm_config.available:
            raise Exception(
                "LLM API key not configured. Set OPENAI_API_KEY, DEEPSEEK_API_KEY, or LLM_API_KEY "
                "environment variable."
            )

        self.api_key = llm_config.api_key
        self.provider = llm_config.provider
        self.model = model or llm_config.model

        if self.provider == "openai":
            self.base_url = "https://api.openai.com/v1"
        else:
            self.base_url = "https://api.deepseek.com/v1"

    def chat(self, messages: List[Dict[str, str]], **kwargs) -> str:
        """
        Send chat request to LLM API.

        Args:
            messages: List of message dicts with 'role' and 'content'
            **kwargs: Additional parameters (temperature, max_tokens, etc.)

        Returns:
            Assistant's response text

        Raises:
            Exception: If API call fails
        """
        if self.provider == "openai":
            return self._chat_openai(messages, **kwargs)
        else:
            return self._chat_deepseek(messages, **kwargs)

    def _chat_openai(self, messages: List[Dict[str, str]], **kwargs) -> str:
        """Call OpenAI API."""
        try:
            from openai import OpenAI
            client = OpenAI(api_key=self.api_key)
            response = client.chat.completions.create(
                model=self.model,
                messages=messages,
                **kwargs,
            )
            return response.choices[0].message.content
        except Exception as e:
            logger.error(f"OpenAI API error: {e}")
            raise Exception(f"OpenAI API error: {e}")

    def _chat_deepseek(self, messages: List[Dict[str, str]], **kwargs) -> str:
        """Call DeepSeek API."""
        import urllib.request
        import urllib.error

        url = f"{self.base_url}/chat/completions"
        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self.api_key}",
        }

        payload = {
            "model": self.model,
            "messages": messages,
            **kwargs,
        }

        req = urllib.request.Request(
            url,
            data=json.dumps(payload).encode("utf-8"),
            headers=headers,
            method="POST",
        )

        try:
            with urllib.request.urlopen(req, timeout=60) as response:
                result = json.loads(response.read().decode("utf-8"))
                return result["choices"][0]["message"]["content"]
        except urllib.error.HTTPError as e:
            error_body = e.read().decode("utf-8")
            logger.error(f"DeepSeek API error: {e.code} - {error_body}")
            raise Exception(f"DeepSeek API error: {e.code}")
        except Exception as e:
            logger.error(f"DeepSeek API request failed: {e}")
            raise


# Backwards compatibility alias
DeepSeekAPI = LLMAPI


class ClinGenRAGRetriever:
    """
    RAG retriever for ClinGen gene-specific rules.

    Retrieves relevant ClinGen rule text chunks for LLM context.
    """

    def __init__(self, db_path: Path):
        self.db_path = db_path
        self._conn = None

    def _get_connection(self):
        if self._conn is None:
            import sqlite3
            self._conn = sqlite3.connect(str(self.db_path))
            self._conn.row_factory = sqlite3.Row
        return self._conn

    def close(self):
        if self._conn:
            self._conn.close()
            self._conn = None

    def retrieve_rules_for_gene(
        self,
        gene: str,
        rule_codes: Optional[List[str]] = None,
    ) -> str:
        """
        Retrieve ClinGen rules for a gene as formatted text.

        Args:
            gene: Gene symbol
            rule_codes: Optional filter for specific rule codes

        Returns:
            Formatted text of ClinGen rules
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        if rule_codes:
            placeholders = ",".join("?" * len(rule_codes))
            cursor.execute(
                f"""SELECT * FROM clingen_rules
                   WHERE gene = ? AND rule_code IN ({placeholders})
                   ORDER BY rule_code, strength_order DESC""",
                [gene] + rule_codes,
            )
        else:
            cursor.execute(
                """SELECT * FROM clingen_rules
                   WHERE gene = ?
                   ORDER BY rule_code, strength_order DESC""",
                (gene,),
            )

        lines = [f"ClinGen ACMG Rules for {gene}", "=" * 50]

        current_code = None
        for row in cursor.fetchall():
            rc = row["rule_code"]
            if rc != current_code:
                lines.append(f"\n[{rc}]")
                current_code = rc

            lines.append(f"  Strength: {row['strength']}")
            if row['description']:
                lines.append(f"  Description: {row['description'][:150]}...")
            if row['application_criteria']:
                lines.append(f"  Criteria: {row['application_criteria'][:150]}...")

        return "\n".join(lines)

    def retrieve_pm1_domains(self, gene: str) -> str:
        """Retrieve PM1 domain specifications for a gene."""
        conn = self._get_connection()
        cursor = conn.cursor()

        cursor.execute(
            "SELECT * FROM pm1_details WHERE gene = ?",
            (gene,)
        )

        lines = [f"ClinGen PM1 Domains for {gene}", "=" * 50]

        for row in cursor.fetchall():
            lines.append(f"\n{row['strength']}: {row['protein_positions']}")
            if row['domain']:
                lines.append(f"  Domain: {row['domain']}")
            if row['notes']:
                lines.append(f"  Notes: {row['notes'][:200]}...")

        if len(lines) == 1:
            return f"No ClinGen PM1 domains defined for {gene}."

        return "\n".join(lines)

    def retrieve_frequency_cutoffs(self, gene: str) -> str:
        """Retrieve frequency cutoffs for a gene."""
        conn = self._get_connection()
        cursor = conn.cursor()

        cursor.execute(
            "SELECT * FROM frequency_cutoffs WHERE gene = ?",
            (gene,)
        )

        lines = [f"ClinGen Frequency Cutoffs for {gene}", "=" * 50]

        for row in cursor.fetchall():
            lines.append(
                f"\n{row['rule_code']} ({row['strength']}): "
                f"{row['cutoff_operator']} {row['cutoff_value']} ({row['frequency_type']})"
            )
            if row['notes']:
                lines.append(f"  Notes: {row['notes'][:200]}...")

        if len(lines) == 1:
            return f"No ClinGen frequency cutoffs defined for {gene}."

        return "\n".join(lines)


class EvidenceAdjuster:
    """
    RAG + LLM evidence adjuster using ClinGen gene-specific rules.

    Workflow:
    1. Retrieve ClinGen rules via RAG
    2. Construct prompt with variant info and rule results
    3. Send to LLM for adjustment decision (auto-selects OpenAI or DeepSeek)
    4. Parse and return adjustments
    """

    @staticmethod
    def _get_rag_llm_prompt() -> str:
        return get_prompt("rag_llm_system")

    @property
    def SYSTEM_PROMPT(self) -> str:
        return self._get_rag_llm_prompt()

    def __init__(
        self,
        clingen_db_path: Path,
        llm_api_key: Optional[str] = None,
        provider: str = None,
        model: str = None,
    ):
        """
        Initialize EvidenceAdjuster.

        Args:
            clingen_db_path: Path to ClinGen rules database
            llm_api_key: Optional LLM API key override
            provider: "openai" or "deepseek" (auto-detected if not provided)
            model: Model name (auto-selected based on provider if not provided)
        """
        self.clingen_db_path = clingen_db_path
        self.llm_api = LLMAPI(api_key=llm_api_key, provider=provider, model=model)
        self.rag = ClinGenRAGRetriever(clingen_db_path)

    def adjust_evidence(
        self,
        gene: str,
        rule_results: List[Dict[str, Any]],
        variant_info: Optional[Dict[str, Any]] = None,
    ) -> AdjustmentResult:
        """
        Adjust ACMG evidence based on ClinGen gene-specific rules.

        Args:
            gene: Gene symbol
            rule_results: List of rule assessment results
                Each dict should have: rule_code, strength, status, comment
            variant_info: Optional variant information dict
                (protein_position, frequency, is_splice, etc.)

        Returns:
            AdjustmentResult with recommended adjustments
        """
        # Step 1: RAG retrieval
        clin_gen_rules = self.rag.retrieve_rules_for_gene(gene)
        pm1_domains = self.rag.retrieve_pm1_domains(gene)
        freq_cutoffs = self.rag.retrieve_frequency_cutoffs(gene)

        # Step 2: Construct user prompt
        user_prompt = self._construct_prompt(
            gene=gene,
            rule_results=rule_results,
            variant_info=variant_info,
            clin_gen_rules=clin_gen_rules,
            pm1_domains=pm1_domains,
            freq_cutoffs=freq_cutoffs,
        )

        # Step 3: Call DeepSeek LLM
        messages = [
            {"role": "system", "content": self.SYSTEM_PROMPT},
            {"role": "user", "content": user_prompt},
        ]

        try:
            response = self.llm_api.chat(
                messages,
                temperature=0.1,
                max_tokens=2000,
            )
        except Exception as e:
            logger.error(f"LLM call failed: {e}")
            return AdjustmentResult(
                gene=gene,
                adjustments=[],
                summary=f"LLM call failed: {e}",
                raw_llm_response=None,
            )

        # Step 4: Parse response
        adjustments = self._parse_llm_response(response)

        return AdjustmentResult(
            gene=gene,
            adjustments=adjustments,
            summary=self._summarize_adjustments(adjustments),
            raw_llm_response=response,
        )

    def _construct_prompt(
        self,
        gene: str,
        rule_results: List[Dict[str, Any]],
        variant_info: Optional[Dict[str, Any]],
        clin_gen_rules: str,
        pm1_domains: str,
        freq_cutoffs: str,
    ) -> str:
        """Construct user prompt for LLM."""
        lines = [
            f"# Gene: {gene}",
            "",
            "## Variant Information:",
        ]

        if variant_info:
            for key, value in variant_info.items():
                lines.append(f"- {key}: {value}")
        else:
            lines.append("(Not provided)")

        lines.extend([
            "",
            "## Current ACMG Evidence:",
        ])

        for result in rule_results:
            lines.append(
                f"- {result.get('rule_code', 'UNKNOWN')}: "
                f"{result.get('strength', 'N/A')} "
                f"(status: {result.get('status', 'N/A')}) "
                f"- {result.get('comment', '')}"
            )

        lines.extend([
            "",
            "## ClinGen Gene-Specific Rules:",
            clin_gen_rules,
            "",
            "## ClinGen PM1 Domains:",
            pm1_domains,
            "",
            "## ClinGen Frequency Cutoffs:",
            freq_cutoffs,
            "",
            "Based on the ClinGen rules above, provide your evidence adjustments in JSON format.",
        ])

        return "\n".join(lines)

    def _parse_llm_response(self, response: str) -> List[EvidenceAdjustment]:
        """Parse LLM JSON response into EvidenceAdjustment objects."""
        try:
            # Try to extract JSON from response
            json_str = response
            if "```json" in response:
                json_str = response.split("```json")[1].split("```")[0]
            elif "```" in response:
                json_str = response.split("```")[1].split("```")[0]

            data = json.loads(json_str.strip())
            adjustments = []

            for adj in data.get("adjustments", []):
                action = AdjustmentAction(adj.get("action", "keep"))

                adjustments.append(EvidenceAdjustment(
                    rule_code=adj.get("rule_code", "UNKNOWN"),
                    original_strength=adj.get("original_strength", "N/A"),
                    adjusted_strength=adj.get("adjusted_strength"),
                    action=action,
                    reason=adj.get("reason", ""),
                    clingen_criteria=adj.get("clingen_criteria", ""),
                ))

            return adjustments

        except json.JSONDecodeError as e:
            logger.warning(f"Failed to parse LLM JSON response: {e}")
            return []

    def _summarize_adjustments(self, adjustments: List[EvidenceAdjustment]) -> str:
        """Create a summary of adjustments."""
        if not adjustments:
            return "No adjustments recommended."

        actions_count = {}
        for adj in adjustments:
            actions_count[adj.action.value] = actions_count.get(adj.action.value, 0) + 1

        parts = []
        for action, count in actions_count.items():
            parts.append(f"{count} {action}(s)")

        return " + ".join(parts)


def adjust_evidence(
    gene: str,
    rule_results: List[Dict[str, Any]],
    clingen_db_path: Path,
    llm_api_key: Optional[str] = None,
    variant_info: Optional[Dict[str, Any]] = None,
    provider: str = None,
) -> AdjustmentResult:
    """
    Convenience function to adjust evidence using ClinGen rules.

    Args:
        gene: Gene symbol
        rule_results: List of rule assessment results
        clingen_db_path: Path to clingen_rules.db
        llm_api_key: LLM API key (optional, uses OPENAI_API_KEY/DEEPSEEK_API_KEY/LLM_API_KEY env vars)
        variant_info: Optional variant information
        provider: "openai" or "deepseek" (auto-detected if not provided)

    Returns:
        AdjustmentResult with recommended adjustments
    """
    adjuster = EvidenceAdjuster(
        clingen_db_path=clingen_db_path,
        llm_api_key=llm_api_key,
        provider=provider,
    )
    return adjuster.adjust_evidence(
        gene=gene,
        rule_results=rule_results,
        variant_info=variant_info,
    )
