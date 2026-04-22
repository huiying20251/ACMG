#!/usr/bin/env python3
"""
文献相关性分类与信息提取模块

在文献检索后使用 LLM 判断：
1. 文献是否与变异相关
2. 文献类型：病例报道 / 队列研究 / 功能研究 / 综述
3. 根据类型触发不同的证据提取流程

API Key Configuration (unified via llm_config):
    - OPENAI_API_KEY: OpenAI API key (preferred)
    - DEEPSEEK_API_KEY: DeepSeek API key (fallback)
    - LLM_API_KEY: Generic fallback key

Prompts Configuration:
    - Prompts are loaded from config/llm_prompts.yaml
    - Fallback prompts are available if config file is not found
"""

import logging
from typing import Optional, List, Dict, Any
from dataclasses import dataclass
from enum import Enum

from llm_config import get_llm_config
from llm_prompts import get_prompt

logger = logging.getLogger("GenOtoScope_Classify.literature_classifier")


class LiteratureType(Enum):
    """文献类型"""
    CASE_REPORT = "case_report"           # 病例报道
    COHORT_STUDY = "cohort_study"        # 队列研究
    FUNCTIONAL_STUDY = "functional_study"  # 功能研究
    REVIEW = "review"                    # 综述
    OTHER = "other"                      # 其他


class FunctionalResult(Enum):
    """功能研究结果"""
    LOSS_OF_FUNCTION = "lof"             # 功能丧失 (影响功能 → PS3)
    GAIN_OF_FUNCTION = "gof"            # 功能获得 (影响功能 → PS3)
    NORMAL_FUNCTION = "normal"           # 功能正常 (BS3)
    UNCERTAIN = "uncertain"              # 不确定 → 不使用


@dataclass
class LiteratureClassification:
    """LLM 文献分类结果"""
    pmid: str
    is_variant_related: bool             # 是否与变异相关
    literature_type: LiteratureType       # 文献类型
    confidence: float                     # 置信度 0-1
    reasoning: str                        # 判断理由


@dataclass
class CaseReportEvidence:
    """病例报道证据（用于 PS2/PM3/PP1）"""
    pmid: str
    case_count: int                      # 病例数
    is_de_novo: bool                     # 是否 de novo
    de_novo_confirmed: bool               # 是否经 parental testing 确认
    inheritance_pattern: str              # 遗传模式
    is_compound_het: bool                # 是否 compound het
    trans_variant: Optional[str] = None  # 反式位置的另一变异
    segregation_data: Optional[Dict] = None  # 共分离数据
    phenotype: Optional[str] = None       # 表型
    hpo_terms: Optional[List[str]] = None  # HPO terms
    statistical_significance: Optional[Dict] = None  # 统计数据 (OR, p-value)
    # ClinVar fields
    clinvar_status: str = ""              # e.g. "Pathogenic", "VUS", etc.
    clinvar_significance: str = ""        # e.g. "single hit", "double hit"
    individual_id: str = "patient"       # Patient identifier

    def to_individual_observation(self, gene: str, variant_description: str) -> 'IndividualVariantObservation':
        """
        Convert CaseReportEvidence to IndividualVariantObservation.

        Args:
            gene: Gene symbol
            variant_description: Variant description

        Returns:
            IndividualVariantObservation object
        """
        from literature_retrieval.literature_utils import IndividualVariantObservation

        # Determine zygosity based on compound het status
        zygosity = "unknown"
        if self.is_compound_het:
            zygosity = "compound_heterozygous"

        # Determine variant inheritance string
        variant_inheritance = "unknown"
        if self.is_de_novo:
            variant_inheritance = "de novo"

        return IndividualVariantObservation(
            pmid=self.pmid,
            individual_id=self.individual_id,
            variant_description=variant_description,
            gene=gene,
            variant_inheritance=variant_inheritance,
            parental_testing=self.de_novo_confirmed,
            zygosity=zygosity,
            clinvar_status=self.clinvar_status,
            clinvar_significance=self.clinvar_significance,
            phenotype=self.phenotype or "",
            hpo_terms=self.hpo_terms or [],
            segregation_data=self.segregation_data,
            confidence="medium",
        )


@dataclass
class FunctionalStudyEvidence:
    """功能研究证据（用于 PS3/BS3）"""
    pmid: str
    functional_result: FunctionalResult  # 功能结果
    technique: str                       # 实验技术
    percentage_effect: Optional[float] = None  # 功能影响百分比
    confidence: str = "medium"           # confidence: high/medium/low
    reasoning: str = ""


# ========== LLM Prompt Templates (loaded from config) ==========

# Functions to get prompts from config file
def _get_literature_classification_prompt() -> str:
    return get_prompt("literature_classification")

def _get_functional_study_extraction_prompt() -> str:
    return get_prompt("functional_study_extraction")

def _get_case_report_extraction_prompt() -> str:
    return get_prompt("case_report_extraction")


# ========== LLM Classifier Class ==========

class LiteratureClassifier:
    """
    LLM-based literature classifier.

    Classifies retrieved articles and routes to appropriate evidence extraction.
    """

    def __init__(self, llm_api_key: Optional[str] = None):
        self.llm_api_key = llm_api_key

    def classify_article(
        self,
        pmid: str,
        title: str,
        abstract: str,
        gene: str,
        variant: str,
    ) -> LiteratureClassification:
        """
        Classify a single article using LLM.

        Args:
            pmid: PubMed ID
            title: Article title
            abstract: Article abstract
            gene: Gene symbol
            variant: Variant description

        Returns:
            LiteratureClassification with type and relevance
        """
        prompt = _get_literature_classification_prompt().format(
            gene=gene,
            variant=variant,
            title=title or "",
            abstract=abstract or "",
        )

        result = self._call_llm(prompt)

        return LiteratureClassification(
            pmid=pmid,
            is_variant_related=result.get("is_variant_related", False),
            literature_type=LiteratureType(result.get("literature_type", "other")),
            confidence=result.get("confidence", 0.0),
            reasoning=result.get("reasoning", ""),
        )

    def extract_case_report_evidence(
        self,
        pmid: str,
        title: str,
        abstract: str,
        gene: str,
        variant: str,
        inheritance_pattern: str = "UNKNOWN",
    ) -> Optional[CaseReportEvidence]:
        """
        Extract case report evidence (PS2/PM3/PP1/PS4).

        Returns:
            CaseReportEvidence if applicable, None otherwise
        """
        prompt = _get_case_report_extraction_prompt().format(
            gene=gene,
            variant=variant,
            inheritance_pattern=inheritance_pattern,
            title=title or "",
            abstract=abstract or "",
        )

        result = self._call_llm(prompt)

        if not result or not result.get("case_count", 0) > 0:
            return None

        return CaseReportEvidence(
            pmid=pmid,
            case_count=result.get("case_count", 1),
            is_de_novo=result.get("is_de_novo", False),
            de_novo_confirmed=result.get("de_novo_confirmed", False),
            inheritance_pattern=result.get("inheritance_pattern", "UNKNOWN"),
            is_compound_het=result.get("is_compound_het", False),
            trans_variant=result.get("trans_variant"),
            segregation_data=result.get("segregation_data"),
            phenotype=result.get("phenotype"),
            hpo_terms=result.get("hpo_terms", []),
            statistical_significance={
                "case_count": result.get("case_count"),
                "control_count": result.get("control_count"),
                "odds_ratio": result.get("odds_ratio"),
                "p_value": result.get("p_value"),
            } if result.get("case_count") else None,
        )

    def extract_functional_evidence(
        self,
        pmid: str,
        title: str,
        abstract: str,
        gene: str,
        variant: str,
    ) -> Optional[FunctionalStudyEvidence]:
        """
        Extract functional study evidence (PS3/BS3).

        Returns:
            FunctionalStudyEvidence if the variant affects function (LOF or GOF),
            None if uncertain or normal function
        """
        prompt = _get_functional_study_extraction_prompt().format(
            gene=gene,
            variant=variant,
            title=title or "",
            abstract=abstract or "",
        )

        result = self._call_llm(prompt)

        if not result:
            return None

        func_result = result.get("functional_result", "uncertain")

        # Only return if LOF or GOF (both indicate affected function → PS3)
        # Don't return for normal or uncertain
        if func_result not in ("lof", "gof"):
            return None

        return FunctionalStudyEvidence(
            pmid=pmid,
            functional_result=FunctionalResult(func_result),
            technique=result.get("technique", "unspecified"),
            percentage_effect=result.get("percentage_effect"),
            confidence=result.get("confidence", "medium"),
            reasoning=result.get("reasoning", ""),
        )

    def _call_llm(self, prompt: str) -> Dict[str, Any]:
        """
        Call LLM API with prompt.

        Uses OpenAI GPT-4o API for JSON responses by default.
        Falls back to DeepSeek if OpenAI key not available.
        """
        import json

        # Get API config (prefer OpenAI)
        llm_config = get_llm_config(provider="openai")
        if not llm_config.available:
            logger.warning("LLM API key not configured - set OPENAI_API_KEY, DEEPSEEK_API_KEY, or LLM_API_KEY")
            return {"is_variant_related": False, "literature_type": "other", "confidence": 0.0}

        try:
            if llm_config.provider == "openai":
                return self._call_openai(prompt, llm_config.api_key)
            else:
                return self._call_deepseek(prompt, llm_config.api_key, llm_config.model)

        except Exception as e:
            logger.warning(f"LLM API call failed: {e}")
            return {"is_variant_related": False, "literature_type": "other", "confidence": 0.0}

    def _call_openai(self, prompt: str, api_key: str) -> Dict[str, Any]:
        """Call OpenAI API."""
        import json
        from openai import OpenAI

        client = OpenAI(api_key=api_key)
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {
                    "role": "system",
                    "content": "You are a helpful assistant that responds in JSON format only."
                },
                {
                    "role": "user",
                    "content": prompt + "\n\nRespond ONLY with valid JSON, no additional text."
                }
            ],
            temperature=0.3,
            max_tokens=1024,
            response_format={"type": "json_object"}
        )

        content = response.choices[0].message.content
        if content:
            return json.loads(content)
        else:
            logger.warning("Empty response from OpenAI API")
            return {"is_variant_related": False, "literature_type": "other", "confidence": 0.0}

    def _call_deepseek(self, prompt: str, api_key: str, model: str = "deepseek-chat") -> Dict[str, Any]:
        """Call DeepSeek API."""
        import json
        import urllib.request

        url = "https://api.deepseek.com/v1/chat/completions"
        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {api_key}",
        }

        payload = {
            "model": model,
            "messages": [
                {
                    "role": "system",
                    "content": "You are a helpful assistant that responds in JSON format only."
                },
                {
                    "role": "user",
                    "content": prompt + "\n\nRespond ONLY with valid JSON, no additional text."
                }
            ],
            "temperature": 0.3,
            "max_tokens": 1024,
        }

        req = urllib.request.Request(
            url,
            data=json.dumps(payload).encode("utf-8"),
            headers=headers,
            method="POST",
        )

        with urllib.request.urlopen(req, timeout=60) as response:
            result = json.loads(response.read().decode("utf-8"))
            content = result["choices"][0]["message"]["content"]
            if content:
                return json.loads(content)
            else:
                logger.warning("Empty response from DeepSeek API")
                return {"is_variant_related": False, "literature_type": "other", "confidence": 0.0}


# ========== Convenience Functions ==========

def classify_and_route_literature(
    articles: List[Any],
    gene: str,
    variant: str,
    inheritance_pattern: str = "UNKNOWN",
    llm_api_key: Optional[str] = None,
) -> tuple[List[CaseReportEvidence], List[FunctionalStudyEvidence]]:
    """
    分类文献并路由到相应证据提取

    Args:
        articles: List of Article objects from literature retrieval
        gene: Gene symbol
        variant: Variant description
        inheritance_pattern: Expected inheritance pattern
        llm_api_key: LLM API key

    Returns:
        tuple: (case_report_evidences, functional_study_evidences)
    """
    classifier = LiteratureClassifier(llm_api_key)

    case_reports = []
    functional_studies = []

    for article in articles:
        # Step 1: Classify article
        classification = classifier.classify_article(
            pmid=article.pmid,
            title=article.title,
            abstract=article.abstract,
            gene=gene,
            variant=variant,
        )

        if not classification.is_variant_related:
            logger.info(f"PMID {article.pmid}: Not variant-related, skipping")
            continue

        logger.info(
            f"PMID {article.pmid}: {classification.literature_type.value} "
            f"(confidence: {classification.confidence:.2f})"
        )

        # Step 2: Route based on type
        if classification.literature_type == LiteratureType.CASE_REPORT:
            evidence = classifier.extract_case_report_evidence(
                pmid=article.pmid,
                title=article.title,
                abstract=article.abstract,
                gene=gene,
                variant=variant,
                inheritance_pattern=inheritance_pattern,
            )
            if evidence:
                case_reports.append(evidence)

        elif classification.literature_type == LiteratureType.COHORT_STUDY:
            evidence = classifier.extract_case_report_evidence(
                pmid=article.pmid,
                title=article.title,
                abstract=article.abstract,
                gene=gene,
                variant=variant,
                inheritance_pattern=inheritance_pattern,
            )
            if evidence:
                case_reports.append(evidence)

        elif classification.literature_type == LiteratureType.FUNCTIONAL_STUDY:
            evidence = classifier.extract_functional_evidence(
                pmid=article.pmid,
                title=article.title,
                abstract=article.abstract,
                gene=gene,
                variant=variant,
            )
            if evidence:
                functional_studies.append(evidence)

    return case_reports, functional_studies
