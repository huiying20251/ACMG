#!/usr/bin/env python3
"""
文献检索触发模块

在 variant normalize 后自动触发多策略文献检索：
- rsID
- 基因名称
- 基因 + c.HGVS
- 基因 + p.HGVS
- 染色体位置 (VCF 格式)

扩展检索（如果一级检索无结果）：
- 基因 + c.核苷酸位置 (如 "BRCA1 c.68")
- 基因 + 密码子位置 (如 "BRCA1 Arg1699")

检索结果用于 ACMG 证据评估：
- PS2: de novo (显性遗传)
- PM3: biallelic (隐性遗传)
- PP1: co-segregation (家系共分离)
- BP2: 反式位置 (隐性遗传)
"""

import logging
import re
from typing import Optional, List, Any, Tuple

from normalizer import VariantNormalizer, VariantInfo
from literature_retrieval import (
    NCBILiteratureRetriever,
    LiteratureRetriever,
    LiteratureConfig,
    VariantLiterature,
    assess_pm3_from_literature,
    assess_ps2_from_literature,
    assess_pp1_from_pedigree,
    assess_bp2_from_partners,
    EvidenceCacheDB,
    cache_evidence,
    get_cached_evidence,
    has_cached_evidence,
)
from literature_retrieval.literature_utils import InheritancePattern

logger = logging.getLogger("GenOtoScope_Classify.literature_trigger")


def build_expanded_queries(normalized: VariantInfo) -> List[str]:
    """
    构建扩展检索词（当一级检索无结果时使用）

    扩展策略：
    - 基因 + c.核苷酸位置 (如 "BRCA1 c.68")
    - 基因 + 密码子位置 (如 "BRCA1 Arg1699" 或 "BRCA1 1699")
    - 染色体位置 (如 "chr17:41244938" 或 "17:41244938")

    Args:
        normalized: VariantInfo from normalizer.normalize()

    Returns:
        List of expanded search queries
    """
    expanded_queries = []
    gene = normalized.gene

    if not gene:
        return expanded_queries

    # 1. 基因 + c.核苷酸位置
    if normalized.hgvs_c:
        # 提取 c. 后面的位置数字
        c_match = re.match(r'^c\.(\d+)', normalized.hgvs_c, re.IGNORECASE)
        if c_match:
            c_position = c_match.group(1)
            # 提取核苷酸变化部分（如果有）
            c_change_match = re.search(r'\d+([A-Z>]+)$', normalized.hgvs_c)
            if c_change_match:
                c_change = c_change_match.group(1)
                expanded_queries.append(f"{gene} c.{c_position}{c_change}")
            else:
                # 只有位置没有变化
                expanded_queries.append(f"{gene} c.{c_position}")
            # 纯数字位置
            expanded_queries.append(f"{gene} c.{c_position}")

    # 2. 基因 + 密码子位置 (氨基酸位置)
    if normalized.hgvs_p:
        # 提取 p. 后面的部分，如 "Gly1699Arg" -> "Gly1699" 或 "Arg1699"
        p_match = re.match(r'^p\.([A-Za-z]+)(\d+)([A-Za-z]+)?$', normalized.hgvs_p)
        if p_match:
            ref_aa = p_match.group(1)
            position = p_match.group(2)
            alt_aa = p_match.group(3)

            # 多种格式检索
            if alt_aa:
                # 完整氨基酸改变: "BRCA1 Gly1699Arg"
                expanded_queries.append(f"{gene} {ref_aa}{position}{alt_aa}")
            # 只用位置: "BRCA1 1699"
            expanded_queries.append(f"{gene} {position}")
            # 用参考氨基酸: "BRCA1 Gly1699"
            expanded_queries.append(f"{gene} {ref_aa}{position}")

    # 3. 染色体位置（带和不带 chr 前缀）
    if normalized.chromosome and normalized.position:
        chrom = normalized.chromosome
        position = normalized.position

        # chr17:41244938 和 17:41244938
        if chrom.startswith("chr"):
            expanded_queries.append(f"chr{chrom.replace('chr', '')}:{position}")
            expanded_queries.append(f"{chrom}:{position}")
        else:
            expanded_queries.append(f"chr{chrom}:{position}")
            expanded_queries.append(f"{chrom}:{position}")

    # 去重
    unique_queries = list(dict.fromkeys(expanded_queries))
    logger.info(f"Built {len(unique_queries)} expanded queries")
    return unique_queries


def search_literature_for_variant(
    normalized: VariantInfo,
    inheritance_pattern: InheritancePattern = InheritancePattern.UNKNOWN,
    max_results_per_query: int = 20,
    max_total_articles: int = 100,
    enable_expanded_search: bool = True,
) -> Tuple[VariantLiterature, bool]:
    """
    为标准化后的变异触发文献检索

    检索策略（两级）：
    1. 一级检索：使用完整检索词
       - rs号、基因、基因+c.HGVS、基因+p.HGVS、VCF格式
    2. 扩展检索（如果一级无结果）：
       - 基因 + c.核苷酸位置
       - 基因 + 密码子位置

    Args:
        normalized: VariantInfo from normalizer.normalize()
        inheritance_pattern: 遗传模式 (AD, AR, XLD, etc.)
        max_results_per_query: 每个检索词最大返回数
        max_total_articles: 总计最大文章数
        enable_expanded_search: 一级无结果时是否启用扩展检索

    Returns:
        Tuple[VariantLiterature, bool]: (文献结果, 是否使用了扩展检索)
    """
    normalizer = VariantNormalizer()

    # ========== 一级检索 ==========
    queries = normalizer.build_search_queries(normalized)

    logger.info(f"Generated {len(queries)} primary search queries")
    for i, query in enumerate(queries, 1):
        logger.info(f"  Primary Query {i}: {query}")

    # 配置文献检索
    config = LiteratureConfig(
        max_results_per_query=max_results_per_query,
        max_total_articles=max_total_articles,
        default_inheritance=inheritance_pattern,
        include_pubtator=True,
    )

    # 初始化检索器
    retriever = LiteratureRetriever(config)

    # 执行一级检索
    literature = retriever.search(normalized)

    logger.info(
        f"Primary search complete: {literature.total_articles} articles found, "
        f"{literature.num_case_reports} case reports"
    )

    # ========== 扩展检索（如果一级无结果） ==========
    used_expanded_search = False
    if enable_expanded_search and literature.total_articles == 0:
        logger.info("No results from primary search, attempting expanded search...")

        expanded_queries = build_expanded_queries(normalized)
        if expanded_queries:
            logger.info(f"Generated {len(expanded_queries)} expanded queries")
            for i, query in enumerate(expanded_queries, 1):
                logger.info(f"  Expanded Query {i}: {query}")

            # 使用扩展检索词重新检索
            # 注意：这里需要修改 NCBILiteratureRetriever 来支持自定义查询
            # 或者直接使用扩展查询作为新的检索
            literature_expanded = _search_with_queries(
                retriever, normalized, expanded_queries
            )

            if literature_expanded.total_articles > 0:
                literature = literature_expanded
                used_expanded_search = True
                logger.info(
                    f"Expanded search successful: {literature.total_articles} articles found"
                )
            else:
                logger.info("Expanded search also returned no results")
        else:
            logger.info("No expanded queries available for this variant")

    return literature, used_expanded_search


def _search_with_queries(
    retriever: LiteratureRetriever,
    variant_info: VariantInfo,
    queries: List[str],
) -> VariantLiterature:
    """
    使用自定义查询词执行文献检索

    Args:
        retriever: LiteratureRetriever instance
        variant_info: Variant information
        queries: List of custom search queries

    Returns:
        VariantLiterature with combined results
    """
    from literature_retrieval.literature_utils import Article

    all_articles = []
    seen_pmids = set()

    for query in queries:
        try:
            # 使用 NCBILiteratureRetriever 的底层方法
            if hasattr(retriever, 'ncbi_retriever'):
                articles = retriever.ncbi_retriever.search(query, max_results=20)
                for article in articles:
                    if article.pmid not in seen_pmids:
                        seen_pmids.add(article.pmid)
                        all_articles.append(article)

            # 如果有 PubTator，也用扩展查询提取
            if hasattr(retriever, 'variant_extractor'):
                mentions = retriever.variant_extractor.extract_variant_mentions(
                    all_articles,
                    variant_rsids=[variant_info.rs_id] if variant_info.rs_id else [],
                    variant_hgvs=[variant_info.hgvs_c, variant_info.hgvs_p],
                )
        except Exception as e:
            logger.warning(f"Expanded query '{query}' failed: {e}")
            continue

    # 创建合并的 literature 对象
    literature = retriever.evidence_router.assess_evidence(
        VariantLiterature(
            variant_id=variant_info.rs_id or variant_info.hgvs_c or "",
            gene=variant_info.gene or "",
            articles=all_articles,
            total_articles=len(all_articles),
        ),
        InheritancePattern.UNKNOWN,
    )

    return literature


def _format_comment_with_pmids(comment: str, pmids: list) -> str:
    """
    Add PMID list to comment.

    Args:
        comment: Original comment
        pmids: List of PMIDs

    Returns:
        Comment with PMIDs appended
    """
    if not pmids:
        return comment
    pmid_str = ", ".join([f"PMID:{p}" for p in pmids])
    if comment:
        return f"{comment} ({pmid_str})"
    return f"Evidence from {pmid_str}"


def assess_literature_evidence(
    literature: VariantLiterature,
    inheritance_pattern: InheritancePattern = InheritancePattern.UNKNOWN,
    gene: Optional[str] = None,
    variant: Optional[str] = None,
    llm_api_key: Optional[str] = None,
) -> dict[str, Any]:
    """
    从文献检索结果中评估 ACMG 证据

    流程:
    1. 使用 LLM 分类文献类型（病例报道/队列/功能研究）
    2. 功能研究文献 → PS3_Literature 证据
    3. 病例/队列文献 → PS2/PM3/PP1/PS4 证据

    Args:
        literature: VariantLiterature from search_literature_for_variant()
        inheritance_pattern: 遗传模式
        gene: Gene symbol (for LLM classification)
        variant: Variant description (for LLM classification)
        llm_api_key: LLM API key (optional)

    Returns:
        dict: 包含各证据评估结果
            - ps2: PS2 (de novo) 评估结果
            - pm3: PM3 (biallelic) 评估结果
            - pp1: PP1 (co-segregation) 评估结果
            - bp2: BP2 (反式位置) 评估结果
            - ps3_literature: PS3 from functional literature
    """
    results = {}

    # Step 1: LLM 文献分类与证据提取
    if gene and variant and llm_api_key:
        try:
            from literature_classifier import classify_and_route_literature, LiteratureClassifier

            classifier = LiteratureClassifier(llm_api_key)

            # 对每篇文献进行 LLM 分类
            case_report_evidences = []
            functional_study_evidences = []

            for article in literature.articles:
                # 分类文献
                classification = classifier.classify_article(
                    pmid=article.pmid,
                    title=article.title,
                    abstract=article.abstract,
                    gene=gene,
                    variant=variant,
                )

                if not classification.is_variant_related:
                    continue

                logger.info(
                    f"PMID {article.pmid}: {classification.literature_type.value} "
                    f"(confidence: {classification.confidence:.2f})"
                )

                # 根据类型提取证据
                if classification.literature_type.value in ("case_report", "cohort_study"):
                    evidence = classifier.extract_case_report_evidence(
                        pmid=article.pmid,
                        title=article.title,
                        abstract=article.abstract,
                        gene=gene,
                        variant=variant,
                        inheritance_pattern=inheritance_pattern.value,
                    )
                    if evidence:
                        case_report_evidences.append(evidence)

                elif classification.literature_type.value == "functional_study":
                    evidence = classifier.extract_functional_evidence(
                        pmid=article.pmid,
                        title=article.title,
                        abstract=article.abstract,
                        gene=gene,
                        variant=variant,
                    )
                    if evidence:
                        functional_study_evidences.append(evidence)

            # 存储文献证据
            results["case_report_evidences"] = case_report_evidences
            results["functional_study_evidences"] = functional_study_evidences

            # Convert to individual-level observations and generate strings for LLM
            individual_observation_strings = []
            for i, evidence in enumerate(case_report_evidences):
                # Create individual observation for each case
                evidence.individual_id = f"case_{i+1}"
                obs = evidence.to_individual_observation(gene=gene, variant_description=variant)
                individual_observation_strings.append(obs.to_string())

            if individual_observation_strings:
                results["individual_observation_strings"] = individual_observation_strings
                logger.info(f"Generated {len(individual_observation_strings)} individual observation strings for LLM analysis")

            # PS3 from literature
            if functional_study_evidences:
                results["ps3_literature"] = _assess_ps3_from_literature(functional_study_evidences)
                logger.info(f"PS3 from literature: {len(functional_study_evidences)} functional studies support pathogenicity")

        except ImportError as e:
            logger.warning(f"LiteratureClassifier not available: {e}")
        except Exception as e:
            logger.warning(f"LLM literature classification failed: {e}")

    # Step 2: 传统证据评估 (基于关键词回退)
    # PS2: de novo 评估
    try:
        ps2_result = assess_ps2_from_literature(literature, inheritance_pattern)
        if ps2_result and ps2_result.comment:
            # Add PMIDs to comment
            ps2_result.comment = _format_comment_with_pmids(ps2_result.comment, ps2_result.pmids)
        results["ps2"] = ps2_result
        if ps2_result:
            logger.info(f"PS2 assessment: {ps2_result.strength.value if ps2_result.strength else 'Not applicable'}")
    except Exception as e:
        logger.warning(f"PS2 assessment failed: {e}")
        results["ps2"] = None

    # PM3: biallelic 评估
    try:
        pm3_result = assess_pm3_from_literature(literature, inheritance_pattern)
        if pm3_result and pm3_result.comment:
            # Add PMIDs to comment
            pm3_result.comment = _format_comment_with_pmids(pm3_result.comment, pm3_result.pmids)
        results["pm3"] = pm3_result
        if pm3_result:
            logger.info(f"PM3 assessment: {pm3_result.strength.value if pm3_result.strength else 'Not applicable'}")
    except Exception as e:
        logger.warning(f"PM3 assessment failed: {e}")
        results["pm3"] = None

    # PP1: co-segregation 评估 (需要家系数据)
    try:
        pp1_result = assess_pp1_from_pedigree(literature, inheritance_pattern)
        if pp1_result and hasattr(pp1_result, 'comment') and pp1_result.comment:
            # Add PMIDs to comment
            ppids = pp1_result.pmids if hasattr(pp1_result, 'pmids') else []
            pp1_result.comment = _format_comment_with_pmids(pp1_result.comment, ppids or [])
        results["pp1"] = pp1_result
    except Exception as e:
        logger.warning(f"PP1 assessment failed: {e}")
        results["pp1"] = None

    # BP2: 反式位置评估 (需要 partner variant 信息)
    try:
        bp2_result = assess_bp2_from_partners(literature)
        if bp2_result and hasattr(bp2_result, 'comment') and bp2_result.comment:
            # Add PMIDs to comment
            bpmids = bp2_result.pmids if hasattr(bp2_result, 'pmids') else []
            bp2_result.comment = _format_comment_with_pmids(bp2_result.comment, bpmids or [])
        results["bp2"] = bp2_result
    except Exception as e:
        logger.warning(f"BP2 assessment failed: {e}")
        results["bp2"] = None

    # PS3 from literature - ensure PMIDs in comment
    if "ps3_literature" in results and results["ps3_literature"]:
        ps3_lit = results["ps3_literature"]
        if "reason" in ps3_lit and "evidences" in ps3_lit:
            # Extract PMIDs from functional study evidences
            lit_pmids = [e.get("pmid", "") for e in ps3_lit.get("evidences", []) if e.get("pmid")]
            if lit_pmids:
                ps3_lit["reason"] = _format_comment_with_pmids(ps3_lit.get("reason", ""), lit_pmids)

    return results


def _assess_ps3_from_literature(
    functional_evidences: list,
) -> dict[str, Any]:
    """
    从功能研究文献证据评估 PS3

    取所有功能文献中最高证据强度

    Args:
        functional_evidences: List of FunctionalStudyEvidence

    Returns:
        dict with PS3 assessment result
    """
    if not functional_evidences:
        return {
            "applicable": False,
            "strength": None,
            "score": 0,
            "reason": "No functional literature"
        }

    # 统计功能研究结果 (LOF + GOF = affected function)
    affected_count = sum(1 for e in functional_evidences
                          if e.functional_result.value in ("lof", "gof"))
    lof_count = sum(1 for e in functional_evidences if e.functional_result.value == "lof")
    gof_count = sum(1 for e in functional_evidences if e.functional_result.value == "gof")
    total = len(functional_evidences)

    if affected_count == 0:
        return {
            "applicable": False,
            "strength": None,
            "score": 0,
            "reason": f"No functional literature supporting affected function ({total} studies)"
        }

    # 取最高置信度
    highest_confidence = max(e.confidence for e in functional_evidences)

    # 证据强度：high confidence = Moderate (2分), medium/low = Supporting (1分)
    if highest_confidence == "high":
        strength = "MODERATE"
        score = 2
    else:
        strength = "SUPPORTING"
        score = 1

    return {
        "applicable": True,
        "strength": strength,
        "score": score,
        "lof_count": lof_count,
        "gof_count": gof_count,
        "affected_count": affected_count,
        "total_studies": total,
        "confidence": highest_confidence,
        "reason": f"{affected_count}/{total} functional studies support affected function (LOF: {lof_count}, GOF: {gof_count})",
        "evidences": [
            {"pmid": e.pmid, "technique": e.technique, "confidence": e.confidence,
             "result": e.functional_result.value}
            for e in functional_evidences
        ]
    }


def get_final_ps3(
    ps3_db: Optional[dict],
    ps3_literature: Optional[dict],
) -> dict[str, Any]:
    """
    综合 PS3 最终结果

    取功能数据库结果和功能文献结果的最高分值

    Args:
        ps3_db: PS3 from SpliceVarDB/functional_evidence.db
        ps3_literature: PS3 from functional literature

    Returns:
        dict with final PS3 assessment including score
    """
    # Score hierarchy: VERY_STRONG(8) > STRONG(4) > MODERATE(2) > SUPPORTING(1)
    # 如果 ps3_db 没有返回 score，默认按 STRONG=4
    def get_score(result: dict) -> int:
        if not result or not result.get("applicable"):
            return 0
        return result.get("score", 4)  # 默认 STRONG=4

    best_result = {
        "applicable": False,
        "strength": None,
        "score": 0,
        "source": "none",
        "reason": ""
    }

    # Check database result
    db_score = get_score(ps3_db)
    if db_score > best_result["score"]:
        best_result = (ps3_db or {}).copy()
        best_result["source"] = "database"

    # Check literature result
    lit_score = get_score(ps3_literature)
    if lit_score > best_result["score"]:
        best_result = (ps3_literature or {}).copy()
        best_result["source"] = "literature"

    if best_result["score"] > 0:
        best_result["applicable"] = True

    return best_result


def retrieve_and_assess_literature(
    normalized: VariantInfo,
    inheritance_pattern: InheritancePattern = InheritancePattern.UNKNOWN,
    max_results_per_query: int = 20,
    max_total_articles: int = 100,
    use_cache: bool = True,
) -> tuple[VariantLiterature, dict[str, Any], bool]:
    """
    检索并评估文献证据的完整流程

    Args:
        normalized: VariantInfo from normalizer.normalize()
        inheritance_pattern: 遗传模式
        max_results_per_query: 每个检索词最大返回数
        max_total_articles: 总计最大文章数
        use_cache: 是否使用缓存 (默认 True)

    Returns:
        tuple: (VariantLiterature, evidence_assessment_dict, used_expanded_search)
    """
    # 获取 variant_id 用于缓存
    variant_id = normalized.rs_id or normalized.hgvs_c or ""
    gene = normalized.gene or ""

    # Step 0: 检查缓存
    if use_cache and variant_id and gene:
        cached_evidence = get_cached_evidence(
            variant_id, gene,
            rules=["PS2", "PM3", "PP1", "PS4", "BP2"]
        )
        if cached_evidence:
            logger.info(f"Using cached literature evidence for {variant_id}/{gene}")
            # 构建返回格式 (模拟 literature 和 evidence)
            literature = VariantLiterature(
                variant_id=variant_id,
                gene=gene,
                articles=[],
            )
            evidence = cached_evidence
            evidence["used_expanded_search"] = False
            return literature, evidence, False

    # Step 1: 检索文献
    literature, used_expanded_search = search_literature_for_variant(
        normalized,
        inheritance_pattern=inheritance_pattern,
        max_results_per_query=max_results_per_query,
        max_total_articles=max_total_articles,
    )

    # Step 2: 评估证据
    evidence = assess_literature_evidence(literature, inheritance_pattern)

    # 添加扩展检索标记到 evidence
    evidence["used_expanded_search"] = used_expanded_search

    # Step 3: 缓存证据结果
    if use_cache and variant_id and gene:
        try:
            cache_evidence(variant_id, gene, {
                "articles": literature.articles if hasattr(literature, 'articles') else [],
                "individual_observations": [],
                "evidence_results": evidence,
                "summary": {
                    "total_articles": len(literature.articles) if hasattr(literature, 'articles') else 0,
                    "has_de_novo_evidence": evidence.get("ps2", {}).get("applicable", False),
                    "has_cosegregation_evidence": evidence.get("pp1", {}).get("applicable", False),
                    "has_case_control_evidence": evidence.get("ps4", {}).get("applicable", False),
                }
            })
            logger.info(f"Cached literature evidence for {variant_id}/{gene}")
        except Exception as e:
            logger.warning(f"Failed to cache literature evidence: {e}")

    return literature, evidence, used_expanded_search


# ========== 便捷函数 ==========

def quick_literature_search(
    query_type: str,
    input_string: str,
    gene_symbol: Optional[str] = None,
    inheritance_pattern: InheritancePattern = InheritancePattern.UNKNOWN,
    use_cache: bool = True,
) -> tuple[VariantInfo, VariantLiterature, dict[str, Any], bool]:
    """
    便捷函数：一步完成 normalize + 文献检索 + 证据评估

    Args:
        query_type: 输入类型 (rsid, vcf, position)
        input_string: 输入字符串
        gene_symbol: 基因符号（用于辅助解析 HGVS p. 格式）
        inheritance_pattern: 遗传模式
        use_cache: 是否使用缓存 (默认 True)

    Returns:
        tuple: (normalized_variant, literature, evidence_assessment, used_expanded_search)
    """
    # Normalize
    normalizer = VariantNormalizer()
    normalized = normalizer.normalize(
        query_type=query_type,
        input_string=input_string,
        gene_symbol=gene_symbol,
    )

    # Retrieve and assess literature
    literature, evidence, used_expanded = retrieve_and_assess_literature(
        normalized,
        inheritance_pattern=inheritance_pattern,
        use_cache=use_cache,
    )

    return normalized, literature, evidence, used_expanded
