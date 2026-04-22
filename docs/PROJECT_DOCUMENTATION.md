# ACMG Variant Classification System - 项目文档

## 目录

1. [项目概述](#项目概述)
2. [ACMG证据规则详解](#acmg证据规则详解)
3. [数据流框架](#数据流框架)
4. [项目架构](#项目架构)
5. [文献检索模块](#文献检索模块)
6. [配置说明](#配置说明)
7. [预测工具集成](#预测工具集成)
8. [数据库资源](#数据库资源)

---

## 项目概述

本项目是基于 **ACMG (American College of Medical Genetics and Genomics)** 指南的变异分类系统，采用**贝叶斯评分**和多规则综合评估方法。

### 核心特性

- 支持 VEP/MyVariant API 自动化注释
- 支持 ClinVar 数据库检索
- 支持多种预测工具 (REVEL, SpliceAI, SIFT, PolyPhen)
- 支持基因特异性 ACMG 规则 (ClinGen)
- 集成 RAG + LLM 证据调整
- 可配置的规则引擎

---

## ACMG证据规则详解

### 致病性证据 (Pathogenic Evidence)

| 规则 | 名称 | 权重 | 数据来源 | 实现文件 |
|------|------|------|----------|----------|
| **PVS1** | 非常强 | +8 | LOF变异 (NMD/关键区域) | `pvs1*.py`, `pvs1_general.py` |
| **PS1** | 强 | +4 | ClinVar 同氨基酸改变 | `ps1*.py` |
| **PS2** | 强 | +4 | de novo 变异确认 | `ps2.py` |
| **PS3** | 强 | +4 | 功能学研究 | `ps3*.py` |
| **PS4** | 强 | +4 | 病例对照研究 | 文献检索 |
| **PM1** | 中等 | +2 | 热点/关键区域 | `pm1*.py` |
| **PM2** | 中等 | +2 | gnomAD 缺失 | `pm2*.py` |
| **PM3** | 中等 | +2 | 反式位置致病 | `pm3*.py` |
| **PM4** | 中等 | +2 | 蛋白长度改变 | `pm4*.py` |
| **PM5** | 中等 | +2 | 同位置不同错义 | `pm5*.py` |
| **PP1** | 支持 | +1 | 共分离分析 | `pp1.py` |
| **PP2** | 支持 | +1 | gnomAD missense z-score | `pp2.py` |
| **PP3** | 支持 | +1 | **预测工具 (REVEL等)** | `pp3*.py` |
| **PP4** | 支持 | +1 | 表型匹配 | 文献检索 |

### 良性证据 (Benign Evidence)

| 规则 | 名称 | 权重 | 数据来源 | 实现文件 |
|------|------|------|----------|----------|
| **BA1** | 独立良性 | -8 | gnomAD 高频 | `ba1*.py` |
| **BS1** | 强 | -4 | gnomAD 频率>病例 | `bs1*.py` |
| **BS2** | 强 | -4 | **FLOSSIES** 健康对照 | `bs2*.py` |
| **BS3** | 强 | -4 | 功能研究无损伤 | `bs3.py` |
| **BS4** | 强 | -4 | 非分离证据 | `bs4.py` |
| **BP1** | 支持 | -1 | 良性基因错义变异 | `bp1.py` |
| **BP2** | 支持 | -1 | 反式位置良性 | `bp2.py` |
| **BP3** | 支持 | -1 | 非关键区域框移 | `bp3.py` |
| **BP4** | 支持 | -1 | **预测工具 (良性)** | `bp4*.py` |
| **BP5** | 支持 | -1 | 另一解释变异 | `bp5.py` |
| **BP6** | 支持 | -1 | **ClinGen 可靠来源** | `bp6.py` |
| **BP7** | 支持 | -1 | 同义变异无剪接影响 | `bp7.py` |

### 分类阈值 (贝叶斯评分)

| 得分范围 | 分类 | 英文 |
|----------|------|------|
| ≥ 10 | 致病 | Pathogenic (Class 5) |
| 6-9 | 可能致病 | Likely Pathogenic (Class 4) |
| 0-5 | 意义不确定 | Uncertain Significance (Class 3) |
| -6 to -1 | 可能良性 | Likely Benign (Class 2) |
| ≤ -7 | 良性 | Benign (Class 1) |

---

## PVS1 证据规则详解

### PVS1 文件说明

项目中有两个 PVS1 实现文件:

| 文件 | 用途 | 特点 |
|------|------|------|
| `pvs1.py` | 基因特异性PVS1 | 依赖基因特异性配置 (disease relevant region BED文件) |
| `pvs1_general.py` | 通用型PVS1 | 不依赖基因特异性配置，适用于任何基因 |

### PVS1 强度判定逻辑

**核心判断因素:**
1. **变异类型**: 必须是 NULL/LOF 变异 (stop_gained, frameshift, splice_donor/acceptor, start_lost)
2. **NMD (无义介导的mRNA降解)**: 是否导致NMD
3. **关键区域**: 截断区域是否在疾病相关区域内
4. **蛋白长度变化**: 非NMD情况下，评估蛋白长度变化百分比

**强度判定矩阵:**

| 情况 | PVS1 | PVS1_general |
|------|------|--------------|
| NULL变异 + NMD | Very Strong | Very Strong |
| NULL变异 + 非NMD + 疾病相关区域 | Strong | Strong |
| NULL变异 + 非NMD + 非疾病相关区域 + 蛋白长度大 | Strong | Moderate |
| NULL变异 + 非NMD + 非疾病相关区域 + 蛋白长度小 | Moderate | Moderate |

### 蛋白长度变化判断

**计算方式:**
- `diff_len_protein_percent`: 蛋白长度变化百分比
- 计算公式: `变化后长度 / 原始长度 × 100%`

**判断阈值:**
- 阈值配置位置: `config.yaml` → `functional_thresholds.threshold_diff_len_prot_percent`
- 默认阈值: 需参考具体基因配置
- 当 `diff_len_protein_percent > threshold` 时，认为蛋白长度变化显著

**示例:**
```
原始蛋白长度: 1000 aa
变异后长度: 950 aa
蛋白长度变化: 50 aa
diff_len_protein_percent = 50/1000 × 100% = 5%

如果阈值 threshold_diff_len_prot_percent = 10%
则 5% < 10%，认为蛋白长度变化较小
```

**蛋白长度数据来源:**
- 参考蛋白序列长度来自 `pyensembl` (ENSEMBL基因组注释)
- 每个转录本有独立的蛋白序列长度
- 项目使用 MANE (Matched Annotation from NCBI and EMBL-EBI) 转录本作为标准
- 蛋白长度计算在预处理阶段完成，存储在 `TranscriptInfo.diff_len_protein_percent` 属性中

### 使用建议

- **已有基因特异性配置** (BRCA1, BRCA2, ATM, PTEN, PALB2, CDH1等): 使用 `pvs1`, `pvs1_brca1`, `pvs1_brca2` 等
- **新基因或无特异性配置**: 使用 `pvs1_general`

---

## 数据流框架

### 完整数据流图

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              输入层 (Input Layer)                            │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   ┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐  │
│   │    rsID     │    │  VCF格式    │    │   HGVS     │    │  直接JSON   │  │
│   │ rs123456    │    │ 17:41245:G:A│    │c.68G>A     │    │  预定义数据  │  │
│   └──────┬──────┘    └──────┬──────┘    └──────┬──────┘    └──────┬──────┘  │
│          │                  │                  │                  │         │
└──────────┼──────────────────┼──────────────────┼──────────────────┼─────────┘
           │                  │                  │                  │
           ▼                  ▼                  ▼                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                          normalizer.py                                      │
│                    (VEP/MyVariant API 统一标准化)                            │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   1. 解析输入格式                                                            │
│   2. 调用 VEP REST API 或 MyVariant.info                                    │
│   3. 获取:                                                                  │
│      - 染色体坐标                                                           │
│      - HGVS 命名 (c., p.)                                                  │
│      - 基因和转录本                                                          │
│      - gnomAD 频率                                                          │
│      - 预测工具分数 (REVEL, SpliceAI, SIFT, PolyPhen) ◄── 新增支持           │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                       variant_converter.py                                   │
│                      (转换为内部 JSON 格式)                                   │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   convert_normalizer_to_variant_json()                                       │
│                                                                             │
│   输出结构:                                                                  │
│   {                                                                         │
│     "chr": "17",                                                            │
│     "pos": 43045678,                                                        │
│     "ref": "G",                                                             │
│     "alt": "A",                                                             │
│     "gene": "BRCA1",                                                        │
│     "variant_type": ["missense_variant"],                                    │
│     "pathogenicity_prediction_tools": {"REVEL": 0.85},  ◄── 新增            │
│     "splicing_prediction_tools": {"SpliceAI": 0.92},  ◄── 新增             │
│     "gnomAD": {...},                                                         │
│     ...                                                                     │
│   }                                                                         │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         load_variant.py                                      │
│                        (创建 Variant 对象)                                   │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   create_variant(variant_json) → Variant 对象                                 │
│                                                                             │
│   ┌─────────────────────────────────────────────────────────────────────┐   │
│   │                        Variant 对象                                  │   │
│   ├─────────────────────────────────────────────────────────────────────┤   │
│   │  variant_info: VariantInfo        - 变异基本信息                      │   │
│   │  transcript_info: [TranscriptInfo]- 转录本列表                       │   │
│   │  prediction_tools: dict          - 预测工具分数                      │   │
│   │  gnomad_popmax: PopulationDatabases_gnomAD                          │   │
│   │  gnomad_faf: PopulationDatabases_gnomAD                              │   │
│   │  flossies: PopulationDatabases     - FLOSSIES 数据库                 │   │
│   │  cancerhotspots: PopulationDatabases - 癌症热点                       │   │
│   │  functional_assay: [FunctionalData]  - 功能数据                       │   │
│   │  splicing_assay: [RNAData]         - RNA 数据                        │   │
│   └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      config_annotation.py                                    │
│                       (注解与规则配置)                                       │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   步骤1: 根据规则列表确定需要的注解 (Annotations)                              │
│   ┌─────────────────────────────────────────────────────────────────────┐   │
│   │  RULE_DICTIONARY = {                                                │   │
│   │    "pp3_protein": Pp3_protein,                                      │   │
│   │    "bp4_protein": Bp4_protein,                                      │   │
│   │    "bp6_clingen": Bp6_clingen,  ◄── 新增                            │   │
│   │    ...                                                              │   │
│   │  }                                                                  │   │
│   └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
│   步骤2: 获取注解数据                                                        │
│   ┌─────────────────────────────────────────────────────────────────────┐   │
│   │  注解类型:                                                           │   │
│   │  1. Variant内置数据 (prediction_tools, gnomad, flossies)            │   │
│   │  2. 文件注解 (ClinVar, SpliceAI, Hotspot)                           │   │
│   │  3. 配置阈值 (thresholds)                                           │   │
│   │  4. 路径配置 (BED文件, 数据库路径)                                    │   │
│   └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        acmg_rules/*                                         │
│                          (证据规则评估)                                      │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   ┌──────────────────┐  ┌──────────────────┐  ┌──────────────────┐        │
│   │  计算证据规则     │  │  频率证据规则     │  │  文献证据规则     │        │
│   ├──────────────────┤  ├──────────────────┤  ├──────────────────┤        │
│   │  PP3 (预测→致病) │  │  PM2 (gnomAD缺失)│  │  PS4 (病例对照)  │        │
│   │  BP4 (预测→良性) │  │  BA1 (gnomAD高频)│  │  PP1 (共分离)    │        │
│   │  BP7 (同义+剪接) │  │  BS2 (FLOSSIES)  │  │  PP4 (表型匹配)  │        │
│   └──────────────────┘  └──────────────────┘  └──────────────────┘        │
│                                                                             │
│   ┌──────────────────┐  ┌──────────────────┐                              │
│   │  ClinVar证据规则  │  │  基因特异性规则   │                              │
│   ├──────────────────┤  ├──────────────────┤                              │
│   │  PS1 (同AA改变) │  │  PVS1 (基因特异)  │                              │
│   │  PM5 (不同AA)   │  │  PM1 (热点区域)  │                              │
│   │  BP6 (ClinGen)  │  │  ...             │                              │
│   └──────────────────┘  └──────────────────┘                              │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                   classification_schemata/                                    │
│                         (贝叶斯评分)                                         │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   rule_results: [RuleResult, RuleResult, ...]                               │
│                              │                                              │
│                              ▼                                              │
│   ┌─────────────────────────────────────────────────────────────────────┐  │
│   │                     Bayes Score Calculator                          │  │
│   │                                                                     │  │
│   │   for rule in rule_results:                                        │  │
│   │       if rule.evidence_type == PATHOGENIC:                          │  │
│   │           score += rule.strength.value                              │  │
│   │       else:  # BENIGN                                               │  │
│   │           score -= rule.strength.value                              │  │
│   │                                                                     │  │
│   │   classification = get_classification(score)                        │  │
│   └─────────────────────────────────────────────────────────────────────┘  │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                           输出层                                             │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   {                                                                         │
│     "classification": "Likely Pathogenic",                                 │
│     "class": 4,                                                            │
│     "bayes_score": 7,                                                      │
│     "rule_results": [                                                      │
│       {"rule": "PP3", "result": true, "strength": "SUPPORTING"},          │
│       {"rule": "BP4", "result": false, "strength": "SUPPORTING"},         │
│       ...                                                                  │
│     ],                                                                     │
│     "evidence_summary": {...}                                              │
│   }                                                                         │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 项目架构

### 目录结构

```
variant_classification/
├── normalizer.py                    # VEP/MyVariant API 统一接口
├── variant_converter.py              # 格式转换 (normalizer → 内部JSON)
├── load_variant.py                  # Variant 对象创建
├── variant.py                        # 核心数据类定义
├── var_type.py                       # 变异类型枚举
├── information.py                     # Info/Classification_Info 数据类
│
├── config_annotation.py              # 注解配置与规则引擎
├── config_schema.json                # 配置文件 schema
│
├── acmg_rules/                       # ACMG 规则实现
│   ├── utils.py                      # RuleResult, evidence_strength 等
│   ├── computation_evidence_utils.py  # 阈值评估工具
│   ├── pp3.py                        # PP3: 预测工具 (致病)
│   ├── bp4.py                        # BP4: 预测工具 (良性)
│   ├── bp6.py                        # BP6: ClinGen 可靠来源 ◄── 新增
│   ├── bp7.py                        # BP7: 同义变异
│   ├── pm1.py                        # PM1: 热点区域
│   ├── pm2.py                        # PM2: gnomAD 缺失
│   ├── pm5.py                        # PM5: 同位置不同错义
│   ├── ps1.py                        # PS1: 相同氨基酸改变
│   ├── ps3.py                        # PS3: 功能学研究
│   ├── pvs1*.py                      # PVS1: 基因特异性 (依赖BED文件)
│   ├── pvs1_general.py                # PVS1: 通用型 (不依赖基因特异性配置)
│   ├── pp2.py                        # PP2: Missense 评估
│   ├── ba1.py, bs1.py, bs2.py        # BA1/BS1/BS2: 频率证据
│   ├── clinvar_annot.py              # ClinVar 注解 (支持 BP6)
│   ├── clinvar_annot_spliceai.py     # ClinVar+SpliceAI 联合注解
│   ├── clinvar_missense.py           # ClinVar 错义变异处理
│   ├── clinvar_splicing.py           # ClinVar 剪接变异处理
│   └── clinvar_utils.py              # ClinVar 数据类 (ClinVarReviewDetails)
│
├── literature_retrieval/            # 文献检索模块 ◄── 新增
│   ├── __init__.py                  # 模块导出
│   ├── literature_utils.py          # 数据类定义
│   ├── ncbi_retriever.py            # NCBI PubMed 检索
│   ├── pubtator_api.py              # PubTator API 客户端
│   ├── evidence_router.py           # 证据路由
│   ├── pm3_evaluator.py            # PM3 LLM 评估器
│   ├── ps2_evaluator.py            # PS2 LLM 评估器 (表型异质性增强)
│   ├── segregation_evaluator.py    # PP1/BS4 共分离评估器 ◄── 新增
│   ├── bp2_evaluator.py             # BP2 评估器 ◄── 新增
│   └── literature_annot.py          # 主注解入口
│
├── classification_schemata/          # 分类方案
│   ├── bayes_scores.py                # 贝叶斯评分计算
│   └── utils.py                      # 分类工具
│
├── ClinGen/                          # ClinGen 规则数据
│   └── clingen_rules.db              # 基因特异性 ACMG 规则
│
├── cloud_api/                        # 云 API 服务
│   └── vep_client.py                 # VEP API 客户端
│
├── functional_evidence.db            # 功能证据数据库
├── clingen_rules.db                  # ClinGen 规则数据库
│
└── docs/                             # 文档
    └── PROJECT_DOCUMENTATION.md      # 本文档
```

### 核心数据流

```
输入 → normalizer.py → variant_converter.py → load_variant.py → Variant对象
                                                              │
                                                              ▼
                                                    config_annotation.py
                                                              │
                                                    根据规则配置选择需要的注解
                                                              │
                                                              ▼
                                            ┌─────────────────┴─────────────────┐
                                            │                                   │
                                   ┌─────────▼─────────┐            ┌───────────▼──────────┐
                                   │  Variant 内置数据   │            │    外部注解函数       │
                                   │  - prediction_tools │            │  - ClinVar 注解      │
                                   │  - gnomad          │            │  - SpliceAI 注解     │
                                   │  - flossies        │            │  - Hotspot 注解      │
                                   └───────────────────┘            └──────────────────────┘
                                                              │
                                                              ▼
                                                    acmg_rules/* (规则评估)
                                                              │
                                                              ▼
                                                    classification_schemata (贝叶斯)
                                                              │
                                                              ▼
                                                          最终分类
```

---

## 文献检索模块

### 概述

文献检索模块 (`literature_retrieval/`) 从多个数据源检索文献证据，支持 ACMG 分类规则 PS4, PP1, PP4, PS2 等。

### 模块架构

```
literature_retrieval/
├── __init__.py                 # 模块导出
├── literature_utils.py         # 数据类定义
├── ncbi_retriever.py           # NCBI PubMed 检索
├── pubtator_api.py            # PubTator 变异注释 API
├── evidence_router.py          # 证据路由
├── pm3_evaluator.py            # PM3 LLM 评估器
├── ps2_evaluator.py            # PS2 LLM 评估器
└── literature_annot.py          # 主注解入口
```

### 数据源

| 数据源 | 功能 | 用途规则 |
|--------|------|----------|
| **NCBI PubMed** | 文献检索 | PS4, PP1, PP4, PS2 |
| **PubTator 3 API** | 变异注释 | PS1, PM5, BP6 |

### 使用方式

```python
from literature_retrieval import (
    NCBILiteratureRetriever,
    EvidenceRouter,
    LiteratureConfig,
    InheritancePattern,
)
from normalizer import VariantNormalizer

# 1. 配置检索
config = LiteratureConfig(
    ncbi_api_key="your_api_key",
    max_total_articles=100,
)

# 2. 检索文献
normalizer = VariantNormalizer()
variant_info = normalizer.normalize("rsid", "rs123456")

retriever = NCBILiteratureRetriever(config=config)
articles = retriever.search_variant_literature(variant_info)

# 3. 路由证据
router = EvidenceRouter()
evidence = router.route_evidence(literature, InheritancePattern.AUTOSOMAL_DOMINANT)

# 4. 快速检索
from literature_retrieval import quick_search
lit = quick_search("rs123456", gene="BRCA1")
```

### 证据路由

| 遗传模式 | 适用规则 |
|----------|----------|
| **AD (常显)** | PS2, PS4, PP1, PP4 |
| **AR (常隐)** | PM3, PP1, PP4 |
| **X连锁** | PS2, PP1, PP4, PM3 |
| **未知** | PP1, PP4 |

### PM3 LLM 评估器 (`pm3_evaluator.py`)

基于 ClinGen SVI PM3 细则，使用 LLM 评估 PM3 证据强度。

#### PM3 证据强度矩阵

**基础强度**（由变异配置类型决定）：

| 配置类型 | 基础强度 |
|----------|----------|
| 反式位置 + 已知致病 + 功能研究支持 LoF | Very Strong (VS) |
| 反式位置 + 已知致病（经分子测试确认） | Strong (S) |
| 反式位置 + 可能致病（经分子测试确认） | Moderate (M) |
| 纯合子/半合子，符合疾病遗传模式 | Supporting (P) |

**病例数升级**（可提升证据强度）：

| 独立先证者数 | 升级强度 |
|--------------|----------|
| 1 个先证者 | Supporting (P) |
| 2 个先证者 | Moderate (M) |
| 3 个先证者 | Strong (S) |
| 4+ 个先证者 | Very Strong (VS) |

**重要**：最终强度取基础强度和升级强度的**较高值**

#### 使用示例

```python
from literature_retrieval import (
    PM3LLMEvaluator,
    InheritancePattern,
)

# 评估 PM3
evaluator = PM3LLMEvaluator(llm_api_key="your_api_key")
assessment = evaluator.assess(
    gene="BRCA2",
    variant_description="c.1000C>T",
    case_reports=lit.case_reports,
    inheritance_pattern=InheritancePattern.AUTOSOMAL_RECESSIVE,
    known_pathogenic_variants=[
        {"description": "c.5000G>A", "significance": "Pathogenic"}
    ],
)

print(f"PM3 applicable: {assessment.applicable}")
print(f"Final strength: {assessment.strength.value}")
print(f"Base strength: {assessment.base_strength}")
print(f"Independent probands: {assessment.num_independent_probands}")
print(f"Confidence: {assessment.confidence}")
print(f"Reasoning: {assessment.reasoning}")
```

#### 评估结果

```python
@dataclass
class PM3Assessment:
    applicable: bool
    strength: PM3Strength  # 最终强度 (VS/S/M/P/NA)
    base_strength: PM3Strength  # 基础强度

    # 配置类型
    is_trans_config: bool  # 复合杂合
    is_homozygous: bool    # 纯合子
    is_hemizygous: bool    # 半合子

    # 证据详情
    num_independent_probands: int  # 独立先证者数
    proband_details: List[Dict]   # 各先证者详情
    known_pathogenic_variant: str  # 反式位置的已知致病变异
    has_lof_functional_study: bool  # 功能研究支持

    confidence: str  # high, medium, low
    reasoning: str   # LLM 推理过程
```

### PS2 LLM 评估器 (`ps2_evaluator.py`)

基于 ClinGen SVI PS2/PM6 de novo 细则，使用 LLM 评估 de novo 证据强度。

#### PS2/PM6 证据规则

| 规则 | 说明 | 基础证据 |
|------|------|----------|
| **PS2** | 已确认 de novo (父母均阴性) | 1 point |
| **PM6** | 推定 de novo (未确认) | 0.5 points |

#### 先证者数升级

| 独立先证者数 | 升级强度 |
|--------------|----------|
| 1 个先证者 | Supporting (P) |
| 2 个先证者 | Moderate (M) |
| 3 个先证者 | Strong (S) |
| 4+ 个先证者 | Very Strong (VS) |

#### 表型异质性调整

表型异质性指病例呈现相似还是不同的表型特征，对PS2强度调整至关重要。

**表型一致性级别:**

| 一致性级别 | 说明 | PS2调整 |
|------------|------|---------|
| **HIGH (高一致)** | 病例来自相同疾病谱，核心特征一致 | 强度可维持或升级 |
| **MEDIUM (中等)** | 存在表型变异，但可能属于同一疾病谱的变异性表达 | 强度维持 |
| **LOW (低一致)** | 病例来自不同疾病系统，特征不一致 | 可能降低强度或不适用 |

**表型异质性评估方法:**
1. 识别每个病例的主要表型（导致诊断的主要特征）
2. 跨病例比较：主要表型是否属于同一疾病谱？
3. 考虑疾病特异性模式：某些疾病的变异性表达仍属同一疾病
4. 检查表型是否能用单一疾病解释

**高一致性示例:**
| 病例1表型 | 病例2表型 | 病例3表型 | 评估 |
|-----------|-----------|-----------|------|
| 发育迟缓,癫痫 | 发育迟缓,癫痫 | 发育迟缓 | HIGH - 同疾病 |
| 长QT综合征,晕厥 | 长QT综合征,心脏骤停 | 长QT综合征 | HIGH - 同疾病 |

**低异质性示例 (不同疾病):**
| 病例1表型 | 病例2表型 | 病例3表型 | 评估 |
|-----------|-----------|-----------|------|
| 癫痫 | 心肌病 | 骨骼发育不良 | LOW - 不同系统 |
| 神经发育障碍 | 结缔组织疾病 | 肾脏疾病 | LOW - 不同疾病 |

**重要提示:**
- 相同基因、不同表型：可能仍是同一疾病的变异性表达（如RASopathies）
- 评估表型时需考虑特定基因及其已知疾病谱

#### 使用示例

```python
from literature_retrieval import (
    PS2LLMEvaluator,
    InheritancePattern,
)

# 评估 PS2
evaluator = PS2LLMEvaluator(llm_api_key="your_api_key")
assessment = evaluator.assess(
    gene="TP53",
    variant_description="c.1000C>T",
    case_reports=lit.case_reports,
    inheritance_pattern=InheritancePattern.AUTOSOMAL_DOMINANT,
)

print(f"PS2 applicable: {assessment.applicable}")
print(f"Rule: {assessment.rule}")  # PS2 or PM6
print(f"De novo status: {assessment.de_novo_status.value}")
print(f"Strength: {assessment.strength.value}")
print(f"Phenotypic consistency: {assessment.phenotype_consistency}")
print(f"Phenotype analysis: {assessment.phenotype_analysis}")
print(f"Confidence: {assessment.confidence}")
```

#### 评估结果

```python
@dataclass
class PS2Assessment:
    applicable: bool
    rule: str  # PS2, PM6, or NA
    strength: PS2Strength  # 最终强度 (VS/S/M/P/NA)

    # de novo 状态
    de_novo_status: DeNovoStatus  # confirmed, assumed, not_de_novo
    has_parental_testing: bool

    # 先证者数
    num_independent_probands: int
    num_families: int
    probands: List[Dict[str, Any]]  # 各先证者详情

    # 表型信息
    phenotype_consistency: str  # high, medium, low
    phenotype_details: str  # 表型描述
    phenotype_analysis: str  # 异质性评估解释

    confidence: str
    reasoning: str
comment: str  # 强度调整说明
```

### PP1/BS4 共分离评估器 (`segregation_evaluator.py`)

基于家系共分离数据，使用 LLM 或启发式方法评估 PP1（共分离）或 BS4（非分离）证据。

#### PP1/BS4 证据规则

**PP1 (共分离):**
- 变异在家系中与疾病共分离
- 证据强度由 LOD Score 决定

**BS4 (非分离):**
- 变异在家系中不与疾病共分离
- 常见于健康携带者

#### LOD Score 阈值 (PP1)

| LOD Score | 证据强度 |
|-----------|----------|
| LOD >= 2.0 | Very Strong (VS) |
| LOD >= 1.5 | Strong (S) |
| LOD >= 1.0 | Moderate (M) |
| LOD >= 0.5 | Supporting (P) |
| LOD < 0.5 | Not applicable |

#### LOD Score 阈值 (BS4)

| LOD Score | 证据强度 |
|-----------|----------|
| LOD <= -2.0 | Very Strong (VS) |
| LOD <= -1.5 | Strong (S) |
| LOD <= -1.0 | Moderate (M) |
| LOD <= -0.5 | Supporting (P) |
| LOD > -0.5 | Not applicable |

#### 使用示例

```python
from literature_retrieval import (
    SegregationEvaluator,
    FamilyMember,
    PedigreeData,
)

# 构建家系数据
members = [
    FamilyMember(individual_id="I-1", affected=False, genotype="non-carrier", relationship="parent"),
    FamilyMember(individual_id="I-2", affected=False, genotype="carrier", relationship="parent"),
    FamilyMember(individual_id="II-1", affected=True, genotype="carrier", relationship="proband"),
    FamilyMember(individual_id="II-2", affected=True, genotype="carrier", relationship="sibling"),
]

pedigree = PedigreeData.from_family_members(
    family_id="FAM001",
    proband_id="II-1",
    members=members
)

# 评估 PP1
evaluator = SegregationEvaluator()
assessment = evaluator.assess_pp1(
    gene="BRCA2",
    variant_description="c.1000C>T",
    pedigree_data=[pedigree],
    inheritance_pattern="autosomal dominant",
)

print(f"PP1 applicable: {assessment.applicable}")
print(f"LOD Score: {assessment.lod_score}")
print(f"Strength: {assessment.strength.value}")
print(f"Segregation status: {assessment.segregation_status.value}")
```

#### 评估结果

```python
@dataclass
class SegregationAssessment:
    applicable: bool
    rule: str  # PP1 or BS4
    segregation_status: SegregationStatus  # segregation/no_segregation/inconclusive
    lod_score: float
    strength: SegregationStrength

    num_affected_carriers: int
    num_affected_non_carriers: int
    num_unaffected_carriers: int
    num_unaffected_non_carriers: int

    families: List[str]
    confidence: str
    reasoning: str
```

### BP2 评估器 (`bp2_evaluator.py`)

BP2 用于常染色体隐性遗传疾病，当变异与已知良性变异呈反式（复合杂合）排列时。

#### BP2 证据规则

| 配置 | 证据强度 |
|------|---------|
| 已知良性变异反式排列 + 相位已确认 | Supporting (P) |
| 已知良性变异反式排列 + 多个病例 | Supporting (P) |
| 经分子检测确认 | Supporting (P) |

#### 使用示例

```python
from literature_retrieval import (
    BP2Evaluator,
    TransPartnerVariant,
)

# 定义反式位置的已知良性变异
partner = TransPartnerVariant(
    variant_description="c.5000G>A",
    clinical_significance="Benign",
    evidence_source="ClinVar",
    pmid=None
)

# 评估 BP2
evaluator = BP2Evaluator()
assessment = evaluator.assess(
    gene="BRCA2",
    variant_description="c.1000C>T",
    partner_variants=[partner],
    inheritance_pattern="autosomal recessive",
)

print(f"BP2 applicable: {assessment.applicable}")
print(f"Partner variant: {assessment.partner_variant.variant_description}")
print(f"Strength: {assessment.strength.value}")
```

### 进一步拓展建议

#### 1. PS4 (病例对照研究) - 优先级中

**现状**: 仅有权重赋值，无评估程序

**问题**:
- 需要病例对照数据来源
- 需计算病例组vs对照组频率差异
- 需考虑统计显著性

**建议实现**:
```python
# ps4_evaluator.py
class Ps4_case_control:
    """
    PS4: Prevalence in affected individuals significantly increased vs controls

    评估方法:
    - 从文献/数据库获取病例组和对照组变异频率
    - 计算比值比(OR)或相对风险
    - 根据置信区间和p值确定强度
    """
    # 需集成病例对照数据库或文献数据
```

**数据源依赖**:
- 疾病特异性变异频率数据库
- 公开GWAS数据
- 病例对照研究文献

#### 2. PP4 (表型匹配) - 优先级中

**现状**: 基于多因素似然数据

**问题**:
- 需要标准化表型评估
- 需HPO术语结构化
- 疾病特异表型谱定义

**建议实现**:
```python
# pp4_evaluator.py
class Pp4_phenotype_matcher:
    """
    PP4: Patient's phenotype or family history is highly specific for a disease

    评估方法:
    - 输入: 患者HPO表型术语
    - 比对: 疾病-表型数据库
    - 计算: 表型相似度评分
    - 输出: PP4证据强度
    """
```

**数据源依赖**:
- HPO (Human Phenotype Ontology) 数据库
- 疾病-表型映射 (如 OMIM, Orphan

```bash
# 设置 NCBI API Key (可选，提升速率限制)
export NCBI_API_KEY="your_key"
export NCBI_EMAIL="your_email@example.com"

# Python 使用
python3 -c "
from literature_retrieval import quick_search
results = quick_search('rs123456', gene='TP53')
print(f'Found {results.total_articles} articles')
print(f'Case reports: {results.num_case_reports}')
"
```

---

## 配置说明

### 配置文件结构

```yaml
# config.yaml
name: "BRCA1_Classification"

rules:
  - "pvs1_general"
  - "ps1_protein"
  - "pp3_protein"
  - "bp4_protein"
  - "bp6_clingen"      # 新增 BP6
  - "pm2"
  - "bs2"

prediction_tool_threshold:
  pathogenicity_prediction:
    name: "REVEL"
    benign:
      direction: "less"
      supporting: 0.1
      moderate: 0.2
      strong: 0.4
      very_strong: 0.6
    pathogenic:
      direction: "greater"
      supporting: 0.5
      moderate: 0.6
      strong: 0.7
      very_strong: 0.8
  splicing_prediction:
    name: "SpliceAI"
    benign:
      direction: "less"
      supporting: 0.2
      moderate: 0.3
      strong: 0.5
    pathogenic:
      direction: "greater"
      supporting: 0.5
      moderate: 0.6
      strong: 0.7

allele_frequency_thresholds:
  threshold_ba1: 0.001
  threshold_bs1: 0.001
  threshold_pm2: 0.0001
  threshold_bs2: 10

annotation_files:
  root: "/data/annotations"
  clinvar:
    root: "/data/clinvar"
    clinvar_snv: "clinvar_20240303.vcf.gz"
  critical_regions:
    root: "/data/regions"
    hotspot_region: "hotspots.bed"
    coldspot_region: "coldspots.bed"
```

---

## 预测工具集成

### 支持的预测工具

| 工具 | 类型 | 用途 | 数据来源 |
|------|------|------|----------|
| **REVEL** | 蛋白质致病性 | PP3/BP4 | VEP MyVariant (dbNSFP) |
| **SpliceAI** | 剪接效应 | PP3/BP4/BP7 | VEP MyVariant |
| **SIFT** | 蛋白质功能 | PP3/BP4 | VEP MyVariant (dbNSFP) |
| **PolyPhen** | 蛋白质功能 | PP3/BP4 | VEP MyVariant (dbNSFP) |

### normalizer.py 中的处理

```python
# VEP 参数配置
"dbnsfp": 1,      # 获取 SIFT, PolyPhen, REVEL
"spliceai": 1,    # 获取 SpliceAI

# 解析预测分数
if "dbnsfp" in transcript:
    revel = transcript["dbnsfp"].get("revel_score")
    sift = transcript["dbnsfp"].get("sift_score")
    polyphen = transcript["dbnsfp"].get("polyphen_score")

if "SpliceAI" in variant:
    spliceai = parse_spliceai(variant["SpliceAI"])
```

### 规则应用

```python
# PP3: 预测致病
if revel > 0.5:  # 或配置的阈值
    result = True  # PP3 阳性

# BP4: 预测良性
if revel < 0.1:  # 或配置的阈值
    result = True  # BP4 阳性
```

---

## 数据库资源

### 功能证据数据库 (functional_evidence.db)

| 表名 | 内容 | 用途规则 |
|------|------|----------|
| splicevardb | SpliceVarDB 剪接变异 | PS3, BP7 |
| functional_studies | 功能研究证据 | PS3, BS3 |

### ClinGen 规则数据库 (clingen_rules.db)

| 表名 | 内容 | 用途规则 |
|------|------|----------|
| gene_rules | 基因特异性 ACMG 规则 | PVS1, PM1, PP1 |
| expert_panels | Expert Panel 评审状态 | BP6 |

### ClinVar VCF 文件

| 文件 | 内容 | 用途规则 |
|------|------|----------|
| clinvar_snv.vcf.gz | SNV 临床注释 | PS1, PM5, BP6 |
| clinvar_indel.vcf.gz | InDel 临床注释 | PS1, PM5, BP6 |

### 基因组注释数据库

项目使用 `pyensembl` (ENSEMBL) 作为基因组注释来源。

#### ENSEMBL vs NCBI RefSeq 区别

| 特性 | ENSEMBL (pyensembl) | NCBI RefSeq |
|------|---------------------|-------------|
| **维护机构** | EMBL-EBI (欧洲) | NCBI (美国) |
| **基因ID格式** | ENSG00000xxxxx | GeneID (如 672) |
| **转录本ID** | ENST00000xxxxx | NM_xxxxxx |
| **蛋白ID** | ENSP00000xxxxx | NP_xxxxxx |
| **注释策略** | 包含更多转录本亚型 | 相对保守 |
| **更新频率** | 定期与NCBI同步 | 基准来源 |

#### 具体差异示例 (BRCA1基因)

| 属性 | ENSEMBL | NCBI |
|------|---------|------|
| 基因ID | ENSG00000012048 | 672 |
| 转录本 | ENST00000352993 | NM_007294 |
| 蛋白 | ENSP00000384933 | NP_009225 |

#### 注释差异说明

1. **转录本数量**: ENSEMBL通常包含更多转录本亚型，NCBI相对保守
2. **基因边界**: 两者偶有细微差异
3. **剪接位点**: 可能有不同的5'/3'端
4. **更新同步**: 两者会互相参考，但有延迟

#### 临床应用建议

| 场景 | 推荐来源 |
|------|----------|
| 临床报告 | NCBI RefSeq (更保守) |
| 研究分析 | 两者皆可，注意ID转换 |
| 变异注释 | 使用两者的并集或交叉验证 |

#### 数据使用说明

| 数据库 | 用途 | 说明 |
|--------|------|------|
| pyensembl | 获取转录本蛋白序列 | 用于计算 `diff_len_protein_percent` |
| ClinVar | 临床变异注释 | 用于 PS1, PM5, BP6 等规则 |
| gnomAD | 人群频率数据 | 使用ENSEMBL格式的基因ID |
| MANE | 官方推荐转录本 | 选择标准转录本进行最终评估 |

---

## BP6 特殊说明

### BP6 规则 (新增)

**BP6: Reputable source (benign classification)**

触发条件:
1. ClinGen/Expert Panel 已分类为 B (Benign) 或 LB (Likely Benign)
2. 多家独立实验室 (>2家) 提交了仅 B/LB 的分类

```python
@dataclass
class ClinVarReviewDetails:
    variation_id: str
    review_status: str          # "reviewed by expert panel"
    clinical_significance: str  # "Benign" 或 "Likely benign"
    is_clingen_classified: bool # ClinGen 是否已分类
    submissions: [ClinVarSubmission]  # 提交详情

    @property
    def meets_bp6_criteria(self) -> bool:
        # 条件1: ClinGen Expert Panel 分类
        if self.is_clingen_classified and self.is_benign:
            return True
        # 条件2: >2家实验室独立提交
        if len(self.benign_submissions) >= 3:
            return True
        return False
```

---

## 快速开始

### 1. 安装依赖

```bash
pip install -r requirements.txt
```

### 2. 配置

```bash
cp config_example.yaml config.yaml
# 编辑 config.yaml 设置文件路径和阈值
```

### 3. 使用示例

```python
from load_variant import from_normalizer
from config_annotation import get_annotations_needed_from_rules, apply_rules
from classification_schemata import calculate_classification

# 方式1: 通过 rsID
variant = from_normalizer("rsid", "rs123456")

# 方式2: 通过 VCF 格式
variant = from_normalizer("vcf", "17:43045678:G:A")

# 获取需要的注解并评估规则
rule_list = ["pp3_protein", "bp4_protein", "bp6_clingen", "pm2"]
annotations = get_annotations_needed_from_rules(rule_list, class_info)

# 应用规则
results = apply_rules(rule_info_dict)

# 计算最终分类
classification = calculate_classification(results)
```

---

## 版本历史

| 版本 | 日期 | 变更 |
|------|------|------|
| 0.2.0 | 2026-04-15 | 新增 BP6 规则, REVEL/SpliceAI 支持 |
| 0.1.0 | 2026-04-10 | 初始版本 |
