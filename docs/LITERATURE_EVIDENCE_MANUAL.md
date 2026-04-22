# 文献检索与 LLM 证据评估系统使用手册

## 版本信息
- **更新日期**: 2026-04-19
- **版本**: 1.0

---

## 目录

1. [系统概述](#1-系统概述)
2. [安装与配置](#2-安装与配置)
3. [文献检索流程](#3-文献检索流程)
4. [LLM 证据评估](#4-llm-证据评估)
5. [证据规则详解](#5-证据规则详解)
6. [配置示例](#6-配置示例)
7. [故障排除](#7-故障排除)

---

## 1. 系统概述

### 1.1 功能介绍

本系统对标准化后的变异自动触发多策略文献检索，并使用 LLM 对文献进行分类和证据提取，评估以下 ACMG 证据：

| 证据 | 规则 | 说明 |
|------|------|------|
| PS2 | de novo | 新发变异（需父母验证） |
| PM3 | biallelic | 双等位基因变异（隐性遗传） |
| PP1 | co-segregation | 家系共分离 |
| BP2 | trans | 反式位置（良性证据） |
| PS3 | functional | 功能研究证据 |

### 1.2 整体架构

```
输入: VariantInfo (normalize 后的变异信息)
     ↓
build_search_queries() 生成检索词
     ↓
NCBILiteratureRetriever.search() 检索 PubMed
     ↓
LiteratureClassifier.classify_article()  # LLM 分类文献
     ↓  分流
┌────────────────────────────────────────────────────────────┐
│ 功能研究文献 → extract_functional_evidence() → PS3_Literature│
│ 病例/队列文献 → extract_case_report_evidence() → PS2/PM3/PP1│
└────────────────────────────────────────────────────────────┘
     ↓
EvidenceRouter.route_evidence() 基于遗传模式路由证据
     ↓
输出: VariantLiterature + evidence assessment
     ↓
ACMG 规则使用 evidence dict 评估
```

---

## 2. 安装与配置

### 2.1 安装依赖

```bash
pip install openai>=1.0.0,<2.0.0
```

### 2.2 配置 API Key

**方式一：环境变量**
```bash
export OPENAI_API_KEY="sk-..."
```

**方式二：代码中传入**
```python
from literature_classifier import LiteratureClassifier

classifier = LiteratureClassifier(llm_api_key="sk-...")
```

### 2.3 配置文件

在 YAML 配置文件的 `rules` 列表中添加所需规则：

```yaml
name: "example_config"
rules:
  # 文献证据规则
  - "ps3_combined"  # 综合 PS3 (SplcieVarDB + functional_evidence.db + 文献)
  - "ps2"           # de novo 变异
  - "pm3"           # 双等位基因变异
  - "bp2"           # 反式位置（良性）
  - "pp1"           # 共分离（优先使用文献，回退到 ClinVar）
  # ... 其他规则
```

---

## 3. 文献检索流程

### 3.1 两级检索策略

#### 一级检索（完整检索词）

系统自动生成多种检索词：

| 检索类型 | 示例 |
|----------|------|
| rs号 | `rs12345` |
| 基因 | `BRCA1` |
| 基因+c.HGVS | `BRCA1 c.68_69delAG` |
| 基因+p.HGVS | `BRCA1 p.Arg1699Ter` |
| VCF格式 | `chr17:43045678:G:A`, `17:43045678:G:A` |

#### 扩展检索（一级无结果时）

当一级检索无结果时，自动使用扩展检索词：

| 检索类型 | 示例 |
|----------|------|
| 基因 + c.核苷酸位置 | `"BRCA1 c.68"`, `"BRCA1 c.68_69del"` |
| 基因 + 密码子位置 | `"BRCA1 Gly1699"`, `"BRCA1 1699"` |
| 染色体位置 | `"chr17:41244938"`, `"17:41244938"` |

### 3.2 文献检索调用

**便捷函数（推荐）**

```python
from literature_trigger import quick_literature_search
from literature_retrieval.literature_utils import InheritancePattern

normalized, literature, evidence, used_expanded = quick_literature_search(
    query_type="rsid",
    input_string="rs123456",
    gene_symbol="BRCA1",
    inheritance_pattern=InheritancePattern.AUTOSOMAL_DOMINANT,
)
```

**分步调用**

```python
from normalizer import VariantNormalizer, VariantInfo
from literature_trigger import search_literature_for_variant, assess_literature_evidence
from literature_retrieval.literature_utils import InheritancePattern

# Step 1: 检索文献
normalizer = VariantNormalizer()
normalized = normalizer.normalize(query_type="rsid", input_string="rs123456")

literature, used_expanded = search_literature_for_variant(
    normalized,
    inheritance_pattern=InheritancePattern.AUTOSOMAL_DOMINANT,
)

# Step 2: 评估证据（需要 LLM API）
evidence = assess_literature_evidence(
    literature,
    inheritance_pattern=InheritancePattern.AUTOSOMAL_DOMINANT,
    gene="BRCA1",
    variant="p.Arg1699Ter",
    llm_api_key="sk-...",
)
```

---

## 4. LLM 证据评估

### 4.1 LLM 分类流程

对每篇文献，LLM 执行两步判断：

```
Step 1: classify_article() - 判断文献是否相关
        ↓
        输出: is_variant_related, literature_type, confidence
        ↓
Step 2a: 如果是病例/队列文献 → extract_case_report_evidence()
         → 输出: CaseReportEvidence (病例数、de_novo状态、遗传模式等)
         ↓
Step 2b: 如果是功能研究文献 → extract_functional_evidence()
         → 输出: FunctionalStudyEvidence (功能结果、技术方法等)
```

### 4.2 Prompt 模板

#### LITERATURE_CLASSIFICATION_PROMPT
用于判断文献是否与变异相关及文献类型。

**输入**: 基因、变异、标题、摘要
**输出**:
```json
{
    "is_variant_related": true/false,
    "literature_type": "case_report/cohort_study/functional_study/review/other",
    "confidence": 0.0-1.0,
    "reasoning": "..."
}
```

#### FUNCTIONAL_STUDY_EXTRACTION_PROMPT
用于提取功能研究证据（PS3）。

**输入**: 基因、变异、标题、摘要
**输出**:
```json
{
    "functional_result": "lof/normal/uncertain",
    "technique": "实验技术描述",
    "percentage_effect": 0.0-100.0,
    "confidence": "high/medium/low",
    "reasoning": "..."
}
```

#### CASE_REPORT_EXTRACTION_PROMPT
用于提取病例报道证据（PS2/PM3/PP1）。

**输入**: 基因、变异、遗传模式、标题、摘要
**输出**:
```json
{
    "case_count": 1,
    "is_de_novo": true/false,
    "de_novo_confirmed": true/false,
    "inheritance_pattern": "AD/AR/XLD/UNKNOWN",
    "is_compound_het": true/false,
    "trans_variant": "description or null",
    "segregation_data": {...} or null,
    "phenotype": "...",
    "hpo_terms": ["HP:xxxxx", ...],
    "odds_ratio": number or null,
    "p_value": number or null
}
```

---

## 5. 证据规则详解

### 5.1 PS2 - de novo 变异证据

#### 证据定义
PS2: 变异经父母验证为 de novo（双方父母均阴性）。

PM6: 假设为 de novo 但未经父母验证。

#### 证据强度

| 条件 | 强度 | 分数 |
|------|------|------|
| 1 个确认 de novo | Supporting (P) | 1 分 |
| 2 个确认 de novo | Moderate (M) | 2 分 |
| 3 个确认 de novo | Strong (S) | 4 分 |
| 4+ 个确认 de novo | Very Strong (VS) | 8 分 |

#### LLM 评估逻辑

```
输入: CaseReportEvidence 列表
     ↓
1. 检查每篇文献的 is_de_novo 和 de_novo_confirmed
2. 统计独立先证者数量
3. 检查表型一致性 (phenotype_consistency)
     - HIGH: 病例表型一致 → 维持/升级强度
     - LOW: 表型异质性高 → 降低强度或不适用
4. 输出: PS2Assessment
```

#### 使用方式

```python
# 配置文件添加规则
rules:
  - "ps2"
```

#### 代码评估

```python
variant.variant_literature["evidence"]["ps2"]
# 返回: EvidenceResult 或 None
```

---

### 5.2 PM3 - 双等位基因变异证据

#### 证据定义
PM3: 变异在反式位置（复合杂合）或纯合/半合子状态下被检出，符合隐性遗传模式。

#### 证据强度

| 配置类型 | 基础强度 |
|----------|----------|
| 反式 + 致病性变异 + 功能研究 | Very Strong |
| 反式 + 致病性变异 | Strong |
| 反式 + 可能致病变异 | Moderate |
| 纯合/半合子 | Supporting |

#### 基于先证者数量的升级

| 先证者数量 | 最终强度 |
|------------|----------|
| 1 | Supporting (P) |
| 2 | Moderate (M) |
| 3 | Strong (S) |
| 4+ | Very Strong (VS) |

**注意**: 最终强度取基础强度和先证者升级强度的较高者。

#### LLM 评估逻辑

```
输入: CaseReportEvidence 列表
     ↓
1. 检查每篇文献的 is_compound_het, is_homozygous, is_hemizygous
2. 确定配置类型 (trans/homozygous/hemizygous)
3. 查找已知致病变异 (known_pathogenic_variant)
4. 统计独立先证者数量
5. 输出: PM3Assessment
```

#### 使用方式

```python
# 配置文件添加规则
rules:
  - "pm3"
```

#### 代码评估

```python
variant.variant_literature["evidence"]["pm3"]
# 返回: EvidenceResult 或 None
```

---

### 5.3 PP1 - 共分离证据

#### 证据定义
PP1: 变异在家系中与疾病共分离（患病家系成员携带变异，健康家系成员不携带）。

#### 证据强度

| 条件 | 强度 |
|------|------|
| 基于 LOD score | Supporting/Moderate/Strong |
| <10 个信息减数分裂 | Supporting |
| ≥10 个信息减数分裂 | Moderate |
| ≥20 个信息减数分裂 | Strong |

#### 证据来源优先级

1. **优先**: 文献检索证据 (`variant.variant_literature["evidence"]["pp1"]`)
2. **回退**: ClinVar 多因素证据 (`variant.multifactorial_likelihood.co_segregation`)

#### 使用方式

```python
# 配置文件添加规则
rules:
  - "pp1"
```

#### 代码评估

```python
variant.variant_literature["evidence"]["pp1"]
# 如果文献无证据，自动回退到 ClinVar
```

---

### 5.4 BP2 - 反式位置（良性证据）

#### 证据定义
BP2: 在隐性遗传病基因中，候选良性变异与已知致病变异处于反式位置。

#### 证据强度

| 条件 | 强度 |
|------|------|
| 1 个先证者 | Supporting |
| 2+ 先证者 | Moderate |
| 多名家系 | Strong |

#### 使用方式

```python
# 配置文件添加规则
rules:
  - "bp2"
```

#### 代码评估

```python
variant.variant_literature["evidence"]["bp2"]
# 返回: EvidenceResult 或 None
```

---

### 5.5 PS3 - 功能研究证据

#### 证据来源（综合）

| 来源 | 证据强度 | 分数 |
|------|----------|------|
| SpliceVarDB | Supporting | 1 分 |
| functional_evidence.db (≥6 篇) | Moderate | 2 分 |
| functional_evidence.db (<6 篇) | Supporting | 1 分 |
| Literature (high confidence) | Moderate | 2 分 |
| Literature (medium/low) | Supporting | 1 分 |

#### PS3_combined 评估逻辑

```
输入: 三种来源的 PS3 结果
     ↓
1. SpliceVarDB 评估
2. functional_evidence.db 评估
3. Literature 评估 (如果可用)
     ↓
取最高分数的结果作为最终 PS3
     ↓
如果有冲突证据 (PS3 vs BS3)，不应用
```

#### 使用方式

```python
# 配置文件添加规则
rules:
  - "ps3_combined"
```

#### 代码评估

```python
# variant.variant_literature["evidence"] 中包含
variant.variant_literature["evidence"]["ps3_literature"]
# 返回: dict with applicable, strength, score, reason
```

---

## 6. 配置示例

### 6.1 完整配置示例

```yaml
# config_example.yaml
name: "BRCA1_example"
rules:
  # ========== 文献证据规则 ==========
  - "ps3_combined"    # 综合功能证据 (SplcieVarDB + functional_evidence.db + 文献)
  - "ps2"             # de novo 变异
  - "pm3"             # 双等位基因变异
  - "bp2"             # 反式位置（良性）
  - "pp1"             # 共分离

  # ========== 其他 ACMG 规则 ==========
  - "pvs1"
  - "ps1_protein"
  - "ps1_splicing"
  - "pm1"
  - "pm2"
  - "pm4"
  - "pm5_protein"
  - "pp2"
  - "pp3_protein"
  - "pp3_splicing"
  - "ba1"
  - "bs1"
  - "bs2"
  - "bs3"
  - "bs4"
  - "bp1"
  - "bp3"
  - "bp4_protein"
  - "bp7"
  - "bp6_clingen"

# 预测工具阈值
prediction_tool_threshold:
  pathogenicity_prediction:
    name: "pathogenicity_prediction"
    benign:
      direction: "less_than"
      supporting: 0.001
      moderate: 0.005
      strong: 0.01
      very_strong: 0.05
    pathogenic:
      direction: "greater_than"
      supporting: 0.01
      moderate: 0.05
      strong: 0.1
      very_strong: 0.5

# 注释文件路径
annotation_files:
  root: "/path/to/data"
  clinvar:
    root: "/path/to/data"
    clinvar_snv: "clinvar_snv.txt"
  splicevardb:
    root: "/path/to/data"
    splicevardb_file: "splicevardb_with_pmid.csv"
  functional_evidence_db:
    root: "/path/to/data"
    functional_evidence_db_file: "functional_evidence.db"
```

### 6.2 运行分类

```bash
# 设置 API Key
export OPENAI_API_KEY="sk-..."

# 运行分类
python classify.py -c config_example.yaml -i rs123456 -g BRCA1 --inheritance AD
```

### 6.3 API 调用示例

```python
from literature_trigger import quick_literature_search
from literature_retrieval.literature_utils import InheritancePattern

# 文献检索 + 评估
normalized, literature, evidence, used_expanded = quick_literature_search(
    query_type="rsid",
    input_string="rs123456",
    gene_symbol="BRCA1",
    inheritance_pattern=InheritancePattern.AUTOSOMAL_DOMINANT,
)

# 检查证据
print("PS2:", evidence.get("ps2"))
print("PM3:", evidence.get("pm3"))
print("PP1:", evidence.get("pp1"))
print("BP2:", evidence.get("bp2"))
print("PS3 literature:", evidence.get("ps3_literature"))
```

---

## 7. 故障排除

### 7.1 LLM API 未配置

**症状**: 文献分类返回默认值（is_variant_related=False）

**解决**:
1. 确认设置了 `OPENAI_API_KEY` 环境变量
2. 或在代码中传入 `llm_api_key` 参数
3. 检查网络连接

### 7.2 无文献检索结果

**症状**: literature.total_articles = 0

**解决**:
1. 检查扩展检索是否触发 (`used_expanded_search`)
2. 尝试提供更完整的变异描述（如同时提供基因符号）
3. 检查 PubMed 连接

### 7.3 规则未生效

**症状**: 分类结果中未见预期的证据

**解决**:
1. 确认规则名称在配置文件中正确拼写
2. 确认 `variant.variant_literature["evidence"]` 包含相应证据
3. 检查日志中的评估信息

### 7.4 证据强度不符合预期

**检查项**:
1. PS2: 确认 de_novo_confirmed = True
2. PM3: 确认 is_compound_het = True 或 is_homozygous = True
3. PP1: 确认 segregation_data 包含家系信息
4. 表型一致性是否影响强度

---

## 附录 A: EvidenceResult 数据结构

```python
@dataclass
class EvidenceResult:
    rule: str              # e.g., "PS2", "PP1"
    applicable: bool       # 证据是否适用
    strength: str          # "STRONG", "MODERATE", "SUPPORTING"
    num_cases: int         # 病例数
    num_controls: int     # 对照数（PS4）
    segregation_data: dict  # 共分离数据（PP1）
    phenotype: str        # 表型描述
    pmids: List[str]      # 相关 PMID 列表
    confidence: str        # "high", "medium", "low"
    comment: str          # 评估说明
```

---

## 附录 B: LLM 模型选择

当前实现使用 OpenAI GPT-4o，可通过修改 `literature_classifier.py` 中的 `_call_llm()` 方法切换：

```python
# 切换模型
response = client.chat.completions.create(
    model="gpt-4o",      # 或 "gpt-4-turbo", "gpt-3.5-turbo"
    messages=[...],
    temperature=0.3,
    max_tokens=1024,
    response_format={"type": "json_object"}
)
```

---

## 附录 C: 日志级别

系统使用 Python logging，日志名称前缀为 `GenOtoScope_Classify.literature_*`

```python
import logging

# 设置 DEBUG 级别查看详细信息
logging.getLogger("GenOtoScope_Classify.literature_classifier").setLevel(logging.DEBUG)
logging.getLogger("GenOtoScope_Classify.literature_trigger").setLevel(logging.DEBUG)
```
