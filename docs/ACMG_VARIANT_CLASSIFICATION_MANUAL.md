# ACMG 变异分类系统使用手册

## 版本信息
- **更新日期**: 2026-04-20
- **版本**: 1.3.0

**版本 1.3.0 更新**:
- PS1/PM5 规则改为 ClinVar API 优先查询 + 本地 VCF fallback 策略
- ClinVar 数据文件不再是分类的必要条件，改为可选的 fallback 数据源
- 新增 PS1_splicing 检索窗口详细说明（基于 ClinGen Splice Variant Classification 规则）

---

## 目录

1. [系统概述](#1-系统概述)
2. [安装与配置](#2-安装与配置)
3. [输入格式](#3-输入格式)
4. [ACMG 证据规则详解](#4-acmg-证据规则详解)
5. [文献检索与 LLM 证据](#5-文献检索与-llm-证据)
   - [5.4 RAG + LLM 证据校正](#54-rag--llm-证据校正)
6. [配置指南](#6-配置指南)
7. [运行方式](#7-运行方式)
8. [故障排除](#8-故障排除)

---

## 1. 系统概述

### 1.1 系统功能

ACMG (American College of Medical Genetics and Genomics) 变异分类系统用于对遗传变异进行致病性分类。系统支持：

- **多格式输入**: JSON、rsID、VCF、HGVS
- **自动化评估**: 基于 ACMG/AMP 指南的证据规则
- **文献检索**: 自动 PubMed 文献检索与证据提取
- **LLM 辅助**: 使用 LLM 分类文献和评估证据
- **基因特异性规则**: 支持 ClinGen 基因特异性规则
- **RAG+LLM 校正**: 基于 ClinGen 规则的证据强度自动调整

### 1.2 分类级别

| 分类 | 说明 | Bayes 分数 |
|------|------|------------|
| **Pathogenic** | 致病 | ≥ 10 |
| **Likely Pathogenic** | 可能致病 | 6 - 9 |
| **Uncertain Significance (VUS)** | 意义不明 | 0 - 5 |
| **Likely Benign** | 可能良性 | -6 - -1 |
| **Benign** | 良性 | ≤ -7 |

### 1.3 证据强度分数

| 强度 | 致病证据分数 | 良性证据分数 |
|------|-------------|-------------|
| Stand Alone (SA) | 10 分 | -10 分 |
| Very Strong (VS) | 8 分 | -8 分 |
| Strong (S) | 4 分 | -4 分 |
| Moderate (M) | 2 分 | -2 分 |
| Supporting (P) | 1 分 | -1 分 |

---

## 2. 安装与配置

### 2.1 系统要求

- Python 3.8+
- 依赖包: 见 `requirements.txt`

### 2.2 安装

```bash
pip install -r requirements.txt
```

### 2.3 LLM API 配置（用于文献检索和 RAG+LLM 校正）

系统支持统一的 LLM API 配置，支持 OpenAI 和 DeepSeek 自动选择：

```bash
# 安装依赖
pip install openai>=1.0.0,<2.0.0

# OpenAI 配置（优先）
export OPENAI_API_KEY="sk-..."

# DeepSeek 配置（备用）
export DEEPSEEK_API_KEY="sk-..."

# 通用配置（备用）
export LLM_API_KEY="sk-..."
```

**API Key 优先级**:
1. `OPENAI_API_KEY` → OpenAI GPT-4o
2. `DEEPSEEK_API_KEY` → DeepSeek deepseek-chat
3. `LLM_API_KEY` → OpenAI GPT-4o（默认）

**统一配置模块**: `llm_config.py` 提供跨模块一致性配置

### 2.4 数据文件要求

**查询策略**: ClinVar 采用 **API优先 + 本地VCF fallback** 的混合策略。

| 数据类型 | 文件 | 必要性 | 说明 |
|----------|------|--------|------|
| ClinVar | clinvar_snv.txt | **Fallback备选** | 本地VCF文件，当ClinVar API查询失败时使用 |
| ClinVar API | 在线查询 | **主要查询方式** | 通过NCBI ClinVar E-utilities API实时查询 |
| SpliceVarDB | splicevardb_with_pmid.csv | 推荐 | 剪接变异功能证据（PS3来源之一） |
| functional_evidence.db | functional_evidence.db | 推荐 | MAVE 功能实验数据（PS3来源之一） |
| ClinGen Rules | clingen_rules.db | 推荐 | ClinGen 基因特异性规则（RAG+LLM 校正） |
| 参考数据 | uniprot, critical_regions 等 | 可选 | 详见配置 |

#### ClinVar 查询流程

```
变异输入
    ↓
优先: ClinVar E-utilities API 在线查询
    ↓ (API失败或无结果时)
Fallback: 本地 clinvar_snv.txt VCF文件查询
    ↓
返回 PS1/PM5 证据
```

**注意**: ClinVar 数据文件不再是必要条件，但建议配置以提高查询可靠性。

---

## 3. 输入格式

### 3.1 JSON 格式

```json
{
  "variant_info": {
    "gene_name": "BRCA1",
    "hgvs_protein": "p.Arg1699Ter",
    "hgvs_cdna": "c.5096C>G",
    "rs_id": "rs80357362",
    "var_type": "missense",
    "chr": "17",
    "genomic_start": 43045678,
    "var_ref": "G",
    "var_obs": "A"
  },
  "transcript_info": [{
    "transcript_id": "NM_007294.3",
    "gene_symbol": "BRCA1"
  }]
}
```

### 3.2 rsID 格式

```
rs123456
```

### 3.3 VCF 格式

```
17:43045678:G:A
chr17:43045678:G:A
```

### 3.4 HGVS 格式

```
NM_007294.3:c.68_69delAG
```

---

## 4. ACMG 证据规则详解

### 4.1 致病性证据 (Pathogenic)

#### PVS1 - 极高变异 (Very Strong Pathogenic)

**定义**: 变异为已知机制下的功能丧失 (LoF)，且该基因无良性变异记录。

**适用条件**:
- 无义变异
- 移码变异 (frameshift)
- 剪接位点变异 (±1,2)
- 起始密码子变异
- 单外显子缺失

**PVS1 与 PS3_SpliceVarDB 冲突处理**:
- 对于经典剪接变异 (±1, ±2)，如果 SpliceVarDB 功能证据显示剪接影响，**返回 PVS1 而非 PS3**
- 这是因为经典剪接位点变异已有明确的 PVS1 证据，不应同时叠加 PS3_Supporting
- SpliceVarDB 对经典剪接变异的评估自动升级为 PVS1 (Very Strong)

**基因特异性规则**:
| 基因 | 规则类 | 说明 |
|------|--------|------|
| BRCA1 | Pvs1_brca1 | 5' 端和剪接位点 LoF |
| BRCA2 | Pvs1_brca2 | 5' 端和剪接位点 LoF |
| ATM | Pvs1_atm | 5' 端 LoF |
| PALB2 | Pvs1_palb2 | 5' 端和剪接位点 LoF |
| PTEN | Pvs1_pten | 5' 端 LoF |
| CDH1 | Pvs1_cdh1 | 5' 端 LoF |

#### PS1 - 相同氨基酸改变的已知致病变异

**定义**: 变异导致相同氨基酸改变 (如 c.76A>T → p.K26X)，且该变异已知为致病性。

**ClinVar 查询**: 通过 ClinVar E-utilities API 在线查询优先，本地 VCF 文件作为 fallback。

**评估要点**:
- 相同氨基酸位置的其他变异已知为致病
- 考虑HGVS表达式的等价性

#### PS1_splicing - 剪接变异位点的已知致病变异

**定义**: 剪接变异（±1~±2经典位点或±3~±20非经典位点）与已知致病变异位于同一核苷酸位置或同一剪接基序。

**ClinVar 查询**: 通过 ClinVar E-utilities API 在线查询优先，本地 VCF 文件作为 fallback。

##### 剪接区域分类

| 区域类型 | HGVS示例 | 说明 |
|----------|----------|------|
| CLASSIC_DONOR | c.321+1G>A, c.321+2G>C | 经典donor位点 (+1, +2) |
| CLASSIC_ACCEPTOR | c.123-1G>T, c.123-2C>G | 经典acceptor位点 (-1, -2) |
| NON_CLASSIC_DONOR | c.321+3A>G, c.321+6T>G | 非经典donor区域 (+3~+6) |
| NON_CLASSIC_ACCEPTOR | c.123-5G>A, c.123-20C>T | 非经典acceptor区域 (-3~-20) |
| NOT_SPLICING | p.Gly12Val (VEP=splice_site) | Branch 2: 错义+VEP标注splice |

##### ClinVar 检索窗口

**基于 HGVS 格式自动计算搜索窗口**:

| 变异类型 | 搜索窗口 (cDNA相对位置) | 基因组搜索范围 |
|----------|------------------------|----------------|
| Donor (+1~+6) | base-1 到 base+6 | c.{base-1} ~ c.{base+6} |
| Acceptor (-1~-20) | base-20 到 base | c.{base-20} ~ c.{base} |
| Branch 2 (missense+splice) | base-1 到 base+6 | c.{base-1} ~ c.{base+6} |

**示例说明**:

```
c.321+1G>A (经典donor):
  - base_position = 321
  - search_window: c.320 ~ c.327 (base-1 到 base+6)
  - ClinVar检索: 321位点的所有+/-1~6范围内的变异

c.123-1G>T (经典acceptor):
  - base_position = 123
  - search_window: c.103 ~ c.123 (base-20 到 base)
  - ClinVar检索: 123位点的所有/-1~-20范围内的变异
```

##### PS1_splicing 强度判定

| 条件 | ClinVar证据 | PS1类型 | 证据强度 |
|------|-------------|---------|----------|
| 同一核苷酸位置 | Pathogenic | PS1 | STRONG |
| 同一核苷酸位置 | Likely Pathogenic | PS1 | MODERATE |
| 同一剪接基序，不同位点 | Pathogenic | PS1_Moderate | MODERATE |
| 同一剪接基序，不同位点 | Likely Pathogenic | PS1_Supporting | SUPPORTING |

##### 分支判断逻辑 (Branch 1 vs Branch 2)

```
输入: HGVS字符串 + VEP consequence
    ↓
Branch 1: HGVS格式为剪接变异 (c.数字+/-数字碱基>碱基)
    ↓ (是)
使用splicing_info中的窗口进行ClinVar检索
    ↓
Branch 2: HGVS为错义变异 (p.XXX) 但VEP标注为splice_site
    ↓ (是，但不属于Branch 1)
搜索窗口: c.{base-1} ~ c.{base+6}
```

#### PS2 - de novo 变异 (Very Strong)

**定义**: 变异经父母验证为 de novo（双方父母均阴性），且先证者表型符合疾病特征。

**证据强度**:

| 先证者数量 | 强度 |
|------------|------|
| 1 个确认 de novo | Supporting |
| 2 个确认 de novo | Moderate |
| 3 个确认 de novo | Strong |
| 4+ 个确认 de novo | Very Strong |

**PM6**: 假设 de novo（未经父母验证）

**文献检索来源**: `variant.variant_literature["evidence"]["ps2"]`

#### PS3 - 功能研究证据 (Strong)

**定义**: 功能研究证明变异对基因/蛋白功能有有害影响。

**证据来源及强度**:

| 来源 | 条件 | 强度 |
|------|------|------|
| SpliceVarDB | 非经典剪接变异 | Supporting |
| SpliceVarDB | 经典剪接变异 (±1, ±2) | **PVS1 (Very Strong)** |
| functional_evidence.db | ≥6 篇文献 | Moderate |
| functional_evidence.db | <6 篇文献 | Supporting |
| Literature (功能研究) | high confidence | Moderate |
| Literature (功能研究) | medium/low confidence | Supporting |

**综合评估 (ps3_combined)**:
```
取 SpliceVarDB、functional_evidence.db、Literature 三者中的最高分
如存在冲突证据 (PS3 vs BS3)，不应用
```

**PVS1 与 PS3_Supporting 冲突处理**:
- 对于经典剪接变异 (±1, ±2)，如果 SpliceVarDB 显示剪接影响，**使用 PVS1 而非 PS3**
- 原因: 经典剪接位点变异已有明确的 PVS1 证据，不应同时叠加 PS3_Supporting
- SpliceVarDB 对经典剪接变异的评估会返回 PVS1 (Very Strong)

#### PS4 - 病例对照研究 (VS/S/M/P)

**定义**: 变异在病例组中的频率显著高于对照组，或多个无关病例携带相同变异且表型一致。

**注意**: PS4 不适用于 de novo 变异（应使用 PS2）。所有病例数均指**不相关的个体（非同一家族）**。

**证据强度** (基于先证者数量和疾病罕见程度):

**1. 标准显性遗传疾病**（如家族性高胆固醇血症、马凡综合征、遗传性耳聋）:
| 先证者数（不相关） | 强度 | 分值 |
|-------------------|------|------|
| ≥15 | Very Strong | 8分 |
| ≥10 | Strong | 4分 |
| ≥6 | Moderate | 2分 |
| ≥2 | Supporting | 1分 |

**2. 非常罕见的发育迟缓疾病（显性遗传）**:
| 先证者数（不相关） | 强度 | 分值 |
|-------------------|------|------|
| ≥4 | Strong | 4分 |
| ≥2 | Moderate | 2分 |
| ≥1 | Supporting | 1分 |

**3. 比家族性高胆固醇血症/马凡/遗传性耳聋更罕见，但比 Dravet 综合征更常见（显性遗传）**:
| 先证者数（不相关） | 强度 | 分值 |
|-------------------|------|------|
| ≥7 | Strong | 4分 |
| ≥4 | Moderate | 2分 |
| ≥2 | Supporting | 1分 |

#### PM1 - 热点区域/关键功能域

**定义**: 变异位于已知致病的热点区域或蛋白质关键功能域。

**ClinGen 热点区域数据库**:
- PM1 表格: ClinGen_PM1_Table_v2.xlsx

#### PM2 - 罕见变异 (Absent from Controls)

**定义**: 变异在对照人群中未发现或极其罕见 (≤0.001)。

**强度变体**:
| 条件 | 强度 |
|------|------|
| 完全不存在于对照人群 | PM2 |
| 极其罕见 (<0.0001) | PM2_supporting |
| 罕见 (<0.001) | PM2_supporting_less |
| 仅见于纯合子 | PM2_supporting_no_indel |

#### PM3 - 双等位基因变异 (Strong)

**定义**: 变异在隐性遗传病基因中处于反式位置（复合杂合）或纯合状态。

**配置类型**:

| 配置 | 基础强度 |
|------|----------|
| 反式 + 致病性变异 + 功能研究 | Very Strong |
| 反式 + 致病性变异 | Strong |
| 反式 + 可能致病变异 | Moderate |
| 纯合/半合子 | Supporting |

**先证者数量升级**:

| 数量 | 最终强度 |
|------|----------|
| 1 | Supporting |
| 2 | Moderate |
| 3 | Strong |
| 4+ | Very Strong |

**文献检索来源**: `variant.variant_literature["evidence"]["pm3"]`

#### PM4 - 蛋白质长度改变 (Moderate)

**定义**: 变异导致蛋白质长度改变（非终止密码子变异）。

**适用条件**:
- 剪接变异导致外显子跳跃
- 框内缺失/插入
- 终止密码子丢失

#### PM5 - 相同位置不同氨基酸变异 (Moderate)

**定义**: 变异位于与已知致病变异相同的氨基酸位置，但改变为不同氨基酸。

**ClinVar 查询**: 通过 ClinVar E-utilities API 在线查询优先，本地 VCF 文件作为 fallback。

**示例**: p.Arg1699Gln (已知致病) → p.Arg1699Ter (待评估)

**强度变体**:
| 条件 | 规则 |
|------|------|
| 相同位置已知致病 | PM5_protein |
| 相同位置已知可能致病 | PM5_protein_pathogenic |
| 终止密码子变异 (PTC) | PM5_protein_ptc |
| 剪接变异 | PM5_splicing_ptc |

#### PP1 - 共分离证据 (S/M/P)

**定义**: 变异在家系中与疾病共分离。

**注意**: PS4 不适用于 de novo 变异（应使用 PS2）。所有病例数均指**不相关的个体（非同一家族）**。

**证据强度** (基于共分离次数):

**1. 显性或X连锁遗传**:
| 共分离次数 | 强度 | 分值 |
|-----------|------|------|
| 2-3 | Supporting (PP1) | 1分 |
| 4 | Moderate (PP1_Moderate) | 2分 |
| 5+ | Strong (PP1_Strong) | 4分 |

注: "共分离" = 一次减数分裂中变异与疾病共分离
对于显性遗传：除先证者外，每个受累携带者亲属 = 1次共分离

**2. 隐性遗传**:
| 证据 | 强度 |
|------|------|
| 1个家系，2个患者均携带变异 = 1次共分离 | Supporting (PP1) |
| 携带者父母，3个受累携带者 + 1个正常不携带者 | Moderate (PP1_Moderate) |
| 携带者父母，2个受累携带者 + 2个正常不携带者 | Moderate (PP1_Moderate) |
| 2个受累携带者 + 3个正常不携带者 | Strong (PP1_Strong) |

**证据来源优先级**:
1. 优先: 文献检索 `variant.variant_literature["evidence"]["pp1"]`
2. 回退: ClinVar 多因素证据

#### PP2 - 基因特异性良性阈值 (Supporting)

**定义**: 变异在罕见良性变异稀少的基因中，且为罕见变异。

**评估依据**: gnomAD missense_zscore（错义变异容忍度评分）
- z-score > 2.09 表示基因对错义变异不耐受，错义变异是常见致病机制
- z-score ≤ 2.09 表示基因对错义变异较容忍，PP2 不适用

**数据来源**: gnomAD Gene Constraint API（自动调用）

**评估流程**:
```
1. 检查是否为错义变异
   ↓ 不是 → PP2 不适用
   ↓ 是
2. 通过 gnomAD API 获取基因的 missense_zscore
3. 比较 z-score 与阈值 2.09
   ↓ z-score > 2.09 → PP2 适用 (Supporting)
   ↓ z-score ≤ 2.09 → PP2 不适用
```

#### PP3 - 预测有害 (S/M/P)

**预测方法**:
- **剪接预测**: PP3_splicing, PP3_splicing_cdh1
- **蛋白预测**: PP3_protein

**当前使用的工具**:

| 工具类型 | 工具名称 | PP3 (Pathogenic) | BP4 (Benign) | 数据来源 |
|----------|----------|-----------------|--------------|----------|
| **Missense** | REVEL | ≥0.644 (S), ≥0.773 (M), ≥0.932 (Str) | ≤0.391 | VEP MyVariant (dbNSFP) |
| **Splicing** | SpliceAI | ≥0.2 | ≤0.1 | SpliceVarDB |

**强度等级说明**:
- **S** = Supporting (1分)
- **M** = Moderate (2分)
- **Str** = Strong (4分)

**PP3 + PS3 合计最高 4 分** (避免重复计算)

**评估流程**:
```
获取 variant.prediction_tools dict
       ↓
检查变异类型 (Missense vs Splicing/Intronic)
       ↓
Missense: REVEL 有值按阈值判断
Splicing: SpliceAI 有值按阈值判断
       ↓
根据阈值触发 PP3/BP4 及对应强度
```

#### PP4 - 表型特异性 (Supporting)

**定义**: 患者表型与单一病因疾病高度特异匹配。

**条件**:
- 基因已知导致该特定表型
- 无其他解释

#### PP5 - 信誉良好的来源报告为致病变异

**定义**: 可靠来源报告为致病变异，但无独立验证。

**注意**: 已被 ClinVar 注释替代

---

### 4.2 良性证据 (Benign)

#### BA1 - 等位基因频率过高 (Stand Alone)

**定义**: 变异在对照人群中常见 (>5%)。

**阈值**: 等位基因频率 > 5%

#### BS1 - 频率高于疾病预期 (Strong)

**定义**: 变异频率高于疾病在人群中的预期。

**示例**: 疾病发病率 1:1000，变异频率 1%

#### BS2 - 早发型疾病中的成年对照 (Strong)

**定义**: 在严重早发型疾病中，健康成年人群存在该变异。

**条件**:
- 疾病为严重早发
- 对照为健康成年人
- 纯合子观察

#### BS3 - 功能研究证明无有害影响 (Strong)

**定义**: 功能研究证明变异对蛋白功能无影响。

**来源**: 同 PS3，但证据指向良性

#### BS4 - 不共分离 (Strong)

**定义**: 变异在家系中与疾病不共分离。

**条件**:
- 多个患病家系成员不携带变异
- 或健康家系成员携带变异

#### BP1 - 主要为 LoF 基因中的错义变异 (Supporting)

**定义**: 主要功能为 LoF 的基因（如截短变异）中出现的错义变异。

#### BP2 - 反式位置的致病变异 (Supporting)

**定义**: 候选良性变异与已知致病变异在反式位置。

**文献检索来源**: `variant.variant_literature["evidence"]["bp2"]`

#### BP3 - 剪接预测但无实验验证 (Supporting)

**定义**: 预测影响剪接但缺乏功能验证。

#### BP4 - 多个预测工具一致认为无害

**定义**: 多种计算预测方法均预测为无害。

**当前使用的工具**:

| 工具类型 | 工具名称 | BP4 (Benign) | PP3 (Pathogenic) | 数据来源 |
|----------|----------|--------------|------------------|----------|
| **Missense** | REVEL | ≤0.391 | ≥0.644 | VEP MyVariant (dbNSFP) |
| **Splicing** | SpliceAI | ≤0.1 | ≥0.2 | SpliceVarDB |

**注意**: BP4 和 PP3 互为反向证据，使用相同的工具但阈值方向相反。

#### BP5 - 多种疾病中的变异 (Supporting)

**定义**: 变异在多种不相关疾病中观察到。

#### BP6 - 可靠来源报告为良性

**定义**: 可靠来源报告为良性/可能良性。

#### BP7 - 同义变异不影响剪接

**定义**: 同义变异不影响剪接共识序列。

---

### 4.3 规则组合与互斥

#### 互斥规则

| 规则 | 互斥 |
|------|------|
| PVS1 + PVS1_general | 互斥 |
| BP4 + BP4_mult_strength | 互斥 |

#### 需同时应用的规则

| 规则 | 需配合 |
|------|--------|
| PS2 | 需确认表型符合 |
| PM3 | 需确认隐性遗传 |
| PP1 | 需家系数据 |

---

## 5. 文献检索与 LLM 证据

### 5.1 文献检索流程

```
输入: VariantInfo
     ↓
build_search_queries() 生成检索词
     ↓
NCBILiteratureRetriever.search() 检索 PubMed
     ↓
LiteratureClassifier.classify_article()  # LLM 分类
     ↓
┌──────────────────────────────────────────────┐
│ 功能研究 → extract_functional_evidence() → PS3 │
│ 病例报道 → extract_case_report_evidence() → PS2/PM3/PP1 │
└──────────────────────────────────────────────┘
     ↓
EvidenceRouter.route_evidence()
     ↓
variant.variant_literature["evidence"]
```

### 5.2 两级检索策略

#### 一级检索词

| 类型 | 示例 |
|------|------|
| rsID | rs12345 |
| 基因 | BRCA1 |
| 基因+c.HGVS | BRCA1 c.68_69delAG |
| 基因+p.HGVS | BRCA1 p.Arg1699Ter |
| VCF | chr17:43045678:G:A |

#### 扩展检索（一级无结果时）

| 类型 | 示例 |
|------|------|
| 基因 + c.位置 | "BRCA1 c.68" |
| 基因 + 密码子 | "BRCA1 Gly1699" |

### 5.3 LLM Prompt 模板

#### 文献分类

```
判断文献是否与变异相关，以及文献类型：
- case_report: 病例报道
- cohort_study: 队列研究
- functional_study: 功能研究
- review: 综述
```

#### 功能研究证据提取

```
判断功能研究是否支持 Loss of Function (LOF)：
- lof: 功能丧失 → PS3
- normal: 功能正常 → BS3
- uncertain: 不确定 → 不使用
```

#### 病例报道证据提取

```
提取病例信息：
- is_de_novo: 是否 de novo
- de_novo_confirmed: 是否经父母验证
- is_compound_het: 是否复合杂合
- inheritance_pattern: 遗传模式
- segregation_data: 共分离数据
```

### 5.4 RAG + LLM 证据校正

RAG+LLM 证据校正模块使用 ClinGen 基因特异性规则对 ACMG 证据进行自动调整。

#### 工作流程

```
分类结果规则字典
     ↓
ClinGen RAG 检索（基因特异性规则）
     ↓
LLM 评估（基于 ClinGen 规则调整建议）
     ↓
调整后规则字典
     ↓
最终 Bayes 分类
```

#### 支持的调整操作

| 操作 | 说明 | 示例 |
|------|------|------|
| **remove** | 移除不适用的证据 | ClinGen 规则说明该基因不使用 PM5 |
| **upgrade** | 升级证据强度 | ClinGen 规则支持升级为 Strong |
| **downgrade** | 降级证据强度 | ClinGen 规则建议降级为 Supporting |
| **keep** | 保持原证据 | ClinGen 规则支持原评估 |

#### 启用方式

```bash
# 命令行启用
python classify.py -c config.yaml -i rs123456 -g BRCA1 --enable-rag-llm

# Python API
from classify import classify
result = classify(
    config_path="config.yaml",
    variant_str="rs123456",
    enable_rag_llm=True,
)
```

#### ClinGen 规则数据库

RAG 检索使用 `clingen_rules.db` 数据库，包含：
- ClinGen 基因特异性规则
- PM1 热点域定义
- 频率阈值

数据库默认路径: `./clingen_rules.db`

---

## 6. 配置指南

### 6.1 配置文件结构

```yaml
name: "config_name"
rules:
  - "rule1"
  - "rule2"
  # ...

prediction_tool_threshold:
  pathogenicity_prediction:
    name: "..."
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

annotation_files:
  root: "/path/to/data"
  clinvar:
    root: "/path/to/data"
    clinvar_snv: "clinvar_snv.txt"
  # ...
```

### 6.2 常用配置模板

#### BRCA1/BRCA2

```yaml
name: "brca_config"
rules:
  # 致病性证据
  - "pvs1_brca1"          # 或 pvs1_brca2
  - "ps1_protein"
  - "ps1_splicing"
  - "ps3_combined"         # 功能证据（综合）
  - "ps2"                  # de novo（如适用）
  - "pm1"
  - "pm2"
  - "pm4"
  - "pm5_protein"
  - "pm5_protein_ptc"
  - "pp1"                  # 共分离
  - "pp3_protein"
  - "pp3_splicing"
  - "pp4_enigma"
  # 良性证据
  - "ba1"
  - "bs1"
  - "bs1_faf"
  - "bs2"
  - "bs3"
  - "bp1"
  - "bp3"
  - "bp4_protein"
  - "bp7"
  - "bp6_clingen"
```

#### ATM

```yaml
name: "atm_config"
rules:
  - "pvs1_atm"
  - "ps1_protein"
  - "ps3_combined"
  - "ps2"
  - "pm1"
  - "pm2"
  - "pm4"
  - "pp1"
  - "pp3_protein"
  - "pp4_enigma"
  - "ba1"
  - "bs1"
  - "bs2"
  - "bp1"
  - "bp3"
  - "bp4_protein"
  - "bp7"
  - "bp7_deep_intronic_atm"
```

#### PTEN

```yaml
name: "pten_config"
rules:
  - "pvs1_pten"
  - "ps1_protein"
  - "ps1_splicing_pten"
  - "ps3_combined"
  - "pm1_tp53"
  - "pm1"
  - "pm2"
  - "pm4_pten"
  - "pm5_protein"
  - "pm5_protein_ptc"
  - "pp1"
  - "pp3_protein"
  - "pp4_enigma"
  - "ba1"
  - "bs1"
  - "bs2"
  - "bp1"
  - "bp3"
  - "bp4_protein"
  - "bp7"
```

### 6.3 文献证据规则配置

```yaml
rules:
  # ========== 文献证据规则 ==========
  - "ps3_combined"    # 综合功能证据
  - "ps2"             # de novo
  - "pm3"             # 双等位基因
  - "bp2"             # 反式位置（良性）
  - "pp1"             # 共分离（优先文献）
```

---

## 7. 运行方式

### 7.1 命令行

```bash
# JSON 输入
python classify.py -c config.yaml -i '{"variant_info": {...}}' -o result.json

# rsID 输入
python classify.py -c config.yaml -i rs123456 -g BRCA1 --inheritance AD

# VCF 输入
python classify.py -c config.yaml -i "17:43045678:G:A" -g BRCA1

# HGVS 输入
python classify.py -c config.yaml -i "NM_007294.3:c.68_69delAG" -g BRCA1

# 启用 RAG+LLM 证据校正
python classify.py -c config.yaml -i rs123456 -g BRCA1 --enable-rag-llm
```

### 7.2 Web API

```bash
# 启动服务
python webservice.py

# 分类请求
curl -X POST http://localhost:8000/classify \
  -H "Content-Type: application/json" \
  -d '{"text": "BRCA1 p.Arg1699Ter", "config_path": "/path/to/config.yaml"}'
```

### 7.3 Python API

```python
from classify import classify
import pathlib

# 分类
final_config, result = classify(
    config_path=pathlib.Path("config.yaml"),
    variant_str="rs123456",
    query_type="rsid",
    gene_symbol="BRCA1",
    inheritance_pattern="AD",
)

print(result)
```

### 7.4 文献检索单独使用

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

# 查看证据
print("PS2:", evidence.get("ps2"))
print("PM3:", evidence.get("pm3"))
print("PP1:", evidence.get("pp1"))
print("PS3 literature:", evidence.get("ps3_literature"))
```

---

## 8. 故障排除

### 8.1 常见错误

| 错误 | 原因 | 解决方案 |
|------|------|----------|
| Config file not found | 配置文件路径错误 | 检查路径 |
| No rule matching | 规则名称拼写错误 | 检查 config_annotation.py |
| LLM API not configured | 未设置 API Key | 设置 OPENAI_API_KEY/DEEPSEEK_API_KEY/LLM_API_KEY |
| No literature results | 检索词无结果 | 使用扩展检索 |
| ClinGen DB not found | clingen_rules.db 不存在 | 检查数据库路径 |
| RAG+LLM adjustment failed | LLM 调用失败 | 检查 API Key 配置 |

### 8.2 证据未应用

| 症状 | 检查项 |
|------|--------|
| PS2 未应用 | 确认 de_novo_confirmed=True |
| PM3 未应用 | 确认 is_compound_het=True |
| PP1 未应用 | 检查文献或 ClinVar 共分离数据 |
| PS3 未应用 | 检查 SpliceVarDB 或文献功能研究 |

### 8.3 分类结果异常

| 症状 | 检查项 |
|------|--------|
| 分数过低 | 检查是否缺少证据规则 |
| 分数过高 | 检查是否有冲突证据 |
| VUS 误判 | 检查规则权重配置 |

---

## 附录 A: 规则速查表

### 致病性证据

| 规则 | 强度 | 条件 |
|------|------|------|
| PVS1 | VS | LoF 变异 |
| PS1 | VS | 相同改变的已知致病 |
| PS2 | VS/S/M/P | de novo 变异 |
| PS3 | S/M/P | 功能研究有害 |
| PS4 | VS/S/M/P | 病例对照研究 |
| PM1 | M | 热点区域 |
| PM2 | P | 罕见变异 |
| PM3 | VS/S/M/P | 双等位基因 |
| PM4 | M | 蛋白长度改变 |
| PM5 | M | 相同位置不同氨基酸 |
| PP1 | S/M/P | 共分离 |
| PP2 | P | 基因特异性良性阈值 |
| PP3 | S/M/P | 预测有害 |
| PP4 | P | 表型匹配 |
| PP5 | P | 来源报告致病（已废弃） |

### 良性证据

| 规则 | 强度 | 条件 |
|------|------|------|
| BA1 | Stand Alone | 频率 >5% |
| BS1 | S | 频率高于预期 |
| BS2 | S | 健康对照中存在 |
| BS3 | S | 功能研究正常 |
| BS4 | S | 不共分离 |
| BP1 | P | LoF 基因错义变异 |
| BP2 | P | 反式致病变异 |
| BP3 | P | 剪接预测无验证 |
| BP4 | P | 预测无害 |
| BP5 | P | 多疾病中观察到 |
| BP6 | P | 来源报告良性 |
| BP7 | P | 同义不响剪接 |

---

## 附录 B: 基因特异性规则

| 基因 | 专用规则 |
|------|----------|
| BRCA1 | pvs1_brca1, ps1_splicing, pp3_splicing_cdh1 |
| BRCA2 | pvs1_brca2, pp3_splicing_cdh1 |
| ATM | pvs1_atm, bp7_deep_intronic_atm |
| PTEN | pvs1_pten, ps1_splicing_pten, pm4_pten |
| PALB2 | pvs1_palb2, bp7_deep_intronic_palb2 |
| CDH1 | pvs1_cdh1, pm5_protein_cdh1, pm5_splicing_cdh1 |
| TP53 | pm1_tp53 |
| SCN1A/2A/3A/8A | ClinGen PM1 表格 |

---

## 附录 C: 日志配置

```python
import logging

# 详细日志
logging.getLogger("Classify").setLevel(logging.DEBUG)

# 文献检索日志
logging.getLogger("GenOtoScope_Classify.literature_classifier").setLevel(logging.DEBUG)
logging.getLogger("GenOtoScope_Classify.literature_trigger").setLevel(logging.DEBUG)
```
