# ACMG Variant Classification System

基于 ACMG（American College of Medical Genetics and Genomics）指南的变异分类系统，支持贝叶斯评分和多规则综合评估。

**版本**: 1.2.0

## 主要特性

- **多格式输入**: JSON、rsID、VCF、HGVS
- **自动化评估**: 基于 ACMG/AMP 指南的证据规则
- **文献检索**: 自动 PubMed 文献检索与证据提取
- **LLM 辅助**: 使用 LLM 分类文献和评估证据
- **ClinGen 基因特异性规则**: 支持 ClinGen 基因特异性规则
- **RAG+LLM 校正**: 基于 ClinGen 规则的证据强度自动调整

## 项目结构

```
variant_classification/
├── .vscode/                  # VSCode 配置
│   ├── launch.json          # 调试配置
│   └── settings.json         # IDE 设置
├── acmg_rules/               # ACMG 规则实现
│   ├── ba1.py, bs1.py, pm2.py  # 良性/致病证据规则
│   └── pvs1.py, ps1.py, pp1.py # 致病证据规则
├── classification_schemata/  # 分类方案
├── docs/                     # 使用手册
├── literature_retrieval/     # 文献检索模块
├── webservice.py             # FastAPI 服务
├── classify.py               # 主分类入口
├── llm_config.py             # 统一 LLM 配置
├── rag_llm_evidence_adjuster.py  # RAG+LLM 证据校正
├── functional_evidence.db   # 功能证据数据库
├── clingen_rules.db         # ClinGen 规则数据库
├── config_example.yaml       # 配置示例
├── test_classify.py         # 测试脚本
└── requirements.txt         # Python 依赖
```

## 快速开始

### 1. 安装依赖

```bash
pip install -r requirements.txt
```

### 2. 配置环境变量（可选，用于 LLM 功能）

```bash
# OpenAI（优先）
export OPENAI_API_KEY="sk-..."

# DeepSeek（备用）
export DEEPSEEK_API_KEY="sk-..."
```

### 3. 在 VSCode 中运行

1. 用 VSCode 打开项目目录
2. 按 `F5` 选择调试配置
3. 选择 "Python: Test Classification" 运行测试

### 4. 命令行运行

```bash
# JSON 输入
python classify.py -c config_example.yaml -i '{"variant_info": {...}}'

# rsID 输入
python classify.py -c config_example.yaml -i rs123456 -g BRCA1 --inheritance AD

# VCF 输入
python classify.py -c config_example.yaml -i "17:43045678:G:A" -g BRCA1

# HGVS 输入
python classify.py -c config_example.yaml -i "NM_007294.3:c.68_69delAG" -g BRCA1

# 启用 RAG+LLM 证据校正
python classify.py -c config_example.yaml -i "17:43045678:G:A" -g BRCA1 --enable-rag-llm
```

### 5. FastAPI 服务

```bash
# 启动服务
python webservice.py --port 8080

# API 调用
curl -X POST http://localhost:8080/classify_variant \
  -H "Content-Type: application/json" \
  -d '{"variant_json": "17:43045678:G:A", "config_path": "config_example.yaml", "gene_symbol": "BRCA1", "enable_rag_llm": true}'
```

## VSCode 调试配置

| 配置名称 | 说明 |
|---------|------|
| Python: Test Classification | 运行测试脚本 |
| Python: Classify CLI | 命令行分类 |
| Python: Classify with RAG+LLM | 启用 RAG+LLM 的分类 |
| Python: FastAPI Server | 启动 Web 服务 |
| Python: Current File | 调试当前文件 |

## 分类阈值（Bayes 评分）

| 分类 | 分数范围 |
|------|----------|
| Pathogenic (致病) | ≥ 10 |
| Likely Pathogenic (可能致病) | 6 - 9 |
| VUS (意义不明) | 0 - 5 |
| Likely Benign (可能良性) | -6 - -1 |
| Benign (良性) | ≤ -7 |

## 证据强度权重

| 强度 | 致病证据 | 良性证据 |
|------|---------|---------|
| Stand Alone | 10 | -10 |
| Very Strong | 8 | -8 |
| Strong | 4 | -4 |
| Moderate | 2 | -2 |
| Supporting | 1 | -1 |

## 三个入口的关系

```
classify() ← 核心函数，所有入口共享同一逻辑
    │
    ├── classify.py (CLI)
    │       python classify.py -i rs123456 -g BRCA1
    │
    ├── webservice.py /classify_variant (FastAPI)
    │       POST /classify_variant
    │
    └── webservice.py /chatbot/classify (FastAPI)
            POST /chatbot/classify
```

三者最终调用相同的 `classify()` 函数，分类结果完全一致。

## 统一 LLM 配置

`llm_config.py` 提供跨模块一致的 LLM API 配置：

```python
from llm_config import get_llm_config

config = get_llm_config()  # 自动选择可用 API
# 或指定 provider
config = get_llm_config("openai")
config = get_llm_config("deepseek")
```

## 依赖说明

- **pandas**: 数据处理
- **pyensembl**: Ensembl 数据库访问
- **cyvcf2**: VCF 文件处理
- **pybedtools**: BED 文件操作
- **hgvs**: HGVS 变异命名解析
- **fastapi**: Web API 服务
- **openai**: LLM API 调用

## 使用手册

详细文档见 [docs/ACMG_VARIANT_CLASSIFICATION_MANUAL.md](docs/ACMG_VARIANT_CLASSIFICATION_MANUAL.md)