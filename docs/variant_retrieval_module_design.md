# 变异文献检索与证据专家模块设计方案

## 1. 概述

### 1.1 目标
构建一个变异文献检索与LLM证据专家定级系统，整合多源文献检索结果，根据遗传模式路由至不同的证据专家LLM，最终输出可用于ACMG分类的证据。

### 1.2 模块架构

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         VariantRetrievalModule                                │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐                       │
│  │   NCBI        │  │   ClinVar    │  │  PubTator    │                       │
│  │  Literature   │  │   Database   │  │    API       │                       │
│  │  Retrieval    │  │   Retrieval  │  │  Retrieval   │                       │
│  │  (healthfutures│  │              │  │               │                       │
│  │   -evagg)     │  │              │  │               │                       │
│  └──────┬───────┘  └──────┬───────┘  └──────┬───────┘                       │
│         │                 │                  │                               │
│         └────────────┬────┴──────────────────┘                               │
│                      ▼                                                        │
│              ┌───────────────┐                                                │
│              │  Observation   │                                                │
│              │    Finder      │  ← 病例+变异+表型+遗传模式提取                     │
│              │  (from evagg)  │                                                │
│              └───────┬───────┘                                                │
│                      │                                                        │
│                      ▼                                                        │
│  ┌─────────────────────────────────────────────────────────────┐             │
│  │              Evidence Router (基于遗传模式)                     │             │
│  │                                                              │             │
│  │  ┌─────────────────────┐   ┌─────────────────────┐          │             │
│  │  │  Autosomal Dominant │   │  Autosomal Recessive│          │             │
│  │  │  PS2 / PP1 / PS4 /  │   │  PM3 / PP1 / PP4    │          │             │
│  │  │  PP4 Expert         │   │  Expert             │          │             │
│  │  └─────────────────────┘   └─────────────────────┘          │             │
│  └─────────────────────────────────────────────────────────────┘             │
│                              │                                               │
│                              ▼                                               │
│                    ┌─────────────────┐                                      │
│                    │  Final Review   │                                      │
│                    │     LLM         │                                      │
│                    └─────────────────┘                                      │
│                              │                                               │
│                              ▼                                               │
│                    ┌─────────────────┐                                      │
│                    │  ACMG Evidence  │                                      │
│                    │     Output      │                                      │
│                    └─────────────────┘                                      │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 2. 数据源检索模块

### 2.1 NCBI文献检索 (NCBILiteratureRetriever)

**参考实现**: 
- `healthfutures-evagg/lib/evagg/ref/ncbi.py` - NCBI API调用
- `healthfutures-evagg/lib/evagg/content/observation.py` - 病例提取
- `normalizer.py` - **检索关键词生成 (`build_search_queries` 方法)**

**数据来源**: NCBI E-utilities (PubMed)

**核心功能**: 从基因检索改为变异检索，复用healthfutures-evagg的病例提取逻辑

**检索策略**:

使用 `VariantNormalizer.build_search_queries()` 自动生成多种检索关键词:

```python
# normalizer.py 中的 build_search_queries 方法生成:
# 1. rs编号: "rs12345"
# 2. 基因名称: "BRCA1"
# 3. 染色体位置 (带和不带chr前缀): "chr17:43045678:G:A", "17:43045678:G:A"
# 4. 基因名称+核苷酸改变: "BRCA1 c.68_69delAG"
# 5. 基因名称+氨基酸改变: "BRCA1 p.Gly12Val"
```

**支持的输入格式** (来自normalizer.py):
| 输入格式 | 示例 | 查询类型 |
|---------|------|---------|
| rsID | `rs699`, `rs123456` | rsid |
| VCF/染色体位置 | `17:43045678:G:A`, `chr17:43045678:G>A` | vcf/position |
| HGVS | `NM_007294.3:c.68_69delAG` | vcf (需预处理) |

**关键词生成示例** (输入: `rs123456` 或 VCF `17:43045678:G:A`):
```python
queries = normalizer.build_search_queries(variant_info)
# 输出:
# ["rs123456", "BRCA1", 
#  "chr17:43045678:G:A", "17:43045678:G:A",
#  "BRCA1 c.68_69delAG", "BRCA1 c.68_69delAG",
#  "BRCA1 p.Gly12Val"]
```

**关键代码参考** (healthfutures-evagg ncbi.py):

```python
# NCBI E-utilities 端点
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

class NcbiClientBase:
    EUTILS_HOST = "https://eutils.ncbi.nlm.nih.gov"
    EUTILS_SEARCH_SITE = "/entrez/eutils/esearch.fcgi"
    EUTILS_FETCH_SITE = "/entrez/eutils/efetch.fcgi"

    def _esearch(self, db: str, term: str, sort: str, **extra_params) -> Any:
        # NCBI检索入口

    def _efetch(self, db: str, id: str, retmode: str = None, rettype: str = None) -> Any:
        # NCBI获取详情入口

    # IPaperLookupClient
    def search(self, query: str, **extra_params) -> Sequence[str]:
        root = self._esearch(db="pubmed", term=query, sort="relevance", **extra_params)
        pmids = [id.text for id in root.findall("./IdList/Id") if id.text]
        return pmids

    def fetch(self, paper_id: str, include_fulltext: bool = False) -> Optional[Paper]:
        # 获取论文详情，包括摘要、全文XML等
```

**PMC Open Access全文获取** (支持Licensed论文的病例提取):

```python
PMCOA_GET_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmcid}"
BIOC_GET_URL = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/{pmcid}/ascii"

def _get_license_props(self, pmcid: str) -> Dict[str, str | bool]:
    """检查论文是否可通过PMC访问"""
    # 返回: can_access, license, OA

def _get_full_text(self, props: Dict) -> str:
    """从PMC获取BioC格式的全文XML"""
```

**实现文件**: `literature_retrieval/ncbi_retriever.py`

```python
class NCBILiteratureRetriever:
    """
    NCBI PubMed literature retriever for variants.

    Uses VariantNormalizer.build_search_queries() to generate multiple
    search keywords from a normalized VariantInfo object.
    """

    ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    def __init__(self, normalizer: VariantNormalizer, deepseek_api_key: str,
                 ncbi_api_key: str = None, timeout: int = 30):
        self.normalizer = normalizer
        self.deepseek_api = DeepSeekAPI(api_key=deepseek_api_key)
        self.base_client = BaseCloudClient(base_url=self.EFETCH_URL, timeout=timeout)
        self.ncbi_api_key = ncbi_api_key

    def search_variant_literature(
        self,
        variant_info: VariantInfo,  # 来自normalizer.py的标准化变异信息
        max_results_per_query: int = 20,
        max_total: int = 100,
    ) -> List[Dict[str, Any]]:
        """
        Search PubMed for variant-related literature.

        Uses build_search_queries() to generate multiple query strategies:
        1. rsID (most specific)
        2. Gene + cDNA change
        3. Gene + protein change
        4. Chromosome position (chr:pos:ref:alt)

        Each query returns up to max_results_per_query PMIDs,
        total limited to max_total.
        """
        # Step 1: Generate all search queries from normalized variant
        queries = self.normalizer.build_search_queries(variant_info)

        # Step 2: Execute searches in parallel
        all_pmids = set()
        for query in queries:
            pmids = self._search_pubmed(query, max_results_per_query)
            all_pmids.update(pmids)
            if len(all_pmids) >= max_total:
                break

        # Step 3: Fetch and process articles
        articles = self._fetch_articles(list(all_pmids)[:max_total])

        return articles

    def _search_pubmed(self, query: str, max_results: int) -> List[str]:
        """Execute PubMed search and return PMIDs."""
        # Reference: healthfutures-evagg/lib/evagg/ref/ncbi.py
        params = {
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "retmode": "json",
            "sort": "relevance",
        }
        if self.ncbi_api_key:
            params["api_key"] = self.ncbi_api_key

        response = self._make_request(self.ESEARCH_URL, params=params)
        if response.success:
            return response.data.get("esearchresult", {}).get("idlist", [])
        return []
```

### 2.1.1 病例与变异提取 (ObservationFinder)

**参考实现**: `healthfutures-evagg/lib/evagg/content/observation.py`

这是healthfutures-evagg最核心的模块，负责从文献中提取：
1. **患者识别** (find_patients)
2. **变异描述提取** (find_variants)
3. **实体链接** (link_entities) - 将患者与变异关联
4. **遗传模式判断** (variant_inheritance)

**完整流程**:

```python
class ObservationFinder:
    """
    Extract observations (patient + variant + phenotype + inheritance) from literature.

    Adapted from healthfutures-evagg observation.py
    """

    async def find_observations(self, gene_symbol: str, paper: Paper) -> Sequence[Observation]:
        """
        Main entry point for observation extraction.

        Flow:
        1. Sanity check paper for variant mention
        2. Find variant descriptions in text
        3. Find patient identifiers in text
        4. Link variants to patients
        5. Return observations
        """
        # Step 1: Get full text sections
        paper_text, table_texts = self._get_text_sections(paper)

        # Step 2: Sanity check
        if not await self._sanity_check_paper(paper_text, gene_symbol, metadata):
            return []

        # Step 3: Find variants
        variant_descriptions = await self._find_variant_descriptions(
            paper_text, table_texts, gene_symbol, metadata
        )

        # Step 4: Find patients
        patients = await self._find_patients(paper_text, table_texts, metadata)

        # Step 5: Link entities
        if patients and variant_descriptions:
            descriptions_by_patient = await self._link_entities(
                paper_text, patients, variant_descriptions, metadata
            )
        else:
            descriptions_by_patient = {"unknown": variant_descriptions}

        # Step 6: Build observations
        observations = self._build_observations(...)
        return observations
```

**Prompt模板** (来自healthfutures-evagg):

#### 2.1.1.1 find_patients.txt
```
{{$text}}

Above is text from a paper describing genetic variants and potentially the human
patients or subjects who possessed these variants. Provide for me a list of all
of the human patients described in this text.

Patient identifiers should be short, not complete sentences. Examples:
 - "I-1", "II6", "1:IV-8", "patient", "proband", "mother", "father"

If no specific human patients are identified, provide "unknown" as your response.

Output format:
{
    "patients": ["patient I.4", "patient II.1"]
}
```

#### 2.1.1.2 find_variants.txt
```
{{$text}}

Above is text from a paper. Provide a list of all genetic variants associated
with the gene {{$gene_symbol}} described in this text.

Variant formats include:
 - "c.1234A>T", "p.K34T", "chr1:1234567A>T"
 - "rs123456789", "NM_000123.1:c.1234A>T"
 - "p.Ala55Lys", "Ala55Lys"

Output format:
{
    "variants": ["c.1234A>T", "NM_000123.1:c.2345del"]
}
```

#### 2.1.1.3 link_entities.txt
```
{{$text}}

Above is text describing patients and variants. Provide a mapping between
patient identifiers and their variants.

Patients: {{$patients}}
Variants: {{$variants}}
Gene: {{$gene_symbol}}

Output format:
{
    "patient 1": ["c.123A>T", "p.Ala44Ter"],
    "patient 2": [],
    "unmatched_variants": ["c.234del"]
}
```

#### 2.1.1.4 variant_inheritance.txt
```
{{$passage}}

Above is text describing a patient and variant.

Patient: {{$patient_descriptions}}
Variant: {{$variant_descriptions}}
Gene: {{$gene}}

Determine the mode of inheritance for this patient/variant:
- "inherited": variant was inherited from parent(s)
- "de novo": variant arose de novo (not inherited)
- "unknown": insufficient information

Output format:
{
    "variant_inheritance": "inherited" | "de novo" | "unknown"
}
```

### 2.1.2 变异解析工厂 (HGVSVariantFactory)

**参考实现**: `healthfutures-evagg/lib/evagg/content/variant.py`

```python
class HGVSVariantFactory(ICreateVariants):
    """
    Parse and normalize HGVS variant descriptions from free text.

    Handles various variant formats and normalizes them to standard HGVS.
    """

    def parse(self, text_desc: str, gene_symbol: str | None, refseq: str | None = None) -> HGVSVariant:
        """
        Parse variant description and return validated HGVSVariant.

        Supports:
        - Protein variants: p.Gly12Val, Gly12Val
        - Coding variants: c.68_69delAG, 68_69delAG
        - Genomic variants: g.12345678A>T (requires refseq)
        - rsID: rs123456789
        """

    def _create_variant_from_text(self, variant_str: str, gene_symbol: str, genome_build: str | None) -> HGVSVariant | None:
        """
        Preprocess variant string from free text before parsing.

        Handles:
        - Extracting rsID: "rs123456789"
        - Stripping gene prefix: "BRCA1:c.123A>T" -> "c.123A>T"
        - Adding missing prefixes: "Gly12Val" -> "p.Gly12Val"
        - Fixing common errors in variant descriptions
        """
```

**变异解析示例**:
```
Input: "BRCA1 p.(Gly12Val)"
Processing: Remove gene prefix -> "p.(Gly12Val)" -> "p.Gly12Val"
Output: HGVSVariant(gene_symbol="BRCA1", hgvs_desc="p.Gly12Val", refseq="NP_...")

Input: "rs123456789"
Processing: Extract rsID
Output: HGVSVariant(rsid="rs123456789", gene_symbol="BRCA1")
```

### 2.2 ClinVar数据库检索 (ClinVarRetriever)

**数据来源**: ClinVar E-utilities API (已有 `cloud_api/clinvar_client.py`)

**定位**: 获取变异在ClinVar数据库的收录情况

**核心功能**:
- 输入: 变异信息
- 检索内容:
  - ClinVar Variation ID
  - Clinical significance (临床意义分类)
  - Submission信息 (提交者、提交日期、临床意义描述)
  - Star rating (星级评分) - 0-4星
  - Condition (相关疾病)
  - Inheritance (遗传模式)

**ClinVar检索端点**:
```
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
```

**检索示例**:
```python
# Search term format
term = f'{gene_name}[Gene Name] AND {protein_change}[Protein Change]'
# e.g., "BRCA1[Gene Name] AND p.Gly12Val[Protein Change]"

# 临床意义筛选
term += ' AND (pathogenic[Clinical Significance] OR likely_pathogenic[Clinical Significance])'
```

**扩展字段**:
```python
@dataclass
class ClinVarVariantInfo:
    """Extended ClinVar variant information."""
    variation_id: str
    rs_id: Optional[str]
    clinical_significance: str
    star_rating: int  # 0-4 stars
    submissions: List[SubmissionInfo]
    conditions: List[ConditionInfo]
    inheritance: Optional[str]  # "autosomal dominant", "autosomal recessive", etc.
    review_status: str
    interpreted_conditions: List[str]
```

**实现文件**: `literature_retrieval/clinvar_retriever.py`

```python
class ClinVarRetriever:
    """
    ClinVar retriever for variant clinical information.

    Extends existing cloud_api/clinvar_client.py to include:
    - Detailed submission information
    - Star ratings
    - Inheritance mode determination
    """

    def get_variant_clinvar_info(
        self,
        gene: str,
        protein_change: str,
    ) -> ClinVarVariantInfo:
        """
        Get comprehensive ClinVar information for a variant.
        """
        # Search ClinVar for variant
        variants = self.client.search_variants(
            gene_name=gene,
            protein_change=protein_change,
            max_results=10,
        )

        if not variants:
            return None

        # Get detailed information including submissions
        variant = variants[0]
        detailed_info = self._get_detailed_variant_info(variant["variation_id"])

        return self._parse_clinvar_response(detailed_info)

    def get_inheritance_mode(self, clinvar_info: ClinVarVariantInfo) -> str:
        """
        Determine inheritance mode from ClinVar data.

        Logic:
        1. Check conditions for inheritance patterns
        2. If multiple conditions with different inheritance, default to "unknown"
        3. Otherwise return the detected inheritance pattern
        """
        # Parse conditions and inheritance information
        pass
```

### 2.3 PubTator API检索 (PubTatorRetriever)

**参考实现**: `healthfutures-evagg/lib/evagg/ref/ncbi.py` 中的 `annotate` 方法

**数据来源**: PubTator API (https://pubtator.org/)

**定位**: 从生物医学文献中提取变异相关实体和关联

**PubTator API端点**:
```
PUBTATOR_GET_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/pmc_export/biocjson?pmcids={pmcid}"
BIOC_GET_URL = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/{pmcid}/ascii"
```

**核心功能**:
- 输入: PMC论文ID
- 检索策略:
  - 获取论文的PubTator标注
  - 提取变异、基因、疾病实体
  - 提取实体之间的关系
- 信息提取: 从标注结果中提取结构化实体和关系

**实现文件**: `literature_retrieval/pubtator_retriever.py`

```python
class PubTatorRetriever:
    """
    PubTator API retriever for variant annotations.

    Retrieves variant-disease relationships and case information
    from biomedical literature using PubTator3 annotations.
    """

    PUBTATOR_API = "https://pubtator.org/api/v2"

    def __init__(self, deepseek_api_key: str, timeout: int = 60):
        self.deepseek_api = DeepSeekAPI(api_key=deepseek_api_key)
        self.base_client = BaseCloudClient(base_url=self.PUBTATOR_API, timeout=timeout)

    def search_variant_annotations(
        self,
        pmcid: str,
    ) -> Dict[str, Any]:
        """
        Search PubTator for variant annotations in a paper.

        Args:
            pmcid: PMC identifier for the paper

        Returns:
            Dict containing:
            - entities: List of annotated entities (variants, genes, diseases)
            - relations: List of entity relationships
            - passages: Text passages with annotations
        """
        # Search PubTator
        url = f"{self.PUBTATOR_API}/publications/pmc_annotations"

        response = self.base_client.get(
            endpoint=f"/publications/pmc_annotations",
            params={"pmcids": pmcid, "concepts": "variant,gene,disease"},
        )

        return self._parse_pubtator_response(response.data)

    def _parse_pubtator_response(self, response: Dict) -> Dict[str, Any]:
        """
        Parse PubTator response into structured format.

        Response format:
        {
            "pmc_id": "PMC...",
            "annotations": [
                {"text": "BRCA1", "type": "Gene", "start": 123, "end": 128},
                {"text": "p.Gly12Val", "type": "Variant", "start": 200, "end": 210},
                ...
            ],
            "relations": [
                {"type": "associated", "source": 0, "target": 1},
                ...
            ]
        }
        """
        pass
```

---

## 3. 病例与表型信息提取

### 3.1 ObservationFinder (整合自healthfutures-evagg)

**参考实现**: `healthfutures-evagg/lib/evagg/content/observation.py`

**用途**: 从文献全文中提取结构化的病例-变异-表型-遗传模式信息

**核心Prompt模板**:

#### 3.1.1 查找患者 (find_patients.txt)
```
{{$text}}

Above is text from a paper describing genetic variants and patients.
Provide a list of all patient identifiers described.

Examples of valid patient identifiers:
 - "I-1", "II6", "IV-8" (family numbering)
 - "patient 1", "patient 2" (study numbering)
 - "proband", "mother", "father" (relationships)

Output JSON:
{
    "patients": ["patient I.4", "patient II.1", "proband"]
}
```

#### 3.1.2 查找变异 (find_variants.txt)
```
{{$text}}

Above is text from a paper. List all genetic variants for gene {{$gene_symbol}}.

Variant formats to find:
 - c.1234A>T, c.1234_1235delAG
 - p.Gly12Val, p.Ala55Lys
 - rs123456789
 - chr17:43045678:G:A

Output JSON:
{
    "variants": ["c.68_69delAG", "p.Gly12Val"]
}
```

#### 3.1.3 实体链接 (link_entities.txt)
```
{{$text}}

Link patients to their variants.

Patients: {{$patients}}
Variants: {{$variants}}
Gene: {{$gene_symbol}}

Output JSON:
{
    "proband": ["c.68_69delAG", "p.Gly12Val"],
    "mother": ["c.68_69delAG"],
    "unmatched_variants": []
}
```

#### 3.1.4 遗传模式 (variant_inheritance.txt)
```
{{$passage}}

Determine inheritance for patient {{$patient}} with variant {{$variant}}.

Options: "inherited", "de novo", "unknown"

Output JSON:
{
    "variant_inheritance": "de novo"
}
```

#### 3.1.5 表型提取 (phenotypes_observation.txt)
```
{{$passage}}

Extract phenotypes for this observation.

Previously defined phenotypes: {{$candidates}}

Output JSON:
{
    "phenotypes": ["breast cancer", "ovarian cancer"]
}
```

### 3.2 全文XML处理 (from healthfutures-evagg fulltext.py)

```python
def get_sections(doc: Optional[str], include: List[str] = None, exclude: List[str] = None):
    """
    Parse BioC XML full-text document into sections.

    Section types include:
    - TITLE, ABSTRACT
    - INTRODUCTION, METHODS, RESULTS, DISCUSSION
    - TABLE, FIGURE
    - REFERENCES
    """
    # Parse XML, extract passages with section_type, text_type, offset, text

def get_fulltext(doc: Optional[str], exclude: List[str] = None) -> str:
    """
    Get full text from document, excluding specified sections.
    """
    # exclude=["AUTH_CONT", "ACK_FUND", "COMP_INT", "REF"]
    pass
```

---

## 4. 相关性判断与信息提取LLM

### 4.1 RelevanceAssessmentLLM

**用途**: 判断检索到的文献是否与目标变异相关，并提取结构化信息

**Prompt模板**:
```
You are a biomedical literature expert specializing in variant classification.

Given the following variant information:
- Gene: {gene}
- HGVS (cDNA): {hgvs_c}
- HGVS (Protein): {hgvs_p}
- VCF: {vcf}

And the following article:
- PMID: {pmid}
- Title: {title}
- Abstract: {abstract}

Tasks:
1. Determine if this article is relevant to the specific variant above.
2. If relevant, extract:
   - Case information (number of patients, demographics)
   - Variant information (genotype, allele frequency)
   - Phenotype description
   - Segregation data (if family studies)
   - Functional assay results (if any)

Output in JSON format:
{
    "relevant": true/false,
    "relevance_confidence": "high/medium/low",
    "reason": "brief explanation",
    "extracted_info": {
        "cases": [...],
        "variant_details": {...},
        "phenotype": "...",
        "segregation": {...},
        "functional_evidence": {...}
    }
}
```

**实现文件**: `literature_retrieval/llm_relevance_assessor.py`

---

## 4. 证据专家路由模块

### 4.1 路由逻辑

```python
class EvidenceRouter:
    """
    Routes variant evidence assessment based on inheritance pattern.

    Autosomal Dominant → PS2, PP1, PS4, PP4 experts
    Autosomal Recessive → PM3, PP1, PP4 experts
    """

    def __init__(self, llm_provider: LLMProvider):
        self.llm_provider = llm_provider

    def route_and_assess(
        self,
        variant_info: Dict,
        inheritance_mode: str,  # "autosomal dominant" or "autosomal recessive"
        retrieved_evidence: List[EvidenceRecord],
    ) -> EvidenceAssessmentResult:
        """
        Route to appropriate evidence experts based on inheritance.

        Returns:
            EvidenceAssessmentResult with per-expert assessments
        """
        if inheritance_mode == "autosomal dominant":
            return self._assess_dominant(
                variant_info, retrieved_evidence
            )
        elif inheritance_mode == "autosomal recessive":
            return self._assess_recessive(
                variant_info, retrieved_evidence
            )
        else:
            # Unknown inheritance - return empty
            return EvidenceAssessmentResult(...)

    def _assess_dominant(
        self,
        variant_info: Dict,
        evidence: List[EvidenceRecord],
    ) -> EvidenceAssessmentResult:
        """
        Autosomal dominant routing:
        - PS2: De novo occurrence assessment
        - PP1: Co-segregation analysis
        - PS4: Case-control studies / Case series
        - PP4: Phenotype specific for disease
        """
        ps2_result = self._call_ps2_expert(variant_info, evidence)
        pp1_result = self._call_pp1_expert(variant_info, evidence)
        ps4_result = self._call_ps4_expert(variant_info, evidence)
        pp4_result = self._call_pp4_expert(variant_info, evidence)

        return EvidenceAssessmentResult(
            ps2=ps2_result,
            pp1=pp1_result,
            ps4=ps4_result,
            pp4=pp4_result,
        )

    def _assess_recessive(
        self,
        variant_info: Dict,
        evidence: List[EvidenceRecord],
    ) -> EvidenceAssessmentResult:
        """
        Autosomal recessive routing:
        - PM3: Compound heterozygote / Homozygote assessment
        - PP1: Co-segregation analysis
        - PP4: Phenotype specific for disease
        """
        pm3_result = self._call_pm3_expert(variant_info, evidence)
        pp1_result = self._call_pp1_expert(variant_info, evidence)
        pp4_result = self._call_pp4_expert(variant_info, evidence)

        return EvidenceAssessmentResult(
            pm3=pm3_result,
            pp1=pp1_result,
            pp4=pp4_result,
        )
```

---

## 5. 证据专家LLM定义

### 5.1 PS2证据专家 (De novo变异)

**适用**: 常染色体显性遗传疾病

**职责**: 评估文献中的de novo变异证据

**Prompt模板**:
```
You are an ACMG variant classification expert evaluating PS2 evidence.

PS2 Definition: Variant is a de novo occurrence (confirmed by appropriate studies).

Given variant information:
- Gene: {gene}
- Disease: {disease}
- Inheritance: {inheritance}
- Patient: {patient_info}

And literature evidence:
{evidence_text}

Tasks:
1. Is this evidence from a confirmed de novo variant?
2. Are the parent studies adequate to confirm de novo status?
3. What is the evidence strength?

Output JSON:
{
    "applicable": true/false,
    "evidence_strength": "Very Strong/Strong/Moderate/Supporting",
    "confidence": "high/medium/low",
    "reasoning": "...",
    "case_details": {
        "confirmed_de_novo": true/false,
        "parent_testing": "description",
        "affected_status": "..."
    }
}
```

### 5.2 PP1证据专家 (共分离分析)

**适用**: 常染色体显性/隐性遗传

**职责**: 评估共分离数据证据

**Prompt模板**:
```
You are an ACMG variant classification expert evaluating PP1 evidence.

PP1 Definition: Co-segregation with disease in multiple affected family members.

Given variant information:
- Gene: {gene}
- Disease: {disease}
- Inheritance: {inheritance}

And co-segregation data:
{segregation_data}

Tasks:
1. Calculate the number of meioses supporting segregation
2. Assess if the gene is definitively known to cause the disease
3. Determine PP1 evidence strength based on:
   - Number of affected individuals tested
   - Pathogenicity of variant
   - LOD score if provided

Output JSON:
{
    "applicable": true/false,
    "evidence_strength": "Very Strong/Strong/Moderate/Supporting",
    "lod_score": float,
    "segregation_ratio": "7/9",
    "reasoning": "..."
}
```

### 5.3 PS4证据专家 (病例对照/病例系列)

**适用**: 常染色体显性遗传

**职责**: 评估病例对照研究或病例系列证据

**Prompt模板**:
```
You are an ACMG variant classification expert evaluating PS4 evidence.

PS4 Definition: Prevalence in affected individuals is significantly increased
compared to controls OR case series with well-defined controls.

Given:
- Gene: {gene}
- Disease: {disease}

And case data:
{case_data}

Tasks:
1. Is this from a case-control study or well-characterized case series?
2. Is the statistical significance adequate (p-value)?
3. What is the odds ratio or relative risk?
4. Determine PS4 strength:
   - OR > 5, p < 0.001 → Strong
   - OR > 5, p < 0.01 → Moderate
   - OR > 3, p < 0.01 → Supporting

Output JSON:
{
    "applicable": true/false,
    "evidence_strength": "Strong/Moderate/Supporting",
    "odds_ratio": float,
    "p_value": float,
    "case_count": int,
    "control_count": int,
    "reasoning": "..."
}
```

### 5.4 PP4证据专家 (表型特异性)

**适用**: 常染色体显性/隐性遗传

**职责**: 评估疾病表型与基因特异性匹配

**Prompt模板**:
```
You are an ACMG variant classification expert evaluating PP4 evidence.

PP4 Definition: Patient's phenotype or family history is highly specific
for a disease with a single genetic etiology.

Given:
- Gene: {gene}
- Reported phenotype: {phenotype}
- Disease associated with gene: {disease}

And evidence:
{evidence}

Tasks:
1. Is the phenotype highly specific for this gene's disease?
2. Are all clinical features explained by the disease?
3. Is this a single-gene etiology disease?

Output JSON:
{
    "applicable": true/false,
    "evidence_strength": "Very Strong/Strong/Moderate/Supporting",
    "phenotype_specificity": "high/medium/low",
    "features_explained": ["..."],
    "features_unexplained": ["..."],
    "reasoning": "..."
}
```

### 5.5 PM3证据专家 (复合杂合/纯合子)

**适用**: 常染色体隐性遗传

**职责**: 评估PM3证据 (反式位置检测)

**Prompt模板**:
```
You are an ACMG variant classification expert evaluating PM3 evidence.

PM3 Definition: For recessive disorders, variant found in trans with a
pathogenic variant.

Given:
- Gene: {gene}
- Disease: {disease}
- Inheritance: autosomal recessive

And evidence of second variant:
{second_variant_evidence}

Tasks:
1. Is there evidence of a second pathogenic variant in trans?
2. Is the phase confirmed (trans vs cis)?
3. Is the second variant independently pathogenic?

Output JSON:
{
    "applicable": true/false,
    "evidence_strength": "Strong/Moderate/Supporting",
    "second_variant_confirmed": true/false,
    "phase_confirmed": true/false,
    "trans_configuration": "description",
    "reasoning": "..."
}
```

---

## 6. 最终审核专家LLM

### 6.1 FinalReviewExpert

**用途**: 整合所有证据专家输出，结合其他ACMG证据，进行最终分类

```python
class FinalReviewExpert:
    """
    Final review LLM expert that synthesizes all evidence.

    Takes:
    1. Evidence expert outputs (PS2, PP1, PS4, PP4 or PM3, PP1, PP4)
    2. Other ACMG evidence from existing rules
    3. Variant information

    Outputs:
    - Final ACMG classification
    - Reasoning and confidence
    """

    SYSTEM_PROMPT = """
You are an ACMG variant classification expert providing final review.

You will receive:
1. Evidence assessment from specialized experts (PS2, PP1, PS4, PP4, PM3)
2. Other ACMG evidence from computational and database sources
3. Variant and disease information

Your task:
1. Review all evidence for consistency
2. Apply ACMG combining rules if needed
3. Provide final classification with reasoning

Return JSON:
{
    "classification": "Pathogenic/Likely Pathogenic/VUS/Likely Benign/Benign",
    "criteria": {
        "very_strong": ["PVS1", ...],
        "strong": ["PS1", ...],
        "moderate": ["PM1", ...],
        "supporting": ["PP3", ...],
        "benign": ["BS1", ...],
        "supporting_benign": ["BP1", ...]
    },
    "bayes_score": int,
    "reasoning": "...",
    "confidence": "high/medium/low"
}
"""
```

---

## 7. 数据结构定义

### 7.1 EvidenceRecord

```python
@dataclass
class EvidenceRecord:
    """Base evidence record from retrieval."""
    source: str  # "ncbi", "clinvar", "pubtator"
    pmid: Optional[str]
    variation_id: Optional[str]
    relevance_score: float
    extracted_info: Dict[str, Any]
    raw_text: str

@dataclass
class CaseInfo:
    """Extracted case information from literature."""
    patient_count: int
    demographics: str
    genotype: str
    phenotype: str
    inheritance_pattern: str
    segregation_data: Optional[SegregationData]
    functional_evidence: Optional[FunctionalData]

@dataclass
class SegregationData:
    """Co-segregation data."""
    lod_score: Optional[float]
    segregation_ratio: str  # e.g., "7/9"
    affected_tested: int
    unaffected_tested: int

@dataclass
class FunctionalData:
    """Functional assay evidence."""
    assay_type: str
    result: str  # "pathogenic", "benign", "neutral"
    p_value: Optional[float]
    description: str
```

### 7.2 EvidenceAssessmentResult

```python
@dataclass
class EvidenceAssessmentResult:
    """Result from evidence expert assessment."""
    applicable_rules: List[str]
    rule_details: Dict[str, RuleAssessment]
    final_classification: str
    bayes_score: int
    reasoning: str
    raw_llm_response: str

@dataclass
class RuleAssessment:
    """Assessment for a single rule."""
    rule_code: str
    applicable: bool
    evidence_strength: Optional[str]
    confidence: str
    reasoning: str
    supporting_evidence: List[str]
```

---

## 8. 实现文件结构

```
literature_retrieval/
├── __init__.py
├── base_retriever.py              # 基础检索类
├── ncbi_retriever.py              # NCBI文献检索 (基于healthfutures-evagg)
├── clinvar_retriever.py           # ClinVar数据库检索
├── pubtator_retriever.py          # PubTator API检索
├── observation_finder.py           # 病例/变异/表型提取 (基于healthfutures-evagg)
├── hgvs_variant_factory.py        # HGVS变异解析 (基于healthfutures-evagg)
├── fulltext_processor.py          # 全文XML处理 (基于healthfutures-evagg)
├── llm_relevance_assessor.py      # LLM相关性判断
├── evidence_router.py             # 证据路由
├── experts/
│   ├── __init__.py
│   ├── base_expert.py             # 专家基类
│   ├── ps2_expert.py              # PS2证据专家 (de novo)
│   ├── pp1_expert.py              # PP1证据专家 (共分离)
│   ├── ps4_expert.py              # PS4证据专家 (病例对照)
│   ├── pp4_expert.py              # PP4证据专家 (表型)
│   ├── pm3_expert.py              # PM3证据专家 (复合杂合)
│   └── final_review_expert.py     # 最终审核专家
├── prompts/                       # Prompt模板 (从healthfutures-evagg适配)
│   ├── system.txt
│   ├── observation/
│   │   ├── find_patients.txt
│   │   ├── find_variants.txt
│   │   ├── link_entities.txt
│   │   ├── variant_inheritance.txt
│   │   └── phenotypes_observation.txt
│   └── evidence/
│       ├── ps2_evaluation.txt
│       ├── pp1_evaluation.txt
│       ├── ps4_evaluation.txt
│       ├── pp4_evaluation.txt
│       └── pm3_evaluation.txt
├── llm_provider.py                # LLM提供者封装 (DeepSeek)
├── config.py                      # 配置管理
└── main.py                        # 主入口
```

**依赖现有模块**:
- `normalizer.py` - 变异标准化 + `build_search_queries()` 生成检索关键词
- `cloud_api/clinvar_client.py` - ClinVar API调用
- `cloud_api/base_client.py` - HTTP基础客户端

**从healthfutures-evagg复用的核心文件**:
- `lib/evagg/ref/ncbi.py` → `literature_retrieval/ncbi_retriever.py`
- `lib/evagg/content/observation.py` → `literature_retrieval/observation_finder.py`
- `lib/evagg/content/variant.py` → `literature_retrieval/hgvs_variant_factory.py`
- `lib/evagg/content/fulltext.py` → `literature_retrieval/fulltext_processor.py`
- `lib/evagg/content/prompts/*` → `literature_retrieval/prompts/*`

---

## 9. API接口设计

### 9.1 主入口函数

```python
from normalizer import VariantNormalizer, VariantInfo

def assess_variant_evidence(
    query_type: str,  # "rsid" or "vcf"
    input_string: str,  # "rs123456" or "17:43045678:G:A"
    disease: str = None,
    inheritance: str = None,  # optional, will infer from ClinVar
    config_path: str = None,
    deepseek_api_key: str = None,
) -> EvidenceAssessmentResult:
    """
    Main entry point for variant evidence assessment.

    Workflow:
    1. Normalize variant using VariantNormalizer (from normalizer.py)
    2. Generate search queries using build_search_queries()
    3. Search NCBI/ClinVar/PubTator
    4. Route to evidence experts based on inheritance mode
    5. Return final assessment

    Args:
        query_type: "rsid" or "vcf" (from normalizer.py)
        input_string: "rs123456" or "17:43045678:G:A"
        disease: Associated disease name (optional)
        inheritance: "autosomal dominant" or "autosomal recessive" (optional, infers from ClinVar)
        config_path: Path to config.yaml
        deepseek_api_key: DeepSeek API key

    Returns:
        EvidenceAssessmentResult with all evidence assessments
    """
    # Step 0: Normalize variant
    normalizer = VariantNormalizer()
    variant_info = normalizer.normalize(query_type, input_string)
    # variant_info is a VariantInfo object with:
    # - rs_id, gene, hgvs_c, hgvs_p, chromosome, position, ref, alt
    # - all_variants (for multi-allelic sites)

    # Step 1: Generate search queries
    search_queries = normalizer.build_search_queries(variant_info)
    # Returns: ["rs123456", "BRCA1", "chr17:43045678:G:A",
    #           "BRCA1 c.68_69delAG", "BRCA1 p.Gly12Val", ...]

    # Step 2: Retrieve from all sources (parallel)
    ncbi_evidence = ncbi_retriever.search_variant_literature(variant_info, max_total=100)
    clinvar_info = clinvar_retriever.get_variant_clinvar_info(variant_info)
    pubtator_evidence = pubtator_retriever.search_variant_annotations(variant_info)

    # Step 3: Determine inheritance mode
    inheritance = clinvar_info.inheritance or inheritance or infer_inheritance(variant_info)

    # Step 4: Route to appropriate experts based on inheritance mode
    router = EvidenceRouter(llm_provider)
    assessment = router.route_and_assess(
        variant_info=variant_info,
        inheritance_mode=inheritance,
        retrieved_evidence=[ncbi_evidence, clinvar_info, pubtator_evidence],
    )

    # Step 4: Final review
    final_result = final_review_expert.review(
        variant_info=variant_info,
        evidence_assessment=assessment,
        other_acmg_evidence=other_evidence,
    )

    return final_result
```

---

## 10. 配置扩展

### 10.1 config_annotation.py 扩展

需要新增配置以支持新模块:

```python
# In config_annotation.py
LITERATURE_RETRIEVAL = {
    "ncbi": {
        "enabled": True,
        "max_results": 50,
    },
    "clinvar": {
        "enabled": True,
        "include_submissions": True,
    },
    "pubtator": {
        "enabled": True,
        "timeout": 60,
    },
}

EVIDENCE_EXPERTS = {
    "model": "deepseek-chat",
    "temperature": 0.1,
    "max_tokens": 2000,
}
```

---

## 11. 与现有系统集成

### 11.1 与现有ACMG Rules的关系

```
┌─────────────────────────────────────────────────────────────┐
│                  Variant Classification Flow                 │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  1. Existing ACMG Rules (acmg_rules/)                       │
│     - PVS1, PS1, PS3, PM1, PM2, PP2, PP3, etc.             │
│     - 这些规则产生基础ACMG证据                               │
│                                                              │
│  2. Literature Retrieval (新增)                              │
│     - 产生: PS2, PP1, PS4, PP4, PM3 证据                    │
│     - 通过 EvidenceRouter 根据遗传模式路由                    │
│                                                              │
│  3. ClinGen RAG Adjuster (rag_llm_evidence_adjuster.py)    │
│     - 调整现有规则证据强度                                    │
│     - 如: PM5 Not should combined with PM1 for BRCA1        │
│                                                              │
│  4. Final Review Expert (新增)                              │
│     - 整合所有证据                                            │
│     - 计算Bayes分值                                           │
│     - 输出最终分类                                            │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

---

## 12. 优化建议

### 12.1 性能优化
1. **并行检索**: 三个数据源可并行检索，减少总等待时间
2. **缓存机制**: 对ClinVar检索结果进行缓存，避免重复查询
3. **批量处理**: PubTator支持批量查询，可合并多个变异

### 12.2 准确性优化
1. **多源交叉验证**: PS4等证据可从多个数据源交叉验证
2. **置信度传递**: 将检索相关性分数传递给证据专家调整阈值
3. **规则冲突检测**: 检测不同专家返回的矛盾结论

### 12.3 成本优化
1. **LLM调用优化**:
   - 对低相关性文献跳过专家评估
   - 批量评估多个相似证据
2. **检索策略优化**:
   - ClinVar已有数据时减少PubTator检索
   - 根据基因知名度调整检索范围

---

## 13. 实施步骤

### Phase 1: 基础复用 (1-2周)
复用healthfutures-evagg的核心代码，适配变异检索:

1. **移植NCBILiteratureRetriever**
   - 从 `healthfutures-evagg/lib/evagg/ref/ncbi.py` 移植
   - 修改为变异-centric检索策略
   - 添加rsID、HGVS等多种检索策略

2. **移植ObservationFinder**
   - 从 `healthfutures-evagg/lib/evagg/content/observation.py` 移植
   - 适配目标基因变异的病例提取
   - 保留患者识别、变异提取、实体链接逻辑

3. **移植HGVSVariantFactory**
   - 从 `healthfutures-evagg/lib/evagg/content/variant.py` 移植
   - 支持多种变异格式解析

4. **移植全文处理模块**
   - 从 `healthfutures-evagg/lib/evagg/content/fulltext.py` 移植
   - BioC XML解析

### Phase 2: ClinVar检索增强 (1周)
扩展现有 `cloud_api/clinvar_client.py`:

1. 添加星级评分获取
2. 添加Submission详情获取
3. 添加遗传模式推断

### Phase 3: 证据专家实现 (2-3周)
实现五个证据专家LLM:

1. **PS2专家** - de novo变异评估
2. **PP1专家** - 共分离分析
3. **PS4专家** - 病例对照评估
4. **PP4专家** - 表型特异性评估
5. **PM3专家** - 复合杂合评估

### Phase 4: 路由与集成 (1-2周)
1. 实现EvidenceRouter
2. 实现FinalReviewExpert
3. 与现有ACMG分类系统集成

### Phase 5: 测试与优化 (1周)
1. 功能测试
2. 性能测试
3. 错误处理完善

---

## 14. 依赖项

```python
# 新增依赖
dependencies = [
    "requests>=2.28.0",
    "urllib3>=1.26.0",
    "defusedxml>=0.7.1",      # 安全XML解析 (来自healthfutures-evagg)
    "aiohttp>=3.8.0",         # 异步HTTP (来自healthfutures-evagg)
]

# 现有依赖
existing = [
    "pybedtools",
    "uvicorn",
    "fastapi",
    "pydantic",
]
```

---

## 15. 健康futures-evagg关键代码参考

### 15.1 类型定义 (healthfutures-evagg/lib/evagg/types/)

```python
class Paper(BaseModel):
    """论文类型"""
    id: str
    pmid: str
    title: str
    abstract: str
    journal: str
    first_author: str
    pub_year: str
    doi: Optional[str]
    pmcid: Optional[str]
    citation: str
    OA: bool
    can_access: bool
    license: str
    link: str
    fulltext_xml: Optional[str] = None
    props: Dict[str, Any] = {}  # 存储所有属性的访问

class HGVSVariant(BaseModel):
    """HGVS变异类型"""
    hgvs_desc: str
    gene_symbol: Optional[str]
    refseq: Optional[str]
    refseq_predicted: bool
    valid: bool
    validation_error: Optional[str]
    protein_consequence: Optional['HGVSVariant'] = None
    coding_equivalents: List['HGVSVariant'] = []

class Observation(BaseModel):
    """观察记录类型 (病例+变异)"""
    variant: HGVSVariant
    individual: str
    variant_descriptions: List[str]
    patient_descriptions: List[str]
    texts: List[TextSection]
    paper_id: str

class TextSection(BaseModel):
    """文本段落"""
    section_type: str
    text_type: str
    offset: int
    text: str
    id: str
```

### 15.2 LLM调用接口 (healthfutures-evagg/lib/evagg/llm/)

```python
class IPromptClient(Protocol):
    """LLM客户端接口"""
    async def prompt_file(
        self,
        user_prompt_file: str,
        system_prompt: str,
        params: Dict[str, str],
        prompt_settings: Dict[str, Any],
    ) -> str:
        """调用LLM处理prompt文件"""
```

### 15.3 系统Prompt (healthfutures-evagg/lib/evagg/content/prompts/system.txt)

```
You are an intelligent assistant to a genetic analyst. Their task is to identify
the genetic variant or variants that are causing a patient's disease. One approach
they use to solve this problem is to seek out evidence from the academic literature
that supports (or refutes) the potential causal role that a given variant is
playing in a patient's disease.

As part of that process, you will assist the analyst in collecting specific details
about genetic variants that have been observed in the literature.

All of your responses should be provided in the form of a JSON object.
```

---

## 16. 与现有系统集成接口

### 16.1 webservice.py 扩展

```python
# 在现有 /chatbot/classify 端点后添加新端点

@app.post("/literature/assess")
async def assess_variant_literature(
    variant_info: VariantInfo,
    config_path: str = "/path/to/config.yaml"
) -> LiteratureAssessmentResult:
    """
    检索文献并评估PS2/PP1/PS4/PP4/PM3证据

    输入:
    - variant_info: 变异信息
    - config_path: 配置文件路径

    输出:
    - LiteratureAssessmentResult: 包含检索结果和证据评估
    """
    # 调用literature_retrieval模块
    result = assess_variant_evidence(
        variant_info=variant_info.dict(),
        config_path=config_path,
        deepseek_api_key=settings.DEEPSEEK_API_KEY,
    )
    return result
```

### 16.2 集成现有ACMG Rules

```
输入变异信息
       │
       ▼
┌─────────────────┐
│  现有ACMG Rules  │ ← PVS1, PS1, PM1, PM2, PP2, PP3等
│  (acmg_rules/)  │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Literature模块   │ ← PS2, PP1, PS4, PP4, PM3 (新增)
│ (literature_     │
│  retrieval/)    │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ ClinGen RAG     │ ← 调整证据强度
│ Adjuster        │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Final Review    │ ← 综合所有证据
│ Expert          │
└────────┬────────┘
         │
         ▼
   最终分类结果
```
