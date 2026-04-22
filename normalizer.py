"""
变异标准化模块
支持多种变异格式输入，统一标准化为基因组坐标，并通过VEP/MyVariant进行注释

工作流程：
1. 格式验证 - 使用正则表达式验证输入格式
2. 坐标转换 - 通过MyVariant/VEP获取基因组坐标和HGVS命名
3. 变异注释 - 通过VEP获取标准转录本（MANE Select）的核苷酸改变和蛋白改变
4. 构建检索关键词 - 生成用于PubMed/PubTator3检索的关键词

优先级：
- MyVariant.info (主方案) → VEP (备用方案)
"""

import re
import asyncio
import aiohttp
import requests
import urllib.parse
from dataclasses import dataclass, field
from typing import Optional, Dict, Any, List, Union
import logging

logger = logging.getLogger(__name__)


# ========== VEP 客户端 (集成) ==========
class VEPClient:
    """VEP注释器 - 用于获取标准转录本的HGVS命名"""
    
    BASE_URL = "https://rest.ensembl.org"
    HEADERS = {"Content-Type": "application/json"}
    
    def __init__(self, species: str = "human"):
        self.species = species
    
    def _get_base_params(self) -> Dict:
        """获取VEP基础参数"""
        return {
            "mane": 1,                # MANE Select转录本
            "gencode_primary": 1,     # GENCODE primary
            "refseq": 1,              # 同时返回RefSeq ID
            "canonical": 1,           # 规范转录本（备选）
            "hgvs": 1,                # 生成c.和p.命名
            "protein": 1,             # 蛋白改变
            "uniprot": 1,             # Uniprot ID
            "af_gnomadg": 1,         # gnomAD基因组频率
            "af_gnomade": 1,         # gnomAD外显子组频率
            "max_af": 1,             # 全局最大AF
            "pick": 1,               # 只选最佳转录本
            "pick_allele": 1,        # 每个等位基因选一个
            # 预测工具参数
            "dbnsfp": 1,             # 返回dbNSFP预测（包含SIFT, PolyPhen, REVEL等）
            "spliceai": 1,           # SpliceAI 剪接预测
        }
    
    def annotate_rsid(self, rsid: str) -> Optional[Dict]:
        """通过rs号注释变异"""
        url = f"{self.BASE_URL}/vep/{self.species}/id/{rsid}"
        params = self._get_base_params()
        
        try:
            response = requests.get(url, params=params, headers=self.HEADERS, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            if isinstance(data, list) and len(data) > 0:
                return self._parse_response(data[0])
            return None
        except Exception as e:
            logger.error(f"VEP rsID注释失败: {e}")
            return None
    
    def annotate_region(self, chrom: str, pos: int, ref: str, alt: str) -> Optional[Dict]:
        """通过chr:pos:ref:alt注释变异"""
        # 构建region格式
        region = f"{chrom}:{pos}:{pos}:1/{alt}"
        url = f"{self.BASE_URL}/vep/{self.species}/region/{region}"
        params = self._get_base_params()
        
        try:
            response = requests.get(url, params=params, headers=self.HEADERS, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            if isinstance(data, list) and len(data) > 0:
                return self._parse_response(data[0])
            return None
        except Exception as e:
            logger.error(f"VEP区域注释失败: {e}")
            return None
    
    def _parse_response(self, data: Dict) -> Dict:
        """解析VEP响应"""
        result = {
            "input": data.get("input"),
            "most_severe_consequence": data.get("most_severe_consequence"),
            "chromosome": None,
            "position": None,
            "ref": None,
            "alt": None,
            "rs_id": None,
            "gene": None,
            "transcript": None,
            "refseq_transcript": None,
            "hgvs_c": None,
            "hgvs_p": None,
            "hgvs_g": None,
            "consequence": None,
            "impact": None,
            "mane_status": None,
            "gnomad_frequencies": {},
            "max_af": None,
            "revel": None,
            "spliceai": None,
            "sift": None,
            "polyphen": None,
        }
        
        # 提取染色体位置信息
        input_str = data.get("input", "")
        if input_str:
            # 尝试解析输入字符串获取位置
            coord_match = re.match(r'chr?(\d+|X|Y|M|MT):(\d+):([ACTG]+):([ACTG]+)', input_str, re.IGNORECASE)
            if coord_match:
                result["chromosome"] = coord_match.group(1)
                result["position"] = int(coord_match.group(2))
                result["ref"] = coord_match.group(3)
                result["alt"] = coord_match.group(4)
            else:
                # 尝试从data中获取seq_region_name和start
                result["chromosome"] = data.get("seq_region_name")
                result["position"] = data.get("start")
        
        # 提取rsID
        colocated = data.get("colocated_variants", [])
        for variant in colocated:
            if "id" in variant and variant["id"].startswith("rs"):
                result["rs_id"] = variant["id"]
                break
        
        # 提取转录本信息（标准临床转录本）
        transcripts = data.get("transcript_consequences", [])
        for transcript in transcripts:
            # 优先选择MANE Select
            if transcript.get("mane_select") or transcript.get("gencode_primary"):
                result["gene"] = transcript.get("gene_symbol")
                result["transcript"] = transcript.get("transcript_id")
                result["refseq_transcript"] = transcript.get("refseq_transcript_id")
                result["hgvs_c"] = transcript.get("hgvsc")
                result["hgvs_p"] = transcript.get("hgvsp")
                ct = transcript.get("consequence_terms")
                result["consequence"] = ct[0] if ct else None
                result["all_consequences"] = ct  # 保存完整列表
                result["impact"] = transcript.get("impact")
                result["mane_status"] = "MANE Select" if transcript.get("mane_select") else "GENCODE Primary"
                break

        # 如果没有MANE，取第一个有hgvs_p的转录本
        if not result["gene"]:
            for transcript in transcripts:
                if transcript.get("hgvsp"):
                    result["gene"] = transcript.get("gene_symbol")
                    result["transcript"] = transcript.get("transcript_id")
                    result["refseq_transcript"] = transcript.get("refseq_transcript_id")
                    result["hgvs_c"] = transcript.get("hgvsc")
                    result["hgvs_p"] = transcript.get("hgvsp")
                    ct = transcript.get("consequence_terms")
                    result["consequence"] = ct[0] if ct else None
                    result["all_consequences"] = ct  # 保存完整列表
                    result["impact"] = transcript.get("impact")
                    result["mane_status"] = "Canonical" if transcript.get("canonical") else "Other"
                    break
        
        # 提取gnomAD频率
        for variant in colocated:
            if "frequencies" in variant:
                freqs = variant["frequencies"]
                result["gnomad_frequencies"] = {
                    "gnomad": freqs.get("gnomad"),
                    "gnomad_afr": freqs.get("gnomad_afr"),
                    "gnomad_amr": freqs.get("gnomad_amr"),
                    "gnomad_eas": freqs.get("gnomad_eas"),
                    "gnomad_nfe": freqs.get("gnomad_nfe"),
                    "gnomad_sas": freqs.get("gnomad_sas"),
                }
            
            if "max_af" in variant:
                result["max_af"] = variant["max_af"]
        
        # 从顶层字段获取max_af
        if "max_af" in data:
            result["max_af"] = data["max_af"]

        # 提取预测工具分数 (从transcript_consequences的dbnsfp字段)
        for transcript in transcripts:
            if transcript.get("hgvsp") or transcript.get("mane_select") or transcript.get("gencode_primary"):
                # 尝试从transcript的dbnsfp字段获取预测分数
                dbnsfp = transcript.get("dbnsfp", {})
                if dbnsfp:
                    result["revel"] = dbnsfp.get("revel_score")
                    result["sift"] = dbnsfp.get("sift_score")
                    result["polyphen"] = dbnsfp.get("polyphen_score")
                break

        # 如果从transcript没找到，尝试从顶层dbnsfp字段获取
        if result["revel"] is None and "dbnsfp" in data:
            dbnsfp = data["dbnsfp"]
            if isinstance(dbnsfp, list) and len(dbnsfp) > 0:
                dbnsfp = dbnsfp[0]
            if isinstance(dbnsfp, dict):
                result["revel"] = dbnsfp.get("revel_score")
                result["sift"] = dbnsfp.get("sift_score")
                result["polyphen"] = dbnsfp.get("polyphen_score")

        # 提取SpliceAI分数 (从colocated_variants)
        for variant in colocated:
            if "SpliceAI" in variant:
                spliceai_str = str(variant.get("SpliceAI", ""))
                if spliceai_str and spliceai_str != ".":
                    # SpliceAI格式: "SpliceAI_score|position|distance|allocated|ds|di"
                    try:
                        parts = spliceai_str.split("|")
                        if len(parts) >= 2:
                            result["spliceai"] = float(parts[0])
                    except (ValueError, TypeError):
                        pass
                break

        return result


# 创建全局VEP客户端实例
_vep_client = VEPClient()


@dataclass
class VariantInfo:
    """标准化变异信息"""
    # 基因坐标
    chromosome: str
    position: int
    reference: str
    alternate: str
    genome_build: str

    # HGVS命名
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    hgvs_g: Optional[str] = None

    # 数据库ID
    rs_id: Optional[str] = None
    clinvar_id: Optional[str] = None

    # 基因信息
    gene: Optional[str] = None
    transcript: Optional[str] = None
    refseq_transcript: Optional[str] = None

    # 其他等位基因（用于多等位基因位点）
    other_alleles: Optional[List[Dict[str, Any]]] = None

    # 临床信息
    consequence: Optional[str] = None      # 第一个/主要consequence
    all_consequences: Optional[List[str]] = None  # 所有consequence（用于PS1_splicing判断）
    impact: Optional[str] = None
    mane_status: Optional[str] = None
    gnomad_frequencies: Optional[Dict] = None
    max_af: Optional[float] = None

    # 预测工具分数
    revel: Optional[float] = None  # REVEL 致病性预测
    spliceai: Optional[float] = None  # SpliceAI 剪接预测
    sift: Optional[float] = None  # SIFT 预测
    polyphen: Optional[float] = None  # PolyPhen 预测

    # 原始输入
    original_input: str = ""
    input_type: str = ""

    # 数据来源
    data_source: str = ""  # "myvariant" or "vep"

    def to_dict(self) -> Dict[str, Any]:
        result = {
            "chromosome": self.chromosome,
            "position": self.position,
            "reference": self.reference,
            "alternate": self.alternate,
            "genome_build": self.genome_build,
            "hgvs_c": self.hgvs_c,
            "hgvs_p": self.hgvs_p,
            "hgvs_g": self.hgvs_g,
            "rs_id": self.rs_id,
            "clinvar_id": self.clinvar_id,
            "gene": self.gene,
            "transcript": self.transcript,
            "refseq_transcript": self.refseq_transcript,
            "consequence": self.consequence,
            "impact": self.impact,
            "mane_status": self.mane_status,
            "gnomad_frequencies": self.gnomad_frequencies,
            "max_af": self.max_af,
            "revel": self.revel,
            "spliceai": self.spliceai,
            "sift": self.sift,
            "polyphen": self.polyphen,
            "original_input": self.original_input,
            "input_type": self.input_type,
            "data_source": self.data_source
        }
        if self.other_alleles:
            result["other_alleles"] = self.other_alleles
        return result

    @property
    def vcf_format(self) -> str:
        """VCF格式表示"""
        return f"{self.chromosome}:{self.position}:{self.reference}:{self.alternate}"

    @property
    def chrom_with_chr(self) -> str:
        """带chr前缀的染色体"""
        chrom = self.chromosome
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        return chrom

    @property
    def is_multiallelic(self) -> bool:
        """是否为多等位基因位点"""
        return self.other_alleles is not None and len(self.other_alleles) > 0

    @property
    def all_variants(self) -> List[Dict[str, Any]]:
        """获取所有等位基因（包括主等位基因）"""
        variants = [{
            "chromosome": self.chromosome,
            "position": self.position,
            "reference": self.reference,
            "alternate": self.alternate,
            "hgvs_c": self.hgvs_c,
            "hgvs_p": self.hgvs_p,
            "hgvs_g": self.hgvs_g,
            "gene": self.gene,
            "transcript": self.transcript
        }]
        if self.other_alleles:
            variants.extend(self.other_alleles)
        return variants


class VariantNormalizer:
    """变异标准化器 - 使用VEP/MyVariant进行注释"""
    
    # ========== 正则表达式模式 ==========
    
    # rsID正则表达式
    RSID_PATTERN = re.compile(r'^rs(\d+)$', re.IGNORECASE)
    
    # VCF/染色体位置模式 (支持 chr1:12345678:G>A 或 1:12345678:G>A)
    VCF_PATTERN = re.compile(
        r'^(?:chr)?(\d+|X|Y|M|MT):(\d+):([ACTG]+):([ACTG]+)$',
        re.IGNORECASE
    )
    
    # ========== API配置 ==========
    
    # MyVariant.info API
    MYVARIANT_API = "https://myvariant.info/v1"
    
    # VEP API
    VEP_API = "https://rest.ensembl.org"

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        self.config = config or {}
        self.genome_build = self.config.get("genome", {}).get("default_build", "GRCh38")
        self.vep_client = VEPClient()

    def normalize(self, query_type: str, input_string: str,
                 gene_symbol: Optional[str] = None,
                 genome_build: Optional[str] = None) -> VariantInfo:
        """
        标准化变异输入
        
        Args:
            query_type: 输入类型 (rsid, vcf, position)
            input_string: 输入字符串
            gene_symbol: 基因符号（可选，用于辅助解析）
            genome_build: 基因组版本
            
        Returns:
            VariantInfo: 标准化后的变异信息
        """
        if genome_build:
            self.genome_build = genome_build

        logger.info(f"Normalizing variant: type={query_type}, input={input_string}")

        # 根据输入类型选择解析方法
        if query_type == "rsid":
            return self._parse_rsid(input_string)
        elif query_type in ("vcf", "position"):
            return self._parse_vcf(input_string)
        else:
            raise ValueError(f"Unsupported query type: {query_type}")

    def _parse_rsid(self, rsid: str) -> VariantInfo:
        """
        解析rsID
        流程: MyVariant (主) → VEP (备用)
        """
        match = self.RSID_PATTERN.match(rsid.strip())
        if not match:
            raise ValueError(f"Invalid rsID format: {rsid}")

        rs_number = match.group(1)
        rs_id = f"rs{rs_number}"

        logger.info(f"Looking up rsID: {rs_id}")

        # ========== 步骤1: 尝试MyVariant ==========
        myvariant_result = self._query_myvariant_by_rsid(rs_id)
        
        if myvariant_result and myvariant_result.get("hgvs_c") and myvariant_result.get("hgvs_p"):
            # MyVariant成功获取HGVS信息
            logger.info(f"MyVariant成功获取HGVS: c.{myvariant_result.get('hgvs_c')}, p.{myvariant_result.get('hgvs_p')}")

            return VariantInfo(
                chromosome=myvariant_result.get("chromosome", ""),
                position=myvariant_result.get("position", 0),
                reference=myvariant_result.get("reference", ""),
                alternate=myvariant_result.get("alternate", ""),
                genome_build=self.genome_build,
                rs_id=rs_id,
                hgvs_c=myvariant_result.get("hgvs_c"),
                hgvs_p=myvariant_result.get("hgvs_p"),
                hgvs_g=myvariant_result.get("hgvs_g"),
                gene=myvariant_result.get("gene"),
                transcript=myvariant_result.get("transcript"),
                revel=myvariant_result.get("revel"),
                spliceai=myvariant_result.get("spliceai"),
                sift=myvariant_result.get("sift"),
                polyphen=myvariant_result.get("polyphen"),
                original_input=rsid,
                input_type="rsid",
                data_source="myvariant"
            )
        
        # ========== 步骤2: MyVariant失败，使用VEP备用 ==========
        logger.info(f"MyVariant未获取完整HGVS，使用VEP备用方案: {rs_id}")
        
        vep_result = self._vep_client.annotate_rsid(rs_id)
        
        if vep_result:
            logger.info(f"VEP成功获取注释: gene={vep_result.get('gene')}, hgvs_c={vep_result.get('hgvs_c')}")

            return VariantInfo(
                chromosome=vep_result.get("chromosome", ""),
                position=vep_result.get("position", 0),
                reference=vep_result.get("ref", ""),
                alternate=vep_result.get("alt", ""),
                genome_build=self.genome_build,
                rs_id=vep_result.get("rs_id") or rs_id,
                hgvs_c=vep_result.get("hgvs_c"),
                hgvs_p=vep_result.get("hgvs_p"),
                hgvs_g=vep_result.get("hgvs_g"),
                gene=vep_result.get("gene"),
                transcript=vep_result.get("transcript"),
                refseq_transcript=vep_result.get("refseq_transcript"),
                consequence=vep_result.get("consequence"),
                impact=vep_result.get("impact"),
                mane_status=vep_result.get("mane_status"),
                gnomad_frequencies=vep_result.get("gnomad_frequencies"),
                max_af=vep_result.get("max_af"),
                revel=vep_result.get("revel"),
                spliceai=vep_result.get("spliceai"),
                sift=vep_result.get("sift"),
                polyphen=vep_result.get("polyphen"),
                original_input=rsid,
                input_type="rsid",
                data_source="vep"
            )
        
        raise ValueError(f"Cannot resolve rsID: {rs_id}")

    def _parse_vcf(self, vcf_input: str) -> VariantInfo:
        """
        解析VCF/染色体位置格式
        流程: MyVariant (主) → VEP (备用)
        """
        # 移除"chr"前缀进行解析
        vcf_input = vcf_input.strip()
        original_input = vcf_input
        if vcf_input.lower().startswith("chr"):
            vcf_input = vcf_input[3:]

        match = self.VCF_PATTERN.match(vcf_input)
        if not match:
            raise ValueError(f"Invalid VCF/position format: {vcf_input}")

        chrom = match.group(1)
        pos = int(match.group(2))
        ref = match.group(3)
        alt = match.group(4)

        # 验证ref和alt不为空
        if not ref or not alt:
            raise ValueError(f"Invalid VCF/position format: ref or alt is empty")

        # 处理chrM -> MT
        if chrom and chrom.upper() == "M":
            chrom = "MT"

        logger.info(f"Parsed VCF: chr{chrom}, pos={pos}, ref={ref}, alt={alt}")

        # ========== 步骤1: 尝试MyVariant ==========
        myvariant_result = self._query_myvariant_by_coords(chrom, pos, ref, alt)
        
        if myvariant_result and (myvariant_result.get("hgvs_c") or myvariant_result.get("hgvs_p")):
            # MyVariant成功获取HGVS信息
            logger.info(f"MyVariant成功获取HGVS: c.{myvariant_result.get('hgvs_c')}, p.{myvariant_result.get('hgvs_p')}")

            return VariantInfo(
                chromosome=myvariant_result.get("chromosome", chrom),
                position=myvariant_result.get("position", pos),
                reference=myvariant_result.get("reference", ref),
                alternate=myvariant_result.get("alternate", alt),
                genome_build=self.genome_build,
                rs_id=myvariant_result.get("rs_id"),
                hgvs_c=myvariant_result.get("hgvs_c"),
                hgvs_p=myvariant_result.get("hgvs_p"),
                hgvs_g=myvariant_result.get("hgvs_g"),
                gene=myvariant_result.get("gene"),
                transcript=myvariant_result.get("transcript"),
                revel=myvariant_result.get("revel"),
                spliceai=myvariant_result.get("spliceai"),
                sift=myvariant_result.get("sift"),
                polyphen=myvariant_result.get("polyphen"),
                original_input=original_input,
                input_type="vcf",
                data_source="myvariant"
            )
        
        # ========== 步骤2: MyVariant失败，使用VEP备用 ==========
        logger.info(f"MyVariant未获取完整HGVS，使用VEP备用方案")
        
        vep_result = self._vep_client.annotate_region(chrom, pos, ref, alt)
        
        if vep_result:
            logger.info(f"VEP成功获取注释: gene={vep_result.get('gene')}, hgvs_c={vep_result.get('hgvs_c')}")

            return VariantInfo(
                chromosome=vep_result.get("chromosome", chrom),
                position=vep_result.get("position", pos),
                reference=vep_result.get("ref", ref),
                alternate=vep_result.get("alt", alt),
                genome_build=self.genome_build,
                rs_id=vep_result.get("rs_id"),
                hgvs_c=vep_result.get("hgvs_c"),
                hgvs_p=vep_result.get("hgvs_p"),
                hgvs_g=vep_result.get("hgvs_g"),
                gene=vep_result.get("gene"),
                transcript=vep_result.get("transcript"),
                refseq_transcript=vep_result.get("refseq_transcript"),
                consequence=vep_result.get("consequence"),
                impact=vep_result.get("impact"),
                mane_status=vep_result.get("mane_status"),
                gnomad_frequencies=vep_result.get("gnomad_frequencies"),
                max_af=vep_result.get("max_af"),
                revel=vep_result.get("revel"),
                spliceai=vep_result.get("spliceai"),
                sift=vep_result.get("sift"),
                polyphen=vep_result.get("polyphen"),
                original_input=original_input,
                input_type="vcf",
                data_source="vep"
            )
        
        # 如果VEP也失败，返回基本信息
        logger.warning(f"MyVariant和VEP都未获取完整注释，返回基本坐标信息")
        return VariantInfo(
            chromosome=chrom,
            position=pos,
            reference=ref,
            alternate=alt,
            genome_build=self.genome_build,
            original_input=original_input,
            input_type="vcf",
data_source="coordinate_only"
        )

    def _query_myvariant_by_rsid(self, rsid: str) -> Optional[Dict[str, Any]]:
        """通过MyVariant.info API查询rsID"""
        try:
            url = f"{self.MYVARIANT_API}/variant/{rsid}"
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                data = response.json()
                return self._parse_myvariant_response(data)
            else:
                logger.warning(f"MyVariant API failed for {rsid}: {response.status_code}")
                return None
        except Exception as e:
            logger.error(f"MyVariant query error: {e}")
            return None

    def _query_myvariant_by_coords(self, chrom: str, pos: int, ref: str, alt: str) -> Optional[Dict[str, Any]]:
        """通过MyVariant.info API查询染色体位置"""
        try:
            # 尝试通过VCF格式查询
            vcf_str = f"{chrom}:{pos}:{ref}:{alt}"
            url = f"{self.MYVARIANT_API}/variant/{vcf_str}"
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                data = response.json()
                return self._parse_myvariant_response(data)
            else:
                logger.warning(f"MyVariant API failed for {vcf_str}: {response.status_code}")
                return None
        except Exception as e:
            logger.error(f"MyVariant query error: {e}")
            return None

    def _parse_myvariant_response(self, data: Any) -> Optional[Dict[str, Any]]:
        """解析MyVariant响应"""
        try:
            # MyVariant返回格式可能是一个列表
            if isinstance(data, list):
                data = data[0] if data else {}
            
            if not data:
                return None
            
            result = {
                "chromosome": None,
                "position": None,
                "reference": None,
                "alternate": None,
                "rs_id": None,
                "gene": None,
                "transcript": None,
                "hgvs_c": None,
                "hgvs_p": None,
                "hgvs_g": None,
                "revel": None,
                "spliceai": None,
                "sift": None,
                "polyphen": None,
            }
            
            # 提取rsID
            if "dbsnp" in data and data["dbsnp"]:
                dbsnp = data["dbsnp"]
                result["rs_id"] = dbsnp.get("rsid") or dbsnp.get("snp")
            
            # 提取染色体和位置
            # 方式1: 从vcf字段获取
            if "vcf" in data and data["vcf"]:
                vcf = data["vcf"]
                if isinstance(vcf, dict):
                    result["chromosome"] = vcf.get("chrom")
                    result["position"] = vcf.get("pos")
                    result["reference"] = vcf.get("ref")
                    result["alternate"] = vcf.get("alt")
                elif isinstance(vcf, str):
                    parts = vcf.split(":")
                    if len(parts) >= 4:
                        result["chromosome"] = parts[0].replace("chr", "")
                        result["position"] = int(parts[1])
                        result["reference"] = parts[2]
                        result["alternate"] = parts[3]
            
            # 方式2: 从顶层字段获取
            if not result["chromosome"] and "chrom" in data:
                result["chromosome"] = str(data["chrom"])
            if not result["position"] and "pos" in data:
                result["position"] = int(data["pos"])
            
            # 清理染色体名称
            if result["chromosome"]:
                result["chromosome"] = str(result["chromosome"]).replace("chr", "")
                if result["chromosome"].upper() in ("X", "Y", "M", "MT"):
                    result["chromosome"] = result["chromosome"].upper()
            
            # 提取HGVS信息
            # 从cadd/snpeff/dbnsfp/clinvar等字段提取
            for source in ["clinvar", "dbnsfp", "snpeff"]:
                if source in data and data[source]:
                    source_data = data[source]
                    
                    # 提取hgvs
                    if "hgvs" in source_data and source_data["hgvs"]:
                        hgvs = source_data["hgvs"]
                        if isinstance(hgvs, dict):
                            result["hgvs_c"] = hgvs.get("coding") or hgvs.get("c")
                            result["hgvs_p"] = hgvs.get("protein") or hgvs.get("p")
                    
                    # 提取gene
                    if "gene" in source_data and source_data["gene"]:
                        gene_data = source_data["gene"]
                        if isinstance(gene_data, dict):
                            result["gene"] = gene_data.get("symbol") or gene_data.get("gene_symbol")
                        elif isinstance(gene_data, str):
                            result["gene"] = gene_data
            
            # 如果从特定源没找到，尝试从顶层_hgvs字段
            if not result["hgvs_c"] and not result["hgvs_p"]:
                if "_hgvs" in data and data["_hgvs"]:
                    hgvs_str = data["_hgvs"]
                    if isinstance(hgvs_str, str):
                        if ":c." in hgvs_str:
                            result["hgvs_c"] = hgvs_str.split(":c.")[-1] if ":" in hgvs_str else hgvs_str
                        if ":p." in hgvs_str:
                            result["hgvs_p"] = hgvs_str.split(":p.")[-1] if ":" in hgvs_str else hgvs_str
                        # 提取基因
                        if ":" in hgvs_str:
                            result["gene"] = hgvs_str.split(":")[0]
            
            # 提取transcript
            if "transcript" in data and data["transcript"]:
                result["transcript"] = data["transcript"]

            # 提取预测工具分数 (从dbnsfp字段)
            if "dbnsfp" in data and data["dbnsfp"]:
                dbnsfp = data["dbnsfp"]
                if isinstance(dbnsfp, dict):
                    result["revel"] = dbnsfp.get("revel_score")
                    result["sift"] = dbnsfp.get("sift_score")
                    result["polyphen"] = dbnsfp.get("polyphen_score")

            # 提取SpliceAI分数 (从cadd或spliceai字段)
            if "cadd" in data and data["cadd"]:
                cadd = data["cadd"]
                if isinstance(cadd, dict):
                    # SpliceAI may be in CADD's raw field or separate field
                    result["spliceai"] = cadd.get("spliceai") or cadd.get("SpliceAI")
                elif isinstance(cadd, (int, float)):
                    result["spliceai"] = cadd

            # 也可以直接从顶层spliceai字段获取
            if result["spliceai"] is None and "spliceai" in data:
                result["spliceai"] = data["spliceai"]

            # 构建hgvs_g
            if result["chromosome"] and result["position"]:
                result["hgvs_g"] = f"g.{result['position']}{result.get('reference', '')}>{result.get('alternate', '')}"

            # 验证有基本的HGVS信息
            if result["hgvs_c"] or result["hgvs_p"]:
                logger.info(f"MyVariant解析成功: hgvs_c={result['hgvs_c']}, hgvs_p={result['hgvs_p']}")
                return result
            
            return None
            
        except Exception as e:
            logger.error(f"Error parsing MyVariant response: {e}")
            return None

    def build_search_queries(self, variant: VariantInfo) -> List[str]:
        """
        构建文献检索查询词
        
        使用OR组合策略，以下关键词都会被用于检索：
        1. 染色体位置 (如 "chr17:43045678:G:A" 或 "17:43045678:G:A")
        2. rs编号 (如 "rs12345")
        3. 基因名称+核苷酸改变 (如 "BRCA1 c.68_69delAG")
        4. 基因名称+氨基酸改变 (如 "BRCA1 p.Arg1699Ter")
        
        对于多等位基因位点，会为每个等位基因生成对应的检索词
        """
        queries = []

        # 1. rs编号（所有等位基因共享）
        if variant.rs_id:
            queries.append(variant.rs_id)

        # 2. 基因名称（所有等位基因共享）
        if variant.gene:
            queries.append(variant.gene)

        # 获取所有等位基因（包括主等位基因和其他等位基因）
        all_alleles = variant.all_variants if variant.is_multiallelic else [{
            "chromosome": variant.chromosome,
            "position": variant.position,
            "reference": variant.reference,
            "alternate": variant.alternate,
            "hgvs_c": variant.hgvs_c,
            "hgvs_p": variant.hgvs_p,
            "gene": variant.gene
        }]

        # 为每个等位基因生成检索词
        for i, allele in enumerate(all_alleles):
            allele_prefix = f"[Allele {i+1}] " if len(all_alleles) > 1 else ""

            # 3. 染色体位置 (带和不带chr前缀)
            if allele.get("chromosome") and allele.get("position"):
                chrom = allele["chromosome"]
                if not chrom.startswith("chr"):
                    chrom = f"chr{chrom}"

                vcf_with_chr = f"{chrom}:{allele['position']}:{allele['reference']}:{allele['alternate']}"
                vcf_without_chr = f"{chrom.replace('chr', '')}:{allele['position']}:{allele['reference']}:{allele['alternate']}"
                queries.append(vcf_with_chr)
                queries.append(vcf_without_chr)

            # 4. 基因名称+核苷酸改变
            gene = allele.get("gene") or variant.gene
            hgvs_c = allele.get("hgvs_c")
            if gene and hgvs_c:
                # 提取c.后面的部分
                c_match = re.search(r'c\.(.+)$', hgvs_c, re.IGNORECASE)
                if c_match:
                    c_change = c_match.group(1)
                    queries.append(f"{gene} {c_change}")
                    queries.append(f"{gene} c.{c_change}")

            # 5. 基因名称+氨基酸改变
            hgvs_p = allele.get("hgvs_p")
            if gene and hgvs_p:
                # 提取p.后面的部分
                p_match = re.search(r'p\.(.+)$', hgvs_p, re.IGNORECASE)
                if p_match:
                    p_change = p_match.group(1)
                    queries.append(f"{gene} {p_change}")
                    queries.append(f"{gene} p.{p_change}")

        # 去重
        unique_queries = list(dict.fromkeys(queries))

        logger.info(f"Built {len(unique_queries)} search queries for {len(all_alleles)} alleles")
        return unique_queries


# 辅助函数
def normalize_variant(query_type: str, input_string: str, 
                     genome_build: str = "GRCh38") -> VariantInfo:
    """
    便捷函数：标准化变异输入
    
    Args:
        query_type: 输入类型 (rsid, vcf, position)
        input_string: 输入字符串 (如 "rs12345" 或 "17:43045678:G:A")
        genome_build: 基因组版本 (默认 GRCh38)
        
    Returns:
        VariantInfo: 标准化后的变异信息
    """
    normalizer = VariantNormalizer({"genome": {"default_build": genome_build}})
    return normalizer.normalize(query_type, input_string)


if __name__ == "__main__":
    # 测试代码
    import logging
    logging.basicConfig(level=logging.INFO)
    
    normalizer = VariantNormalizer()
    
    # 测试rsID
    print("=== 测试rsID: rs699 ===")
    try:
        result = normalizer.normalize("rsid", "rs699")
        print(f"结果: {result.to_dict()}")
        queries = normalizer.build_search_queries(result)
        print(f"检索词: {queries}")
    except Exception as e:
        print(f"错误: {e}")
    
    print("\n=== 测试VCF: 17:43045678:G:A ===")
    try:
        result = normalizer.normalize("vcf", "17:43045678:G:A")
        print(f"结果: {result.to_dict()}")
        queries = normalizer.build_search_queries(result)
        print(f"检索词: {queries}")
    except Exception as e:
        print(f"错误: {e}")
