#!/usr/bin/env python3
"""
PS1_splicing 判断工具模块

根据 ClinGen Splice Variant Classification 文档实现变异分类逻辑：

Branch 1: HGVS格式符合剪接变异 (±1~+6, -1~-20)
Branch 2: VEP注释为splice site但实际是missense变异

PS1 强度判定规则:
- 经典剪接变异 (±1, ±2) + 同一核苷酸 + P → PS1
- 经典剪接变异 (±1, ±2) + 同一核苷酸 + LP → PS1_Moderate
- 同一剪接基序的非经典剪接变异 + 同一核苷酸 + P → PS1_Moderate
- 同一剪接基序的非经典剪接变异 + 同一核苷酸 + LP → PS1_Supporting
"""

import re
from dataclasses import dataclass
from enum import Enum
from typing import Optional
import logging

logger = logging.getLogger("GenOtoScope_Classify.ps1_splicing_utils")


class SplicingRegion(Enum):
    """剪接区域分类"""
    CLASSIC_DONOR = "classic_donor"      # +1, +2 (donor site)
    CLASSIC_ACCEPTOR = "classic_acceptor"  # -1, -2 (acceptor site)
    NON_CLASSIC_DONOR = "non_classic_donor"  # +3~+6 (donor region)
    NON_CLASSIC_ACCEPTOR = "non_classic_acceptor"  # -3~-20 (acceptor region)
    NOT_SPLICING = "not_splicing"


@dataclass
class SplicingVariantInfo:
    """剪接变异信息"""
    is_splicing_variant: bool
    region: SplicingRegion
    base_position: int  # c.321 中的 321
    offset: int        # +3 中的 3 (或 -5 中的 -5)
    is_donor: bool     # True for +, False for -
    hgvs_string: str   # 原始HGVS字符串
    search_window_start: int  # 检索起始位置
    search_window_end: int    # 检索结束位置


# HGVS剪接变异正则: c.数字+/-数字核苷酸>核苷酸
SPLICING_HGVS_PATTERN = re.compile(
    r'^c\.(\d+)([+-])(\d+)([A-Za-z])>([A-Za-z])$'
)


def parse_splicing_hgvs(hgvs_string: str) -> Optional[SplicingVariantInfo]:
    """
    解析剪接变异HGVS格式

    Args:
        hgvs_string: HGVS格式，如 "c.321+2G>C" 或 "c.123-5G>A"

    Returns:
        SplicingVariantInfo 或 None (如果不是剪接变异)
    """
    match = SPLICING_HGVS_PATTERN.match(hgvs_string.strip())
    if not match:
        return None

    base_position = int(match.group(1))
    sign = match.group(2)
    offset = int(match.group(3))
    ref_base = match.group(4)
    alt_base = match.group(5)

    is_donor = (sign == '+')
    is_classic = offset <= 2

    # 确定搜索窗口
    if is_donor:
        # Donor site: +1, +2 (经典) 或 +3~+6 (非经典)
        # 搜索范围: base-1 到 base+6
        search_window_start = base_position - 1
        search_window_end = base_position + 6
    else:
        # Acceptor site: -1, -2 (经典) 或 -3~-20 (非经典)
        # 搜索范围: base-20 到 base
        search_window_start = base_position - 20
        search_window_end = base_position

    # 确定区域分类
    if is_donor:
        if is_classic:
            region = SplicingRegion.CLASSIC_DONOR
        else:
            region = SplicingRegion.NON_CLASSIC_DONOR
    else:
        if is_classic:
            region = SplicingRegion.CLASSIC_ACCEPTOR
        else:
            region = SplicingRegion.NON_CLASSIC_ACCEPTOR

    return SplicingVariantInfo(
        is_splicing_variant=True,
        region=region,
        base_position=base_position,
        offset=offset,
        is_donor=is_donor,
        hgvs_string=hgvs_string,
        search_window_start=search_window_start,
        search_window_end=search_window_end,
    )


def is_classic_splicing_region(region: SplicingRegion) -> bool:
    """判断是否为经典剪接区域 (±1, ±2)"""
    return region in (SplicingRegion.CLASSIC_DONOR, SplicingRegion.CLASSIC_ACCEPTOR)


def get_ps1_strength(
    is_same_nucleotide: bool,
    is_classic_region: bool,
    is_pathogenic: bool,  # True for P, False for LP
) -> tuple[str, str]:
    """
    根据ClinGen规则确定PS1强度

    Args:
        is_same_nucleotide: ClinVar条目是否与待评估变异在同一核苷酸位置
        is_classic_region: 待评估变异是否在经典剪接区域 (±1, ±2)
        is_pathogenic: ClinVar条目是否为Pathogenic (vs Likely Pathogenic)

    Returns:
        (PS1类型, 证据强度)
    """
    if is_same_nucleotide:
        # 相同位点：P → PS1, LP → PS1_Moderate
        if is_pathogenic:
            return ("PS1", "STRONG")
        else:
            return ("PS1_Moderate", "MODERATE")
    else:
        # 不同核苷酸但同一剪接基序：PS1_Moderate 或 PS1_Supporting
        if is_classic_region and is_pathogenic:
            # 经典剪接区域的不同位点 + P → PS1_Moderate
            return ("PS1_Moderate", "MODERATE")
        else:
            # 非经典剪接区域的不同位点，或 LP → PS1_Supporting
            return ("PS1_Supporting", "SUPPORTING")

    return ("PS1", "STRONG")  # 默认


def should_use_ps1_splicing(
    hgvs_string: str,
    vep_consequence: list[str],
) -> tuple[bool, Optional[SplicingVariantInfo]]:
    """
    判断变异是否应该走PS1_splicing流程

    Args:
        hgvs_string: HGVS字符串，如 "c.321+2G>C" 或 "p.Gly12Val"
        vep_consequence: VEP注释的consequence terms列表

    Returns:
        (是否走PS1_splicing, SplicingVariantInfo)
    """
    # Branch 1: 直接HGVS格式判断 (剪接变异 c.数字+/-数字核苷酸>核苷酸)
    splicing_info = parse_splicing_hgvs(hgvs_string)
    if splicing_info and splicing_info.is_splicing_variant:
        return (True, splicing_info)

    # Branch 2: VEP注释为splice site但实际是missense (p.XXX)
    # 错义变异 + VEP=splice_donor/acceptor_site → 走PS1_splicing
    splice_terms = {
        'splice_donor_variant',
        'splice_acceptor_variant',
        'splice_donor_5th_base_variant',
        'splice_region_variant',
    }
    is_splice_vep = any(term in splice_terms for term in vep_consequence)

    if is_splice_vep:
        # 尝试从hgvs_string解析base_position (用于missense+splice_site)
        # p.Gly12Val → 提取数字12作为base_position
        missense_pattern = re.compile(r'p\.[A-Za-z]+(\d+)')
        match = missense_pattern.match(hgvs_string.strip())
        if match:
            base_position = int(match.group(1))
            # 错义+splice_site检索窗口: base-1 到 base+6
            return (True, SplicingVariantInfo(
                is_splicing_variant=True,
                region=SplicingRegion.NOT_SPLICING,  # 特殊处理标志
                base_position=base_position,
                offset=0,
                is_donor=True,  # 默认donor方向
                hgvs_string=hgvs_string,
                search_window_start=base_position - 1,
                search_window_end=base_position + 6,
            ))

    return (False, None)


def build_splicing_search_positions(splicing_info: SplicingVariantInfo) -> list[int]:
    """
    构建剪接变异检索位置列表

    对于c.321+1G>A (Donor):
    - 搜索范围: c.320, c.321, c.321+1 ~ c.321+6
    - 即: base-1, base, base+1, ..., base+6

    对于c.320-1G>A (Acceptor):
    - 搜索范围: c.300, c.301, ..., c.320
    - 即: base-20, base-19, ..., base

    对于错义+splice_site (Branch 2):
    - 搜索范围: c.位点-1 ~ c.位点+6

    Returns:
        相对位置列表 (相对于base_position)
    """
    positions = []

    if splicing_info.region == SplicingRegion.NOT_SPLICING:
        # Branch 2: missense + splice site
        # 搜索范围: base-1 到 base+6
        for offset in range(-1, 7):  # -1, 0, 1, 2, 3, 4, 5, 6
            positions.append(splicing_info.base_position + offset)
    elif splicing_info.is_donor:
        # Donor site: +1, +2 (经典) 或 +3~+6 (非经典)
        # 搜索范围: base-1 到 base+6
        for offset in range(-1, 7):  # -1, 0, 1, 2, 3, 4, 5, 6
            positions.append(splicing_info.base_position + offset)
    else:
        # Acceptor site: -1, -2 (经典) 或 -3~-20 (非经典)
        # 搜索范围: base-20 到 base
        for offset in range(-20, 1):  # -20, -19, ..., -1, 0
            positions.append(splicing_info.base_position + offset)

    return positions


# ========== 测试代码 ==========
if __name__ == "__main__":
    test_cases = [
        "c.321+1G>A",   # 经典donor
        "c.321+2G>C",   # 经典donor
        "c.321+3A>G",   # 非经典donor
        "c.321+6T>G",   # 非经典donor
        "c.123-1G>T",   # 经典acceptor
        "c.123-2C>G",   # 经典acceptor
        "c.123-5G>A",   # 非经典acceptor
        "c.123-20C>T",  # 非经典acceptor
        "c.100A>G",     # 不是剪接变异
    ]

    for hgvs in test_cases:
        result = parse_splicing_hgvs(hgvs)
        if result:
            print(f"{hgvs}: {result.region.value}, donor={result.is_donor}, "
                  f"base={result.base_position}, offset={result.offset}, "
                  f"search=[{result.search_window_start}, {result.search_window_end}]")
        else:
            print(f"{hgvs}: Not a splicing variant")