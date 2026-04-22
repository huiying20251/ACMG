#!/usr/bin/env python3
"""
Null variant 判断工具模块

使用 HGVS 正则表达式判断变异是否为 null/LOF 类型：
- 经典 splice site (±1, ±2)
- Stop gained (*, X, Ter 结尾)
- Frameshift (fs)
- Start lost (Met1)

不包含 transcript_ablation（大片段缺失走另一套流程）
"""

import re
from typing import Optional


# ========== 正则表达式定义 ==========

# 经典 splice site: ±1, ±2
# 格式: c.数字+1核苷酸>核苷酸 或 c.数字-2核苷酸>核苷酸
CLASSIC_SPLICE_PATTERN = re.compile(
    r'^c\.(\d+)([+-])([12])([A-Za-z])>([A-Za-z])$'
)

# Stop gained: p.氨基酸* 或 p.氨基酸X 或 p.氨基酸Ter
# 例如: p.Trp123*, p.Trp123X, p.Trp123Ter
STOP_GAINED_PATTERN = re.compile(
    r'^p\.\w+\d+(\*|X|Ter)$'
)

# Frameshift: 包含 fs
# 例如: p.Gly123fs*12, p.Gly123fs
FRAMESHIFT_PATTERN = re.compile(
    r'^p\.\w+\d+fs'
)

# Stop lost: 终止密码子丢失，变为氨基酸
# 格式: p.*数字氨基酸 或 p.*数字
# 例如: p.*123Gly, p.*123
STOP_LOST_PATTERN = re.compile(
    r'^p\.\*?\d+\w+$'
)

# Start lost: p.Met1 开头
# 例如: p.Met1?, p.Met1X, p.Met1fs
START_LOST_PATTERN = re.compile(
    r'^p\.Met1[^\s]*$'
)


# ========== 判断函数 ==========

def is_classic_splice(hgvs_string: str) -> bool:
    """
    判断是否为经典 splice site 变异 (±1, ±2)

    Args:
        hgvs_string: HGVS 字符串，如 "c.321+1G>C" 或 "c.123-2G>A"

    Returns:
        True if 经典 splice site
    """
    return CLASSIC_SPLICE_PATTERN.match(hgvs_string.strip()) is not None


def is_stop_gained(hgvs_string: str) -> bool:
    """
    判断是否为 stop gained (nonsense) 变异

    Args:
        hgvs_string: HGVS 字符串，如 "p.Trp123*" 或 "p.Trp123X"

    Returns:
        True if stop gained
    """
    return STOP_GAINED_PATTERN.match(hgvs_string.strip()) is not None


def is_frameshift(hgvs_string: str) -> bool:
    """
    判断是否为 frameshift 变异

    Args:
        hgvs_string: HGVS 字符串，如 "p.Gly123fs*12"

    Returns:
        True if frameshift
    """
    return FRAMESHIFT_PATTERN.match(hgvs_string.strip()) is not None


def is_stop_lost(hgvs_string: str) -> bool:
    """
    判断是否为 stop lost (终止密码子丢失) 变异

    Args:
        hgvs_string: HGVS 字符串，如 "p.*123Gly"

    Returns:
        True if stop lost
    """
    return STOP_LOST_PATTERN.match(hgvs_string.strip()) is not None


def is_start_lost(hgvs_string: str) -> bool:
    """
    判断是否为 start lost 变异

    Args:
        hgvs_string: HGVS 字符串，如 "p.Met1?" 或 "p.Met1X"

    Returns:
        True if start lost
    """
    return START_LOST_PATTERN.match(hgvs_string.strip()) is not None


def is_null_variant_by_hgvs(hgvs_string: str) -> bool:
    """
    综合判断是否为 null variant (仅通过 HGVS 字符串)

    Args:
        hgvs_string: HGVS 字符串

    Returns:
        True if any null variant type
    """
    if not hgvs_string:
        return False
    return (
        is_classic_splice(hgvs_string)
        or is_start_lost(hgvs_string)
        or is_stop_gained(hgvs_string)
        or is_frameshift(hgvs_string)
        or is_stop_lost(hgvs_string)
    )


def get_null_variant_type(hgvs_string: str) -> Optional[str]:
    """
    获取 null variant 的具体类型

    Args:
        hgvs_string: HGVS 字符串

    Returns:
        类型字符串: "classic_splice", "stop_gained", "frameshift", "stop_lost", "start_lost", 或 None

    注意: start_lost 优先检查，因为 p.Met1X, p.Met1fs 等形式虽然符合
          其他模式，但更可能是 start lost
    """
    if is_classic_splice(hgvs_string):
        return "classic_splice"
    if is_start_lost(hgvs_string):
        # start_lost 优先于 stop_gained/frameshift，因为 Met1X, Met1fs 等更可能是 start lost
        return "start_lost"
    if is_stop_gained(hgvs_string):
        return "stop_gained"
    if is_frameshift(hgvs_string):
        return "frameshift"
    if is_stop_lost(hgvs_string):
        return "stop_lost"
    return None


def get_null_variant_description(hgvs_string: str) -> str:
    """
    获取 null variant 的人类可读描述

    Args:
        hgvs_string: HGVS 字符串

    Returns:
        描述字符串: "frameshift", "stop_gained", "splice site", "stop lost", "start loss", 或 "non-NULL"
    """
    type_map = {
        "classic_splice": "splice site",
        "stop_gained": "stop_gained",
        "frameshift": "frameshift",
        "stop_lost": "stop lost",
        "start_lost": "start loss",
    }
    var_type = get_null_variant_type(hgvs_string)
    return type_map.get(var_type, "non-NULL")


# ========== 测试代码 ==========
if __name__ == "__main__":
    test_cases = [
        # 经典 splice (±1, ±2)
        ("c.321+1G>A", "classic_splice"),
        ("c.321+2G>C", "classic_splice"),
        ("c.123-1G>T", "classic_splice"),
        ("c.123-2C>G", "classic_splice"),
        # 非经典 splice (不走 null variant 流程)
        ("c.321+3A>G", None),
        ("c.321+6T>G", None),
        ("c.123-5G>A", None),
        # Stop gained
        ("p.Trp123*", "stop_gained"),
        ("p.Trp123X", "stop_gained"),
        ("p.Trp123Ter", "stop_gained"),
        ("p.Gly123*", "stop_gained"),
        # Frameshift
        ("p.Gly123fs*12", "frameshift"),
        ("p.Gly123fs", "frameshift"),
        ("p.Ser456fs*", "frameshift"),
        # Stop lost (终止密码子丢失)
        ("p.*123Gly", "stop_lost"),
        ("p.*123", "stop_lost"),
        ("p.*456Ter", "stop_lost"),
        # Start lost
        ("p.Met1?", "start_lost"),
        ("p.Met1X", "start_lost"),
        ("p.Met1fs", "start_lost"),
        # 非 null variant
        ("p.Gly123Val", None),
        ("c.321G>A", None),
        ("p.Ser456=", None),
    ]

    print("Testing null variant detection:")
    print("-" * 60)
    for hgvs, expected in test_cases:
        result = get_null_variant_type(hgvs)
        status = "✓" if result == expected else "✗"
        print(f"{status} {hgvs:20} -> {str(result):15} (expected: {expected})")
