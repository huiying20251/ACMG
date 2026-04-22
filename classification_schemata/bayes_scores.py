#!/usr/bin/env python3
"""
ACMG Bayes scoring for variant classification.

Evidence strength weights:
- Very Strong: 8分
- Strong: 4分
- Moderate: 2分
- Supporting: 1分

Pathogenic evidence = positive, Benign evidence = negative

Classification thresholds:
- ≥10: Pathogenic (5)
- 6-9: Likely Pathogenic (4)
- 0-5: Uncertain Significance (3)
- -6 to -1: Likely Benign (2)
- ≤-7: Benign (1)
"""

from typing import Dict, List, Optional, Tuple


# Base weights for each ACMG rule (per ClinGen SVI)
BASE_WEIGHTS = {
    # Very Strong Pathogenic
    "PVS1": 8,

    # Strong Pathogenic
    "PS1": 4,   # Same amino acid change
    "PS2": 4,   # de novo (confirmed)
    "PS3": 4,   # Functional studies (damaging)
    "PS4": 4,   # Prevalence in cases

    # Moderate Pathogenic
    "PM1": 2,   # Hotspot/critical domain
    "PM2": 2,   # Absent from controls
    "PM3": 2,   # In trans with known pathogenic
    "PM4": 2,   # Protein length change
    "PM5": 2,   # Different missense at same position
    "PM6": 2,   # de novo (unconfirmed)

    # Supporting Pathogenic
    "PP1": 1,    # Segregation with disease
    "PP2": 1,    # Missense in gene with low benign variation
    "PP3": 1,    # Computational evidence (damaging)
    "PP4": 1,    # Phenotype specific
    "PP5": 1,    # Reputable source (pathogenic)

    # Stand Alone Benign
    "BA1": -8,

    # Strong Benign
    "BS1": -4,  # Allele frequency in cases > controls
    "BS2": -4,  # Observed in healthy adult
    "BS3": -4,  # Functional studies (no damage)
    "BS4": -4,  # Non-segregation with disease

    # Supporting Benign
    "BP1": -1,   # Missense in gene with high benign variation
    "BP2": -1,   # In trans with known benign
    "BP3": -1,   # Frameshift in non-critical region
    "BP4": -1,   # Computational evidence (benign)
    "BP5": -1,   # Found in case with another explained variant
    "BP6": -1,   # Reputable source (benign)
    "BP7": -1,   # Silent variant, splicing not affected
}


# Evidence strength to weight multiplier
STRENGTH_WEIGHT = {
    "very_strong": 8,
    "strong": 4,
    "moderate": 2,
    "supporting": 1,
    "stand_alone": 10,
}


# PP3/BP4 专用权重 (REVEL/SpliceAI 等计算预测工具)
# PP3 + PS3 合计最高 4 分，因此 PP3 最高为 Strong (4分)
# very_strong 配置会 cap 为 strong (4分)
PP3_STRENGTH_WEIGHT = {
    "strong": 4,
    "moderate": 2,
    "supporting": 1,
    "very_strong": 4,  # capped to strong (4 points)
}


# PS4 专用权重 (病例对照研究证据)
# PS4 可基于先证者数量达到多强度 (VS/S/M/P)
# PS4 无 combined max 限制，可以到 Very Strong (8分)
PS4_STRENGTH_WEIGHT = {
    "very_strong": 8,
    "strong": 4,
    "moderate": 2,
    "supporting": 1,
}


# PP1 专用权重 (共分离证据)
# PP1 可基于共分离次数达到多强度 (S/M/P)
# 无 combined max 限制，可以到 Strong (4分)
PP1_STRENGTH_WEIGHT = {
    "strong": 4,
    "moderate": 2,
    "supporting": 1,
}


# Classification thresholds
SCORE_THRESHOLDS = {
    "pathogenic": 10,          # ≥10 → Pathogenic (Class 5)
    "likely_pathogenic": 6,    # 6-9 → Likely Pathogenic (Class 4)
    "uncertain_min": 0,        # 0-5 → Uncertain Significance (Class 3)
    "uncertain_max": 5,
    "likely_benign_min": -6,   # -6 to -1 → Likely Benign (Class 2)
    "likely_benign_max": -1,
    "benign": -7,              # ≤-7 → Benign (Class 1)
}


# Classification descriptions
CLASSIFICATION_DESCRIPTIONS = {
    5: "Pathogenic",
    4: "Likely Pathogenic",
    3: "Uncertain Significance (VUS)",
    2: "Likely Benign",
    1: "Benign",
}


def get_bayes_weight(rule_name: str, strength: str) -> float:
    """
    Get the Bayes weight for a rule with strength adjustment.

    Per ClinGen SVI recommendations, evidence strength tiers:
    - Very Strong (VS) = 8 points
    - Strong (S) = 4 points
    - Moderate (M) = 2 points
    - Supporting (P) = 1 point

    The strength parameter determines the weight directly.

    Special case for PP3/BP4 (computational prediction tools):
    - PP3 + PS3 combined max is 4 points
    - PP3 uses PP3_STRENGTH_WEIGHT (max 4, no very_strong)
    - BP4 uses same logic (benign direction)

    Args:
        rule_name: Rule name (e.g., "PM1", "BS3", "PP3", "BP4")
        strength: Evidence strength (e.g., "supporting", "moderate", "strong", "very_strong")

    Returns:
        Adjusted Bayes weight (positive for pathogenic, negative for benign)
    """
    base_weight = BASE_WEIGHTS.get(rule_name, 0)

    # Determine sign based on base weight (positive = pathogenic, negative = benign)
    sign = 1 if base_weight > 0 else -1

    # Special handling for PP3 and BP4 (computational prediction rules)
    # PP3 + PS3 combined max is 4 points, so PP3/BP4 max is strong (4 points)
    if rule_name in ("PP3", "BP4"):
        strength_weight = PP3_STRENGTH_WEIGHT.get(strength, 1)
    # Special handling for PS4 (case-control evidence with multiple strength levels)
    # PS4 can reach Very Strong (8 points) based on proband count
    elif rule_name == "PS4":
        strength_weight = PS4_STRENGTH_WEIGHT.get(strength, 1)
    # Special handling for PP1 (co-segregation evidence with multiple strength levels)
    # PP1 can reach Strong (4 points) based on segregation count
    elif rule_name == "PP1":
        strength_weight = PP1_STRENGTH_WEIGHT.get(strength, 1)
    else:
        strength_weight = STRENGTH_WEIGHT.get(strength, 1)

    return sign * strength_weight


def calculate_total_score(
    rule_results: List[Tuple[str, str, bool, str]],
) -> float:
    """
    Calculate total Bayes score from rule results.

    Args:
        rule_results: List of (rule_name, strength, status, evidence_type) tuples
                    status=True means rule applies, False means does not apply

    Returns:
        Total Bayes score
    """
    total_score = 0.0

    for rule_name, strength, status, evidence_type in rule_results:
        if status:
            # Only count rules that apply
            weight = get_bayes_weight(rule_name, strength)

            # Sign is already determined by base weight (positive=pathogenic, negative=benign)
            total_score += weight

    return total_score


def get_classification_from_score(score: float) -> int:
    """
    Get ACMG classification (1-5) from Bayes score.

    Args:
        score: Total Bayes score

    Returns:
        Classification: 1=Benign, 2=Likely Benign, 3=VUS, 4=Likely Pathogenic, 5=Pathogenic
    """
    if score >= SCORE_THRESHOLDS["pathogenic"]:
        return 5  # Pathogenic
    elif score >= SCORE_THRESHOLDS["likely_pathogenic"]:
        return 4  # Likely Pathogenic
    elif score >= SCORE_THRESHOLDS["uncertain_min"] and score <= SCORE_THRESHOLDS["uncertain_max"]:
        return 3  # Uncertain Significance
    elif score >= SCORE_THRESHOLDS["likely_benign_min"] and score <= SCORE_THRESHOLDS["likely_benign_max"]:
        return 2  # Likely Benign
    else:
        return 1  # Benign


def get_classification_description(classification: int) -> str:
    """
    Get human-readable description of classification.

    Args:
        classification: Classification number (1-5)

    Returns:
        Description string
    """
    return CLASSIFICATION_DESCRIPTIONS.get(classification, "Unknown")


def classify_with_bayes(rule_results: List[Tuple[str, str, bool, str]]) -> Tuple[int, float, str]:
    """
    Classify variant using ACMG Bayes scoring.

    Args:
        rule_results: List of (rule_name, strength, status, evidence_type) tuples

    Returns:
        Tuple of (classification, score, description)
    """
    score = calculate_total_score(rule_results)
    classification = get_classification_from_score(score)
    description = get_classification_description(classification)

    return classification, score, description


def get_score_breakdown(
    rule_results: List[Tuple[str, str, bool, str]]
) -> Dict[str, float]:
    """
    Get detailed score breakdown by evidence type.

    Args:
        rule_results: List of (rule_name, strength, status, evidence_type) tuples

    Returns:
        Dict with breakdown of pathogenic and benign scores
    """
    pathogenic_score = 0.0
    benign_score = 0.0
    applied_rules = []

    for rule_name, strength, status, evidence_type in rule_results:
        if status:
            weight = get_bayes_weight(rule_name, strength)
            if evidence_type == "pathogenic" or weight > 0:
                pathogenic_score += abs(weight)
            elif evidence_type == "benign" or weight < 0:
                benign_score -= abs(weight)

            applied_rules.append({
                "rule": rule_name,
                "strength": strength,
                "weight": weight,
                "evidence_type": evidence_type,
            })

    return {
        "total_score": pathogenic_score + benign_score,
        "pathogenic_score": pathogenic_score,
        "benign_score": benign_score,
        "applied_rules": applied_rules,
    }