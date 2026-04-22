#!/usr/bin/env python3

import pandas as pd

import classification_schemata.schemata as Class_schema
from classification_schemata.utils import (
    get_classifications_from_rule_combinations,
    get_final_classification_from_possible_classes,
)
from classification_schemata.bayes_scores import (
    calculate_total_score,
    get_classification_from_score,
    get_classification_description,
    classify_with_bayes,
)


def get_final_classifications(rules: dict, config: dict) -> dict:
    """
    Get final classification for variants
    """
    rules_df = pd.DataFrame(rules).transpose()
    applicable_rules = rules_df[rules_df.status == True]
    # Get final classification splice evidence
    rules_splicing = applicable_rules[
        applicable_rules.rule_type.isin(["splicing", "general"])
    ]
    class_splicing = get_classification(
        rules_splicing, config["name"], config["version"]
    )
    # Get final classification protein evidence
    rules_protein = applicable_rules[
        applicable_rules.rule_type.isin(["protein", "general"])
    ]
    class_protein = get_classification(rules_protein, config["name"], config["version"])
    # Add results to dictionary
    rules["classification_protein"] = class_protein
    rules["classification_splicing"] = class_splicing
    return rules


def get_classification(rule_results: pd.DataFrame, config: str, version: str) -> int:
    """
    Execute final classification
    """
    schema = VERSION_CLASS_SCHEMATA.get(config, {}).get(version, None)
    if schema is None:
        raise ValueError(
            f"No final classification schemata defined for configuration {config} version {version}. Please check."
        )
    abs_rule_results = create_evidence_strength_count(rule_results)
    abs_rule_results["applicable_rules"] = rule_results.index
    possible_classes = get_classifications_from_rule_combinations(
        schema, abs_rule_results
    )
    final_class = get_final_classification_from_possible_classes(possible_classes)
    return final_class


def get_classification_bayes(rules: dict) -> dict:
    """
    Get final classification using Bayes score method.

    Args:
        rules: Dictionary of rule results with status, strength, evidence_type

    Returns:
        Dict with classification, score, and description
    """
    # Convert rules dict to list of tuples for Bayes scoring
    rule_results = []
    for rule_name, rule_data in rules.items():
        if isinstance(rule_data, dict):
            status = rule_data.get("status", False)
            strength = rule_data.get("strength", "supporting")
            evidence_type = rule_data.get("evidence_type", "pathogenic")
        else:
            continue

        # Extract base rule name (remove any suffix)
        base_rule_name = rule_name.upper()
        if "_" in base_rule_name:
            base_rule_name = base_rule_name.split("_")[0]

        rule_results.append((base_rule_name, strength, status, evidence_type))

    classification, score, description = classify_with_bayes(rule_results)

    return {
        "classification": classification,
        "score": score,
        "description": description,
    }


def create_evidence_strength_count(rules: pd.DataFrame) -> dict:
    """
    From rules table create dictionary counting how often each evidence strenght is applicable
    """
    rules["evidence"] = rules.evidence_type + "_" + rules.strength
    counts = rules.groupby("evidence").count()
    return counts.rule_type.to_dict()


VERSION_CLASS_SCHEMATA = {
    "ACMG standard + SVI": {"1.0.0": Class_schema.schema_acmg},
    "ACMG ATM": {"1.3.0": Class_schema.schema_atm},
    "ACMG BRCA1": {"1.1.0": Class_schema.schema_brca1},
    "ACMG BRCA2": {"1.1.0": Class_schema.schema_brca2},
    "ACMG CDH1": {"3.1.0": Class_schema.schema_cdh1},
    "ACMG PALB2": {"1.1.0": Class_schema.schema_palb2},
    "ACMG PTEN": {"3.1.0": Class_schema.schema_pten},
    "ACMG TP53": {"1.4.0": Class_schema.schema_tp53},
}
