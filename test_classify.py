#!/usr/bin/env python3
"""
ACMG Variant Classification - Test Script

Usage:
    python test_classify.py
"""

import json
import pathlib
import sys

# Add current directory to path
sys.path.insert(0, str(pathlib.Path(__file__).parent))

from classify import classify


def test_json_input():
    """Test with JSON input (legacy mode)."""
    print("\n" + "="*60)
    print("Test 1: JSON Input (Legacy Mode)")
    print("="*60)

    config_path = pathlib.Path("config_example.yaml")
    if not config_path.exists():
        print("ERROR: config_example.yaml not found")
        return None

    # JSON input format
    json_input = json.dumps({
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
    })

    final_config, result = classify(
        config_path=config_path,
        variant_str=json_input,
        query_type="json",
        inheritance_pattern="AD",
    )

    result_json = json.loads(result)
    print(f"\nClassification: {result_json.get('classification', 'N/A')}")
    print(f"Bayes Score: {result_json.get('bayes_score', 'N/A')}")
    print(f"\nApplied Rules:")
    for rule in result_json.get('rules', []):
        print(f"  - {rule['rule_code']}: {rule['strength']} ({rule['status']})")

    return result_json


def test_rsid_input():
    """Test with rsID input."""
    print("\n" + "="*60)
    print("Test 2: rsID Input (with literature search)")
    print("="*60)

    config_path = pathlib.Path("config_example.yaml")
    if not config_path.exists():
        print("ERROR: config_example.yaml not found")
        return None

    final_config, result = classify(
        config_path=config_path,
        variant_str="rs123456",  # Example rsID
        query_type="rsid",
        gene_symbol="BRCA1",
        inheritance_pattern="AD",
    )

    result_json = json.loads(result)
    print(f"\nClassification: {result_json.get('classification', 'N/A')}")
    print(f"Bayes Score: {result_json.get('bayes_score', 'N/A')}")
    print(f"\nApplied Rules:")
    for rule in result_json.get('rules', []):
        print(f"  - {rule['rule_code']}: {rule['strength']} ({rule['status']})")

    return result_json


def test_vcf_input():
    """Test with VCF input."""
    print("\n" + "="*60)
    print("Test 3: VCF Input")
    print("="*60)

    config_path = pathlib.Path("config_example.yaml")
    if not config_path.exists():
        print("ERROR: config_example.yaml not found")
        return None

    final_config, result = classify(
        config_path=config_path,
        variant_str="17:43045678:G:A",  # Example VCF
        query_type="vcf",
        gene_symbol="BRCA1",
        inheritance_pattern="AD",
    )

    result_json = json.loads(result)
    print(f"\nClassification: {result_json.get('classification', 'N/A')}")
    print(f"Bayes Score: {result_json.get('bayes_score', 'N/A')}")
    print(f"\nApplied Rules:")
    for rule in result_json.get('rules', []):
        print(f"  - {rule['rule_code']}: {rule['strength']} ({rule['status']})")

    return result_json


def test_with_rag_llm():
    """Test with RAG+LLM evidence adjustment enabled."""
    print("\n" + "="*60)
    print("Test 4: VCF Input with RAG+LLM Adjustment")
    print("="*60)

    config_path = pathlib.Path("config_example.yaml")
    if not config_path.exists():
        print("ERROR: config_example.yaml not found")
        return None

    final_config, result = classify(
        config_path=config_path,
        variant_str="17:43045678:G:A",
        query_type="vcf",
        gene_symbol="BRCA1",
        inheritance_pattern="AD",
        enable_rag_llm=True,  # Enable RAG+LLM
    )

    result_json = json.loads(result)
    print(f"\nClassification: {result_json.get('classification', 'N/A')}")
    print(f"Bayes Score: {result_json.get('bayes_score', 'N/A')}")
    print(f"\nApplied Rules:")
    for rule in result_json.get('rules', []):
        print(f"  - {rule['rule_code']}: {rule['strength']} ({rule['status']})")

    # Print RAG+LLM adjustment info if present
    if 'rag_llm_adjustment' in result_json:
        print(f"\nRAG+LLM Adjustment:")
        for adj in result_json['rag_llm_adjustment'].get('adjustments', []):
            print(f"  - {adj['rule']}: {adj['action']} - {adj['reason']}")

    return result_json


def main():
    print("="*60)
    print("ACMG Variant Classification System - Test Suite")
    print("="*60)

    # Check if config exists
    config_path = pathlib.Path("config_example.yaml")
    if not config_path.exists():
        print("ERROR: config_example.yaml not found!")
        print("Please ensure you are running from the project directory.")
        sys.exit(1)

    # Run tests
    try:
        test_json_input()
    except Exception as e:
        print(f"Test 1 failed: {e}")

    try:
        test_vcf_input()
    except Exception as e:
        print(f"Test 3 failed: {e}")

    try:
        test_with_rag_llm()
    except Exception as e:
        print(f"Test 4 failed: {e}")

    print("\n" + "="*60)
    print("Tests completed")
    print("="*60)


if __name__ == "__main__":
    main()