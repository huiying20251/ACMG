#!/usr/bin/env python3

import pathlib
import argparse
import re
import logging
import os
from typing import Optional, List, Dict, Any

from ensembl import ensembl
from load_config import load_config, get_gene_specific_config
from load_variant import load_variant, from_normalizer
from check_disease_relevant_transcript import check_disease_relevant_transcript
from information import Classification_Info
from config_annotation import (
    get_annotations_needed_from_rules,
    get_annotation_functions,
    get_unique_annotations_needed,
    execute_annotation,
    remove_rules_with_missing_annotation,
    apply_rules,
)
from create_output import create_output, create_rules_dict
from check_incompatible_rules import check_incompatible_rules
from final_classification import get_final_classifications
from variant import Variant

from os import path
import json
import sys

logger = logging.getLogger("Classify")


# Default ClinGen rules database path
DEFAULT_CLINGEN_DB_PATH = pathlib.Path(__file__).parent / "clingen_rules.db"


def _apply_rag_llm_adjustment(
    rule_dict: dict,
    gene: str,
    variant_info: dict,
    clingen_db_path: Optional[pathlib.Path] = None,
    llm_api_key: Optional[str] = None,
    provider: str = None,
) -> tuple[dict, dict]:
    """
    Apply RAG+LLM evidence adjustment based on ClinGen gene-specific rules.

    Args:
        rule_dict: Dictionary of rule results
        gene: Gene symbol
        variant_info: Variant information dict
        clingen_db_path: Path to ClinGen rules database
        llm_api_key: LLM API key (uses unified config if not provided)
        provider: "openai" or "deepseek" (auto-detected if not provided)

    Returns:
        Tuple of (adjusted_rule_dict, adjustment_result)
    """
    from llm_config import get_llm_config

    if clingen_db_path is None:
        clingen_db_path = DEFAULT_CLINGEN_DB_PATH

    if not clingen_db_path.exists():
        logger.warning(f"ClinGen DB not found at {clingen_db_path}, skipping RAG+LLM adjustment")
        return rule_dict, {"adjustments": [], "summary": "ClinGen DB not available"}

    # Check if API key is available using unified config
    llm_config = get_llm_config(provider)
    if not llm_config.available:
        logger.warning("LLM API key not configured, skipping RAG+LLM adjustment. Set OPENAI_API_KEY, DEEPSEEK_API_KEY, or LLM_API_KEY environment variable.")
        return rule_dict, {"adjustments": [], "summary": "API key not configured"}

    try:
        from rag_llm_evidence_adjuster import adjust_evidence

        # Convert rule_dict to list format expected by adjust_evidence
        rule_results = []
        for rule_name, rule_data in rule_dict.items():
            if isinstance(rule_data, dict) and "status" in rule_data:
                rule_results.append({
                    "rule_code": rule_name,
                    "strength": rule_data.get("strength", "supporting"),
                    "status": rule_data.get("status", False),
                    "comment": rule_data.get("comment", ""),
                    "evidence_type": rule_data.get("evidence_type", "pathogenic"),
                })

        if not rule_results:
            return rule_dict, {"adjustments": [], "summary": "No applicable rules"}

        # Call RAG+LLM adjuster
        result = adjust_evidence(
            gene=gene,
            rule_results=rule_results,
            clingen_db_path=clingen_db_path,
            llm_api_key=llm_api_key,
            variant_info=variant_info,
            provider=provider,
        )

        # Apply adjustments to rule_dict
        adjusted_dict = rule_dict.copy()
        for adj in result.adjustments:
            rule_code = adj.rule_code.upper()
            if rule_code in adjusted_dict:
                if adj.action.value == "remove":
                    adjusted_dict[rule_code]["status"] = False
                    adjusted_dict[rule_code]["comment"] = f"[RAG+LLM REMOVED] {adjusted_dict[rule_code].get('comment', '')}"
                elif adj.action.value in ("upgrade", "downgrade"):
                    adjusted_dict[rule_code]["strength"] = adj.adjusted_strength
                    adjusted_dict[rule_code]["comment"] = f"[RAG+LLM {adj.action.value.upper()}: {adj.original_strength}->{adj.adjusted_strength}] {adjusted_dict[rule_code].get('comment', '')}"

        adjustment_info = {
            "adjustments": [
                {"rule": a.rule_code, "action": a.action.value, "reason": a.reason}
                for a in result.adjustments
            ],
            "summary": result.summary,
        }

        logger.info(f"RAG+LLM adjustment: {len(result.adjustments)} adjustment(s) - {result.summary}")
        return adjusted_dict, adjustment_info

    except ImportError as e:
        logger.warning(f"RAG+LLM module not available: {e}")
        return rule_dict, {"adjustments": [], "summary": "Module not available"}
    except Exception as e:
        logger.warning(f"RAG+LLM adjustment failed: {e}")
        return rule_dict, {"adjustments": [], "summary": f"Error: {e}"}


def _format_clinvar_for_output(clinvar_results: dict) -> dict:
    """
    Format ClinVar results for output.

    Args:
        clinvar_results: dict with ClinVar_Type keys and ClinVar values

    Returns:
        dict formatted for JSON output
    """
    from clinvar_utils import ClinVar_Type

    formatted = {}
    for clinvar_type, clinvar_obj in clinvar_results.items():
        type_name = clinvar_type.value if isinstance(clinvar_type, ClinVar_Type) else str(clinvar_type)
        entry = {
            "pathogenic": clinvar_obj.pathogenic,
            "highest_classification": clinvar_obj.highest_classification.value if clinvar_obj.highest_classification else None,
            "ids": clinvar_obj.ids,
            "associated_ids": clinvar_obj.associated_ids,
            # Extended fields for frontend
            "review_status": clinvar_obj.review_status,
            "conditions": clinvar_obj.conditions,
            "submission_count": clinvar_obj.submission_count,
            "pathogenic_submission_count": clinvar_obj.pathogenic_submission_count,
            "is_clingen_reviewed": clinvar_obj.is_clingen_reviewed,
        }
        formatted[type_name] = entry
    return formatted


def _detect_query_type(variant_str: str) -> tuple[str, str]:
    """
    Detect input type from variant string.

    Args:
        variant_str: Input string (JSON, rsID, VCF, or HGVS)

    Returns:
        Tuple of (query_type, processed_input)
        - query_type: "json", "rsid", "vcf", or "position"
        - processed_input: The cleaned input string
    """
    variant_str = variant_str.strip()

    # Try JSON first
    if variant_str.startswith('{'):
        return "json", variant_str

    # rsID pattern: rs123456
    rs_match = re.match(r'^rs(\d+)$', variant_str, re.IGNORECASE)
    if rs_match:
        return "rsid", variant_str

    # VCF pattern: chr17:43045678:G:A or 17:43045678:G:A
    vcf_pattern = r'^(chr)?(\d+|X|Y):(\d+):([A-Z]+):([A-Z]+)$'
    vcf_match = re.match(vcf_pattern, variant_str, re.IGNORECASE)
    if vcf_match:
        chrom = vcf_match.group(2)
        pos = vcf_match.group(3)
        ref = vcf_match.group(4)
        alt = vcf_match.group(5)
        return "vcf", f"{chrom}:{pos}:{ref}:{alt}"

    # HGVS pattern: NM_007294.3:c.68_69delAG
    if variant_str.startswith('NM_') or variant_str.startswith('ENST'):
        return "position", variant_str

    # Default: treat as JSON (may fail if not valid JSON)
    return "json", variant_str


def _add_literature_search(variant: Variant, inheritance_pattern: str) -> Variant:
    """
    Add literature search to a variant loaded from JSON.

    For JSON input, we still need literature evidence (PS3/PS2/PM3/PP1)
    which is independent of VEP annotation.

    Args:
        variant: Variant object from load_variant()
        inheritance_pattern: Inheritance pattern (AD/AR/XLD/UNKNOWN)

    Returns:
        Variant with variant_literature populated
    """
    try:
        from literature_trigger import retrieve_and_assess_literature
        from literature_retrieval.literature_utils import InheritancePattern
        from normalizer import VariantInfo as NormalizerVariantInfo

        # Convert InheritancePattern string to enum
        try:
            inh_pattern = InheritancePattern[inheritance_pattern.upper()]
        except KeyError:
            inh_pattern = InheritancePattern.UNKNOWN

        # Build a minimal VariantInfo from the Variant object for literature retrieval
        var_info = variant.variant_info
        transcript_info = variant.transcript_info[0] if variant.transcript_info else None

        # Get hgvs_c from transcript if available
        hgvs_c = None
        if transcript_info and transcript_info.var_hgvs:
            # var_hgvs is a ParsedCPosedit object, convert to string
            hgvs_c = str(transcript_info.var_hgvs)

        # Create a NormalizerVariantInfo-like object for literature retrieval
        normalized = NormalizerVariantInfo(
            chromosome=var_info.chr,
            position=var_info.genomic_start,
            reference=var_info.var_ref,
            alternate=var_info.var_obs,
            genome_build="GRCh37",  # Default build, could be extracted from JSON if needed
            gene=var_info.gene_name,
            rs_id=var_info.rs_id,
            hgvs_c=hgvs_c,
            hgvs_p=var_info.hgvs_protein,
            transcript=transcript_info.transcript_id if transcript_info else None,
            consequence=None,  # Not needed for literature search
            all_consequences=[],
        )

        # Trigger literature retrieval
        literature, evidence, used_expanded = retrieve_and_assess_literature(
            normalized,
            inheritance_pattern=inh_pattern,
        )

        # Store in variant
        variant.variant_literature = {
            "literature": literature,
            "evidence": evidence,
            "used_expanded_search": used_expanded,
        }

        logger.info(
            f"Literature search completed (JSON mode): {literature.total_articles} articles, "
            f"expanded_search_used={used_expanded}"
        )

    except ImportError as e:
        logger.warning(f"Literature retrieval module not available: {e}")
    except Exception as e:
        logger.warning(f"Literature search failed: {e}")

    return variant


def classify(
    config_path: pathlib.Path,
    variant_str: str,
    query_type: str = None,
    gene_symbol: str = None,
    inheritance_pattern: str = "UNKNOWN",
    enable_rag_llm: bool = False,
    clingen_db_path: Optional[pathlib.Path] = None,
    llm_api_key: Optional[str] = None,
    provider: str = None,
) -> tuple[dict, str]:
    """
    Perform variant classification.

    Args:
        config_path: Path to configuration file
        variant_str: Variant input (JSON string, rsID, VCF format, or HGVS)
        query_type: Input type override - "json", "rsid", "vcf", "position"
                    If None, auto-detects from variant_str
        gene_symbol: Gene symbol for disambiguation (used with rsid/vcf/position)
        inheritance_pattern: Inheritance pattern for literature retrieval (AD, AR, XLD, etc.)
        enable_rag_llm: Enable RAG+LLM evidence adjustment
        clingen_db_path: Path to ClinGen rules database (uses default if not provided)
        llm_api_key: LLM API key (uses unified config if not provided)
        provider: "openai" or "deepseek" (auto-detected if not provided)

    Returns:
        Tuple of (final_config, classification_result_json)
    """
    config = load_config(config_path)

    # Auto-detect query type if not specified
    if query_type is None:
        detected_type, processed_input = _detect_query_type(variant_str)
        query_type = detected_type
        variant_str = processed_input
    else:
        detected_type, processed_input = query_type, variant_str

    # Load variant based on input type
    if query_type == "json":
        # JSON mode: use load_variant
        variant = load_variant(variant_str)
        # Fetch gnomAD missense_zscore for PP2 if gene available
        gene_symbol = variant.variant_info.gene_name
        if gene_symbol:
            try:
                from cloud_api import GnomADConstraintClient
                gnomad_client = GnomADConstraintClient()
                missense_zscore = gnomad_client.get_missense_zscore(gene_symbol)
                if missense_zscore is not None:
                    variant.gnomad_popmax.missense_zscore = missense_zscore
                    variant.gnomad_faf.missense_zscore = missense_zscore
                    logger.info(f"Fetched gnomAD missense_zscore for {gene_symbol}: {missense_zscore}")
            except ImportError:
                logger.warning("GnomADConstraintClient not available")
            except Exception as e:
                logger.warning(f"Failed to fetch gnomAD missense_zscore: {e}")
        # Trigger literature search even for JSON input
        # (literature evidence is independent of VEP annotation)
        variant = _add_literature_search(variant, inheritance_pattern)
    else:
        # VEP Normalizer mode: use from_normalizer (triggers literature retrieval)
        variant = from_normalizer(
            query_type=query_type,
            input_string=variant_str,
            gene_symbol=gene_symbol,
            inheritance_pattern=inheritance_pattern,
            enable_literature_search=True,
        )

    final_config = get_gene_specific_config(config, variant.variant_info.gene_name)
    variant_disease_relevant = check_disease_relevant_transcript(variant, final_config)
    class_info = Classification_Info()
    annotations_needed_by_rules = get_annotations_needed_from_rules(
        final_config["rules"], class_info
    )
    annotations_needed = get_unique_annotations_needed(annotations_needed_by_rules)
    annotations_set_up = get_annotation_functions(
        annotations_needed, variant_disease_relevant, final_config, class_info
    )
    annotation = execute_annotation(annotations_set_up)

    # Extract ClinVar results from annotation Info objects and store in variant for output
    for annot in annotations_set_up:
        if annot.name == 'variant_clinvar' and annot.value is not None:
            variant.clinvar_results = annot.value
            break

    annotations_needed_by_rules_filtered = remove_rules_with_missing_annotation(
        annotations_needed_by_rules
    )
    rule_results = apply_rules(annotations_needed_by_rules_filtered)
    rule_dict = create_rules_dict(rule_results)
    rule_dict_checked = check_incompatible_rules(
        rule_dict, final_config["name"], final_config["rules"]
    )

    # Apply RAG+LLM evidence adjustment if enabled
    rag_llm_adjustment = {"adjustments": [], "summary": "RAG+LLM not enabled"}
    if enable_rag_llm:
        variant_info_dict = {
            "gene": variant.variant_info.gene_name,
            "hgvs_p": variant.variant_info.hgvs_protein,
            "hgvs_c": str(variant.transcript_info[0].var_hgvs) if variant.transcript_info else None,
            "rs_id": variant.variant_info.rs_id,
            "gnomad_popmax": variant.gnomad_popmax.subpopulation_frequency if variant.gnomad_popmax else None,
        }
        rule_dict_checked, rag_llm_adjustment = _apply_rag_llm_adjustment(
            rule_dict_checked,
            gene=variant.variant_info.gene_name,
            variant_info=variant_info_dict,
            clingen_db_path=clingen_db_path,
            llm_api_key=llm_api_key,
            provider=provider,
        )

    rule_final_class = get_final_classifications(rule_dict_checked, final_config)

    # Add RAG+LLM adjustment info to output
    out_result = create_output(rule_final_class)
    out_dict = json.loads(out_result)

    if rag_llm_adjustment.get("adjustments"):
        out_dict["rag_llm_adjustment"] = rag_llm_adjustment

    # Add ClinVar info to output if available
    if hasattr(variant, 'clinvar_results') and variant.clinvar_results:
        out_dict["clinvar"] = _format_clinvar_for_output(variant.clinvar_results)

    out_result = json.dumps(out_dict)

    ensembl.clear_cache()
    return final_config, out_result


if __name__ == "__main__":
    # define CLI arguments
    parser = argparse.ArgumentParser(description="Variant ACMG Classification")

    parser.add_argument(
        "-i", "--input", default="",
        help="Variant input: JSON string, rsID (rs123456), VCF (17:43045678:G:A), or HGVS (NM_007294.3:c.68_69delAG)",
        type=str
    )
    parser.add_argument(
        "-c", "--config", default="",
        help="path to configuration",
        type=str,
    )
    parser.add_argument(
        "-o", "--output", default="",
        help="path to output file",
        type=str,
    )
    parser.add_argument(
        "-q", "--query-type", default=None,
        help="Input type override: json, rsid, vcf, position. If not specified, auto-detected",
        type=str,
        choices=["json", "rsid", "vcf", "position"],
    )
    parser.add_argument(
        "-g", "--gene", default=None,
        help="Gene symbol for disambiguation (used with rsid/vcf/position inputs)",
        type=str,
    )
    parser.add_argument(
        "--inheritance", default="UNKNOWN",
        help="Inheritance pattern for literature retrieval: AD, AR, XLD, etc. (default: UNKNOWN)",
        type=str,
    )
    parser.add_argument(
        "--enable-rag-llm",
        action="store_true",
        help="Enable RAG+LLM evidence adjustment (uses OPENAI_API_KEY, DEEPSEEK_API_KEY, or LLM_API_KEY)",
    )

    # read passed CLI arguments
    args = parser.parse_args()

    # check if arguments were given
    if args.input == "":
        input_file = sys.stdin
    if args.config == "":
        raise ValueError("No config file provided.")

    # Execute classification
    path_config = pathlib.Path(args.config)
    input = args.input
    if path.exists(input):
        with open(input) as infile:
            input = infile.read()

    final_config, result = classify(
        path_config,
        input,
        query_type=args.query_type,
        gene_symbol=args.gene,
        inheritance_pattern=args.inheritance,
        enable_rag_llm=args.enable_rag_llm,
    )

    # write classification to sout or to file
    if args.output != "":
        sys.stdout = open(args.output, "w")  # overwrite print with sout
    result_json = json.loads(result)
    final_result = json.dumps(result_json, indent=4)
    print(final_result)
