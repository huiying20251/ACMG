#!/usr/bin/env python3
"""
Variant format converter.
Converts between normalizer VariantInfo and internal Variant format.
"""

import json
import logging
from typing import Optional, Dict, Any

from normalizer import VariantNormalizer, VariantInfo

logger = logging.getLogger(__name__)


def convert_normalizer_to_variant_json(normalized: VariantInfo) -> str:
    """
    Convert normalizer VariantInfo to internal variant JSON format.

    Args:
        normalized: VariantInfo from normalizer.normalize()

    Returns:
        JSON string in internal variant format
    """
    # Parse HGVS to extract cDNA position and variant type
    hgvs_c = normalized.hgvs_c or ""
    hgvs_p = normalized.hgvs_p or ""

    # Determine variant type from consequences
    var_type = determine_variant_type(normalized.consequence)

    # Build prediction tools dictionaries
    # Pathogenicity prediction tools (REVEL, SIFT, PolyPhen)
    pathogenicity_prediction_tools = {}
    if normalized.revel is not None:
        pathogenicity_prediction_tools["REVEL"] = normalized.revel
    if normalized.sift is not None:
        pathogenicity_prediction_tools["SIFT"] = normalized.sift
    if normalized.polyphen is not None:
        pathogenicity_prediction_tools["PolyPhen"] = normalized.polyphen

    # Splicing prediction tools (SpliceAI)
    splicing_prediction_tools = {}
    if normalized.spliceai is not None:
        splicing_prediction_tools["SpliceAI"] = normalized.spliceai

    # Build the variant JSON structure expected by load_variant
    variant_json = {
        "chr": normalized.chromosome or "",
        "pos": normalized.position or 0,
        "ref": normalized.reference or "",
        "alt": normalized.alternate or "",
        "variant_type": [var_type] if var_type else ["missense_variant"],
        "gene": normalized.gene or "",
        "hgvs_c": hgvs_c,
        "hgvs_p": hgvs_p,
        "rs_id": normalized.rs_id or "",
        "mane_select": normalized.mane_status or "",
        "transcript": normalized.transcript or "",
        "refseq_transcript": normalized.refseq_transcript or "",
        "variant_effect": build_variant_effect(normalized),
        "gnomAD": build_gnomad_entry(normalized),
        "pathogenicity_prediction_tools": pathogenicity_prediction_tools,
        "splicing_prediction_tools": splicing_prediction_tools,
        "cancer_hotspots": {},
        "functional_data": [],
        "mRNA_analysis": [],
    }

    return json.dumps(variant_json)


def determine_variant_type(consequence: Optional[str]) -> str:
    """
    Determine variant type from consequence string.

    Args:
        consequence: VEP consequence term

    Returns:
        Variant type string
    """
    if not consequence:
        return "missense_variant"

    consequence_lower = consequence.lower()

    if "frameshift" in consequence_lower:
        return "frameshift_variant"
    elif "stop" in consequence_lower and "gain" in consequence_lower:
        return "stop_gained"
    elif "start_lost" in consequence_lower:
        return "start_lost"
    elif "splice" in consequence_lower:
        return "splice_region_variant"
    elif "synonymous" in consequence_lower or "silent" in consequence_lower:
        return "synonymous_variant"
    elif "missense" in consequence_lower:
        return "missense_variant"
    elif "nonsense" in consequence_lower:
        return "stop_gained"
    elif "inframe" in consequence_lower:
        return "inframe_indel"
    elif "insertion" in consequence_lower:
        return "insertion"
    elif "deletion" in consequence_lower:
        return "deletion"
    elif "dup" in consequence_lower:
        return "duplication"
    else:
        return "missense_variant"


def build_variant_effect(normalized: VariantInfo) -> list[Dict[str, Any]]:
    """
    Build variant_effect list for transcript annotation.

    Args:
        normalized: VariantInfo from normalizer

    Returns:
        List of transcript effect dicts
    """
    effects = []

    # Main transcript
    if normalized.transcript:
        effect = {
            "transcript": normalized.transcript,
            "hgvs_c": normalized.hgvs_c or "",
            "hgvs_p": normalized.hgvs_p or "",
            "variant_type": [determine_variant_type(normalized.consequence)],
            "exon": None,
            "intron": None,
        }
        effects.append(effect)

    # RefSeq transcript if different
    if normalized.refseq_transcript and normalized.refseq_transcript != normalized.transcript:
        effect = {
            "transcript": normalized.refseq_transcript,
            "hgvs_c": normalized.hgvs_c or "",
            "hgvs_p": normalized.hgvs_p or "",
            "variant_type": [determine_variant_type(normalized.consequence)],
            "exon": None,
            "intron": None,
        }
        effects.append(effect)

    return effects


def build_gnomad_entry(normalized: VariantInfo) -> Dict[str, Any]:
    """
    Build gnomAD entry from normalized data.

    Args:
        normalized: VariantInfo from normalizer

    Returns:
        gnomAD dict
    """
    gnomad = {}

    if normalized.gnomad_frequencies:
        if "gnomadg" in normalized.gnomad_frequencies:
            gnomad["AF"] = normalized.gnomad_frequencies["gnomadg"].get("AF", 0)
            gnomad["AC"] = normalized.gnomad_frequencies["gnomadg"].get("AC", 0)
        elif "gnomade" in normalized.gnomad_frequencies:
            gnomad["AF"] = normalized.gnomad_frequencies["gnomade"].get("AF", 0)
            gnomad["AC"] = normalized.gnomad_frequencies["gnomade"].get("AC", 0)

    if normalized.max_af is not None:
        gnomad["popmax_AF"] = normalized.max_af

    # Note: missense_zscore is gene-level and not returned by VEP.
    # It should be fetched from gnomAD gene constraint API or precomputed file.
    # The field is included here for downstream processing.
    gnomad["missense_zscore"] = None

    return gnomad


def load_from_normalizer(
    query_type: str,
    input_string: str,
    gene_symbol: Optional[str] = None,
) -> str:
    """
    Load variant from normalizer and convert to internal JSON format.

    Args:
        query_type: Input type (rsid, vcf, position)
        input_string: Input string
        gene_symbol: Optional gene symbol for disambiguation

    Returns:
        JSON string in internal variant format

    Example:
        >>> variant_json = load_from_normalizer("rsid", "rs123456")
        >>> variant_json = load_from_normalizer("vcf", "17:43045678:G:A")
    """
    normalizer = VariantNormalizer()

    normalized = normalizer.normalize(
        query_type=query_type,
        input_string=input_string,
        gene_symbol=gene_symbol,
    )

    if not normalized:
        raise ValueError(f"Could not normalize variant: {input_string}")

    return convert_normalizer_to_variant_json(normalized)


def load_variant_from_text(variant_text: str) -> str:
    """
    Load variant from natural language text containing variant description.

    Supports formats:
    - "17:43045678:G:A" or "chr17:43045678:G:A"
    - "rs123456"
    - "NM_007294.3:c.68_69delAG"

    Args:
        variant_text: Text containing variant description

    Returns:
        JSON string in internal variant format
    """
    import re

    variant_text = variant_text.strip()

    # Try rsID first
    rs_match = re.match(r'^rs(\d+)$', variant_text, re.IGNORECASE)
    if rs_match:
        return load_from_normalizer("rsid", variant_text)

    # Try VCF format
    vcf_pattern = r'^(chr)?(\d+|X|Y):(\d+):([A-Z]+):([A-Z]+)$'
    vcf_match = re.match(vcf_pattern, variant_text, re.IGNORECASE)
    if vcf_match:
        chrom = vcf_match.group(2)
        pos = vcf_match.group(3)
        ref = vcf_match.group(4)
        alt = vcf_match.group(5)
        vcf_str = f"{chrom}:{pos}:{ref}:{alt}"
        return load_from_normalizer("vcf", vcf_str)

    # Try HGVS format
    if variant_text.startswith("NM_") or variant_text.startswith("ENST"):
        return load_from_normalizer("hgvs", variant_text)

    raise ValueError(f"Could not parse variant from text: {variant_text}")
