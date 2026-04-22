#!/usr/bin/env python3

import pathlib
import datetime
import argparse
import re

import uvicorn
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

from classify import classify
from _version import __version__

from fastapi import Request, status
from fastapi.responses import JSONResponse
import traceback

import pybedtools

from normalizer import VariantNormalizer
from variant_converter import convert_normalizer_to_variant_json

app = FastAPI()
@app.exception_handler(Exception)
async def my_exception_handler(request: Request, exc: Exception):
    return JSONResponse(
        status_code=500,
        content={"message": ''.join(traceback.format_exc()).replace('\n', '')}
    )

class Input(BaseModel):
    config_path: str
    variant_json: str
    query_type: str = None  # json/rsid/vcf/position, auto-detected if None
    gene_symbol: str = None
    inheritance_pattern: str = "UNKNOWN"
    enable_rag_llm: bool = False  # Enable RAG+LLM evidence adjustment


class Result(BaseModel):
    result: str
    config_file: str
    scheme_name: str
    scheme_version: str
    date: str
    tool_version: str


class ChatbotInput(BaseModel):
    text: str
    config_path: str = "/path/to/config.yaml"
    inheritance_pattern: str = "UNKNOWN"
    enable_rag_llm: bool = False


class ChatbotResponse(BaseModel):
    triggered: bool
    result: str = None
    error: str = None
    normalized_variant: str = None


class Variant(BaseModel):
    variant: str


TRIGGER_KEYWORD = "变异分类"


def extract_variant_from_text(text: str) -> str:
    """
    Extract variant string from chatbot text.

    Supports formats:
    - "变异分类 17:43045678:G:A"
    - "变异分类 chr17:43045678:G:A"
    - "变异分类 rs123456"
    - "变异分类 NM_007294.3:c.68_69delAG"

    Args:
        text: Input text

    Returns:
        Extracted variant string or empty string
    """
    # Remove trigger keyword
    text_without_trigger = text.replace(TRIGGER_KEYWORD, "").strip()

    # Try to extract VCF-style variant (chr:pos:ref:alt)
    vcf_pattern = r'(chr)?(\d+|X|Y):(\d+):([A-Z]+):([A-Z]+)'
    vcf_match = re.search(vcf_pattern, text_without_trigger)
    if vcf_match:
        chrom, pos, ref, alt = vcf_match.groups()
        return f"{chrom}{pos}:{ref}:{alt}"

    # Try rsID
    rs_pattern = r'rs(\d+)'
    rs_match = re.search(rs_pattern, text_without_trigger)
    if rs_match:
        return f"rs{rs_match.group(1)}"

    # Try HGVS
    hgvs_pattern = r'(NM_\d+\.\d+:c\.\d+[A-Z>]+)'
    hgvs_match = re.search(hgvs_pattern, text_without_trigger)
    if hgvs_match:
        return hgvs_match.group(1)

    # Return whatever is left after removing trigger
    return text_without_trigger.strip()


def contains_trigger(text: str) -> bool:
    """Check if text contains the trigger keyword."""
    return TRIGGER_KEYWORD in text


@app.get("_ping")
async def _ping():
    pass


@app.post("/chatbot/classify", response_model=ChatbotResponse)
async def chatbot_classify(input: ChatbotInput) -> ChatbotResponse:
    """
    Chatbot endpoint for variant classification.

    Triggered by keyword "变异分类" in the input text.

    Example:
        POST /chatbot/classify
        {"text": "变异分类 17:43045678:G:A", "config_path": "/path/to/config.yaml"}

    Supports literature retrieval when variant is rsID, VCF, or HGVS format.
    """
    if not contains_trigger(input.text):
        return ChatbotResponse(triggered=False)

    try:
        # Extract variant from text
        variant_str = extract_variant_from_text(input.text)
        if not variant_str:
            return ChatbotResponse(
                triggered=True,
                error="Could not extract variant from text",
            )

        # Classify directly with the variant string
        # classify() will auto-detect type and trigger literature search if applicable
        config_path = pathlib.Path(input.config_path)
        if not config_path.exists():
            return ChatbotResponse(
                triggered=True,
                error=f"Config file not found: {input.config_path}",
            )

        final_config, classification_result = classify(
            config_path,
            variant_str,
            inheritance_pattern=input.inheritance_pattern,
            enable_rag_llm=input.enable_rag_llm,
        )

        date = datetime.date.today().isoformat()
        pybedtools.cleanup()

        return ChatbotResponse(
            triggered=True,
            result=classification_result,
            normalized_variant=variant_str,
        )

    except Exception as e:
        return ChatbotResponse(
            triggered=True,
            error=str(e),
        )


def build_variant_json(normalized, vep_annotation=None) -> str:
    """
    Build variant JSON from normalized data and VEP annotation.

    Uses convert_normalizer_to_variant_json to ensure consistency
    with local execution via load_variant.py.

    Args:
        normalized: VariantInfo from normalizer
        vep_annotation: VEP annotation dict (unused, kept for compatibility)

    Returns:
        JSON string in internal variant format
    """
    # Use the same conversion function as load_variant.py for consistency
    variant_json = convert_normalizer_to_variant_json(normalized)
    return variant_json


@app.get("_ping")
async def _ping():
    pass


summary = "Classify variant"


@app.post(
    "/classify_variant",
    summary=summary,
    description=f"""
          {summary}

          Supports multiple input formats (auto-detected):
          - JSON string: legacy mode, no literature search
          - rsID: e.g., "rs123456" - triggers literature search
          - VCF format: e.g., "17:43045678:G:A" - triggers literature search
          - HGVS: e.g., "NM_007294.3:c.68_69delAG" - triggers literature search

          Parameters
          ---------
          variant_json: str
              Variant input (JSON, rsID, VCF, or HGVS)
          config_path: str
              Path to classification config
          query_type: str (optional)
              Input type override: json/rsid/vcf/position
          gene_symbol: str (optional)
              Gene symbol for disambiguation
          inheritance_pattern: str (optional)
              Inheritance pattern: AD/AR/XLD/UNKNOWN (default: UNKNOWN)
          enable_rag_llm: bool (optional)
              Enable RAG+LLM evidence adjustment based on ClinGen rules (default: False)


          Returns
          ---------
          result: str
              Json string containing the variant classification result
          config: str
              Name of classification configuration file
          date: str
              Date of automated classification
          version: str
              The version of the classification tool

          """,
)
async def classify_variant(input: Input) -> Result:
    """
    Execute classification of variant.

    Supports multiple input formats:
    - JSON string (legacy mode, no literature search)
    - rsID (e.g., "rs123456") - triggers literature search
    - VCF format (e.g., "17:43045678:G:A") - triggers literature search
    - HGVS (e.g., "NM_007294.3:c.68_69delAG") - triggers literature search
    """
    variant_str = input.variant_json
    config_path = pathlib.Path(input.config_path)
    if not config_path.exists():
        raise HTTPException(
            status_code=404,
            detail=f"The config path {input.config_path} does not exist.",
        )
    final_config, classification_result = classify(
        config_path,
        variant_str,
        query_type=input.query_type,
        gene_symbol=input.gene_symbol,
        inheritance_pattern=input.inheritance_pattern,
        enable_rag_llm=input.enable_rag_llm,
    )
    date = datetime.date.today().isoformat()
    pybedtools.cleanup()
    return Result(
        result=classification_result,
        config_file=config_path.name,
        scheme_name=final_config["name"],
        scheme_version=final_config["version"],
        date=date,
        tool_version=__version__
    )


def main():
    # define CLI arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--port", action="store", default=8080, help="Port to listen on", type=int
    )
    parser.add_argument(
        "--host", action="store", default="0.0.0.0", help="Hosts to listen on", type=str
    )

    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )

    # read passed CLI arguments
    args = parser.parse_args()

    # create and run the web service
    uvicorn.run(app, host=args.host, port=args.port, reload=False)


if __name__ == "__main__":
    main()
