#!/usr/bin/env python3
"""
LLM Prompts Configuration Loader

Loads prompts from YAML configuration file.
All prompts used for ACMG variant classification LLM analysis are centralized here.
"""

import pathlib
import logging
from typing import Dict, Any

import yaml

logger = logging.getLogger(__name__)

# Default prompts path (relative to project root)
DEFAULT_PROMPTS_PATH = pathlib.Path(__file__).parent / "config" / "llm_prompts.yaml"

# Cache for loaded prompts
_prompts_cache: Dict[str, Any] = {}


def load_prompts(prompts_path: pathlib.Path = None) -> Dict[str, Any]:
    """
    Load prompts from YAML configuration file.

    Args:
        prompts_path: Path to prompts YAML file. Uses default if not provided.

    Returns:
        Dictionary of prompts
    """
    global _prompts_cache

    if prompts_path is None:
        prompts_path = DEFAULT_PROMPTS_PATH

    cache_key = str(prompts_path)

    # Return cached prompts if available
    if cache_key in _prompts_cache:
        return _prompts_cache[cache_key]

    if not prompts_path.exists():
        logger.warning(f"Prompts file not found: {prompts_path}. Using fallback prompts.")
        return _get_fallback_prompts()

    try:
        with open(prompts_path, 'r', encoding='utf-8') as f:
            prompts = yaml.safe_load(f)

        if prompts is None or 'prompts' not in prompts:
            logger.warning(f"Invalid prompts file format: {prompts_path}")
            return _get_fallback_prompts()

        _prompts_cache[cache_key] = prompts['prompts']
        logger.info(f"Loaded {len(prompts['prompts'])} prompts from {prompts_path}")
        return prompts['prompts']

    except Exception as e:
        logger.error(f"Error loading prompts from {prompts_path}: {e}")
        return _get_fallback_prompts()


def get_prompt(prompt_name: str, prompts_path: pathlib.Path = None) -> str:
    """
    Get a specific prompt by name.

    Args:
        prompt_name: Name of the prompt (e.g., 'literature_classification')
        prompts_path: Optional path to prompts file

    Returns:
        Prompt string
    """
    prompts = load_prompts(prompts_path)

    if prompt_name not in prompts:
        logger.warning(f"Prompt '{prompt_name}' not found. Using fallback.")
        fallback = _get_fallback_prompts()
        if prompt_name in fallback:
            return fallback[prompt_name]
        raise ValueError(f"Prompt '{prompt_name}' not found in fallback prompts either")

    return prompts[prompt_name]


def _get_fallback_prompts() -> Dict[str, str]:
    """
    Return fallback prompts (same as original hardcoded prompts).
    Used when config file is not available.
    """
    return {
        "literature_classification": """You are a genetic variant literature classifier.

## Task 1: Determine if article is variant-related
The article MUST be directly about this specific variant or its functional impact.

## Task 2: Classify literature type

### CASE REPORT (病例报道):
- Single patient or small number of patients
- Usually describes individual clinical presentation
- Keywords: "patient", "case", "proband", "diagnosed with"
- Example: "A 5-year-old girl with developmental delay..."

### COHORT STUDY (队列研究):
- Systematic study of multiple patients
- Usually has statistical analysis
- Keywords: "cohort", "case-control", "odds ratio", "p-value", "statistically significant"
- Example: "We studied 150 patients with..."

### FUNCTIONAL STUDY (功能研究):
- In vitro or in vivo functional experiments
- Measures protein function, splicing, enzymatic activity, etc.
- Keywords: "functional assay", "minigene", "RNA analysis", "patient RNA", "splicing assay"
- Example: "Minigene assay showed 80% exon skipping..."

### REVIEW (综述):
- Summarizes multiple studies
- No new patient data
- Keywords: "review", "meta-analysis", "systematic review"

## Input Information:
- Gene: {gene}
- Variant: {variant}
- Title: {title}
- Abstract: {abstract}

## Output Format (JSON):
{{
    "is_variant_related": true/false,
    "literature_type": "case_report/cohort_study/functional_study/review/other",
    "confidence": 0.0-1.0,
    "reasoning": "brief explanation"
}}""",

        "functional_study_extraction": """You are extracting functional study evidence from literature.

## Variant Information:
- Gene: {gene}
- Variant: {variant}

## Article:
Title: {title}
Abstract: {abstract}

## Task:
Determine if this functional study shows the variant AFFECTS function compared to wild-type.

### Evidence of Loss of Function (功能丧失 → PS3):
- Significant splicing alteration (>50% abnormal)
- Protein truncation/destabilization
- Loss of enzymatic activity
- Defective DNA binding
- Impaired protein-protein interaction
- Reduced protein expression/localization

### Evidence of Gain of Function (功能获得 → PS3):
- Increased enzymatic activity
- Enhanced protein function
- Constitutive activation
- Novel pathological function
- Increased expression/localization beyond wild-type

### Evidence of Normal Function (功能正常 → BS3):
- Normal splicing pattern
- Normal protein function
- Wild-type like activity
- No significant difference from control

### Uncertain (不确定 → 不使用):
- Insufficient data
- Conflicting results
- Technical limitations mentioned
- Results not statistically significant

## Output Format (JSON):
{{
    "functional_result": "lof/gof/normal/uncertain",
    "technique": "实验技术描述",
    "percentage_effect": 0.0-100.0 (if applicable),
    "confidence": "high/medium/low",
    "reasoning": "explanation"
}}""",

        "case_report_extraction": """You are extracting case report evidence from literature.

## Variant Information:
- Gene: {gene}
- Variant: {variant}
- Inheritance pattern: {inheritance_pattern}

## Article:
Title: {title}
Abstract: {abstract}

## Task:
Extract evidence for ACMG criteria:

### PS2 (de novo):
- Confirmed de novo = parental testing shows variant not in either parent
- Assumed de novo = variant present but parents not tested

### PM3 (biallelic/compound het):
- Variant in trans with another pathogenic variant
- Homozygous for this variant
- Hemizygous (X-linked)

### PP1 (co-segregation):
- Family members tested
- LOD score if reported
- Affected carriers vs unaffected non-carriers

### PS4 (case-control):
- Number of cases with variant
- Number of controls
- Odds ratio, p-value

## Output Format (JSON):
{{
    "case_count": number,
    "is_de_novo": true/false,
    "de_novo_confirmed": true/false (parental testing),
    "inheritance_pattern": "AD/AR/XLD/UNKNOWN",
    "is_compound_het": true/false,
    "trans_variant": "description or null",
    "segregation_data": {{
        "families_tested": number,
        "affected_carriers": number,
        "unaffected_carriers": number,
        "lod_score": number or null
    }} or null,
    "phenotype": "description or null",
    "hpo_terms": ["HP:xxxxx", ...] or [],
    "case_count": number,
    "control_count": number or null,
    "odds_ratio": number or null,
    "p_value": number or null
}}""",

        "ps2_assessment": """You are evaluating a genetic variant for ACMG/AMP classification evidence PS2 (de novo).

## Variant Information:
- Gene: {gene}
- Variant: {variant}
- Inheritance: {inheritance_pattern}

## Evidence Criteria:
PS2: De novo variant (confirmed by parental testing)
PM6: Assumed de novo (without parental confirmation)

## Input Data:
{input_data}

## Task:
Assess whether the variant meets PS2 or PM6 criteria based on literature evidence.

## Output Format (JSON):
{{
    "applicable": true/false,
    "strength": "STRONG/MODERATE/SUPPORTING",
    "confirmed_de_novo_count": number,
    "assumed_de_novo_count": number,
    "phenotype_consistency": "HIGH/LOW/UNKNOWN",
    "reason": "explanation"
}}""",

        "pm3_assessment": """You are evaluating a genetic variant for ACMG/AMP classification evidence PM3.

## Variant Information:
- Gene: {gene}
- Variant: {variant}
- Inheritance: {inheritance_pattern}

## Evidence Criteria:
PM3: Variant in trans with a known pathogenic variant (biallelic/compound het)

## Input Data:
{input_data}

## Task:
Assess whether the variant meets PM3 criteria based on literature evidence.

## Output Format (JSON):
{{
    "applicable": true/false,
    "strength": "VERY_STRONG/STRONG/MODERATE/SUPPORTING",
    "biallelic_count": number,
    "compound_het_count": number,
    "trans_variant": "description if available",
    "reason": "explanation"
}}""",

        "bp2_assessment": """You are evaluating a genetic variant for ACMG/AMP classification evidence BP2 (benign variant in trans with a known benign variant).

## Variant Information:
- Gene: {gene}
- Variant: {variant}
- Inheritance: {inheritance_pattern}

## Evidence Criteria:
BP2: Benign variant in trans with a known pathogenic variant in a recessive disease gene

## Input Data:
{input_data}

## Task:
Assess whether the variant meets BP2 criteria based on literature evidence.

## Output Format (JSON):
{{
    "applicable": true/false,
    "strength": "STRONG/MODERATE/SUPPORTING",
    "trans_with_benign_count": number,
    "reason": "explanation"
}}""",

        "segregation_assessment": """You are evaluating a genetic variant for ACMG/AMP classification evidence of co-segregation (PP1) or lack of co-segregation (BS4).

## Variant Information:
- Gene: {gene}
- Variant: {variant}
- Inheritance: {inheritance_pattern}

## Evidence Criteria:
PP1: Co-segregation with disease in family (LOD score if available)
BS4: Lack of co-segregation (affected non-carriers or unaffected carriers)

## Input Data:
{input_data}

## Task:
Assess co-segregation or lack thereof based on family study data.

## Output Format (JSON):
{{
    "applicable": true/false,
    "direction": "supporting_pathogenic/supporting_benign",
    "lod_score": number or null,
    "families_tested": number,
    "affected_carriers": number,
    "unaffected_carriers": number,
    "reason": "explanation"
}}""",

        "rag_llm_system": """You are an ACMG variant classification expert. Your task is to adjust ACMG evidence rules based on ClinGen gene-specific recommendations.

## Gene: {gene}
## Variant: {variant_description}

You will receive:
1. Current rule assessment results
2. ClinGen gene-specific rules (if available)
3. Variant-specific information

## Task:
Review and adjust evidence rules according to:
- ClinGen recommendations
- Gene-specific classification guidelines
- Published expert panel recommendations

## Important Guidelines:
- Be conservative when evidence is uncertain
- Follow ClinGen approved criteria when available
- Document any adjustments with reason

## Output Format:
{{
    "adjustments": [
        {{
            "rule_code": "PS2",
            "action": "apply/downgrade/upgrade/remove",
            "reason": "explanation based on ClinGen/evidence"
        }}
    ],
    "summary": "overall summary of changes"
}}""",
    }


def clear_cache():
    """Clear the prompts cache."""
    global _prompts_cache
    _prompts_cache.clear()
