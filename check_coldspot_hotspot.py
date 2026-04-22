#!/usr/bin/env python3

import pathlib
import csv
from typing import Optional, Tuple

from collections.abc import Callable

import pyensembl
import logging

from pybedtools import BedTool

from ensembl import ensembl
from information import Classification_Info, Info
from variant import TranscriptInfo, VariantInfo
from utils import create_bed_line
from cloud_api import ClinVarEutilsClient

logger = logging.getLogger("GenOtoScope_Classify.check_coldspot_hotspot")


def check_variant_intersection_with_bed(
    variant_hotspot_annotation_path: pathlib.Path,
    variant: VariantInfo,
    transcript: list[TranscriptInfo],
) -> bool:
    """
    Check hotspot annotation file for location of variant
    """
    # The reference transcript is only needed for the strand, therefore the specific transcript does not matter
    strand = get_variant_strand(transcript, variant)
    var_start = variant.genomic_start
    var_end = variant.genomic_end
    is_in_hotspot = check_intersection_with_bed_no_strand(
        variant, var_start, var_end, strand, variant_hotspot_annotation_path
    )
    return is_in_hotspot


def get_check_hotspot(
    class_info: Classification_Info,
) -> tuple[Callable, tuple[Info, ...]]:
    """
    Get function for hotspot annotation and needed classification_information objects
    """
    return (
        check_variant_intersection_with_bed,
        (
            class_info.VARIANT_HOTSPOT_ANNOTATION_PATH,
            class_info.VARIANT,
            class_info.TRANSCRIPT,
        ),
    )


def get_check_coldspot(
    class_info: Classification_Info,
) -> tuple[Callable, tuple[Info, ...]]:
    """
    Get function for hotspot annotation and needed classification_information objects
    """
    return (
        check_variant_intersection_with_bed,
        (
            class_info.VARIANT_COLDSPOT_ANNOTATION_PATH,
            class_info.VARIANT,
            class_info.TRANSCRIPT,
        ),
    )


def check_intersection_with_bed_no_strand(
    variant: VariantInfo,
    gen_start: int,
    gen_end: int,
    strand: str,
    path_bed: pathlib.Path,
) -> bool:
    """
    Check if variant overlaps region in given bed file without checking for strand
    """
    variant_interval = BedTool(
        create_bed_line(variant, gen_start, gen_end, strand),
        from_string=True,
    )[0]
    bed = BedTool(path_bed).sort()
    annotation_hits = bed.all_hits(variant_interval)
    if len(annotation_hits) > 0:
        return True
    return False


def get_variant_strand(transcripts: list[TranscriptInfo], variant: VariantInfo) -> str:
    """
    Try to get strand for gene
    """
    try:
        # Try getting strand from transcript
        ref_transcript = ensembl.transcript_by_id(transcripts[0].transcript_id)
        strand = ref_transcript.strand
        return strand
    except Exception:
        # In case getting strand from transcript fails, try getting strand from gene
        genes = ensembl.genes_by_name(variant.gene_name)
        # If there is only one gene matching, assume it's the correct one and return strand
        if len(genes) == 1:
            return genes[0].strand
        # If more than one gene was found, check if any of the genes is on the same chromosome as the variant
        genes_same_chr = []
        for gene in genes:
            if gene.contig in variant.chr:
                genes_same_chr.append(gene)
        if not len(genes_same_chr):
            raise ValueError(
                "The reconstruction of the varaint strand failed. None of the matching genes is located on the same strand as the variant."
            )
        elif len(genes_same_chr) == 1:
            return genes_same_chr[0].strand
        else:
            # If more than one gene is on the same chromosome, check the location, pick the closest gene
            closest_gene = find_gene_closest_to_variant(variant, genes_same_chr)
            return closest_gene.strand


def find_gene_closest_to_variant(
    variant: VariantInfo, genes: list[pyensembl.gene.Gene]
) -> pyensembl.gene.Gene:
    """
    Find the gene located clostest to the variant
    """
    distance_dict = {}
    for gene in genes:
        dist_to_start = min(
            abs(variant.genomic_start - gene.start),
            abs(variant.genomic_end - gene.start),
        )
        dist_to_end = min(
            abs(variant.genomic_start - gene.end),
            abs(variant.genomic_end - gene.end),
        )
        dist = min(dist_to_start, dist_to_end)
        distance_dict[gene.id] = dist
    gene_id_closest = min(distance_dict, key=distance_dict.get)
    for gene in genes:
        if gene.id == gene_id_closest:
            return gene


def extract_protein_position_from_hgvs(hgvs_p: Optional[str]) -> Optional[int]:
    """
    Extract protein position from HGVS p. notation.

    Args:
        hgvs_p: HGVS protein change (e.g., "p.Gly12Val")

    Returns:
        Protein position as integer or None
    """
    if not hgvs_p:
        return None

    import re

    # Match patterns like p.Gly12Val, p.Trp12*, p.Arg12_Arg13del
    pattern = r'p\.[A-Za-z]+(\d+)[A-Za-z*]*'
    match = re.search(pattern, hgvs_p)
    if match:
        return int(match.group(1))

    return None


def get_uniprot_domain_from_bed(
    gene_name: str,
    protein_position: int,
    uniprot_domain_path: pathlib.Path,
) -> Optional[dict]:
    """
    Get UniProt protein domain for a given protein position from BED file.

    Args:
        gene_name: Gene name (e.g., "BRCA1")
        protein_position: Protein position (amino acid number)
        uniprot_domain_path: Path to UniProt domain BED file

    Returns:
        Dict with domain info or None
    """
    if not uniprot_domain_path or not uniprot_domain_path.exists():
        return None

    with open(uniprot_domain_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row.get("gene_name") or row.get("gene") or row.get("Gene")
            if gene and gene.upper() == gene_name.upper():
                try:
                    domain_start = int(row["protein_start"])
                    domain_end = int(row["protein_end"])
                    if domain_start <= protein_position <= domain_end:
                        return {
                            "domain_name": row.get("domain_name", row.get("name", "unknown")),
                            "start": domain_start,
                            "end": domain_end,
                            "source": row.get("source", "UniProt"),
                        }
                except (KeyError, ValueError):
                    continue
    return None


def is_pathogenic_in_clinvar(variant: dict) -> bool:
    """
    Check if a ClinVar variant is pathogenic or likely pathogenic.

    Args:
        variant: ClinVar variant dict

    Returns:
        True if pathogenic/likely pathogenic
    """
    significance = variant.get("clinical_significance", "")
    if isinstance(significance, list) and len(significance) > 0:
        significance = significance[0]
    if isinstance(significance, dict):
        significance = significance.get("description", "")
    significance_lower = significance.lower() if significance else ""
    return "pathogenic" in significance_lower or "likely pathogenic" in significance_lower


def check_domain_and_clinvar_density(
    variant: VariantInfo,
    transcript: TranscriptInfo,
    uniprot_domain_path: pathlib.Path,
    clinvar_client: Optional[ClinVarEutilsClient] = None,
    domain_min_pathogenic: int = 3,
    window_min_pathogenic: int = 4,
    window_bp: int = 25,
) -> Tuple[bool, str]:
    """
    Check domain importance based on ClinVar variant density.

    This implements the PM1 enhanced logic:
    1. Check if variant is in predefined hotspot BED (done by caller)
    2. Query UniProt for protein domain
    3. Check if domain has >= domain_min_pathogenic pathogenic variants in ClinVar
    4. Check if ±window_bp genomic region has >= window_min_pathogenic pathogenic variants

    Args:
        variant: VariantInfo object
        transcript: TranscriptInfo object
        uniprot_domain_path: Path to UniProt domain BED file
        clinvar_client: ClinVar E-utilities API client
        domain_min_pathogenic: Min pathogenic in domain to apply PM1 (default 3)
        window_min_pathogenic: Min pathogenic in window to apply PM1 (default 4)
        window_bp: Genomic window size in bp (default 25)

    Returns:
        Tuple of (pm1_applies, reason)
    """
    # Step 1: Extract protein position from HGVS
    protein_position = None
    if transcript and hasattr(transcript, 'var_protein') and transcript.var_protein:
        protein_position = extract_protein_position_from_hgvs(transcript.var_protein)

    if not protein_position:
        return False, "Could not extract protein position from HGVS"

    # Step 2: Get UniProt domain
    gene_name = variant.gene_name
    if not gene_name:
        return False, "No gene name available"

    domain_info = get_uniprot_domain_from_bed(
        gene_name, protein_position, uniprot_domain_path
    )

    if not domain_info:
        return False, f"No UniProt domain found for {gene_name} at position {protein_position}"

    # Step 3: Check ClinVar density in domain
    if clinvar_client:
        domain_count = clinvar_client.count_pathogenic_in_gene(
            gene_name=gene_name,
            domain_start=domain_info["start"],
            domain_end=domain_info["end"],
        )
        if domain_count >= domain_min_pathogenic:
            return True, f"Domain {domain_info['domain_name']} has {domain_count} pathogenic variants (≥{domain_min_pathogenic})"

    # Step 4: Check ±window_bp genomic range density
    if clinvar_client:
        variants = clinvar_client.get_pathogenic_variants_in_region(
            gene_name=gene_name,
            chr=variant.chr,
            start=variant.genomic_start,
            end=variant.genomic_start,
            window_bp=window_bp,
        )
        window_count = len([v for v in variants if is_pathogenic_in_clinvar(v)])
        if window_count >= window_min_pathogenic:
            return True, f"{window_count} pathogenic variants within ±{window_bp}bp (≥{window_min_pathogenic})"

    return False, f"No significant ClinVar density found for domain {domain_info['domain_name']}"
