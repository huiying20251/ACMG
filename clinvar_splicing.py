#!/usr/bin/env python3

import logging
import pathlib
from collections.abc import Iterable
from typing import Optional
import pandas as pd

import pyensembl

from variant import VariantInfo, TranscriptInfo
from var_type import VARTYPE_GROUPS
from ps1_splicing_utils import SplicingVariantInfo, SplicingRegion
from clinvar_utils import (
    ClinVar,
    ClinVar_Type,
    create_ClinVar,
    get_affected_transcript,
)
from genotoscope_exon_skipping import (
    parse_variant_intron_pos,
    find_exon_by_ref_pos,
)
from custom_exceptions import No_transcript_with_var_type_found
from cloud_api import ClinVarEutilsClient
from clinvar_utils import convert_vcf_gen_to_df

logger = logging.getLogger("GenOtoScope_Classify.clinvar.splicing")


def check_clinvar_splicing(
    variant: VariantInfo,
    transcripts: Iterable[TranscriptInfo],
    path_clinvar: Optional[pathlib.Path] = None,
    splicing_info: Optional[SplicingVariantInfo] = None,
) -> tuple[ClinVar, ClinVar]:
    """
    Check ClinVar for entries supporting pathogenicity of splice site.

    Strategy: API-first with fallback to local VCF file.

    Args:
        variant: The variant to check
        transcripts: List of transcript infos
        path_clinvar: Path to local ClinVar VCF file (optional, fallback)
        splicing_info: SplicingVariantInfo from ps1_splicing_utils for determining search window

    Returns:
        Tuple of (ClinVar_same_nucleotide, ClinVar_same_splice_site)
    """
    if len(variant.var_obs) != 1 or len(variant.var_ref) != 1:
        logger.warning(
            "Variant is not a SNV. PS1/PM5 currently not implemented for delins."
        )
        ClinVar_same_nucleotide = create_ClinVar(
            pd.DataFrame(), ClinVar_Type.SAME_NUCLEOTIDE
        )
        ClinVar_same_splice_site = create_ClinVar(
            pd.DataFrame(), ClinVar_Type.SAME_SPLICE_SITE
        )
        return (ClinVar_same_nucleotide, ClinVar_same_splice_site)

    # Check if we have splicing_info to determine search window
    # If not provided, fall back to old behavior (VARTYPE_GROUPS.INTRONIC)
    if splicing_info is None or splicing_info.region == SplicingRegion.NOT_SPLICING:
        # Branch 2: missense + splice site variant
        try:
            affected_transcript, ref_transcript = get_affected_transcript(
                transcripts, VARTYPE_GROUPS.INTRONIC
            )
        except No_transcript_with_var_type_found:
            logger.warning("No transcript with variant type splicing found.")
            ClinVar_same_nucleotide = create_ClinVar(
                pd.DataFrame(), ClinVar_Type.SAME_NUCLEOTIDE
            )
            ClinVar_same_splice_site = create_ClinVar(
                pd.DataFrame(), ClinVar_Type.SAME_SPLICE_SITE
            )
            return (ClinVar_same_nucleotide, ClinVar_same_splice_site)

        clinvar_client = ClinVarEutilsClient()

        # ========== Step 1: Try ClinVar API first ==========
        clinvar_same_pos = []
        clinvar_splice_site = []
        api_success = False

        try:
            # Original logic: exact position + splice site ±5bp
            clinvar_same_pos = clinvar_client.search_splicing_by_genomic_position(
                gene_name=variant.gene_name,
                chromosome=variant.chr,
                genomic_start=variant.genomic_start,
                genomic_end=variant.genomic_end,
                window_bp=0,
            )

            (start_splice_site, end_splice_site) = find_corresponding_splice_site(
                affected_transcript, ref_transcript, variant
            )

            clinvar_splice_site = clinvar_client.search_splicing_by_genomic_position(
                gene_name=variant.gene_name,
                chromosome=variant.chr,
                genomic_start=start_splice_site,
                genomic_end=end_splice_site,
                window_bp=5,
            )

            if clinvar_same_pos or clinvar_splice_site:
                api_success = True
                logger.info(
                    f"ClinVar API returned {len(clinvar_same_pos)} same-pos + "
                    f"{len(clinvar_splice_site)} splice-site entries"
                )

        except Exception as e:
            logger.warning(f"ClinVar API query failed: {e}")

        # ========== Step 2: Fallback to local VCF if API failed or returned empty ==========
        if not api_success or (not clinvar_same_pos and not clinvar_splice_site and path_clinvar and path_clinvar.exists()):
            if path_clinvar and path_clinvar.exists():
                logger.info(f"Falling back to local ClinVar VCF: {path_clinvar}")
                try:
                    clinvar_same_pos, clinvar_splice_site = _search_splicing_vcf(
                        path_clinvar, variant, affected_transcript, ref_transcript
                    )
                except Exception as e:
                    logger.warning(f"Local ClinVar VCF fallback failed: {e}")
                    clinvar_same_pos = []
                    clinvar_splice_site = []
            else:
                logger.info("No local ClinVar VCF file provided for fallback")

        # Process results
        clinvar_same_pos_df = process_clinvar_splicing_api_results(
            clinvar_same_pos, variant, exclude_same_nucleotide=True
        )
        ClinVar_same_pos = create_ClinVar(
            clinvar_same_pos_df, ClinVar_Type.SAME_NUCLEOTIDE
        )

        clinvar_splice_site_df = process_clinvar_splicing_api_results(
            clinvar_splice_site, variant, exclude_same_nucleotide=True
        )
        ClinVar_splice_site = create_ClinVar(
            clinvar_splice_site_df, ClinVar_Type.SAME_SPLICE_SITE
        )

        return (ClinVar_same_pos, ClinVar_splice_site)

    # ========== Branch 1 & 2: New logic based on splicing_info ==========
    clinvar_client = ClinVarEutilsClient()

    # 计算相对偏移量
    base_to_start_offset = splicing_info.search_window_start - splicing_info.base_position
    base_to_end_offset = splicing_info.search_window_end - splicing_info.base_position

    # 计算基因组搜索范围
    search_start = variant.genomic_start + base_to_start_offset
    search_end = variant.genomic_start + base_to_end_offset

    logger.info(f"PS1_splicing search: {splicing_info.hgvs_string}, "
                f"region={splicing_info.region.value}, "
                f"search window: {search_start}-{search_end}")

    # ========== Step 1: Try ClinVar API first ==========
    clinvar_region = []
    api_success = False

    try:
        # Search for all variants in the splicing region
        clinvar_region = clinvar_client.search_splicing_by_genomic_position(
            gene_name=variant.gene_name,
            chromosome=variant.chr,
            genomic_start=search_start,
            genomic_end=search_end,
            window_bp=0,
        )

        if clinvar_region:
            api_success = True
            logger.info(f"ClinVar API returned {len(clinvar_region)} entries")
        else:
            logger.info(f"ClinVar API returned no entries for splicing region")

    except Exception as e:
        logger.warning(f"ClinVar API query failed: {e}")

    # ========== Step 2: Fallback to local VCF if API failed or returned empty ==========
    if not api_success or (not clinvar_region and path_clinvar and path_clinvar.exists()):
        if path_clinvar and path_clinvar.exists():
            logger.info(f"Falling back to local ClinVar VCF: {path_clinvar}")
            try:
                clinvar_region = _search_splicing_region_vcf(
                    path_clinvar, variant.chr, search_start, search_end
                )
            except Exception as e:
                logger.warning(f"Local ClinVar VCF fallback failed: {e}")
                clinvar_region = []
        else:
            logger.info("No local ClinVar VCF file provided for fallback")

    # ========== Step 3: Process results ==========
    clinvar_region_df = process_clinvar_splicing_api_results(
        clinvar_region, variant, exclude_same_nucleotide=True
    )

    # Split into SAME_NUCLEOTIDE (same position, different alt) and SAME_SPLICE_SITE (different position, same motif)
    df_same_pos = pd.DataFrame()
    df_same_motif = pd.DataFrame()

    if not clinvar_region_df.empty and "is_same_position" in clinvar_region_df.columns:
        df_same_pos = clinvar_region_df[clinvar_region_df["is_same_position"] == True]
        df_same_motif = clinvar_region_df[clinvar_region_df["is_same_position"] == False]

    ClinVar_same_pos = create_ClinVar(
        df_same_pos, ClinVar_Type.SAME_NUCLEOTIDE
    )
    ClinVar_splice_site = create_ClinVar(
        df_same_motif, ClinVar_Type.SAME_SPLICE_SITE
    )

    return (ClinVar_same_pos, ClinVar_splice_site)


def _search_splicing_vcf(
    path_clinvar: pathlib.Path,
    variant: VariantInfo,
    affected_transcript: TranscriptInfo,
    ref_transcript: pyensembl.transcript.Transcript,
) -> tuple[list, list]:
    """
    Search local ClinVar VCF for splicing variants.

    Args:
        path_clinvar: Path to ClinVar VCF file
        variant: The variant to check
        affected_transcript: Affected transcript info
        ref_transcript: Reference transcript

    Returns:
        Tuple of (clinvar_same_pos list, clinvar_splice_site list)
    """
    from cyvcf2 import VCF

    clinvar_vcf = VCF(str(path_clinvar))

    # Search for same position
    clinvar_same_pos = []
    try:
        for rec in clinvar_vcf(f"{variant.chr}:{variant.genomic_start}-{variant.genomic_start}"):
            rec_dict = _parse_vcf_record(rec, variant)
            if rec_dict:
                clinvar_same_pos.append(rec_dict)
    except Exception as e:
        logger.warning(f"Error searching VCF for same position: {e}")

    # Search for splice site region
    clinvar_splice_site = []
    try:
        (start_splice_site, end_splice_site) = find_corresponding_splice_site(
            affected_transcript, ref_transcript, variant
        )
        window_start = start_splice_site - 5
        window_end = end_splice_site + 5

        for rec in clinvar_vcf(f"{variant.chr}:{window_start}-{window_end}"):
            rec_dict = _parse_vcf_record(rec, variant)
            if rec_dict:
                # Check if this is at a different position
                if rec_dict.get("pos") != str(variant.genomic_start):
                    clinvar_splice_site.append(rec_dict)
    except Exception as e:
        logger.warning(f"Error searching VCF for splice site: {e}")

    return clinvar_same_pos, clinvar_splice_site


def _search_splicing_region_vcf(
    path_clinvar: pathlib.Path,
    chromosome: str,
    start: int,
    end: int,
) -> list:
    """
    Search local ClinVar VCF for a genomic region.

    Args:
        path_clinvar: Path to ClinVar VCF file
        chromosome: Chromosome
        start: Start position
        end: End position

    Returns:
        List of ClinVar records in the region
    """
    from cyvcf2 import VCF

    clinvar_vcf = VCF(str(path_clinvar))
    records = []

    try:
        for rec in clinvar_vcf(f"{chromosome}:{start}-{end}"):
            rec_dict = _parse_vcf_record(rec, None)
            if rec_dict:
                records.append(rec_dict)
    except Exception as e:
        logger.warning(f"Error searching VCF region {chromosome}:{start}-{end}: {e}")

    return records


def _parse_vcf_record(rec, variant: Optional[VariantInfo]) -> Optional[dict]:
    """
    Parse a VCF record into API-style format.

    Args:
        rec: cyvcf2 VCF record
        variant: The variant to check (for filtering)

    Returns:
        Dictionary in API-style format or None
    """
    try:
        # Get INFO field
        info = rec.INFO

        # Extract significance
        clnsig = info.get("CLNSIG", "unknown")
        if isinstance(clnsig, list):
            clnsig = clnsig[0] if clnsig else "unknown"

        # Extract gene
        gene_info = info.get("GENEINFO", "")
        gene_symbol = gene_info.split(":")[0] if gene_info else ""

        record = {
            "gene": gene_symbol,
            "significance": str(clnsig),
            "variation_id": rec.ID if rec.ID else "",
            "pos": str(rec.POS),
            "alt": rec.ALT[0] if rec.ALT else "",
        }

        # Get variant type if available
        if hasattr(rec, "var_type"):
            record["variant_type"] = str(rec.var_type)

        return record

    except Exception as e:
        logger.debug(f"Error parsing VCF record: {e}")
        return None


def process_clinvar_splicing_api_results(
    clinvar_entries: list[dict],
    variant: VariantInfo,
    exclude_same_nucleotide: bool = True,
) -> pd.DataFrame:
    """
    Process ClinVar API results for splicing variants.

    Args:
        clinvar_entries: List of ClinVar variant records from API
        variant: The input variant
        exclude_same_nucleotide: Whether to exclude the same nucleotide change

    Returns:
        DataFrame of ClinVar records with 'is_same_position' column
        - is_same_position=True: 相同位点不同碱基
        - is_same_position=False: 同基序不同位点
    """
    if not clinvar_entries:
        return pd.DataFrame()

    records = []
    for entry in clinvar_entries:
        record = parse_clinvar_splicing_api_entry(entry, variant)
        if record:
            records.append(record)

    if not records:
        return pd.DataFrame()

    df = pd.DataFrame(records)

    # 重命名列以匹配 create_ClinVar 的期望
    # API 返回: significance, variation_id
    # create_ClinVar 期望: CLNSIG, id
    if "significance" in df.columns:
        df = df.rename(columns={"significance": "CLNSIG"})
    if "variation_id" in df.columns:
        df = df.rename(columns={"variation_id": "id"})

    # 标记是否为同一位置（不同碱基）
    if "pos" in df.columns:
        df["is_same_position"] = df["pos"] == str(variant.genomic_start)

    # Exclude the same nucleotide change (same position AND same alt base)
    if exclude_same_nucleotide and not df.empty:
        if "pos" in df.columns and "alt" in df.columns:
            df = df[
                ~(
                    (df["alt"] == variant.var_obs)
                    & (df["pos"] == str(variant.genomic_start))
                )
            ]

    return df


def parse_clinvar_splicing_api_entry(entry: dict, variant: VariantInfo) -> dict:
    """
    Parse a ClinVar API entry for splicing variants.

    Args:
        entry: ClinVar API result entry
        variant: The input variant

    Returns:
        Parsed record dictionary or None
    """
    try:
        # Extract variant info from ClinVar API response
        rcvaccession = entry.get("rcvaccession", [])
        if not rcvaccession:
            return None

        # Get the first RCV accession
        rcv = rcvaccession[0] if isinstance(rcvaccession, list) else rcvaccession

        # Extract clinical significance
        clin_sig = rcv.get("clinical_significance", {})
        if isinstance(clin_sig, dict):
            significance = clin_sig.get("description", "unknown")
        else:
            significance = str(clin_sig)

        # Extract gene info
        gene_info = entry.get("gene", {})
        if isinstance(gene_info, dict):
            gene_symbol = gene_info.get("symbol", variant.gene_name)
        else:
            gene_symbol = variant.gene_name

        # Extract variant coordinates if available
        record = {
            "gene": gene_symbol,
            "significance": significance,
            "variation_id": entry.get("variation_id", ""),
        }

        # Position and allele info if available
        if "position" in entry:
            record["pos"] = str(entry["position"])
        if "allele" in entry:
            record["alt"] = entry["allele"]

        return record

    except Exception as e:
        logger.debug(f"Error parsing ClinVar splicing API entry: {e}")
        return None


def find_corresponding_splice_site(
    transcript: TranscriptInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    variant: VariantInfo,
) -> tuple[int, int]:
    """
    Reconstruct splice site
    Splice site is defined as +/- 1,2 as only for these locations varinat is clearly defines as a splice variant
    """
    if "+" in str(transcript.var_hgvs) or "-" in str(transcript.var_hgvs):
        splice_site_start, splice_site_stop = get_splice_site_for_intronic_variant(
            variant, transcript, ref_transcript
        )
    else:
        splice_site_start, splice_site_stop = get_splice_site_for_exonic_variant(
            variant, ref_transcript
        )
    return (splice_site_start, splice_site_stop)


def get_splice_site_for_intronic_variant(
    variant: VariantInfo,
    transcript: TranscriptInfo,
    ref_transcript: pyensembl.transcript.Transcript,
) -> tuple[int, int]:
    """
    Get splice site for intronic variant
    """
    (
        split_symbol,
        distance_to_splice_site,
        direction_to_splice_site,
    ) = parse_variant_intron_pos(transcript.var_hgvs)
    distance_to_splice_site_start = distance_to_splice_site - 1
    distance_to_splice_site_stop = distance_to_splice_site - 2
    if ref_transcript.strand == "-":
        # Reverse direction_to_splice_site, to point into genomic direction of splice site
        direction_to_splice_site = direction_to_splice_site * -1
    splice_site_start = (
        variant.genomic_start + direction_to_splice_site * distance_to_splice_site_start
    )
    splice_site_stop = (
        variant.genomic_start + direction_to_splice_site * distance_to_splice_site_stop
    )
    if splice_site_start > splice_site_stop:
        return splice_site_stop, splice_site_start
    else:
        return splice_site_start, splice_site_stop


def get_splice_site_for_exonic_variant(
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
) -> tuple[int, int]:
    """
    Get splice site for exonic variant
    """
    exon_index, pos_in_exon = find_exon_by_ref_pos(
        ref_transcript, variant.genomic_start, is_genomic=True
    )
    affected_exon = ref_transcript.exons[exon_index]
    if abs(variant.genomic_start - affected_exon.start) < abs(
        variant.genomic_start - affected_exon.end
    ):
        splice_site_start = affected_exon.start - 2
        splice_site_end = affected_exon.start - 1
    else:
        splice_site_start = affected_exon.end + 1
        splice_site_end = affected_exon.end + 2
    return splice_site_start, splice_site_end
