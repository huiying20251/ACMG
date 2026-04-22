#!/usr/bin/env python3
"""
预处理功能证据CSV文件为SQLite数据库。

功能：
1. 读取nm_to_chr.tsv建立转录本→染色体映射
2. 读取Export_20260402085651.csv
3. 按(转录本,cDNA位置)去重合并：
   - 同一变异多篇文献classification一致 → 保留该classification
   - 同一变异多篇文献classification不一致 → classification_review标记"冲突"
4. 输出SQLite数据库

Usage:
    python preprocess_functional_evidence.py [--input-csv EXPORT_CSV] [--nm2chr NM2CHR_TSV] [--output-db OUTPUT_DB]
"""

import argparse
import csv
import re
import sqlite3
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


@dataclass
class FunctionalEvidenceEntry:
    """功能证据条目"""
    identifier: str          # 完整Identifier: NM_005157.6(ABL1):c.80-4C>T
    transcript_id: str      # 转录本ID: NM_005157.6
    gene_name: str          # 基因名: ABL1
    chrom: Optional[str]    # 染色体号: 9
    position: int           # 基因组位置: 130854060
    cdna_pos: str           # cDNA位置: c.80-4
    ref_alt: str            # Ref/Alt: C/T
    molecular_consequence: str  # 分子后果: Splice site
    classification: str     # 功能分类: Functionally normal
    cross_assay_hits: int   # 交叉验证次数
    mave_technique: str    # MAVE技术: CRISPR-Based Genome Editing
    publication: str        # PMID: 33606977
    score: float            # 功能分值: -1.276
    phenotype: str          # 表型: ABL1-mediated cellular proliferation
    functional_description: str  # 功能描述
    clinvar_info: str       # ClinVar信息
    pop_frequency: str      # 群体频率

    # 去重合并后的字段
    classification_review: str = ""  # 冲突/一致
    unique_publications: str = ""     # 所有PMID，用|分隔


class NM2ChrMapper:
    """转录本ID到染色体号的映射"""

    def __init__(self, mapping_file: Path):
        self.mapping: Dict[str, str] = {}  # transcript_id -> chrom
        self._load_mapping(mapping_file)

    def _load_mapping(self, mapping_file: Path):
        """加载转录本→染色体映射"""
        with open(mapping_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # 清理转录本ID中的双引号和版本号
                refseq = row['refseq_mrna'].strip().strip('"')
                chrom = row['chromosome_name'].strip().strip('"')
                self.mapping[refseq] = chrom
        print(f"Loaded {len(self.mapping)} transcript→chromosome mappings")

    def get_chrom(self, transcript_id: str) -> Optional[str]:
        """获取转录本对应的染色体号"""
        # 去除版本号获取基准ID
        base_id = self._get_base_transcript_id(transcript_id)
        return self.mapping.get(base_id)

    @staticmethod
    def _get_base_transcript_id(transcript_id: str) -> str:
        """从完整转录本ID提取基准ID（去除版本号）"""
        # NM_005157.6 -> NM_005157
        match = re.match(r'(NM_\d+)', transcript_id)
        if match:
            return match.group(1)
        return transcript_id


class FunctionalEvidenceProcessor:
    """功能证据数据处理器"""

    def __init__(self, nm2chr_mapper: NM2ChrMapper):
        self.nm2chr = nm2chr_mapper
        # 按(转录本,cDNA位置)聚合条目
        self.aggregated: Dict[str, List[FunctionalEvidenceEntry]] = defaultdict(list)
        # 统计
        self.stats = {
            'total_rows': 0,
            'unique_variants': 0,
            'conflicting': 0,
            'consistent': 0,
            'missing_chrom': 0,
        }

    def parse_identifier(self, identifier: str) -> Tuple[str, str, str, int]:
        """
        解析Identifier，提取转录本ID、基因名、cDNA位置

        Args:
            identifier: NM_005157.6(ABL1):c.80-4C>T

        Returns:
            (transcript_id, gene_name, cdna_pos, position)
        """
        # 提取转录本ID和基因名: NM_005157.6(ABL1)
        match = re.match(r'(NM_\d+\.\d+)\(([^)]+)\):(.+)', identifier)
        if not match:
            raise ValueError(f"Cannot parse identifier: {identifier}")

        transcript_id = match.group(1)
        gene_name = match.group(2)
        remainder = match.group(3)  # c.80-4C>T

        # 提取cDNA位置和Ref/Alt
        # 格式可能是: c.80-4C>T 或 c.80-4C>T (p.Arg31Gln)
        cdna_match = re.match(r'(c\.\d+[-+]\d+|\bc\.\d+)([A-Z])>([A-Z])', remainder)
        if cdna_match:
            cdna_pos = cdna_match.group(1)  # c.80-4
            ref = cdna_match.group(2)  # C
            alt = cdna_match.group(3)  # T
        else:
            # 尝试其他格式
            parts = remainder.split()
            cdna_pos = parts[0]

        # 提取基因组位置（从Identifier中获取）
        return transcript_id, gene_name, cdna_pos

    def parse_row(self, row: dict, line_num: int) -> Optional[FunctionalEvidenceEntry]:
        """解析CSV行"""
        try:
            identifier = row['Identifier'].strip().strip('"')

            # 解析Identifier
            transcript_id, gene_name, cdna_pos = self.parse_identifier(identifier)

            # 获取染色体号
            chrom = self.nm2chr.get_chrom(transcript_id)

            # 解析位置
            position = int(row['Position'].strip().strip('"'))

            # 解析Ref/Alt
            ref_alt = row['Ref/Alt'].strip().strip('"')

            # 解析Score
            score_str = row['Score'].strip().strip('"')
            score = float(score_str) if score_str else 0.0

            # 解析交叉验证次数
            cross_assay_hits = int(row['Cross-assay hits'].strip().strip('"')) if row.get('Cross-assay hits') else 1

            entry = FunctionalEvidenceEntry(
                identifier=identifier,
                transcript_id=transcript_id,
                gene_name=gene_name,
                chrom=chrom,
                position=position,
                cdna_pos=cdna_pos,
                ref_alt=ref_alt,
                molecular_consequence=row['Molecular consequence'].strip().strip('"'),
                classification=row['Functional classification'].strip().strip('"'),
                cross_assay_hits=cross_assay_hits,
                mave_technique=row['MAVE technique'].strip().strip('"'),
                publication=row['Publication'].strip().strip('"'),
                score=score,
                phenotype=row['Phenotype'].strip().strip('"'),
                functional_description=row['Functional description'].strip().strip('"'),
                clinvar_info=row['ClinVar information'].strip().strip('"'),
                pop_frequency=row['Population frequency'].strip().strip('"'),
            )

            if chrom is None:
                self.stats['missing_chrom'] += 1

            return entry

        except (KeyError, ValueError) as e:
            print(f"Warning: Error parsing line {line_num}: {e}", file=sys.stderr)
            return None

    def create_variant_key(self, entry: FunctionalEvidenceEntry) -> str:
        """创建变异位点的唯一键"""
        # 使用(转录本, cDNA位置)作为key，忽略Ref/Alt（同一位置可能不同变异）
        return f"{entry.transcript_id}:{entry.cdna_pos}"

    def aggregate_entries(self, entries: List[FunctionalEvidenceEntry]) -> FunctionalEvidenceEntry:
        """
        聚合同一变异位点的多个条目

        去重合并逻辑：
        1. 如果多篇文献的classification都一致 → classification_review = ""
        2. 如果多篇文献的classification不一致 → classification_review = "冲突"
        """
        if not entries:
            raise ValueError("No entries to aggregate")

        if len(entries) == 1:
            # 只有一个条目，无需合并
            master = entries[0]
            master.classification_review = ""
            master.unique_publications = master.publication
            return master

        # 收集所有不同的classification和publications
        classifications: Set[str] = set()
        publications: List[str] = []
        master = entries[0]

        for entry in entries:
            classifications.add(entry.classification)
            if entry.publication and entry.publication not in publications:
                publications.append(entry.publication)

        # 判断是否有冲突
        if len(classifications) > 1:
            classification_review = "冲突"
            self.stats['conflicting'] += 1
        else:
            classification_review = ""
            self.stats['consistent'] += 1

        # 使用第一个条目作为主条目，合并信息
        master.classification_review = classification_review
        master.unique_publications = "|".join(publications)
        master.cross_assay_hits = len(entries)

        # 如果有冲突，记录所有不同的classification
        if classification_review == "冲突":
            master.functional_description = f"CONFLICT: {', '.join(classifications)} | {master.functional_description}"

        return master

    def process_csv(self, csv_path: Path) -> List[FunctionalEvidenceEntry]:
        """处理CSV文件"""
        print(f"Processing CSV file: {csv_path}")

        # 移除BOM并读取文件
        with open(csv_path, 'r', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)

            for line_num, row in enumerate(reader, start=2):  # start=2因为第一行是表头
                self.stats['total_rows'] += 1

                entry = self.parse_row(row, line_num)
                if entry is None:
                    continue

                key = self.create_variant_key(entry)
                self.aggregated[key].append(entry)

                if line_num % 100000 == 0:
                    print(f"  Processed {line_num:,} rows...")

        print(f"  Total rows: {self.stats['total_rows']:,}")
        print(f"  Unique variants: {len(self.aggregated):,}")

        # 合并每个variant的条目
        merged_entries: List[FunctionalEvidenceEntry] = []
        for key, entries in self.aggregated.items():
            merged = self.aggregate_entries(entries)
            merged_entries.append(merged)

        print(f"  Consistent: {self.stats['consistent']:,}")
        print(f"  Conflicting: {self.stats['conflicting']:,}")
        print(f"  Missing chromosome: {self.stats['missing_chrom']:,}")

        return merged_entries


def create_sqlite_db(entries: List[FunctionalEvidenceEntry], db_path: Path):
    """创建SQLite数据库"""
    print(f"Creating SQLite database: {db_path}")

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # 创建表
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS functional_evidence (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifier TEXT NOT NULL,
            transcript_id TEXT NOT NULL,
            gene_name TEXT NOT NULL,
            chrom TEXT,
            position INTEGER NOT NULL,
            cdna_pos TEXT NOT NULL,
            ref_alt TEXT,
            molecular_consequence TEXT,
            classification TEXT NOT NULL,
            classification_review TEXT DEFAULT '',
            cross_assay_hits INTEGER,
            mave_technique TEXT,
            publication TEXT,
            unique_publications TEXT,
            score REAL,
            phenotype TEXT,
            functional_description TEXT,
            clinvar_info TEXT,
            pop_frequency TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')

    # 创建索引
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_transcript_cdna ON functional_evidence(transcript_id, cdna_pos)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_gene_name ON functional_evidence(gene_name)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_chrom_position ON functional_evidence(chrom, position)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_classification ON functional_evidence(classification)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_classification_review ON functional_evidence(classification_review)')

    # 插入数据
    for entry in entries:
        cursor.execute('''
            INSERT INTO functional_evidence (
                identifier, transcript_id, gene_name, chrom, position, cdna_pos,
                ref_alt, molecular_consequence, classification, classification_review,
                cross_assay_hits, mave_technique, publication, unique_publications,
                score, phenotype, functional_description, clinvar_info, pop_frequency
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            entry.identifier,
            entry.transcript_id,
            entry.gene_name,
            entry.chrom,
            entry.position,
            entry.cdna_pos,
            entry.ref_alt,
            entry.molecular_consequence,
            entry.classification,
            entry.classification_review,
            entry.cross_assay_hits,
            entry.mave_technique,
            entry.publication,
            entry.unique_publications,
            entry.score,
            entry.phenotype,
            entry.functional_description,
            entry.clinvar_info,
            entry.pop_frequency,
        ))

    conn.commit()

    # 验证数据
    cursor.execute('SELECT COUNT(*) FROM functional_evidence')
    count = cursor.fetchone()[0]
    print(f"  Inserted {count:,} records")

    # 显示分类统计
    cursor.execute('''
        SELECT classification, classification_review, COUNT(*) as cnt
        FROM functional_evidence
        GROUP BY classification, classification_review
        ORDER BY cnt DESC
    ''')
    print("\n  Classification distribution:")
    for row in cursor.fetchall():
        review = f" [{row[1]}]" if row[1] else ""
        print(f"    {row[0]}{review}: {row[2]:,}")

    # 显示染色体缺失统计
    cursor.execute('SELECT COUNT(*) FROM functional_evidence WHERE chrom IS NULL')
    missing_chrom = cursor.fetchone()[0]
    if missing_chrom > 0:
        print(f"\n  WARNING: {missing_chrom:,} entries missing chromosome mapping")

    conn.close()
    print("Done!")


def main():
    parser = argparse.ArgumentParser(
        description='预处庺功能证据CSV为SQLite数据库'
    )
    parser.add_argument(
        '--input-csv',
        type=Path,
        default=Path('/home/huiying/variant_classification/Export_20260402085651.csv'),
        help='输入CSV文件路径'
    )
    parser.add_argument(
        '--nm2chr',
        type=Path,
        default=Path('/home/huiying/variant_classification/nm_to_chr.tsv'),
        help='转录本→染色体映射文件'
    )
    parser.add_argument(
        '--output-db',
        type=Path,
        default=Path('/home/huiying/variant_classification/functional_evidence.db'),
        help='输出SQLite数据库路径'
    )

    args = parser.parse_args()

    # 检查输入文件
    if not args.input_csv.exists():
        print(f"Error: Input CSV not found: {args.input_csv}", file=sys.stderr)
        sys.exit(1)
    if not args.nm2chr.exists():
        print(f"Error: NM2Chr mapping file not found: {args.nm2chr}", file=sys.stderr)
        sys.exit(1)

    # 创建mapper
    nm2chr_mapper = NM2ChrMapper(args.nm2chr)

    # 处理CSV
    processor = FunctionalEvidenceProcessor(nm2chr_mapper)
    entries = processor.process_csv(args.input_csv)

    # 创建SQLite数据库
    create_sqlite_db(entries, args.output_db)

    return 0


if __name__ == '__main__':
    sys.exit(main())