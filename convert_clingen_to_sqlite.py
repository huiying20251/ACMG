#!/usr/bin/env python3
"""
Convert ClinGen Excel files to SQLite database for ACMG rule configuration.

Files:
1. ClinGen_rules_frequency_cutoffs.xlsx - BA1/BS1/PM2 frequency cutoffs
2. ClinGen_rules_pm1_details.xlsx - PM1 hotspot domains
3. ClinGen_rules_rules.xlsx - Gene-specific ACMG rules

Usage:
    python convert_clingen_to_sqlite.py [--output-db clingen_rules.db]
"""

import argparse
import csv
import sqlite3
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Dict, Any, Optional


def parse_xlsx(zip_file) -> List[List[str]]:
    """Parse xlsx file and return rows as list of lists."""
    shared_strings = []
    if 'xl/sharedStrings.xml' in zip_file.namelist():
        with zip_file.open('xl/sharedStrings.xml') as f:
            tree = ET.parse(f)
            root = tree.getroot()
            for si in root.iter('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}si'):
                text = ''
                for t in si.iter('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t'):
                    if t.text:
                        text += t.text
                shared_strings.append(text)

    with zip_file.open('xl/worksheets/sheet1.xml') as f:
        tree = ET.parse(f)
        root = tree.getroot()

        rows = []
        for row in root.iter('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}row'):
            row_data = []
            for cell in row.iter('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}c'):
                cell_type = cell.get('t', '')
                v = cell.find('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}v')
                if v is not None and v.text is not None:
                    if cell_type == 's':
                        row_data.append(shared_strings[int(v.text)])
                    else:
                        row_data.append(v.text)
                else:
                    row_data.append('')
            rows.append(row_data)
    return rows


def create_frequency_cutoffs_table(cursor: sqlite3.Cursor):
    """Create table for frequency cutoffs."""
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS frequency_cutoffs (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene TEXT NOT NULL,
            rule_code TEXT NOT NULL,
            strength TEXT NOT NULL,
            cutoff_operator TEXT NOT NULL,
            cutoff_value REAL NOT NULL,
            frequency_type TEXT,
            database TEXT,
            population_group TEXT,
            source TEXT,
            notes TEXT,
            UNIQUE(gene, rule_code, strength)
        )
    ''')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_freq_gene ON frequency_cutoffs(gene)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_freq_rule ON frequency_cutoffs(rule_code)')


def create_pm1_details_table(cursor: sqlite3.Cursor):
    """Create table for PM1 domain details."""
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS pm1_details (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene TEXT NOT NULL,
            strength TEXT NOT NULL,
            domain TEXT,
            cdna_positions TEXT,
            protein_positions TEXT,
            exon_numbers TEXT,
            transcript TEXT,
            notes TEXT,
            no_applicable INTEGER DEFAULT 0,
            same_as_original INTEGER DEFAULT 0,
            pm1_table_ref TEXT,
            application_criteria TEXT,
            UNIQUE(gene, strength, protein_positions)
        )
    ''')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_pm1_gene ON pm1_details(gene)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_pm1_protein ON pm1_details(gene, protein_positions)')


def create_clingen_rules_table(cursor: sqlite3.Cursor):
    """Create table for ClinGen gene-specific rules."""
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS clingen_rules (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene TEXT NOT NULL,
            hgnc_id TEXT,
            rule_code TEXT NOT NULL,
            strength TEXT NOT NULL,
            strength_order INTEGER,
            applicability INTEGER DEFAULT 0,
            description TEXT,
            application_criteria TEXT,
            source TEXT,
            version TEXT,
            last_updated TEXT,
            UNIQUE(gene, rule_code, strength)
        )
    ''')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_rules_gene ON clingen_rules(gene)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_rules_rule_code ON clingen_rules(rule_code)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_rules_gene_code ON clingen_rules(gene, rule_code)')


def import_frequency_cutoffs(cursor: sqlite3.Cursor, xlsx_path: Path):
    """Import frequency cutoffs from xlsx file."""
    print(f"Importing frequency cutoffs from {xlsx_path}...")

    with zipfile.ZipFile(xlsx_path, 'r') as z:
        rows = parse_xlsx(z)

    header = rows[0]
    data_rows = rows[1:]

    for row in data_rows:
        if len(row) < 5:
            continue

        # Parse row into dict
        row_dict = dict(zip(header, row + [''] * (len(header) - len(row))))

        # Parse cutoff value (handle scientific notation)
        try:
            cutoff_value = float(row_dict.get('cutoff_value', 0))
        except ValueError:
            cutoff_value = 0.0

        cursor.execute('''
            INSERT OR REPLACE INTO frequency_cutoffs
            (gene, rule_code, strength, cutoff_operator, cutoff_value,
             frequency_type, database, population_group, source, notes)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            row_dict.get('gene', '').strip(),
            row_dict.get('rule_code', '').strip(),
            row_dict.get('strength', '').strip(),
            row_dict.get('cutoff_operator', '').strip(),
            cutoff_value,
            row_dict.get('frequency_type', '').strip(),
            row_dict.get('database', '').strip(),
            row_dict.get('population_group', '').strip(),
            row_dict.get('source', '').strip(),
            row_dict.get('notes', '').strip(),
        ))

    print(f"  Imported {len(data_rows)} frequency cutoff records")


def import_pm1_details(cursor: sqlite3.Cursor, xlsx_path: Path):
    """Import PM1 details from xlsx file."""
    print(f"Importing PM1 details from {xlsx_path}...")

    with zipfile.ZipFile(xlsx_path, 'r') as z:
        rows = parse_xlsx(z)

    header = rows[0]
    data_rows = rows[1:]

    for row in data_rows:
        if len(row) < 5:
            continue

        row_dict = dict(zip(header, row + [''] * (len(header) - len(row))))

        cursor.execute('''
            INSERT OR REPLACE INTO pm1_details
            (gene, strength, domain, cdna_positions, protein_positions,
             exon_numbers, transcript, notes, no_applicable, same_as_original,
             pm1_table_ref, application_criteria)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            row_dict.get('gene', '').strip(),
            row_dict.get('strength', '').strip(),
            row_dict.get('domain', '').strip(),
            row_dict.get('cdna_positions', '').strip(),
            row_dict.get('protein_positions', '').strip(),
            row_dict.get('exon_numbers', '').strip(),
            row_dict.get('transcript', '').strip(),
            row_dict.get('notes', '').strip(),
            int(row_dict.get('no_applicable', 0)) if row_dict.get('no_applicable', '0').isdigit() else 0,
            int(row_dict.get('same_as_original', 0)) if row_dict.get('same_as_original', '0').isdigit() else 0,
            row_dict.get('pm1_table_ref', '').strip(),
            row_dict.get('application_criteria', '').strip(),
        ))

    print(f"  Imported {len(data_rows)} PM1 detail records")


def import_clingen_rules(cursor: sqlite3.Cursor, xlsx_path: Path):
    """Import ClinGen rules from xlsx file."""
    print(f"Importing ClinGen rules from {xlsx_path}...")

    with zipfile.ZipFile(xlsx_path, 'r') as z:
        rows = parse_xlsx(z)

    header = rows[0]
    data_rows = rows[1:]

    for row in data_rows:
        if len(row) < 5:
            continue

        row_dict = dict(zip(header, row + [''] * (len(header) - len(row))))

        # Parse strength_order
        try:
            strength_order = int(row_dict.get('strength_order', 0))
        except ValueError:
            strength_order = 0

        cursor.execute('''
            INSERT OR REPLACE INTO clingen_rules
            (gene, hgnc_id, rule_code, strength, strength_order,
             applicability, description, application_criteria, source, version, last_updated)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            row_dict.get('gene', '').strip(),
            row_dict.get('hgnc_id', '').strip(),
            row_dict.get('rule_code', '').strip(),
            row_dict.get('strength', '').strip(),
            strength_order,
            int(row_dict.get('applicability', 0)) if str(row_dict.get('applicability', 0)).isdigit() else 0,
            row_dict.get('description', '').strip(),
            row_dict.get('application_criteria', '').strip(),
            row_dict.get('source', '').strip(),
            row_dict.get('version', '').strip(),
            row_dict.get('last_updated', '').strip(),
        ))

    print(f"  Imported {len(data_rows)} ClinGen rule records")


def create_rag_table(cursor: sqlite3.Cursor):
    """Create table for RAG text chunks."""
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS rag_chunks (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            chunk_text TEXT NOT NULL,
            source_table TEXT NOT NULL,
            source_id INTEGER,
            gene TEXT,
            rule_code TEXT,
            chunk_embedding BLOB,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_rag_gene ON rag_chunks(gene)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_rag_rule_code ON rag_chunks(rule_code)')


def main():
    parser = argparse.ArgumentParser(description='Convert ClinGen Excel files to SQLite')
    parser.add_argument(
        '--output-db',
        type=Path,
        default=Path('clingen_rules.db'),
        help='Output SQLite database path'
    )
    parser.add_argument(
        '--clingen-dir',
        type=Path,
        default=Path('ClinGen'),
        help='ClinGen directory containing Excel files'
    )

    args = parser.parse_args()

    # Check input files
    freq_file = args.clingen_dir / 'ClinGen_rules_frequency_cutoffs.xlsx'
    pm1_file = args.clingen_dir / 'ClinGen_rules_pm1_details-已校正--2026-3-30.xlsx'
    rules_file = args.clingen_dir / 'ClinGen_rules_rules.xlsx'

    if not freq_file.exists():
        print(f"Error: Frequency cutoffs file not found: {freq_file}", file=sys.stderr)
        sys.exit(1)
    if not pm1_file.exists():
        print(f"Error: PM1 details file not found: {pm1_file}", file=sys.stderr)
        sys.exit(1)
    if not rules_file.exists():
        print(f"Error: Rules file not found: {rules_file}", file=sys.stderr)
        sys.exit(1)

    # Create database
    print(f"Creating database: {args.output_db}")
    conn = sqlite3.connect(args.output_db)
    cursor = conn.cursor()

    # Create tables
    create_frequency_cutoffs_table(cursor)
    create_pm1_details_table(cursor)
    create_clingen_rules_table(cursor)
    create_rag_table(cursor)

    # Import data
    import_frequency_cutoffs(cursor, freq_file)
    import_pm1_details(cursor, pm1_file)
    import_clingen_rules(cursor, rules_file)

    conn.commit()

    # Print statistics
    print("\n=== Database Statistics ===")

    cursor.execute('SELECT COUNT(*) FROM frequency_cutoffs')
    print(f"Frequency cutoffs: {cursor.fetchone()[0]:,}")

    cursor.execute('SELECT COUNT(*) FROM pm1_details')
    print(f"PM1 details: {cursor.fetchone()[0]:,}")

    cursor.execute('SELECT COUNT(*) FROM clingen_rules')
    print(f"ClinGen rules: {cursor.fetchone()[0]:,}")

    # Print unique genes per table
    cursor.execute('SELECT COUNT(DISTINCT gene) FROM frequency_cutoffs')
    print(f"Genes in frequency_cutoffs: {cursor.fetchone()[0]}")

    cursor.execute('SELECT COUNT(DISTINCT gene) FROM pm1_details')
    print(f"Genes in pm1_details: {cursor.fetchone()[0]}")

    cursor.execute('SELECT COUNT(DISTINCT gene) FROM clingen_rules')
    print(f"Genes in clingen_rules: {cursor.fetchone()[0]}")

    # Print rule distribution
    cursor.execute('''
        SELECT rule_code, COUNT(*) as cnt
        FROM clingen_rules
        GROUP BY rule_code
        ORDER BY cnt DESC
    ''')
    print("\nRule distribution in ClinGen rules:")
    for row in cursor.fetchall():
        print(f"  {row[0]}: {row[1]:,}")

    conn.close()
    print("\nDone!")


if __name__ == '__main__':
    sys.exit(main())
