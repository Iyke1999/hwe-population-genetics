"""
Load the SNP panel and tidy genotype counts into SQLite.

Reads:  data/panel.csv, data/processed/genotype_counts.csv, sql/schema.sql
Writes: db/hwe.sqlite
"""

import sqlite3
from pathlib import Path

import pandas as pd

PANEL_PATH = Path("data/panel.csv")
GENOTYPES_PATH = Path("data/processed/genotype_counts.csv")
SCHEMA_PATH = Path("sql/schema.sql")
DB_PATH = Path("db/hwe.sqlite")


def main():
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)

    with open(SCHEMA_PATH) as f:
        conn.executescript(f.read())

    panel = pd.read_csv(PANEL_PATH)
    panel.to_sql("snps", conn, if_exists="replace", index=False)

    genotypes = pd.read_csv(GENOTYPES_PATH)
    genotypes.to_sql("genotype_counts", conn, if_exists="replace", index=False)

    conn.commit()

    snp_count = conn.execute("SELECT COUNT(*) FROM snps").fetchone()[0]
    row_count = conn.execute("SELECT COUNT(*) FROM genotype_counts").fetchone()[0]
    print(f"Loaded {snp_count} SNPs and {row_count} genotype-count rows into {DB_PATH}")

    conn.close()


if __name__ == "__main__":
    main()