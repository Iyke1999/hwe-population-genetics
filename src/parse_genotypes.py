"""
Parse raw Ensembl JSON responses into a tidy genotype-count table.

Reads:  data/raw/<rsid>.json, data/panel.csv
Writes: data/processed/genotype_counts.csv
"""

import json
from pathlib import Path

import pandas as pd

RAW_DIR = Path("data/raw")
PANEL_PATH = Path("data/panel.csv")
OUT_PATH = Path("data/processed/genotype_counts.csv")

# The 26 real 1000 Genomes Phase 3 sub-populations.
# Excludes "ALL" and the 5 superpopulation aggregates (AFR, AMR, EAS, EUR, SAS),
# which are sums of these and would double-count if included.
VALID_POPULATIONS = {
    "ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN",
    "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL",
    "PEL", "PJL", "PUR", "STU", "TSI", "YRI",
}


def parse_one(rsid: str, gene: str) -> list[dict]:
    """Parse a single raw JSON file into tidy rows."""
    path = RAW_DIR / f"{rsid}.json"
    if not path.exists():
        print(f"  WARNING: no raw file for {rsid}, skipping")
        return []

    with open(path) as f:
        data = json.load(f)

    rows = []
    for entry in data.get("population_genotypes", []):
        pop_field = entry["population"]

        if not pop_field.startswith("1000GENOMES:phase_3:"):
            continue

        pop_code = pop_field.split(":")[-1]
        if pop_code not in VALID_POPULATIONS:
            continue

        rows.append({
            "rsid": rsid,
            "gene": gene,
            "population": pop_code,
            "genotype": entry["genotype"],
            "count": entry["count"],
            "frequency": entry["frequency"],
        })

    return rows


def main():
    panel = pd.read_csv(PANEL_PATH)
    all_rows = []

    for _, row in panel.iterrows():
        print(f"Parsing {row['rsid']} ({row['gene']})...")
        all_rows.extend(parse_one(row["rsid"], row["gene"]))

    tidy = pd.DataFrame(all_rows)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    tidy.to_csv(OUT_PATH, index=False)

    print(f"\nWrote {len(tidy)} rows to {OUT_PATH}")
    print(f"Covering {tidy['rsid'].nunique()} SNPs across {tidy['population'].nunique()} populations")


if __name__ == "__main__":
    main()