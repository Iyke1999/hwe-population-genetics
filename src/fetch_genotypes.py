"""
Fetch population genotype counts for a panel of SNPs from the Ensembl REST API.

Reads:  data/panel.csv (columns: gene, rsid, category, trait, rationale)
Writes: data/raw/<rsid>.json (one raw API response per SNP)
"""

import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

ENSEMBL_SERVER = "https://rest.ensembl.org"
PANEL_PATH = Path("data/panel.csv")
RAW_DIR = Path("data/raw")
REQUEST_DELAY = 1  # seconds between requests — be polite to Ensembl's rate limits


def fetch_variant(rsid: str) -> dict:
    """Fetch population genotype data for a single rsID from Ensembl."""
    url = f"{ENSEMBL_SERVER}/variation/human/{rsid}"
    params = {"population_genotypes": 1, "content-type": "application/json"}
    response = requests.get(url, params=params, headers={"Content-Type": "application/json"})
    response.raise_for_status()
    return response.json()


def main():
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    panel = pd.read_csv(PANEL_PATH)

    for _, row in panel.iterrows():
        rsid = row["rsid"]
        out_path = RAW_DIR / f"{rsid}.json"

        if out_path.exists():
            print(f"Skipping {rsid} (already fetched)")
            continue

        print(f"Fetching {rsid} ({row['gene']})...")
        try:
            data = fetch_variant(rsid)
        except requests.exceptions.RequestException as e:
            print(f"  FAILED: {rsid} — {e}", file=sys.stderr)
            continue

        with open(out_path, "w") as f:
            json.dump(data, f, indent=2)

        time.sleep(REQUEST_DELAY)

    print("Done.")


if __name__ == "__main__":
    main()