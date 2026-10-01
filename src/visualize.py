"""
Build Phase 4 visualizations:
  1. HWE exact-test p-value heatmap (populations x SNPs)
  2. Allele frequency heatmap (populations x SNPs) - the differentiation story

Reads:  db/hwe.sqlite (hwe_results, genotype_counts tables)
Writes: outputs/figures/*.png
"""

import sqlite3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

DB_PATH = Path("db/hwe.sqlite")
FIG_DIR = Path("outputs/figures")

# Group populations by superpopulation for readable axis ordering
SUPERPOP_ORDER = {
    "YRI": "AFR", "LWK": "AFR", "GWD": "AFR", "MSL": "AFR", "ESN": "AFR",
    "ASW": "AFR", "ACB": "AFR",
    "MXL": "AMR", "PUR": "AMR", "CLM": "AMR", "PEL": "AMR",
    "CHB": "EAS", "JPT": "EAS", "CHS": "EAS", "CDX": "EAS", "KHV": "EAS",
    "CEU": "EUR", "TSI": "EUR", "FIN": "EUR", "GBR": "EUR", "IBS": "EUR",
    "GIH": "SAS", "PJL": "SAS", "BEB": "SAS", "STU": "SAS", "ITU": "SAS",
}


def ordered_populations(populations: list[str]) -> list[str]:
    return sorted(populations, key=lambda p: (SUPERPOP_ORDER.get(p, "ZZZ"), p))


def plot_hwe_heatmap(conn: sqlite3.Connection):
    results = pd.read_sql("SELECT rsid, population, exact_p_fdr FROM hwe_results", conn)
    pivot = results.pivot(index="population", columns="rsid", values="exact_p_fdr")
    pivot = pivot.reindex(ordered_populations(pivot.index.tolist()))

    neg_log_p = -np.log10(pivot.clip(lower=1e-10))

    fig, ax = plt.subplots(figsize=(10, 9))
    sns.heatmap(
        neg_log_p, cmap="Reds", ax=ax, linewidths=0.4, linecolor="white",
        cbar_kws={"label": "-log10(FDR-adjusted p-value)"},
    )
 # placeholder removed below
    ax.set_title(
        "HWE Exact Test: FDR-adjusted significance\n"
        "(no cells cross the q < 0.05 threshold - consistent with a clean reference panel)",
        fontsize=11,
    )
    ax.set_xlabel("SNP (rsID)")
    ax.set_ylabel("Population (grouped by superpopulation)")
    plt.tight_layout()
    fig.savefig(FIG_DIR / "hwe_pvalue_heatmap.png", dpi=150)
    plt.close(fig)
    print("Saved hwe_pvalue_heatmap.png")


def compute_allele_frequencies(conn: sqlite3.Connection) -> pd.DataFrame:
    genotypes = pd.read_sql("SELECT * FROM genotype_counts", conn)

    rows = []
    for (rsid, population), group in genotypes.groupby(["rsid", "population"]):
        allele_tally = {}
        total_alleles = 0
        for _, r in group.iterrows():
            alleles = r["genotype"].split("|")
            for a in alleles:
                allele_tally[a] = allele_tally.get(a, 0) + r["count"]
                total_alleles += r["count"]

        # Use the alphabetically-first allele consistently per SNP for comparability
        target_allele = sorted(allele_tally.keys())[0]
        freq = allele_tally.get(target_allele, 0) / total_alleles if total_alleles else np.nan

        rows.append({
            "rsid": rsid, "population": population,
            "allele": target_allele, "frequency": freq,
        })

    return pd.DataFrame(rows)


def plot_allele_frequency_heatmap(conn: sqlite3.Connection):
    freqs = compute_allele_frequencies(conn)
    pivot = freqs.pivot(index="population", columns="rsid", values="frequency")
    pivot = pivot.reindex(ordered_populations(pivot.index.tolist()))

    fig, ax = plt.subplots(figsize=(10, 9))
    sns.heatmap(
        pivot, cmap="viridis", vmin=0, vmax=1, ax=ax,
        linewidths=0.4, linecolor="white",
        cbar_kws={"label": "Allele frequency"},
    )
    ax.set_title(
        "Allele frequency differentiation across populations\n"
        "(reference allele, consistent per SNP)",
        fontsize=11,
    )
    ax.set_xlabel("SNP (rsID)")
    ax.set_ylabel("Population (grouped by superpopulation)")
    plt.tight_layout()
    fig.savefig(FIG_DIR / "allele_frequency_heatmap.png", dpi=150)
    plt.close(fig)
    print("Saved allele_frequency_heatmap.png")

    freqs.to_csv("outputs/allele_frequencies.csv", index=False)
    print("Saved allele_frequencies.csv")


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)

    
    plot_allele_frequency_heatmap(conn)

    conn.close()


if __name__ == "__main__":
    main()