"""
Run Hardy-Weinberg Equilibrium tests (chi-square + exact) for every
SNP x population pair, with multiple-testing correction.

Reads:  db/hwe.sqlite (genotype_counts table)
Writes: db/hwe.sqlite (new hwe_results table), outputs/results.csv
"""

import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

DB_PATH = Path("db/hwe.sqlite")
OUT_PATH = Path("outputs/results.csv")


def wigginton_exact_test(obs_hets: int, obs_hom1: int, obs_hom2: int) -> float:
    """
    Exact test for HWE, following Wigginton, Cutler & Abecasis (2005),
    'A Note on Exact Tests of Hardy-Weinberg Equilibrium', Am J Hum Genet.
    More reliable than chi-square when expected genotype counts are small.
    """
    obs_homc = max(obs_hom1, obs_hom2)
    obs_homr = min(obs_hom1, obs_hom2)
    rare = 2 * obs_homr + obs_hets
    n = obs_hets + obs_homc + obs_homr

    if n == 0:
        return 1.0

    probs = np.zeros(rare + 1)
    mid = rare * (2 * n - rare) // (2 * n)
    if mid % 2 != rare % 2:
        mid += 1

    probs[mid] = 1.0
    total = probs[mid]

    curr_hets, curr_homr, curr_homc = mid, (rare - mid) // 2, n - mid - (rare - mid) // 2
    while curr_hets >= 2:
        probs[curr_hets - 2] = (
            probs[curr_hets] * curr_hets * (curr_hets - 1)
            / (4 * (curr_homr + 1) * (curr_homc + 1))
        )
        total += probs[curr_hets - 2]
        curr_homr += 1
        curr_homc += 1
        curr_hets -= 2

    curr_hets, curr_homr, curr_homc = mid, (rare - mid) // 2, n - mid - (rare - mid) // 2
    while curr_hets <= rare - 2:
        probs[curr_hets + 2] = (
            probs[curr_hets] * 4 * curr_homr * curr_homc
            / ((curr_hets + 2) * (curr_hets + 1))
        )
        total += probs[curr_hets + 2]
        curr_homr -= 1
        curr_homc -= 1
        curr_hets += 2

    target = probs[obs_hets]
    p_value = min(1.0, probs[probs <= target].sum() / total)
    return p_value


def chi_square_test(obs_hom1: int, obs_het: int, obs_hom2: int) -> tuple[float, float]:
    """Standard 1-df HWE chi-square test (allele freq estimated from the sample)."""
    n = obs_hom1 + obs_het + obs_hom2
    if n == 0:
        return np.nan, np.nan

    p = (2 * obs_hom1 + obs_het) / (2 * n)
    q = 1 - p

    exp_hom1 = p ** 2 * n
    exp_het = 2 * p * q * n
    exp_hom2 = q ** 2 * n

    # avoid divide-by-zero on fixed loci (no variation observed)
    stat = 0.0
    for obs, exp in [(obs_hom1, exp_hom1), (obs_het, exp_het), (obs_hom2, exp_hom2)]:
        if exp > 0:
            stat += (obs - exp) ** 2 / exp

    p_value = chi2.sf(stat, df=1)
    return stat, p_value


def genotype_counts_for_group(df: pd.DataFrame) -> tuple[int, int, int]:
    """Collapse a group's genotype rows into (hom1_count, het_count, hom2_count)."""
    counts = {}
    for _, row in df.iterrows():
        alleles = tuple(sorted(row["genotype"].split("|")))
        counts[alleles] = counts.get(alleles, 0) + row["count"]

    homs = [(alleles, c) for alleles, c in counts.items() if alleles[0] == alleles[1]]
    hets = [(alleles, c) for alleles, c in counts.items() if alleles[0] != alleles[1]]

    hom1_count = homs[0][1] if len(homs) > 0 else 0
    hom2_count = homs[1][1] if len(homs) > 1 else 0
    het_count = hets[0][1] if len(hets) > 0 else 0

    return hom1_count, het_count, hom2_count


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """FDR correction. Returns adjusted p-values aligned to the input order."""
    n = len(p_values)
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * n / (np.arange(n) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    out = np.empty(n)
    out[order] = adjusted
    return out


def main():
    conn = sqlite3.connect(DB_PATH)
    genotypes = pd.read_sql("SELECT * FROM genotype_counts", conn)

    results = []
    for (rsid, population), group in genotypes.groupby(["rsid", "population"]):
        hom1, het, hom2 = genotype_counts_for_group(group)
        n = hom1 + het + hom2

        chi2_stat, chi2_p = chi_square_test(hom1, het, hom2)
        exact_p = wigginton_exact_test(het, hom1, hom2)

        results.append({
            "rsid": rsid,
            "population": population,
            "n_samples": n,
            "hom1_count": hom1,
            "het_count": het,
            "hom2_count": hom2,
            "chi2_stat": chi2_stat,
            "chi2_p": chi2_p,
            "exact_p": exact_p,
        })

    results_df = pd.DataFrame(results)
    results_df["chi2_p_fdr"] = benjamini_hochberg(results_df["chi2_p"].values)
    results_df["exact_p_fdr"] = benjamini_hochberg(results_df["exact_p"].values)
    results_df["significant_deviation"] = results_df["exact_p_fdr"] < 0.05

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(OUT_PATH, index=False)
    results_df.to_sql("hwe_results", conn, if_exists="replace", index=False)
    conn.commit()
    conn.close()

    n_sig = results_df["significant_deviation"].sum()
    print(f"Tested {len(results_df)} SNP x population pairs")
    print(f"{n_sig} show significant HWE deviation after FDR correction (q < 0.05)")
    print(f"Results written to {OUT_PATH} and db/hwe.sqlite (hwe_results table)")


if __name__ == "__main__":
    main()