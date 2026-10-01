# Methods

## Research Question
Do single-nucleotide polymorphisms known to differentiate human populations
(pigmentation, metabolic, and forensic ancestry-informative markers) show
departures from Hardy-Weinberg equilibrium within 1000 Genomes Project
reference populations, and how does allele frequency vary across those
populations?

## Motivation
This project extends the Hardy-Weinberg equilibrium testing methodology used
in my published research on morphogenetic traits (fingerprint patterns,
palmar creases, facial and ear morphology in Nigerian populations) into
molecular population genetics. The statistical approach — testing observed
genotype counts against HWE expectations — is the same; here it is applied
to SNP data via a reproducible, API-driven pipeline rather than manually
collected sample data.

## Data Source
Genotype data was pulled from the **1000 Genomes Project Phase 3** reference
panel via the **Ensembl REST API** (`/variation/human/:id?population_genotypes=1`).
Data covers 26 populations across 5 superpopulations (AFR, AMR, EAS, EUR, SAS).

## SNP Panel
13 SNPs were selected across 6 functional categories (pigmentation, hair
morphology, blood group, metabolism, immune/malaria resistance, and athletic
performance), chosen for known strong population differentiation and/or
direct relevance to forensic genetics. Full panel and rationale: `data/panel.csv`.

## Pipeline
1. **Fetch** (`src/fetch_genotypes.py`) — raw genotype counts per SNP pulled from Ensembl
2. **Parse** (`src/parse_genotypes.py`) — filtered to the 26 real 1000 Genomes sub-populations (excluding superpopulation aggregates and non-1000G sources), normalized into a tidy table
3. **Load** (`src/load_to_db.py`) — loaded into a SQLite database (`db/hwe.sqlite`)
4. **Test** (`src/hwe_test.py`) — chi-square and exact (Wigginton et al. 2005) HWE tests per SNP × population pair, with Benjamini-Hochberg FDR correction across all 338 tests
5. **Visualize** (`src/visualize.py`) — HWE significance heatmap and allele-frequency heatmap across populations

## Results Summary
No SNP × population pair showed a statistically significant HWE deviation
after FDR correction (lowest q-value: 0.465). This is an expected and
reassuring result — it is consistent with 1000 Genomes being a well-curated
reference panel, rather than evidence of technical genotyping artifacts.
HWE testing within a population is a quality-control check; it is distinct
from — and does not contradict — the panel's known strength at detecting
**allele frequency differentiation between populations**, which is the more
biologically meaningful signal for this SNP panel (see
`outputs/figures/allele_frequency_heatmap.png`).

## Limitations
- **Reference panel size**: populations range from ~60–110 individuals,
  limiting statistical power to detect subtle HWE deviations
- **Ascertainment bias**: SNPs were selected because they are known to be
  informative, not randomly sampled from the genome — results should not
  be generalized to genome-wide expectations
- **Population labels**: 1000 Genomes population codes reflect sampling
  location, not self-identified ethnicity or genetic ancestry in full
- **Multiple testing**: FDR correction controls the expected proportion of
  false positives across 338 simultaneous tests, but does not eliminate
  the possibility of true effects being missed due to conservative correction