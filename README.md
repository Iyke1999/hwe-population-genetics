# Hardy-Weinberg Equilibrium Testing Across 1000 Genomes Populations

A reproducible Python/SQL pipeline testing Hardy-Weinberg equilibrium and
allele frequency differentiation across a panel of ancestry-informative and
forensically relevant SNPs, using public 1000 Genomes Project reference data.

This project extends the HWE methodology from my published research on
morphogenetic trait inheritance (see [Okolie et al., Google Scholar](https://scholar.google.com/citations?user=xeYI73QAAAAJ))
into molecular population genetics, rebuilt as an automated, API-driven
pipeline rather than manually collected sample data.

## Key Findings
- Tested 13 SNPs across 26 populations (338 SNP × population pairs)
- **No significant HWE deviations after FDR correction** — consistent with
  a clean, well-curated reference dataset
- Clear **allele frequency differentiation** across populations for markers
  like EDAR, the Duffy-null variant, and LCT — visualized in
  `outputs/figures/allele_frequency_heatmap.png`

See [METHODS.md](METHODS.md) for full methodology, panel rationale, and limitations.

## Project Structure
hwe-population-genetics/
├── README.md
├── METHODS.md
├── requirements.txt
├── data/
│ ├── panel.csv # SNP panel with rationale
│ ├── raw/ # cached Ensembl API responses
│ └── processed/ # cleaned genotype count table
├── src/
│ ├── fetch_genotypes.py
│ ├── parse_genotypes.py
│ ├── load_to_db.py
│ ├── hwe_test.py
│ └── visualize.py
├── sql/
│ └── schema.sql
├── outputs/
│ ├── figures/
│ └── results.csv
└── db/
└── hwe.sqlite


## How to Reproduce
```bash
git clone https://github.com/Iyke1999/hwe-population-genetics.git
cd hwe-population-genetics
python3 -m venv venv
source venv/bin/activate  # venv\Scripts\Activate.ps1 on Windows
pip install -r requirements.txt

python3 src/fetch_genotypes.py
python3 src/parse_genotypes.py
python3 src/load_to_db.py
python3 src/hwe_test.py
python3 src/visualize.py
```

## Background
Built by Ikechukwu Emmanuel Okolie — formerly a researcher in forensic genetics and anthropometry at Delta State University, now working as a research assistant. This project bridges that research background with a modern, reproducible data engineering and statistics pipeline.