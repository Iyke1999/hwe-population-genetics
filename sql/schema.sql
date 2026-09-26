CREATE TABLE IF NOT EXISTS snps (
    rsid TEXT PRIMARY KEY,
    gene TEXT NOT NULL,
    category TEXT,
    trait TEXT,
    rationale TEXT
);

CREATE TABLE IF NOT EXISTS genotype_counts (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    rsid TEXT NOT NULL,
    population TEXT NOT NULL,
    genotype TEXT NOT NULL,
    count INTEGER NOT NULL,
    frequency REAL NOT NULL,
    FOREIGN KEY (rsid) REFERENCES snps(rsid)
);

CREATE INDEX IF NOT EXISTS idx_genotype_rsid_pop
    ON genotype_counts (rsid, population);