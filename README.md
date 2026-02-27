# cIMT-EWAS

Collection of core R and shell scripts used in the carotid intima-media thickness
(cIMT) epigenome-wide association study (EWAS). The repository provides a
reproducible workflow from raw Illumina EPIC methylation data through EWAS,
meta-analysis, molecular QTL querying, causal inference (SMR/coloc), and
methylation risk score construction.

## Scripts

- `01_preprocess_methylation.R` – preprocess EPIC array IDAT files and generate methylation betas
- `02_run_ewas.R` – fit EWAS models and generate summary statistics
- `03_metaanalysis_stouffer.sh` – stouffer sample-size weighted meta-analysis 
- `04_molecular_qtl.R` – query/process molecular QTL data (code provided for meQTL; similar parameters used for eQTM and eQTL)
- `05_smr.R` – summary Mendelian randomization analyses (SMR)
- `06_coloc.R` – colocalization analysis between CpGs and traits
- `07_mrs.R` – build and test methylation risk score 
- `08_enrichment_background.R` – generate permutation-based background of cpgs with matched methylation levels and variability for enrichment testing


*Author: Konstanze Tan*
