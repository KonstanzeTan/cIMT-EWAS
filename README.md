# cIMT-EWAS

Lightweight collection of R and shell scripts used in the carotid intima-media thickness
(cIMT) epigenome-wide association study (EWAS). The repository provides a
reproducible workflow from raw Illumina EPIC methylation data through EWAS,
meta-analysis, molecular QTL querying, causal inference (SMR/coloc), and
methylation risk score construction.

## Scripts

- `01_preprocess_methylation.R` – preprocess EPIC array IDAT files
- `02_run_ewas.R` – fit EWAS models and generate summary statistics
- `03_metaanalysis_stouffer.sh` – combine results across cohorts
- `04_molecular_qtl.R` – query/process molecular QTL data
- `05_smr.R` – summary Mendelian randomization analyses (SMR)
- `06_coloc.R` – colocalization analysis between CpGs and traits
- `07_mrs.R` – build and test methylation risk score
- `08_enrichment_background.R` – generate annotation background for enrichment

## Usage

1. Clone repo:
   ```bash
   git clone https://github.com/KonstanzeTan/cIMT-EWAS.git
   cd cIMT-EWAS
   ```
2. Install required R packages (e.g. via `BiocManager`).
3. Run scripts sequentially, adjusting inputs/parameters as needed.

Outputs are saved in the working directory; inspect individual scripts for
requirements and produced files.

## Notes

* Designed for batch execution; adapt to other phenotypes or datasets.
* Scripts are intended as analysis templates rather than polished software.

---

*Author: Konstanze Tan*