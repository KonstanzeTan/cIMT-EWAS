# cIMT-EWAS

This repository contains analysis scripts used in the epigenome-wide association study (EWAS) investigating carotid intima-media thickness (cIMT). The pipeline covers DNA methylation data preprocessing, EWAS execution, meta-analysis, and downstream causal inference using molecular QTLs, Mendelian randomization, colocalization and methylation risk score (MRS) construction.

## 📁 Repository Structure

```
01_preprocess_methylation.R   # methylation array preprocessing pipeline
02_run_ewas.R                 # run EWAS (model fitting & Covariate handling)
03_metaanalysis_stouffer.sh   # meta-analysis using Stouffer's method
04_molecular_qtl.R            # query and process molecular QTL results
05_causal_smr.R               # summary-data based Mendelian randomization (SMR)
06_causal_coloc.R             # colocalization analysis between CpGs and traits
07_mrs.R                      # construct and evaluate methylation risk score
```

## 🚀 Getting Started

1. **Clone the repository**
   ```bash
   git clone https://github.com/KonstanzeTan/cIMT-EWAS.git
   cd cIMT-EWAS
   ```

2. **Install required R packages**
   Each R script specifies its dependencies (e.g. `BiocManager`, `minfi`, etc.). Install them with:
   ```r
   install.packages("BiocManager")
   BiocManager::install(c("minfi", "IlluminaHumanMethylationEPICmanifest",
                          "S4Vectors", "limma"))
   ```
   Additional packages may be required for downstream analyses (e.g. `TwoSampleMR`, `coloc`, etc.). Consult individual scripts for details.

3. **Run the pipeline**
   - Step 1: preprocess raw methylation IDAT files using `01_preprocess_methylation.R`.
   - Step 2: execute EWAS with `02_run_ewas.R`.
   - Step 3: perform meta-analysis via `03_metaanalysis_stouffer.sh`.
   - Steps 4–6: conduct causal inference and QTL lookups.
   - Step 7: build an MRS score.

## 📌 Notes

* The preprocessing script assumes Illumina EPIC array data and requires a manifest file named `manifestb2.csv` in the working directory.
* Shell and R scripts are designed for batch execution and can be adapted to your own data/phenotypes.
* Output files produced by the pipeline are saved in the working directory; each script documents its outputs and required inputs.

## 🧪 Citation

If you use these scripts or adapt this pipeline for your research, please cite the associated manuscript (details to be added).

## 📬 Questions / Contributions

Open an issue or submit a pull request via the GitHub repository.


---

*Developed by Konstanze Tan – data analysis for cIMT EWAS.*
