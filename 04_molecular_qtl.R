################################################################################
#
#  Cis-meQTL Analysis Pipeline (Single Cohort)
#  Software: MatrixEQTL
#
#  Description:
#  This script tests cis associations (within 1 Mb) between genetic variants
#  and DNA methylation levels at sentinel CpG sites, adjusting for
#  demographic and technical covariates. It is designed to be called once
#  per cohort from the companion PBS script.
#
#  Model:
#    CpG ~ SNP + Age + Gender + Ethnicity + WBC (6 cell types) + 30 ctrl PCs
#
#  Notes on covariate choice:
#    - Ethnicity (categorical, recoded to numeric) is used instead of genetic
#      PCs for methodological consistency with the companion cis-eQTL analysis.
#    - HELIOS cohort is excluded from meQTL analysis to prevent issues with
#      sample overlap in downstream MR analysis.
#
#  Inputs:
#    1. Methylation data (RDS):
#         sentinel_cpgs_<cohort>_beta.rds
#         A data frame [CpG probes x samples] of beta values for sentinel
#         CpGs. Row names = CpG IDs, column names = SentrixIDs.
#
#    2. Genetic data (PLINK .traw format):
#         sg10k_maf0.01_hwe1e-03_Rsq0.30_sentinelcpgs1Mb_<cohort>.traw
#         Transposed raw genotype file. Key columns:
#           SNP     : variant ID
#           CHR     : chromosome
#           POS     : base-pair position (1-based)
#           COUNTED : counted (effect) allele
#           ALT     : alternative allele
#           (C)M    : unused column
#           <FID_IID> columns : additive allele dosages per sample
#
#    3. Covariates (RDS):
#         cov_sg10k_wsentrixIDcohorteth.rds
#         A data frame [samples x variables] with columns:
#           SentrixID  : methylation array sample ID
#           FID / IID  : genetic sample IDs (for ID matching)
#           Cohort     : cohort label (used to subset)
#           Ethnicity  : ethnic group code ("C", "M", "I", "O")
#           Age, Gender, WBC proportions, control probe PCs
#
#    4. Position files:
#         sentinel_cpgs.bed : BED file (0-based) with CpG positions
#           Columns: chr, start, end, cpgID
#         SNP positions are extracted from the .traw file directly.
#
#  Output:
#    cis_meqtl_AgeGenderEthWBC30cp_<cohort>
#      Tab-separated MatrixEQTL results for all cis pairs (<=1 Mb), with
#      columns: SNP, gene (CpG), beta, t-stat, p-value, FDR
#
#  Usage:
#    Rscript meqtl.R <cohort_name>
#    e.g. Rscript meqtl.R MEC
#
#  Dependencies:
#    install.packages(c("MatrixEQTL", "dplyr", "tibble", "stringr"))
#    BiocManager::install("Biobase")
#
################################################################################


# ==============================================================================
# 1. Load packages and parse command-line arguments
# ==============================================================================

library(devtools)
library(Biobase)
library(MatrixEQTL)
library(dplyr)
library(tibble)
library(stringr)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a cohort name was provided
if (length(args) == 0) {
  stop("No cohort name provided. Usage: Rscript meqtl.R <cohort_name>")
}

# Set the cohort name from the command-line argument
cohort <- args[1]


# ==============================================================================
# 2. Specify file paths
# ==============================================================================

methylation_data <- paste0("sentinel_cpgs_", tolower(cohort), "_beta.rds")
genetic_data     <- paste0("sg10k_maf0.01_hwe1e-03_Rsq0.30_sentinelcpgs1Mb_", cohort, ".traw")
covariates       <- c("cov_sg10k_wsentrixIDcohorteth.rds")
# genetic_pca    <- paste0("sg10k_maf0.01_hwe1e-03_Rsq0.30_", cohort, ".eigenvec")


# ==============================================================================
# 3. Load and format methylation data
# ==============================================================================

# Load covariates (needed here for sample ID matching)
cov_meqtl_sg10k <- readRDS(covariates)
# genetic_pcs   <- read.delim(genetic_pca, check.names=F)

# Load methylation data and subset to samples present in this cohort's covariates
beta_sentinelcpgs_formatted <- readRDS(methylation_data) %>%
  as.data.frame() %>%
  select(intersect(
    names(.),
    cov_meqtl_sg10k %>%
      filter(Cohort == cohort) %>%
      pull(SentrixID)
  ))


# ==============================================================================
# 4. Format covariates
# ==============================================================================
#
# Performed after methylation data formatting since the covariate file contains
# some samples without methylation data available (e.g. in PRISM).
# Covariates: Age, Gender, Ethnicity, WBC proportions, control probe PCs.
# Ethnicity is recoded to numeric: C=1, M=2, I=3, O=4.
# Final format: [covariates x samples] (transposed for MatrixEQTL).
# ==============================================================================

cov_formatted <- cov_meqtl_sg10k %>%
  filter(SentrixID %in% names(beta_sentinelcpgs_formatted)) %>%
  # left_join(genetic_pcs %>% select(IID, PC1, PC2, PC3, PC4, PC5), by = "IID") %>%
  select(-Cohort, -FID, -IID) %>%
  # Recode ethnicity and ensure all columns are numeric
  mutate(
    Ethnicity = case_when(
      Ethnicity == "C" ~ 1,
      Ethnicity == "M" ~ 2,
      Ethnicity == "I" ~ 3,
      Ethnicity == "O" ~ 4,
    )
  ) %>%
  tibble::column_to_rownames("SentrixID") %>%
  mutate(across(everything(), as.numeric)) %>%
  t() %>%
  as.data.frame()


# ==============================================================================
# 5. Load and format genetic data
# ==============================================================================
#
# Read PLINK .traw file: rows = SNPs, columns = samples (as FID_IID).
# Strip the _IID suffix from column names, then convert FID to SentrixID
# using the covariate lookup table so all matrices share the same IDs.
# ==============================================================================

cis_snps <- read.delim(genetic_data, check.names = FALSE)

cis_snps_formatted <- cis_snps %>%
  {rownames(.) <- .$SNP; .} %>%                       # manually assign row names
  select(-SNP, -CHR, -`(C)M`, -POS, -COUNTED, -ALT) %>%
  rename_with(~ sub("_.*", "", .), everything()) %>%   # change IDs from FID_IID to FID only
  {colnames(.) <- cov_meqtl_sg10k$SentrixID[match(colnames(.), cov_meqtl_sg10k$FID)]; .}  # convert FID to SentrixID


# ==============================================================================
# 6. Align sample order across all matrices
# ==============================================================================
#
# All three matrices (betas, covariates, genotypes) must have identical
# column names in the same order.
# ==============================================================================

cov_formatted <- cov_formatted[, colnames(beta_sentinelcpgs_formatted)]
colnames(cov_formatted) == colnames(beta_sentinelcpgs_formatted)

cis_snps_formatted <- cis_snps_formatted[, colnames(beta_sentinelcpgs_formatted)]
colnames(cis_snps_formatted) == colnames(beta_sentinelcpgs_formatted)


# ==============================================================================
# 7. Load position files for cis-distance calculation
# ==============================================================================

# CpG positions from BED file (0-based coordinates)
# Same files that were used to extract the 1Mb SNPs from the plink files
cpgspos <- read.delim("sentinel_cpgs.bed", header = FALSE)
cpgspos <- cpgspos %>%
  select(4, 1, 2, 3) %>%
  rename(cpgid = 1, chr = 2, pos_start = 3, pos_end = 4) %>%
  mutate(chr = sub("^chr", "", chr))   # remove 'chr' prefix from the chr column

# SNP positions extracted from the .traw file (1-based coordinates)
snpspos <- cis_snps %>%
  select(SNP, CHR, POS) %>%
  rename(snpid = SNP, chr = CHR, pos = POS)


# ==============================================================================
# 8. Create MatrixEQTL SlicedData objects
# ==============================================================================

cvrt <- SlicedData$new()
cvrt$CreateFromMatrix(as.matrix(cov_formatted))

snps <- SlicedData$new()
snps$CreateFromMatrix(as.matrix(cis_snps_formatted))
snps$fileSliceSize <- 2000   # read file in pieces of 2,000 rows

cpgs <- SlicedData$new()
cpgs$CreateFromMatrix(as.matrix(beta_sentinelcpgs_formatted))
cpgs$fileSliceSize <- 2000   # balance: larger = faster compute but slower to load


# ==============================================================================
# 9. Run cis-meQTL analysis
# ==============================================================================
#
# Standard additive linear model: methylation ~ dosage + covariates.
# All cis pairs within 1 Mb are tested; results saved to file.
# pvOutputThreshold.cis = 1 retains all pairs (filter downstream).
# pvOutputThreshold = 0 skips trans analysis entirely.
# ==============================================================================

output_file_name <- paste0("cis_meqtl_AgeGenderEthWBC30cp_", cohort)

me <- Matrix_eQTL_main(
  snps = snps,
  gene = cpgs,
  cvrt = cvrt,
  output_file_name.cis = output_file_name,
  pvOutputThreshold.cis = 1,
  pvOutputThreshold = 0,
  snpspos = snpspos,
  genepos = cpgspos,
  cisDist = 1e6,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)
