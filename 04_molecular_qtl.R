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
#      Tab-separated MatrixEQTL results for all cis pairs (≤1 Mb), with
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
# 0. User-configurable parameters
# ==============================================================================

# -- File path templates (use <COHORT> as placeholder) --
methylation_template <- "sentinel_cpgs_<COHORT>_beta.rds"
genetic_template     <- "sg10k_maf0.01_hwe1e-03_Rsq0.30_sentinelcpgs1Mb_<COHORT>.traw"
covariate_file       <- "cov_sg10k_wsentrixIDcohorteth.rds"
cpg_position_file    <- "sentinel_cpgs.bed"

# -- MatrixEQTL parameters --
cis_distance     <- 1e6   # cis window: 1 Mb
pv_threshold_cis <- 1     # report all cis pairs (filter later)
pv_threshold_trans <- 0   # set to 0 to skip trans analysis
slice_size       <- 2000  # rows per slice (memory vs speed trade-off)


# ==============================================================================
# 1. Parse command-line arguments
# ==============================================================================

library(Biobase)
library(MatrixEQTL)
library(dplyr)
library(tibble)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No cohort name provided. Usage: Rscript meqtl.R <cohort_name>")
}
cohort <- args[1]
message("=== meQTL analysis for cohort: ", cohort, " ===")

# Resolve file paths
methylation_file <- gsub("<COHORT>", tolower(cohort), methylation_template)
genetic_file     <- gsub("<COHORT>", cohort, genetic_template)


# ==============================================================================
# 2. Load and format methylation data
# ==============================================================================

# Beta matrix: [sentinel CpG probes x samples (SentrixIDs)]
cov_meqtl_sg10k <- readRDS(covariate_file)

beta_sentinelcpgs_formatted <- readRDS(methylation_file) %>%
  as.data.frame() %>%
  select(intersect(
    names(.),
    cov_meqtl_sg10k %>%
      filter(Cohort == cohort) %>%
      pull(SentrixID)
  ))

message("Methylation matrix: ", nrow(beta_sentinelcpgs_formatted),
        " CpGs x ", ncol(beta_sentinelcpgs_formatted), " samples")


# ==============================================================================
# 3. Load and format covariates
# ==============================================================================

# Subset to samples present in the methylation data (some covariate-file
# samples may lack methylation data, e.g. in PRISM).
# Ethnicity is recoded to numeric: C=1, M=2, I=3, O=4.
# Final format: [covariates x samples] (transposed for MatrixEQTL).
cov_formatted <- cov_meqtl_sg10k %>%
  filter(SentrixID %in% names(beta_sentinelcpgs_formatted)) %>%
  select(-Cohort, -FID, -IID) %>%
  mutate(Ethnicity = case_when(
    Ethnicity == "C" ~ 1,
    Ethnicity == "M" ~ 2,
    Ethnicity == "I" ~ 3,
    Ethnicity == "O" ~ 4
  )) %>%
  column_to_rownames("SentrixID") %>%
  mutate(across(everything(), as.numeric)) %>%
  t() %>%
  as.data.frame()

message("Covariates: ", nrow(cov_formatted), " variables x ",
        ncol(cov_formatted), " samples")


# ==============================================================================
# 4. Load and format genetic data
# ==============================================================================

# Read PLINK .traw file: rows = SNPs, columns = samples (as FID_IID).
# Strip the _IID suffix from column names, then convert FID to SentrixID
# using the covariate lookup table so all matrices share the same IDs.
cis_snps <- read.delim(genetic_file, check.names = FALSE)

cis_snps_formatted <- cis_snps %>%
  { rownames(.) <- .$SNP; . } %>%
  select(-SNP, -CHR, -`(C)M`, -POS, -COUNTED, -ALT) %>%
  rename_with(~ sub("_.*", "", .), everything()) %>%
  { colnames(.) <- cov_meqtl_sg10k$SentrixID[match(colnames(.), cov_meqtl_sg10k$FID)]; . }

message("Genotype matrix: ", nrow(cis_snps_formatted), " SNPs x ",
        ncol(cis_snps_formatted), " samples")


# ==============================================================================
# 5. Align sample order across all matrices
# ==============================================================================

# All three matrices (betas, covariates, genotypes) must have identical
# column names in the same order.
cov_formatted      <- cov_formatted[, colnames(beta_sentinelcpgs_formatted)]
cis_snps_formatted <- cis_snps_formatted[, colnames(beta_sentinelcpgs_formatted)]

stopifnot("Covariate columns do not match methylation columns" =
            all(colnames(cov_formatted) == colnames(beta_sentinelcpgs_formatted)))
stopifnot("Genotype columns do not match methylation columns" =
            all(colnames(cis_snps_formatted) == colnames(beta_sentinelcpgs_formatted)))


# ==============================================================================
# 6. Load position files for cis-distance calculation
# ==============================================================================

# CpG positions from BED file (0-based coordinates).
# Required columns for MatrixEQTL genepos: cpgid, chr, pos_start, pos_end
cpgspos <- read.delim(cpg_position_file, header = FALSE) %>%
  select(4, 1, 2, 3) %>%
  rename(cpgid = 1, chr = 2, pos_start = 3, pos_end = 4) %>%
  mutate(chr = sub("^chr", "", chr))

# SNP positions extracted from the .traw file.
# Required columns for MatrixEQTL snpspos: snpid, chr, pos (1-based)
snpspos <- cis_snps %>%
  select(SNP, CHR, POS) %>%
  rename(snpid = SNP, chr = CHR, pos = POS)


# ==============================================================================
# 7. Create MatrixEQTL SlicedData objects
# ==============================================================================

cvrt <- SlicedData$new()
cvrt$CreateFromMatrix(as.matrix(cov_formatted))

snps <- SlicedData$new()
snps$CreateFromMatrix(as.matrix(cis_snps_formatted))
snps$fileSliceSize <- slice_size

cpgs <- SlicedData$new()
cpgs$CreateFromMatrix(as.matrix(beta_sentinelcpgs_formatted))
cpgs$fileSliceSize <- slice_size


# ==============================================================================
# 8. Run cis-meQTL analysis
# ==============================================================================

# Standard additive linear model: methylation ~ dosage + covariates
# All cis pairs within 1 Mb are tested; results saved to file.
output_file_name <- paste0("cis_meqtl_AgeGenderEthWBC30cp_", cohort)

me <- Matrix_eQTL_main(
  snps       = snps,
  gene       = cpgs,
  cvrt       = cvrt,
  output_file_name.cis = output_file_name,
  pvOutputThreshold.cis = pv_threshold_cis,
  pvOutputThreshold     = pv_threshold_trans,
  snpspos    = snpspos,
  genepos    = cpgspos,
  cisDist    = cis_distance,
  useModel   = modelLINEAR,
  errorCovariance = numeric(),
  verbose    = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory   = FALSE
)

message("=== meQTL analysis complete for cohort: ", cohort, " ===")
message("Cis-meQTLs tested: ", me$cis$ntests)

sessionInfo()
