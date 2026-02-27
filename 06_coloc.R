################################################################################
#
#  Bayesian Colocalisation Analysis (coloc.abf)
#
#  Description:
#  This script tests whether the meQTL signal at each sentinel CpG and the
#  GWAS signal for coronary artery disease (CAD) share a common causal
#  variant, using approximate Bayes factor colocalisation (coloc.abf). This
#  helps distinguish true causal mediation from linkage disequilibrium-driven
#  confounding in the SMR results.
#
#  Pipeline steps:
#    1. Prepare meQTL data: attach CpG-level standard deviations, format
#       as coloc-compatible lists (one per CpG)
#    2. Prepare GWAS data: format as a single coloc-compatible list
#    3. Run coloc.abf for each CpG against the shared GWAS dataset
#
#  Key coloc hypotheses tested:
#    H0 : no association with either trait at this locus
#    H1 : association with meQTL only
#    H2 : association with GWAS trait only
#    H3 : both traits associated, but with different causal variants
#    H4 : both traits associated, sharing the same causal variant
#    A high PP.H4 supports a shared causal mechanism.
#
#  Inputs:
#    - meqtl_gwas_harmonised.rds :
#        A data frame [SNP-CpG-trait triplets x variables] with harmonised
#        effect alleles. Required columns:
#          snp             : variant identifier
#          cpg             : CpG probe identifier
#          chrpos_hg19     : chr:pos string (hg19); position extracted below
#          beta_meqtl      : meQTL effect size
#          se_meqtl        : meQTL standard error
#          beta_gwas_harm  : GWAS effect size (allele-harmonised)
#          se_gwas         : GWAS standard error
#
#    - sentinel_cpgs_pooled_sd.rds :
#        A data frame [CpGs x 2] with columns:
#          cpg        : CpG probe identifier
#          pooled_sd  : pooled standard deviation of methylation beta
#                       values across the meta-analysed cohorts
#
#  Outputs:
#    - coloc_meqtl_formatted.rds : list of coloc-ready meQTL datasets
#    - coloc_gwas_formatted.rds  : coloc-ready GWAS dataset
#    - all_coloc_results.rds     : data frame [CpGs x 7] with coloc
#        summary statistics (nsnps, PP.H0-PP.H4) for each sentinel CpG
#
#  Dependencies:
#    install.packages(c("coloc", "dplyr", "tidyr", "stringr"))
#
################################################################################


# ==============================================================================
# 0. Load libraries and data
# ==============================================================================

library(coloc)
library(dplyr)
library(tidyr)
library(stringr)

# Load harmonised meQTL-GWAS summary statistics
coloc_df <- readRDS("meqtl_gwas_harmonised.rds")

# Extract numeric position from chr:pos string (hg19 coordinates)
coloc_df <- coloc_df %>%
  mutate(pos_hg19 = as.numeric(str_extract(chrpos_hg19, "(?<=:)\\d+")))

# Load pooled SD of methylation betas for each sentinel CpG
sentinel_cpgs_sd <- readRDS("sentinel_cpgs_pooled_sd.rds")

# Subset into separate meQTL and GWAS data frames
coloc_meqtl <- coloc_df %>%
  select(beta_meqtl, se_meqtl, snp, cpg, pos_hg19)

# A single SNP may be a meQTL for multiple CpGs, so deduplicate to produce
# one GWAS entry per SNP
coloc_gwas <- coloc_df %>%
  select(beta_gwas_harm, se_gwas, snp, pos_hg19) %>%
  distinct()


# ==============================================================================
# 1. Prepare meQTL dataset for coloc
# ==============================================================================
#
# For each CpG:
#   - Append the pooled SD of methylation (sdY), required by coloc for
#     quantitative traits to estimate prior variance on effect sizes
#   - Split into individual coloc-compatible named lists
#   - Set type = "quant" (methylation is a quantitative trait)
# ==============================================================================

# Append pooled SD information to meQTL data
coloc_meqtl_withsd <- coloc_meqtl %>%
  left_join(sentinel_cpgs_sd, by = "cpg") %>%
  rename(sdY = pooled_sd)

# Function to prepare coloc-compatible list for one CpG
prepare_coloc_data_meqtl <- function(df_cpg) {
  list(
    beta     = df_cpg$beta_meqtl,
    varbeta  = df_cpg$se_meqtl^2,       # variance = SE^2
    snp      = df_cpg$snp,
    position = df_cpg$pos_hg19,
    type     = "quant",
    sdY      = unique(df_cpg$sdY)[1]     # same for all SNPs within a CpG
  )
}

# Split by CpG and format; result is a named list of coloc datasets
coloc_meqtl_formatted <- lapply(split(coloc_meqtl_withsd, coloc_meqtl_withsd$cpg), prepare_coloc_data_meqtl)
saveRDS(coloc_meqtl_formatted, "coloc_meqtl_formatted.rds")


# ==============================================================================
# 2. Prepare GWAS dataset for coloc
# ==============================================================================
#
# CAD is a binary trait, so type = "cc" (case-control).
# ==============================================================================

prepare_coloc_data_gwas <- function(df_gwas) {
  list(
    beta     = df_gwas$beta_gwas_harm,   # harmonised effect size
    varbeta  = df_gwas$se_gwas^2,        # variance = SE^2
    snp      = df_gwas$snp,
    position = df_gwas$pos_hg19,
    type     = "cc"                      # case-control (binary trait)
  )
}

coloc_gwas_formatted <- prepare_coloc_data_gwas(coloc_gwas)
saveRDS(coloc_gwas_formatted, "coloc_gwas_formatted.rds")


# ==============================================================================
# 3. Run colocalisation analysis for each CpG
# ==============================================================================

combined_coloc_results <- list()

# Loop through all CpGs in coloc_meqtl_formatted
for (cpg in names(coloc_meqtl_formatted)) {

  # Run coloc.abf for the current CpG
  coloc_res <- coloc.abf(
    dataset1 = coloc_meqtl_formatted[[cpg]],
    dataset2 = coloc_gwas_formatted
  )

  # Extract the summary and add a column for the CpG name
  summary_res     <- as.data.frame(t(coloc_res$summary))
  summary_res$CpG <- cpg

  # Extract detailed SNP-level results and add a column for the CpG name
  snp_res     <- coloc_res$results
  snp_res$CpG <- cpg

  # Store both summary and SNP-level results in the list
  combined_coloc_results[[cpg]] <- list(
    summary     = summary_res,
    snp_results = snp_res
  )
}

# Combine all summaries into a single data frame
all_coloc_results <- do.call(rbind, lapply(combined_coloc_results, `[[`, "summary"))
saveRDS(all_coloc_results, "all_coloc_results.rds")
