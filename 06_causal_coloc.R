################################################################################
#
#  Bayesian Colocalisation Analysis (coloc.abf)
#
#  Description:
#  This script tests whether the meQTL signal at each sentinel CpG and the
#  GWAS signal for the trait of interest share a common causal variant,
#  using approximate Bayes factor colocalisation (coloc.abf). This helps
#  distinguish true causal mediation from linkage disequilibrium-driven
#  confounding in the SMR results.
#
#  Pipeline steps:
#    1. Prepare meQTL data: attach CpG-level standard deviations, format
#       as coloc-compatible lists (one per CpG)
#    2. Prepare GWAS data: format as a single coloc-compatible list
#    3. Run coloc.abf for each CpG against the shared GWAS dataset
#    4. Collect and save summary and SNP-level posterior probabilities
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
#    1. meqtl_gwas_harmonised.rds
#         A data frame [SNP-CpG-trait triplets x variables] with harmonised
#         effect alleles. Required columns:
#           snp             : variant identifier
#           cpg             : CpG probe identifier
#           pos             : base-pair position
#           beta_meqtl      : meQTL effect size
#           se_meqtl        : meQTL standard error
#           beta_gwas_harm  : GWAS effect size (allele-harmonised)
#           se_gwas         : GWAS standard error
#
#    2. sentinel_cpgs_pooled_sd.rds
#         A data frame [CpGs x 2] with columns:
#           cpg        : CpG probe identifier
#           pooled_sd  : pooled standard deviation of methylation beta
#                        values across the meta-analysed cohorts
#
#  Outputs:
#    - coloc_meqtl_formatted.rds : list of coloc-ready meQTL datasets
#    - all_coloc_results.rds     : data frame [CpGs x 7] with coloc
#        summary statistics (nsnps, PP.H0–PP.H4) for each sentinel CpG
#    - combined_coloc_results    : in-memory list containing both summary
#        and SNP-level posterior probabilities per CpG
#
#  Dependencies:
#    install.packages(c("coloc", "dplyr", "tidyr"))
#
################################################################################


# ==============================================================================
# 0. User-configurable parameters
# ==============================================================================

harmonised_file   <- "meqtl_gwas_harmonised.rds"
cpg_sd_file       <- "sentinel_cpgs_pooled_sd.rds"
output_meqtl_fmt  <- "coloc_meqtl_formatted.rds"
output_results    <- "all_coloc_results.rds"

# sdY for the GWAS trait. Set to 1 here because the cIMT phenotype was
# rank inverse-normal transformed prior to GWAS; adjust if using a
# different trait or transformation.
gwas_sdY <- 1


# ==============================================================================
# 1. Load libraries and data
# ==============================================================================

library(coloc)
library(dplyr)
library(tidyr)

coloc_df        <- readRDS(harmonised_file)
sentinel_cpgs_sd <- readRDS(cpg_sd_file)

message("Loaded ", nrow(coloc_df), " harmonised SNP-CpG-trait triplets")


# ==============================================================================
# 2. Prepare meQTL dataset for coloc
# ==============================================================================

# Subset to meQTL-relevant columns
coloc_meqtl <- coloc_df %>%
  select(beta_meqtl, se_meqtl, snp, cpg, pos)

# Append the pooled SD of methylation beta values for each CpG.
# sdY is required by coloc for quantitative traits to estimate
# prior variance on effect sizes.
coloc_meqtl_withsd <- coloc_meqtl %>%
  left_join(sentinel_cpgs_sd, by = "cpg") %>%
  rename(sdY = pooled_sd)

# Convert each CpG's data into a coloc-compatible named list:
#   beta     : vector of effect sizes
#   varbeta  : vector of variances (SE^2)
#   snp      : vector of SNP identifiers
#   position : vector of base-pair positions
#   type     : "quant" for quantitative trait
#   sdY      : standard deviation of the phenotype (methylation)
prepare_coloc_meqtl <- function(df_cpg) {
  list(
    beta     = df_cpg$beta_meqtl,
    varbeta  = df_cpg$se_meqtl^2,
    snp      = df_cpg$snp,
    position = df_cpg$pos,
    type     = "quant",
    sdY      = unique(df_cpg$sdY)[1]
  )
}

# Split by CpG and format; result is a named list of coloc datasets
coloc_meqtl_formatted <- lapply(
  split(coloc_meqtl_withsd, coloc_meqtl_withsd$cpg),
  prepare_coloc_meqtl
)

saveRDS(coloc_meqtl_formatted, output_meqtl_fmt)
message("Formatted meQTL data for ", length(coloc_meqtl_formatted), " CpGs")


# ==============================================================================
# 3. Prepare GWAS dataset for coloc
# ==============================================================================

# A single SNP may be a meQTL for multiple CpGs, so deduplicate to produce
# one GWAS entry per SNP.
coloc_gwas <- coloc_df %>%
  select(beta_gwas_harm, se_gwas, snp, pos) %>%
  distinct()

coloc_gwas_formatted <- list(
  beta     = coloc_gwas$beta_gwas_harm,
  varbeta  = coloc_gwas$se_gwas^2,
  snp      = coloc_gwas$snp,
  position = coloc_gwas$pos,
  type     = "quant",
  sdY      = gwas_sdY
)

message("GWAS dataset: ", length(coloc_gwas_formatted$snp), " unique SNPs")


# ==============================================================================
# 4. Run colocalisation analysis for each CpG
# ==============================================================================

combined_coloc_results <- list()

for (cpg in names(coloc_meqtl_formatted)) {
  message("Running coloc.abf for: ", cpg)

  coloc_res <- coloc.abf(
    dataset1 = coloc_meqtl_formatted[[cpg]],
    dataset2 = coloc_gwas_formatted
  )

  # Summary: one-row data frame with nsnps, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4
  summary_res     <- as.data.frame(t(coloc_res$summary))
  summary_res$CpG <- cpg

  # SNP-level: posterior probabilities for each SNP under each hypothesis
  snp_res     <- coloc_res$results
  snp_res$CpG <- cpg

  combined_coloc_results[[cpg]] <- list(
    summary     = summary_res,
    snp_results = snp_res
  )
}


# ==============================================================================
# 5. Combine and save results
# ==============================================================================

all_coloc_results <- do.call(rbind, lapply(combined_coloc_results, `[[`, "summary"))

saveRDS(all_coloc_results, output_results)
message("Colocalisation complete. Results for ", nrow(all_coloc_results),
        " CpGs saved to: ", output_results)
