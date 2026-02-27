################################################################################
#
#  Summary-based Mendelian Randomisation (SMR) Analysis
#  Method: Wald ratio with delta-method standard errors
#
#  Description:
#  This script estimates the causal effect of DNA methylation on a trait
#  of interest using the Wald ratio (single-instrument MR). For each
#  CpG–SNP–trait triplet, it combines meQTL and GWAS summary statistics
#  to compute an SMR estimate, its standard error, and a p-value.
#
#  Formulae (ref: https://github.com/jianyanglab/gsmr/blob/main/R/gsmr.R):
#    beta_SMR  = beta_GWAS / beta_meQTL
#    se_SMR    = sqrt( (se_GWAS^2 * beta_meQTL^2 +
#                        se_meQTL^2 * beta_GWAS^2) / beta_meQTL^4 )
#    test_stat = beta_SMR^2 / se_SMR^2          (chi-squared, df = 1)
#    p_SMR     = pchisq(test_stat, 1, lower.tail = FALSE)
#
#  Input:
#    meqtl_gwas_harmonised.rds
#      A data frame where each row is a CpG–SNP–trait triplet with
#      harmonised effect alleles. Required columns:
#        beta_meqtl      : meQTL effect size (SNP on CpG methylation)
#        se_meqtl        : meQTL standard error
#        beta_gwas_harm  : GWAS effect size (SNP on trait), harmonised to
#                          the same effect allele as the meQTL
#        se_gwas         : GWAS standard error
#      Additional annotation columns (e.g. SNP, CpG, chr, pos) are
#      retained unchanged.
#
#  Output:
#    smr_cpgcIMT.rds
#      The input data frame with four appended columns:
#        beta_SMR       : Wald ratio causal estimate
#        se_SMR         : delta-method standard error
#        smr_test_stat  : chi-squared test statistic (df = 1)
#        p_smr          : SMR p-value
#
#  Dependencies:
#    install.packages("dplyr")
#
################################################################################


# ==============================================================================
# 0. User-configurable parameters
# ==============================================================================

input_file  <- "meqtl_gwas_harmonised.rds"
output_file <- "smr_cpgcIMT.rds"


# ==============================================================================
# 1. Load data
# ==============================================================================

library(dplyr)

smr_df <- readRDS(input_file)
message("Loaded ", nrow(smr_df), " CpG-SNP-trait triplets")


# ==============================================================================
# 2. Compute SMR statistics (Wald ratio + delta-method SE)
# ==============================================================================

# Remove rows where the meQTL effect is zero (undefined Wald ratio)
smr_df <- smr_df %>%
  filter(beta_meqtl != 0) %>%
  mutate(
    # Wald ratio: causal effect estimate
    beta_SMR = beta_gwas_harm / beta_meqtl,

    # Delta-method standard error for the ratio estimate
    se_SMR = sqrt((se_gwas^2 * beta_meqtl^2 +
                   se_meqtl^2 * beta_gwas_harm^2) / beta_meqtl^4),

    # Chi-squared test statistic (1 df)
    smr_test_stat = beta_SMR^2 / se_SMR^2,

    # Two-sided p-value from chi-squared distribution
    p_smr = pchisq(smr_test_stat, df = 1, lower.tail = FALSE)
  )

message("SMR estimates computed for ", nrow(smr_df), " triplets ",
        "(excluded ", nrow(readRDS(input_file)) - nrow(smr_df),
        " with beta_meqtl = 0)")


# ==============================================================================
# 3. Save results
# ==============================================================================

saveRDS(smr_df, output_file)
message("Results saved to: ", output_file)
