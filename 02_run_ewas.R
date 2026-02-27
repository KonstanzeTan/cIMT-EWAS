################################################################################
#
#  Epigenome-Wide Association Study (EWAS) Regression Pipeline
#
#  Description:
#  This script performs an epigenome-wide association study testing the
#  association between DNA methylation (beta values) at each CpG site and
#  a continuous phenotype of interest, adjusting for covariates. The analysis
#  is run separately for each ethnic group.
#
#  Pipeline steps:
#    1. Load preprocessed data (methylation betas, control probe PCs,
#       white blood cell estimates, phenotype file)
#    2. Prepare phenotype data (recode variables, rename columns)
#    3. Subset all data to the target ethnic group
#    4. Compute ethnic-specific control probe PCs
#    5. Build a master covariate table merging phenotype, WBC proportions,
#       and control probe PCs
#    6. Align methylation beta matrix columns to the covariate table
#    7. Construct the regression formula
#    8. Run a linear regression at each CpG site and extract the
#       methylation coefficient
#
#  Inputs:
#    - beta_QN_rmGSwap_rmDup_HELIOS.RData :
#        Contains `betaRmDup`, a matrix [CpG probes x samples] of quantile-
#        normalised beta values after gender-swap and duplicate removal.
#    - ctrlprobes.RData :
#        Contains `ctrl.all`, a matrix [samples x control probes] of control
#        probe intensities (output of the preprocessing pipeline).
#    - housemanEstimates_QN_detPrelax_betas.RData :
#        Contains `constrainedCoefs`, a matrix [samples x cell types] of
#        estimated white blood cell proportions (Houseman method).
#    - Phenotype file (RDS) :
#        A data frame [samples x variables] with at minimum:
#          - SentrixID      : sample identifier matching beta column names
#          - Ethnic group column (e.g. "C", "M", "I")
#          - Outcome variable (continuous)
#          - Age, Gender, Smoking status
#
#  Outputs:
#    - Ethnic-specific control probe PCs (RDS)
#    - regout_<ethnicity>.RData :
#        A matrix [CpG probes x 4] containing, for each probe, the linear
#        model summary statistics for the methylation term:
#          Column 1: Estimate (regression coefficient)
#          Column 2: Std. Error
#          Column 3: t value
#          Column 4: Pr(>|t|) (p-value)
#
#  Usage:
#    Set the user-configurable parameters in Section 0 below, then run the
#    full script. To analyse a different ethnic group, change `ethnic_code`.
#
#  Dependencies:
#    install.packages(c("dplyr", "forcats"))
#
################################################################################


# ==============================================================================
# 0. User-configurable parameters
# ==============================================================================

# -- File paths (update to match your directory structure) --
beta_file       <- "beta_QN_rmGSwap_rmDup_HELIOS.RData"   # methylation betas
ctrl_file       <- "ctrlprobes.RData"                       # control probe intensities
wbc_file        <- "housemanEstimates_QN_detPrelax_betas.RData" # WBC estimates
phenotype_file  <- "cIMT_samplesheet_wsmoking.rds"          # phenotype data

# -- Analysis parameters --
ethnic_code     <- "C"              # Ethnic group code: "C" (Chinese), "M" (Malay), "I" (Indian)
ethnic_label    <- "CHN"            # Short label used in output filenames
outcome_var     <- "log_CIMT_distal_mean"  # Outcome variable name (will be created below)
sample_id_col   <- "SentrixID"      # Column in phenotype file matching beta column names
ethnic_col      <- "FREG5_Ethnic_Group"    # Column identifying ethnic group

# -- Covariates --
# Fixed covariates specified explicitly in the formula (categorical wrapped in as.factor)
fixed_covariates <- "Age + as.factor(Gender) + as.factor(smoking)"

# Column index ranges for additional covariates (control probe PCs and WBC proportions).
# These are appended programmatically. Adjust after inspecting the merged covariate table.
# Set to NULL to skip automatic covariate detection (and specify all covariates above).
ctrl_pc_col_range <- 21:50    # columns containing control probe PCs (after merge)
wbc_col_range     <- 225:230  # columns containing WBC proportions (after merge)


# ==============================================================================
# 1. Load required libraries and data
# ==============================================================================

library(dplyr)
library(forcats)

load(ctrl_file)          # loads: ctrl.all  [samples x control probes]
load(wbc_file)           # loads: constrainedCoefs  [samples x cell types]
load(beta_file)          # loads: betaRmDup  [CpG probes x samples]

phe <- readRDS(phenotype_file)
message("Phenotype dimensions: ", nrow(phe), " samples x ", ncol(phe), " variables")


# ==============================================================================
# 2. Prepare phenotype data
# ==============================================================================

# Create log-transformed outcome
phe$log_CIMT_distal_mean <- log(as.numeric(phe$CIMT_distal_mean))

# Rename columns for clarity and recode smoking into 3 categories:
#   Original codes: 1 = never, 2 = passive/minimal -> non-smoker
#                   3 = former                      -> ex-smoker
#                   4 = occasional, 5 = daily       -> smoker
phe <- phe %>%
  rename(Age = FREG8_Age, Gender = FREG7_Gender) %>%
  mutate(smoking = case_when(
    smoking %in% c(1, 2) ~ 1,
    smoking == 3          ~ 2,
    smoking %in% c(4, 5) ~ 3,
    is.na(smoking)        ~ NA_real_
  )) %>%
  mutate(smoking = factor(smoking,
                          levels = c(1, 2, 3),
                          labels = c("non-smoker", "ex-smoker", "smoker")))


# ==============================================================================
# 3. Subset to target ethnic group
# ==============================================================================

phe_ethnic <- phe[phe[[ethnic_col]] == ethnic_code, ]
message("Samples in ethnic group '", ethnic_code, "': ", nrow(phe_ethnic))

# Ethnic-specific sample IDs
ethnic_ids <- phe_ethnic[[sample_id_col]]


# ==============================================================================
# 4. Compute ethnic-specific control probe PCs
# ==============================================================================

# Subset control probe matrix to ethnic group and recompute PCA so that the
# PCs reflect within-group variation rather than between-group differences.
ctrl_ethnic <- ctrl.all[rownames(ctrl.all) %in% ethnic_ids, ]
message("Control probe matrix: ", nrow(ctrl_ethnic), " samples x ", ncol(ctrl_ethnic), " probes")

pca_cp <- prcomp(na.omit(ctrl_ethnic))
PC_cp  <- as.data.frame(predict(pca_cp))
PC_cp[[sample_id_col]] <- rownames(PC_cp)

saveRDS(PC_cp, paste0("PC_cp_", ethnic_label, ".rds"))


# ==============================================================================
# 5. Subset methylation betas to ethnic group
# ==============================================================================

# Remove CpG sites that are entirely NA across all samples
betaRmDup <- betaRmDup[rowSums(is.na(betaRmDup)) != ncol(betaRmDup), ]
message("Probes after removing all-NA rows: ", nrow(betaRmDup))

# Subset to ethnic-specific samples
beta_ethnic <- betaRmDup[, colnames(betaRmDup) %in% ethnic_ids]
message("Ethnic-specific beta matrix: ", nrow(beta_ethnic), " probes x ", ncol(beta_ethnic), " samples")


# ==============================================================================
# 6. Build master covariate table
# ==============================================================================

# Merge phenotype, control probe PCs, and WBC proportions by sample ID.
# constrainedCoefs: [samples x cell types] from Houseman estimation
wbc <- as.data.frame(constrainedCoefs)
wbc[[sample_id_col]] <- rownames(wbc)

covar <- merge(phe_ethnic, PC_cp, by = sample_id_col)
covar <- merge(covar, wbc, by = sample_id_col)
message("Master covariate table: ", nrow(covar), " samples x ", ncol(covar), " variables")


# ==============================================================================
# 7. Align beta matrix columns to covariate table rows
# ==============================================================================

# Ensure the sample order in the beta matrix matches the covariate table exactly.
beta_ethnic <- beta_ethnic[, match(covar[[sample_id_col]], colnames(beta_ethnic))]
stopifnot("Sample order mismatch between betas and covariates" =
            all(colnames(beta_ethnic) == covar[[sample_id_col]]))


# ==============================================================================
# 8. Construct regression formula
# ==============================================================================

# Build the covariate side of the formula. Additional covariates (control probe
# PCs and WBC proportions) are identified by column index in the merged table.
additional_covars <- colnames(covar)[c(ctrl_pc_col_range, wbc_col_range)]
message("Number of additional covariates (PCs + WBC): ", length(additional_covars))

rhs <- paste(c("beta_i", fixed_covariates, additional_covars), collapse = " + ")
formula_str <- paste(outcome_var, "~", rhs)
message("Regression formula:\n  ", formula_str)


# ==============================================================================
# 9. Run EWAS: probe-by-probe linear regression
# ==============================================================================

# For each CpG site, fit:
#   outcome ~ beta_at_CpG + Age + Gender + Smoking + control_probe_PCs + WBC
# and extract the coefficient summary for the methylation term (row 2).
#
# Output: regout [CpG probes x 4] with columns:
#   Estimate, Std. Error, t value, Pr(>|t|)

n_probes <- nrow(beta_ethnic)
regout   <- matrix(NA_real_, nrow = n_probes, ncol = 4)

message("Running regression across ", n_probes, " probes...")

for (i in seq_len(n_probes)) {
  if (i %% 10000 == 0) message("  probe ", i, " / ", n_probes)

  # Add current CpG beta values as a column in the covariate frame
  covar$beta_i <- as.numeric(beta_ethnic[i, ])

  reg_formula <- as.formula(formula_str)

  fit <- tryCatch(
    summary(lm(reg_formula, data = covar, na.action = na.exclude)),
    error = function(e) NULL
  )

  if (!is.null(fit) && nrow(fit$coefficients) >= 2) {
    regout[i, ] <- fit$coefficients[2, ]
  }
}

rownames(regout) <- rownames(beta_ethnic)
colnames(regout) <- c("Estimate", "Std.Error", "t.value", "P.value")

message("Regression complete. Results matrix: ",
        nrow(regout), " probes x ", ncol(regout), " statistics")


# ==============================================================================
# 10. Save results
# ==============================================================================

outfile <- paste0("regout_", outcome_var, "_", ethnic_label, "_adjusted.RData")
save(regout, file = outfile)
message("Results saved to: ", outfile)

sessionInfo()
