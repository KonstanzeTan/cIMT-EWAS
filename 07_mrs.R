################################################################################
#
#  Methylation Risk Score (MRS) Pipeline
#
#  Description:
#  This script constructs and evaluates a methylation risk score (MRS) for
#  predicting elevated carotid intima-media thickness (cIMT). The pipeline
#  uses sentinel CpG sites identified from the EWAS, derives weights via
#  ethnic-stratified regressions meta-analysed with METAL (Stouffer method),
#  and evaluates predictive performance in a held-out test set.
#
#  Pipeline steps:
#    1. Stratified train-test split (70/30), preserving ethnicity and
#       risk-status proportions
#    2. Derive MRS weights: ethnic-specific regressions of log(cIMT) on
#       sentinel CpG betas, then meta-analyse across ancestries
#    3. Construct MRS in the test set as a weighted sum of scaled betas
#    4. Define MRS quartiles based on control (non-risk) individuals
#    5. Evaluate discriminative performance (AUC) across models
#    6. Assess MRS-quartile association with cIMT risk (logistic regression)
#
#  This script demonstrates the pipeline for ONE ethnicity (Chinese).
#  To replicate for Malay or Indian, change `ethnic_code` and `ethnic_label`
#  in the relevant sections. The meta-analysis step (METAL) combines all
#  three ancestry groups.
#
#  Inputs:
#    - phe_wcIMTref.rds :
#        Phenotype data frame [samples x variables] with columns:
#          SentrixID      : methylation array sample ID
#          Ethnicity      : ethnic group code ("C", "M", "I")
#          avgcIMT_75c    : cIMT measurement
#          at_risk_75c    : risk classification ("Y" / "N") based on
#                           age- and sex-specific 75th percentile
#          log_CIMT_distal_mean : log-transformed cIMT
#          Age, Gender, smoking
#    - sentinel_cpgs_beta.rds :
#        Matrix [sentinel CpG probes x samples] of QN beta values
#    - ctrlprobes.RData :
#        Contains `ctrl.all` [samples x control probes]
#    - housemanEstimates_QN_detPrelax_betas.RData :
#        Contains `constrainedCoefs` [samples x WBC cell types]
#
#  Outputs:
#    - phe_wcIMTref_traintest.rds        : phenotype with train/test assignment
#    - regout_*_sentinel.rds             : ethnic-specific regression weights
#    - METAL .tbl files                  : meta-analysed MRS weights
#    - sentinel_mrs_*_test_constructed.rds : test-set MRS values
#    - summary_roc_by_ethnicity.rds      : AUC results across models
#    - combined_regout_MRSquartile.xlsx   : quartile-based risk associations
#
#  Dependencies:
#    install.packages(c("dplyr", "caret", "forcats", "pROC", "ggplot2",
#                       "tidyr", "writexl"))
#
################################################################################


# ==============================================================================
# 0. User-configurable parameters
# ==============================================================================

# -- File paths --
phenotype_file         <- "phe_wcIMTref.rds"
sentinel_beta_file     <- "sentinel_cpgs_beta.rds"
ctrl_probe_file        <- "ctrlprobes.RData"
wbc_file               <- "housemanEstimates_QN_detPrelax_betas.RData"

# -- Train-test split --
train_proportion       <- 0.70
random_seed            <- 123

# -- Ethnic group for weight derivation (repeat for "M" and "I") --
ethnic_code            <- "C"
ethnic_label           <- "CHN"

# -- Column references --
sample_id_col          <- "SentrixID"
outcome_var            <- "log_CIMT_distal_mean"
risk_var               <- "at_risk_75c"       # binary risk classification
ethnicity_col          <- "Ethnicity"


# ==============================================================================
# STEP 1: Stratified train-test split
# ==============================================================================

library(dplyr)
library(caret)

set.seed(random_seed)

phe <- readRDS(phenotype_file) %>%
  filter(!is.na(avgcIMT_75c))

# Create a combined stratification variable to preserve joint proportions
# of ethnicity and risk status in both train and test sets.
phe$strata <- interaction(phe[[ethnicity_col]], phe[[risk_var]])
trainIndex  <- createDataPartition(phe$strata, p = train_proportion, list = FALSE)

phe$Set <- "test"
phe$Set[trainIndex] <- "train"

# Verify stratification
message("=== Ethnicity proportions by set ===")
print(prop.table(table(phe$Set, phe[[ethnicity_col]]), margin = 1))
message("=== Risk-status proportions by set ===")
print(prop.table(table(phe$Set, phe[[risk_var]]), margin = 1))

saveRDS(phe, "phe_wcIMTref_traintest.rds")


# ==============================================================================
# STEP 2: Derive MRS weights (ethnic-specific regression in training set)
# ==============================================================================
#
# For each sentinel CpG, regress log(cIMT) ~ CpG_beta + covariates in the
# training set of one ethnic group. Two covariate models are fitted:
#   Model A: Age + Gender                  (minimal)
#   Model B: Age + Gender + Smoking + WBC + 30 control probe PCs (full)
#
# This section is run once per ethnicity. Change `ethnic_code` / `ethnic_label`
# and re-run for Malay ("M"/"MLY") and Indian ("I"/"IND").
# ==============================================================================

library(forcats)

load(ctrl_probe_file)   # loads: ctrl.all
load(wbc_file)          # loads: constrainedCoefs
sentinel_cpgs_beta <- readRDS(sentinel_beta_file)

# -- 2a. Subset phenotype to ethnic-specific training samples --
phe_train <- readRDS("phe_wcIMTref_traintest.rds") %>%
  filter(Set == "train", !!sym(ethnicity_col) == ethnic_code)

message("Training samples for ", ethnic_label, ": ", nrow(phe_train))

# -- 2b. Ethnic-specific control probe PCs --
ctrl_ethnic <- ctrl.all[rownames(ctrl.all) %in% phe_train[[sample_id_col]], ]
pca_cp      <- prcomp(na.omit(ctrl_ethnic))
PC_cp       <- as.data.frame(predict(pca_cp))
PC_cp[[sample_id_col]] <- rownames(PC_cp)

saveRDS(PC_cp, paste0("PC_cp_", ethnic_label, "_train.rds"))

# -- 2c. Subset methylation to ethnic training samples --
beta_ethnic <- sentinel_cpgs_beta[, intersect(colnames(sentinel_cpgs_beta),
                                              phe_train[[sample_id_col]])]

# -- 2d. Build master covariate table --
wbc <- as.data.frame(constrainedCoefs)
wbc[[sample_id_col]] <- rownames(wbc)

covar <- merge(phe_train, PC_cp, by = sample_id_col)
covar <- merge(covar, wbc, by = sample_id_col)
message("Covariate table: ", nrow(covar), " samples x ", ncol(covar), " variables")

# -- 2e. Align beta columns to covariate rows --
beta_ethnic <- beta_ethnic[, match(covar[[sample_id_col]], colnames(beta_ethnic))]
stopifnot(all(colnames(beta_ethnic) == covar[[sample_id_col]]))

# -- 2f. Identify additional covariates by name --
# Control probe PCs (columns starting with "PC") and WBC proportions
pc_cols  <- grep("^PC[0-9]+$", colnames(covar), value = TRUE)
wbc_cols <- colnames(constrainedCoefs)
additional_covars <- intersect(c(pc_cols, wbc_cols), colnames(covar))

# -- 2g. Helper function: run CpG-by-CpG regression --
run_cpg_regression <- function(beta_mat, covar_df, formula_str) {
  n_probes <- nrow(beta_mat)
  regout   <- matrix(NA_real_, nrow = n_probes, ncol = 5)
  colnames(regout) <- c("Estimate", "Std.Error", "t.value", "P.value", "N_eff")

  for (i in seq_len(n_probes)) {
    if (i %% 50 == 0) message("  probe ", i, " / ", n_probes)
    covar_df$beta_i <- as.numeric(beta_mat[i, ])
    fit <- tryCatch(
      lm(as.formula(formula_str), data = covar_df, na.action = na.exclude),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      regout[i, ] <- c(summary(fit)$coefficients[2, ], nobs(fit))
    }
  }
  rownames(regout) <- rownames(beta_mat)
  return(regout)
}

# -- Model A: Age + Gender --
formula_A <- paste0(outcome_var, " ~ beta_i + Age + as.factor(Gender)")
regout_A  <- run_cpg_regression(beta_ethnic, covar, formula_A)
saveRDS(regout_A, paste0("regout_lncIMTmean_", ethnic_label, "_AgeGenderAdj_sentinel.rds"))

# -- Model B: Age + Gender + Smoking + WBC + 30 control probe PCs --
formula_B <- paste0(outcome_var, " ~ beta_i + Age + as.factor(Gender) + as.factor(smoking) + ",
                    paste(additional_covars, collapse = " + "))
regout_B  <- run_cpg_regression(beta_ethnic, covar, formula_B)
saveRDS(regout_B, paste0("regout_lncIMTmean_", ethnic_label, "_AgeGenderSmokWBCPcAdj_sentinel.rds"))

message("Weight derivation complete for: ", ethnic_label)


# ==============================================================================
# STEP 2b: Meta-analyse weights across ethnicities (METAL)
# ==============================================================================
#
# After running Step 2 for all three ethnicities and formatting the output
# files as tab-separated text (with columns: CpG, dir_beta, p, N), run METAL
# with the following script. Save as e.g. "mrs_meta.metal" and execute:
#   metal mrs_meta.metal
#
# ---------- BEGIN METAL SCRIPT (save separately) ----------
# SEPARATOR   TAB
# MARKERLABEL CpG
# PVALUE      p
# EFFECT      dir_beta
# SCHEME      SAMPLESIZE
# WEIGHTLABEL N
#
# PROCESS regout_lncIMTmean_CHN_AgeGenderAdj_sentinel_formatted.txt
# PROCESS regout_lncIMTmean_MLY_AgeGenderAdj_sentinel_formatted.txt
# PROCESS regout_lncIMTmean_IND_AgeGenderAdj_sentinel_formatted.txt
#
# OUTFILE MA_sentinel_ASIAN_HELIOS_CHNMLYIND_stouffer_AgeGender .tbl
# ANALYZE HETEROGENEITY
# QUIT
# ---------- END METAL SCRIPT ----------


# ==============================================================================
# STEP 3: Construct MRS in the test set
# ==============================================================================
#
# MRS = sum( scaled_beta_CpG_i * weight_i )
# Weights are Stouffer Z-scores from the METAL meta-analysis.
# Betas are z-score-scaled per CpG within the test set.
# ==============================================================================

# -- 3a. Load test-set data --
# `sentinel_mrs_df` is a pre-assembled data frame containing:
#   - Phenotype columns (Age, Gender, Ethnicity, at_risk_75c, etc.)
#   - Sentinel CpG beta columns (named "cg...")
#   Rows = test-set samples for one ethnicity.

ethnicity_test <- "CHN"   # change for each ethnicity
sentinel_mrs_df <- readRDS(paste0("sentinel_mrs_", ethnicity_test, "_test_df.rds"))

# -- 3b. Load meta-analysed weights (Z-scores) --
ma_agegender <- read.delim("MA_sentinel_ASIAN_HELIOS_CHNMLYIND_stouffer_AgeGender1.tbl")
ma_fullcov   <- read.delim("MA_sentinel_ASIAN_HELIOS_CHNMLYIND_stouffer_AgeGenderSmokWBCPc1.tbl")

# -- 3c. Identify sentinel CpG columns and scale betas --
sentinel_cpgs <- grep("^cg", colnames(sentinel_mrs_df), value = TRUE)

scaled_betas <- sentinel_mrs_df[, sentinel_cpgs] %>%
  mutate(across(everything(), scale))

# -- 3d. Create named weight vectors --
weights_A <- setNames(ma_agegender$Zscore, ma_agegender$MarkerName)
weights_B <- setNames(ma_fullcov$Zscore, ma_fullcov$MarkerName)

# -- 3e. Compute weighted sum --
sentinel_mrs_df$MRS_age_gender <- apply(scaled_betas, 1, function(row) {
  sum(row * weights_A[names(row)], na.rm = TRUE)
})

sentinel_mrs_df$MRS_fullcov <- apply(scaled_betas, 1, function(row) {
  sum(row * weights_B[names(row)], na.rm = TRUE)
})

# Reorder columns so MRS follows the CpG columns
last_cg <- tail(sentinel_cpgs, 1)
sentinel_mrs_df <- sentinel_mrs_df %>%
  relocate(MRS_age_gender, MRS_fullcov, .after = all_of(last_cg))

saveRDS(sentinel_mrs_df, paste0("sentinel_mrs_", ethnicity_test, "_test_constructed.rds"))
message("MRS constructed for ", nrow(sentinel_mrs_df), " test-set samples")


# ==============================================================================
# STEP 4: Define MRS quartiles (based on non-risk controls)
# ==============================================================================
#
# Quartile boundaries are derived from the non-risk (at_risk_75c == "N")
# subset so that the score distribution reflects healthy individuals.
# ==============================================================================

assign_quartiles <- function(df, mrs_col, risk_col = "at_risk_75c") {
  # Compute quartile cut-points from controls only
  controls <- df[[mrs_col]][df[[risk_col]] == "N"]
  q1 <- quantile(controls, 0.25)
  q2 <- quantile(controls, 0.50)
  q3 <- quantile(controls, 0.75)

  quartile_col <- paste0(mrs_col, "_quartiles")
  df[[quartile_col]] <- case_when(
    df[[mrs_col]] <  q1 ~ "Q1",
    df[[mrs_col]] <  q2 ~ "Q2",
    df[[mrs_col]] <  q3 ~ "Q3",
    df[[mrs_col]] >= q3 ~ "Q4"
  )
  return(df)
}

mrs <- readRDS(paste0("sentinel_mrs_", ethnicity_test, "_test_constructed.rds"))
mrs <- assign_quartiles(mrs, "MRS_age_gender")
mrs <- assign_quartiles(mrs, "MRS_fullcov")

saveRDS(mrs, paste0("sentinel_mrs_quartile_", ethnicity_test, "_test_constructed.rds"))


# ==============================================================================
# STEP 5: Evaluate discriminative performance (AUC)
# ==============================================================================
#
# Models compared (logistic regression predicting at_risk_75c):
#   Model 1 : MRS (age/gender weights) + Age + Gender
#   Model 2 : MRS (full covariate weights) + Age + Gender
#   Model 3 : Model 2 + Smoking + WBC
#   Model 4 : Model 3 + 30 control probe PCs
#   Model 5 : Age + Gender only (baseline, no MRS)
#   Model 6 : Age + Gender + Smoking (no MRS)
#
# AUCs are computed per ethnicity and overall. DeLong tests compare each
# model to Models 1, 5, and 6 as reference.
# ==============================================================================

library(pROC)
library(ggplot2)
library(tidyr)

# Load combined test-set MRS data (all ethnicities merged)
mrs_all <- readRDS("mrs_combined_test.rds")
all_cov <- colnames(mrs_all)

# -- 5a. Define model formulas --
# Ethnic-specific (no Ethnicity term)
models_ethnic <- list(
  Model_1 = "as.factor(MRS_age_gender_quartiles) + Age + as.factor(Gender)",
  Model_2 = "as.factor(MRS_fullcov_quartiles) + Age + as.factor(Gender)",
  Model_5 = "Age + as.factor(Gender)",
  Model_6 = "Age + as.factor(Gender) + as.factor(smoking)"
)

# Models 3 and 4 add WBC and PC covariates programmatically
wbc_cols_test    <- all_cov[38:43]
wbc_pc_cols_test <- all_cov[8:43]

models_ethnic[["Model_3"]] <- paste(
  "as.factor(MRS_fullcov_quartiles) + Age + as.factor(Gender) + as.factor(smoking)",
  paste(wbc_cols_test, collapse = " + "), sep = " + ")

models_ethnic[["Model_4"]] <- paste(
  "as.factor(MRS_fullcov_quartiles) + Age + as.factor(Gender) + as.factor(smoking)",
  paste(wbc_pc_cols_test, collapse = " + "), sep = " + ")

# Overall models (add Ethnicity covariate)
models_overall <- lapply(models_ethnic, function(f) paste(f, "+ as.factor(Ethnicity)"))

# -- 5b. AUC calculation helper --
calculate_auc <- function(data, formula_rhs, model_name, eth_label,
                          roc_ref1 = NULL, roc_ref5 = NULL, roc_ref6 = NULL) {
  tryCatch({
    f   <- as.formula(paste("as.factor(at_risk_75c) ~", formula_rhs))
    fit <- glm(f, data = data, family = "binomial", na.action = na.exclude)
    data$pred <- predict(fit, type = "response")

    roc_obj <- roc(as.factor(data$at_risk_75c) ~ data$pred, quiet = TRUE)
    ci_obj  <- ci.auc(roc_obj, conf.level = 0.95)

    data.frame(
      Model     = model_name,
      Ethnicity = eth_label,
      AUC       = as.numeric(roc_obj$auc),
      CI_Lower  = ci_obj[1],
      CI_Upper  = ci_obj[3],
      p_vs_mod1 = if (!is.null(roc_ref1)) roc.test(roc_obj, roc_ref1, method = "delong")$p.value else NA,
      p_vs_mod5 = if (!is.null(roc_ref5)) roc.test(roc_obj, roc_ref5, method = "delong")$p.value else NA,
      p_vs_mod6 = if (!is.null(roc_ref6)) roc.test(roc_obj, roc_ref6, method = "delong")$p.value else NA,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message("  Error in ", model_name, " for ", eth_label, ": ", e$message)
    data.frame(Model = model_name, Ethnicity = eth_label,
               AUC = NA, CI_Lower = NA, CI_Upper = NA,
               p_vs_mod1 = NA, p_vs_mod5 = NA, p_vs_mod6 = NA,
               stringsAsFactors = FALSE)
  })
}

# -- 5c. Fit reference ROC objects and compute AUCs per ethnicity --
desired_order <- c("C", "M", "I", "Overall")
results <- data.frame()

for (eth in c("C", "M", "I")) {
  sub <- mrs_all %>% filter(Ethnicity == eth)

  # Reference ROCs for DeLong comparison
  ref1 <- roc(as.factor(sub$at_risk_75c) ~ predict(glm(as.formula(paste("as.factor(at_risk_75c) ~", models_ethnic[["Model_1"]])), data = sub, family = "binomial", na.action = na.exclude), type = "response"), quiet = TRUE)
  ref5 <- roc(as.factor(sub$at_risk_75c) ~ predict(glm(as.formula(paste("as.factor(at_risk_75c) ~", models_ethnic[["Model_5"]])), data = sub, family = "binomial", na.action = na.exclude), type = "response"), quiet = TRUE)
  ref6 <- roc(as.factor(sub$at_risk_75c) ~ predict(glm(as.formula(paste("as.factor(at_risk_75c) ~", models_ethnic[["Model_6"]])), data = sub, family = "binomial", na.action = na.exclude), type = "response"), quiet = TRUE)

  for (mod_name in names(models_ethnic)) {
    results <- rbind(results, calculate_auc(sub, models_ethnic[[mod_name]], mod_name, eth, ref1, ref5, ref6))
  }
}

# Overall (all ethnicities pooled, with Ethnicity covariate)
ref1_all <- roc(as.factor(mrs_all$at_risk_75c) ~ predict(glm(as.formula(paste("as.factor(at_risk_75c) ~", models_overall[["Model_1"]])), data = mrs_all, family = "binomial", na.action = na.exclude), type = "response"), quiet = TRUE)
ref5_all <- roc(as.factor(mrs_all$at_risk_75c) ~ predict(glm(as.formula(paste("as.factor(at_risk_75c) ~", models_overall[["Model_5"]])), data = mrs_all, family = "binomial", na.action = na.exclude), type = "response"), quiet = TRUE)
ref6_all <- roc(as.factor(mrs_all$at_risk_75c) ~ predict(glm(as.formula(paste("as.factor(at_risk_75c) ~", models_overall[["Model_6"]])), data = mrs_all, family = "binomial", na.action = na.exclude), type = "response"), quiet = TRUE)

for (mod_name in names(models_overall)) {
  results <- rbind(results, calculate_auc(mrs_all, models_overall[[mod_name]], mod_name, "Overall", ref1_all, ref5_all, ref6_all))
}

results$Ethnicity <- factor(results$Ethnicity, levels = rev(desired_order))
saveRDS(results, "summary_roc_by_ethnicity.rds")

# -- 5d. AUC bar plot --
model_colors <- c(Model_1 = "#FFA500", Model_2 = "#FF0000", Model_3 = "#C80000",
                  Model_4 = "#7B0000", Model_5 = "#C6F5D6", Model_6 = "#1B8000")

p_auc <- ggplot(results, aes(x = AUC, y = Ethnicity, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge2(reverse = TRUE), width = 0.7) +
  geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper),
                position = position_dodge2(width = 2, reverse = TRUE),
                width = 0.7, color = "black") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_fill_manual(values = model_colors) +
  labs(title = "Model Performance by Ethnicity", x = "AUC (95% CI)", y = NULL) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom")

ggsave("plot_auc_bymod_ethnicity.pdf", p_auc, width = 10, height = 8)


# ==============================================================================
# STEP 6: MRS-quartile association with cIMT risk (logistic regression)
# ==============================================================================
#
# For each model, fit a logistic regression of risk status on MRS quartiles
# (Q1 = reference) + covariates. Extract odds ratios and 95% CIs for Q2–Q4.
# ==============================================================================

library(writexl)

mrs_test <- readRDS("mrs_combined_test.rds")

# -- 6a. Helper function: fit logistic model and extract quartile ORs --
extract_quartile_or <- function(data, formula_rhs, mrs_quartile_col, model_label) {
  # Count cases/controls per quartile
  risk_counts <- data %>%
    group_by(!!sym(mrs_quartile_col)) %>%
    summarise(n_risk = sum(at_risk_75c == "Y"),
              n_control = sum(at_risk_75c == "N"), .groups = "drop")

  f   <- as.formula(paste("as.factor(at_risk_75c) ~", formula_rhs))
  fit <- glm(f, data = data, family = "binomial", na.action = na.exclude)

  # Extract Q2, Q3, Q4 coefficients (rows 2–4 of the coefficient table)
  coefs <- summary(fit)$coefficients[2:4, ]

  result <- data.frame(
    MRS_quartile = paste0("Q", 2:4),
    beta         = coefs[, 1],
    se           = coefs[, 2],
    z            = coefs[, 3],
    p            = coefs[, 4]
  ) %>%
    mutate(
      odds_ratio = exp(beta),
      ci_lower   = exp(beta - 1.96 * se),
      ci_upper   = exp(beta + 1.96 * se),
      OR_95CI    = sprintf("%.2f (%.2f–%.2f)", odds_ratio, ci_lower, ci_upper),
      model      = model_label
    ) %>%
    left_join(risk_counts, by = setNames(mrs_quartile_col, "MRS_quartile"))

  return(result)
}

# -- 6b. Define models and extract ORs --
or_mod1 <- extract_quartile_or(
  mrs_test,
  "as.factor(MRS_age_gender_quartiles) + Age + as.factor(Gender) + as.factor(Ethnicity)",
  "MRS_age_gender_quartiles", "Model 1 (Age + Gender)")

or_mod2 <- extract_quartile_or(
  mrs_test,
  "as.factor(MRS_fullcov_quartiles) + Age + as.factor(Gender) + as.factor(Ethnicity)",
  "MRS_fullcov_quartiles", "Model 2 (Full covariates)")

wbc_str    <- paste(all_cov[38:43], collapse = " + ")
wbc_pc_str <- paste(all_cov[8:43], collapse = " + ")

or_mod3 <- extract_quartile_or(
  mrs_test,
  paste("as.factor(MRS_fullcov_quartiles) + Age + as.factor(Gender) + as.factor(Ethnicity) + as.factor(smoking) +", wbc_str),
  "MRS_fullcov_quartiles", "Model 3 (+WBC)")

or_mod4 <- extract_quartile_or(
  mrs_test,
  paste("as.factor(MRS_fullcov_quartiles) + Age + as.factor(Gender) + as.factor(Ethnicity) + as.factor(smoking) +", wbc_pc_str),
  "MRS_fullcov_quartiles", "Model 4 (+WBC + PCs)")

combined_or <- bind_rows(or_mod1, or_mod2, or_mod3, or_mod4)
write_xlsx(combined_or, "combined_regout_MRSquartile.xlsx")
