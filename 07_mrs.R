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
#       sentinel CpG betas in the training set (shown for Chinese; repeat
#       for Malay and Indian by changing the ethnicity filter)
#    3. Meta-analyse weights across ethnicities using METAL (Stouffer method)
#    4. Construct MRS in the test set as a weighted sum of scaled betas
#    5. Define MRS quartiles based on control (non-risk) individuals
#    6. Evaluate discriminative performance (AUC) across models
#    7. Assess MRS-quartile association with cIMT risk (logistic regression)
#
#  Inputs:
#    - phe_wcIMTref.rds :
#        Phenotype data frame [samples x variables] with columns:
#          SentrixID, Ethnicity ("C"/"M"/"I"), avgcIMT_75c, at_risk_75c
#          ("Y"/"N"), log_CIMT_distal_mean, Age, Gender, smoking
#    - sentinel_cpgs_beta.rds :
#        Matrix [sentinel CpG probes x samples] of QN beta values
#    - ctrlprobes.RData :
#        Contains `ctrl.all` [samples x control probes]
#    - housemanEstimates_QN_detPrelax_betas.RData :
#        Contains `constrainedCoefs` [samples x WBC cell types]
#    - METAL output .tbl files :
#        Meta-analysed weights (Z-scores) from Step 3
#    - mrs_combined_test.rds :
#        Combined test-set data with MRS, quartiles, and covariates
#
#  Dependencies:
#    install.packages(c("dplyr", "caret", "forcats", "pROC", "ggplot2",
#                       "tidyr", "writexl"))
#
################################################################################


# ==============================================================================
# STEP 1: Stratified train-test split
# ==============================================================================

set.seed(123)

library(dplyr)
library(caret)

# Load phenotype data and remove samples without risk classification
phe_wcIMTref <- readRDS("phe_wcIMTref.rds") %>%
  filter(!is.na(avgcIMT_75c))

# Create combined stratification variable to preserve joint proportions
# of ethnicity and risk status in both train and test sets
phe_wcIMTref$strata <- interaction(phe_wcIMTref$Ethnicity, phe_wcIMTref$at_risk_75c)
trainIndex <- createDataPartition(phe_wcIMTref$strata, p = 0.7, list = FALSE)

# Assign train/test labels
phe_wcIMTref$Set <- "test"
phe_wcIMTref$Set[trainIndex] <- "train"

# Verify proportions
ethnicity_props <- prop.table(table(phe_wcIMTref$Set, phe_wcIMTref$Ethnicity), margin = 1)
print("Ethnicity proportions:")
print(ethnicity_props)

risk_props <- prop.table(table(phe_wcIMTref$Set, phe_wcIMTref$at_risk_75c), margin = 1)
print("Risk status proportions:")
print(risk_props)

original_risk_props <- prop.table(table(phe_wcIMTref$Ethnicity, phe_wcIMTref$at_risk_75c), margin = 1)
print("Risk status proportions:")

saveRDS(phe_wcIMTref, "phe_wcIMTref_traintest.rds")


# ==============================================================================
# STEP 2: Derive MRS weights (ethnic-specific regression in training set)
# ==============================================================================
#
# Two regression models are fitted per sentinel CpG:
#   Regression 1: log(cIMT) ~ CpG + Age + Gender
#   Regression 2: log(cIMT) ~ CpG + Age + Gender + Smoking + WBC + 30 ctrl PCs
#
# This section shows the Chinese ("C") analysis. To run for Malay or Indian,
# change the Ethnicity filter from "C" to "M" or "I", update variable names
# (e.g. phe_CHN -> phe_MLY), and adjust the WBC/PC column indices as needed
# (these may differ slightly across ethnicities due to merging order).
# ==============================================================================

library(forcats)

# Load shared data
load("ctrlprobes.RData")                              # ctrl.all [samples x control probes]
load("housemanEstimates_QN_detPrelax_betas.RData")    # constrainedCoefs [samples x cell types]
sentinel_cpgs_beta <- readRDS("sentinel_cpgs_beta.rds") # [sentinel CpGs x samples]

# ------------------------------------------------------------------------------
# 2a. Chinese
# ------------------------------------------------------------------------------

# Select ethnic-specific training samples
phe_CHN <- readRDS("phe_wcIMTref_traintest.rds") %>%
  filter(Set == "train") %>%
  filter(Ethnicity == "C")

# Ethnic-specific control probe PCs
ctrl.all_CHN <- ctrl.all[rownames(ctrl.all) %in% phe_CHN[, 4], ]
dim(ctrl.all_CHN)
pca_cp_CHN <- prcomp(na.omit(ctrl.all_CHN))
PC_cp_CHN <- predict(pca_cp_CHN)
saveRDS(PC_cp_CHN, "PC_cp_CHN_train.rds")

# Subset methylation to ethnic training samples
sentinel_cpgs_beta_CHN <- sentinel_cpgs_beta[, intersect(colnames(sentinel_cpgs_beta), phe_CHN$SentrixID)]

# Format control probe PCs for merging
PC_cp_CHN <- as.data.frame(PC_cp_CHN)
PC_cp_CHN$SentrixID <- row.names(PC_cp_CHN)

# Format WBC proportions for merging
constrainedCoefs <- as.data.frame(constrainedCoefs)
constrainedCoefs$SentrixID <- row.names(constrainedCoefs)

# Create master covariate table
phe_CHN_all_variables <- merge(phe_CHN, PC_cp_CHN, by = "SentrixID")
phe_CHN_all_variables <- merge(phe_CHN_all_variables, constrainedCoefs, by = "SentrixID")

# Align beta columns to covariate rows
beta_order <- phe_CHN_all_variables$SentrixID
sentinel_cpgs_beta_CHN <- sentinel_cpgs_beta_CHN[, match(beta_order, colnames(sentinel_cpgs_beta_CHN))]
colnames(sentinel_cpgs_beta_CHN) == phe_CHN_all_variables$SentrixID

# All covariate column names (for index-based selection below)
all_cov <- colnames(phe_CHN_all_variables)

# --- Regression 1: Age + Gender only -----------------------------------------

RHS_equation <- '+Age+as.factor(Gender)'

regout_CHN_lncIMT_AgeGenderAdj <- matrix(NA, nrow(sentinel_cpgs_beta_CHN), 5)
colnames(regout_CHN_lncIMT_AgeGenderAdj) <- c("Estimate", "Std.Error", "t.value", "Pr(>|t|)", "N_eff")

for (i in 1:nrow(sentinel_cpgs_beta_CHN)) {
  reg_formula <- as.formula(paste(
    "log_CIMT_distal_mean~", "sentinel_cpgs_beta_CHN[i,]", RHS_equation, sep = ""
  ))
  fit <- lm(reg_formula, data = phe_CHN_all_variables, na.action = na.exclude)
  regout_CHN_lncIMT_AgeGenderAdj[i, ] <- c(summary(fit)$coefficient[2, ], nobs(fit))
}

rownames(regout_CHN_lncIMT_AgeGenderAdj) <- rownames(sentinel_cpgs_beta_CHN)
saveRDS(regout_CHN_lncIMT_AgeGenderAdj, "regout_lncIMTmean_CHN_AgeGenderAdj_sentinel.rds")

# --- Regression 2: Age + Gender + Smoking + WBC + 30 control probe PCs -------

# NOTE: Column indices for WBC and control probe PCs depend on the merged
# covariate table structure. Inspect `all_cov` and adjust indices accordingly.
# For Chinese: columns 27:56 = control probe PCs, columns 231:236 = WBC
cov_wbc_pc <- all_cov[c(27:56, 231:236)]
print(cov_wbc_pc)

# Build the right-hand side of the regression formula
RHS_equation_cov_wbc_pc <- '+Age+as.factor(Gender)+as.factor(smoking)'
for (element_loop in cov_wbc_pc) {
  RHS_equation_cov_wbc_pc <- paste(RHS_equation_cov_wbc_pc, '+', element_loop, sep = "")
}
print(RHS_equation_cov_wbc_pc)

# Run regression
regout_CHN_lncIMT_AgeGenderSmokWBCPCAdj <- matrix(0, nrow(sentinel_cpgs_beta_CHN), 5)
colnames(regout_CHN_lncIMT_AgeGenderSmokWBCPCAdj) <- c("Estimate", "Std.Error", "t.value", "Pr(>|t|)", "N_eff")

for (i in 1:nrow(sentinel_cpgs_beta_CHN)) {
  reg_formula <- as.formula(paste(
    "log_CIMT_distal_mean~", "sentinel_cpgs_beta_CHN[i,]", RHS_equation_cov_wbc_pc, sep = ""
  ))
  fit <- lm(reg_formula, data = phe_CHN_all_variables, na.action = na.exclude)
  regout_CHN_lncIMT_AgeGenderSmokWBCPCAdj[i, ] <- c(summary(fit)$coefficient[2, ], nobs(fit))
}

rownames(regout_CHN_lncIMT_AgeGenderSmokWBCPCAdj) <- rownames(sentinel_cpgs_beta_CHN)
print(reg_formula)
saveRDS(regout_CHN_lncIMT_AgeGenderSmokWBCPCAdj, "regout_lncIMTmean_CHN_AgeGenderSmokWBCPcAdj_sentinel.rds")

# ------------------------------------------------------------------------------
# 2b. Repeat for Malay (Ethnicity == "M") and Indian (Ethnicity == "I")
# ------------------------------------------------------------------------------
# The code is identical to Section 2a above. Change:
#   - Ethnicity filter: "M" or "I"
#   - Variable names: phe_MLY / phe_IND, ctrl.all_MLY / ctrl.all_IND, etc.
#   - WBC/PC column indices: inspect `all_cov` after merging for each ethnicity
#     e.g. Malay: all_cov[c(27:56, 131:136)], Indian: all_cov[c(27:56, 128:133)]
#   - Output filenames: replace "CHN" with "MLY" or "IND"
# ------------------------------------------------------------------------------


# ==============================================================================
# STEP 3: Meta-analyse weights across ethnicities (METAL)
# ==============================================================================
#
# After running Step 2 for all three ethnicities and formatting the output
# files as tab-separated text (with columns: CpG, dir_beta, p, N), run METAL
# with the script below. Save as e.g. "mrs_meta_agegender.metal" and execute:
#   metal mrs_meta_agegender.metal
#
# Repeat with the full-covariate regression outputs for a second set of weights.
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
# STEP 4: Construct MRS in the test set
# ==============================================================================
#
# MRS = sum( scaled_beta_CpG_i * weight_i )
# Weights are Stouffer Z-scores from the METAL meta-analysis.
# Betas are z-score scaled per CpG within the test set.
#
# Run this section once per ethnicity by changing `ethnicity` below.
# ==============================================================================

library(dplyr)

# Specify ethnicity for calculating MRS
ethnicity <- c("CHN")  # change to "MLY" or "IND" as needed

# File paths
ethnic_variables       <- paste0("sentinel_mrs_", ethnicity, "_test_df.rds")
ethnic_mrs_constructed <- paste0("sentinel_mrs_", ethnicity, "_test_constructed.rds")

# Read ethnic-specific test-set variables
sentinel_mrs_df <- readRDS(ethnic_variables)

# Read meta-analysed weights (Z-scores from METAL)
ma_agegender <- read.delim("MA_sentinel_ASIAN_HELIOS_CHNMLYIND_stouffer_AgeGender1.tbl")
ma_fullcov   <- read.delim("MA_sentinel_ASIAN_HELIOS_CHNMLYIND_stouffer_AgeGenderSmokWBCPc1.tbl")

# Identify columns storing sentinel CpG methylation betas
sentinel_cpgs <- grep("^cg", colnames(sentinel_mrs_df), value = TRUE)

# Scale betas (z-score per CpG)
scaled_sentinel_mrs_df <- sentinel_mrs_df[, sentinel_cpgs] %>%
  mutate(across(everything(), scale))

# Create named weight vectors from meta-analysis Z-scores
weights_age_gender <- setNames(ma_agegender$Zscore, ma_agegender$MarkerName)
weights_fullcov    <- setNames(ma_fullcov$Zscore, ma_fullcov$MarkerName)

# Compute weighted sum for each sample
sentinel_mrs_df$MRS_age_gender <- apply(scaled_sentinel_mrs_df[, sentinel_cpgs], 1, function(row) {
  sum(row * weights_age_gender[names(row)], na.rm = TRUE)
})

sentinel_mrs_df$MRS_fullcov <- apply(scaled_sentinel_mrs_df[, sentinel_cpgs], 1, function(row) {
  sum(row * weights_fullcov[names(row)], na.rm = TRUE)
})

# Relocate MRS columns to after the last CpG column
last_cg_column <- tail(grep("^cg", colnames(sentinel_mrs_df), value = TRUE), 1)

sentinel_mrs_constructed <- sentinel_mrs_df %>%
  relocate(MRS_age_gender, MRS_fullcov, .after = all_of(last_cg_column))

saveRDS(sentinel_mrs_constructed, ethnic_mrs_constructed)


# ==============================================================================
# STEP 5: Define MRS quartiles (based on non-risk controls)
# ==============================================================================
#
# Quartile boundaries are derived from the non-risk (at_risk_75c == "N")
# subset so that the score distribution reflects healthy individuals.
#
# Run this section once per ethnicity by changing `ethnicity` below.
# ==============================================================================

library(dplyr)

ethnicity <- c("CHN")  # change to "MLY" or "IND" as needed

mrs_test_file              <- paste0("sentinel_mrs_", ethnicity, "_test_constructed.rds")
mrs_wquartile_test_outfile <- paste0("sentinel_mrs_quartile_", ethnicity, "_test_constructed.rds")

mrs <- readRDS(mrs_test_file)

# Define quartiles for MRS_age_gender (based on control distribution)
mrs_wquartile <- mrs %>%
  mutate(
    Q1_MRS_age_gender = quantile(MRS_age_gender[at_risk_75c == "N"], 0.25),
    Q2_MRS_age_gender = quantile(MRS_age_gender[at_risk_75c == "N"], 0.50),
    Q3_MRS_age_gender = quantile(MRS_age_gender[at_risk_75c == "N"], 0.75)
  ) %>%
  mutate(MRS_age_gender_quartiles = case_when(
    MRS_age_gender <  Q1_MRS_age_gender ~ "Q1",
    MRS_age_gender >= Q1_MRS_age_gender & MRS_age_gender < Q2_MRS_age_gender ~ "Q2",
    MRS_age_gender >= Q2_MRS_age_gender & MRS_age_gender < Q3_MRS_age_gender ~ "Q3",
    MRS_age_gender >= Q3_MRS_age_gender ~ "Q4"
  )) %>%
  select(-Q1_MRS_age_gender, -Q2_MRS_age_gender, -Q3_MRS_age_gender)

# Define quartiles for MRS_fullcov (based on control distribution)
mrs_wquartile_full <- mrs_wquartile %>%
  mutate(
    Q1_MRS_fullcov = quantile(MRS_fullcov[at_risk_75c == "N"], 0.25),
    Q2_MRS_fullcov = quantile(MRS_fullcov[at_risk_75c == "N"], 0.50),
    Q3_MRS_fullcov = quantile(MRS_fullcov[at_risk_75c == "N"], 0.75)
  ) %>%
  mutate(MRS_fullcov_quartiles = case_when(
    MRS_fullcov <  Q1_MRS_fullcov ~ "Q1",
    MRS_fullcov >= Q1_MRS_fullcov & MRS_fullcov < Q2_MRS_fullcov ~ "Q2",
    MRS_fullcov >= Q2_MRS_fullcov & MRS_fullcov < Q3_MRS_fullcov ~ "Q3",
    MRS_fullcov >= Q3_MRS_fullcov ~ "Q4"
  )) %>%
  select(-Q1_MRS_fullcov, -Q2_MRS_fullcov, -Q3_MRS_fullcov)

saveRDS(mrs_wquartile_full, mrs_wquartile_test_outfile)


# ==============================================================================
# STEP 6: Evaluate discriminative performance (AUC)
# ==============================================================================
#
# Models compared (logistic regression predicting at_risk_75c):
#   Model 1 : MRS (age/gender weights quartiles) + Age + Gender
#   Model 2 : MRS (full-covariate weights quartiles) + Age + Gender
#   Model 3 : Model 2 + Smoking + WBC
#   Model 4 : Model 3 + 30 control probe PCs
#   Model 5 : Age + Gender only (baseline, no MRS)
#   Model 6 : Age + Gender + Smoking (no MRS)
#
# Ethnic-specific models omit Ethnicity; overall models include it.
# DeLong tests compare each model to Models 1, 5, and 6 as reference.
# ==============================================================================

library(pROC)
library(dplyr)
library(ggplot2)
library(tidyr)

# Load combined test-set MRS data (all ethnicities merged)
mrs_all <- readRDS("mrs_combined_test.rds")

all_cov <- colnames(mrs_all)

# --- Define model formulas (ethnic-specific, without Ethnicity term) ----------

mod_1_RHS_equation_base <- 'as.factor(MRS_age_gender_quartiles) + Age + as.factor(Gender)'
mod_2_RHS_equation_base <- 'as.factor(MRS_fullcov_quartiles) + Age + as.factor(Gender)'
mod_5_RHS_equation_base <- 'Age + as.factor(Gender)'
mod_6_RHS_equation_base <- 'Age + as.factor(Gender) + as.factor(smoking)'

# Model 3: + smoking + WBC
cov_wbc <- c(all_cov[38:43])
mod_3_RHS_equation_base <- 'as.factor(MRS_fullcov_quartiles) + Age + as.factor(Gender) + as.factor(smoking)'
for (element_loop in cov_wbc) {
  mod_3_RHS_equation_base <- paste(mod_3_RHS_equation_base, '+', element_loop, sep = "")
}

# Model 4: + smoking + WBC + control probe PCs
cov_wbc_pc <- c(all_cov[8:43])
mod_4_RHS_equation_base <- 'as.factor(MRS_fullcov_quartiles) + Age + as.factor(Gender) + as.factor(smoking)'
for (element_loop in cov_wbc_pc) {
  mod_4_RHS_equation_base <- paste(mod_4_RHS_equation_base, '+', element_loop, sep = "")
}

# --- Overall formulas (add Ethnicity covariate) -------------------------------
mod_1_RHS_equation_overall <- paste(mod_1_RHS_equation_base, "+ as.factor(Ethnicity)")
mod_2_RHS_equation_overall <- paste(mod_2_RHS_equation_base, "+ as.factor(Ethnicity)")
mod_3_RHS_equation_overall <- paste(mod_3_RHS_equation_base, "+ as.factor(Ethnicity)")
mod_4_RHS_equation_overall <- paste(mod_4_RHS_equation_base, "+ as.factor(Ethnicity)")
mod_5_RHS_equation_overall <- paste(mod_5_RHS_equation_base, "+ as.factor(Ethnicity)")
mod_6_RHS_equation_overall <- paste(mod_6_RHS_equation_base, "+ as.factor(Ethnicity)")

# --- Plotting parameters -----------------------------------------------------
model_colors <- c(
  "Model_1" = "#FFA500",  # Orange
  "Model_2" = "#FF0000",  # Red
  "Model_3" = "#C80000",  # Darker Red
  "Model_4" = "#7B0000",  # Very Dark Red
  "Model_5" = "#C6F5D6",  # Light Green
  "Model_6" = "#1B8000"   # Dark Green
)

desired_order <- c("C", "M", "I", "Overall")

# --- Initialise results storage -----------------------------------------------
results <- data.frame(
  Model = character(), Ethnicity = character(),
  AUC = numeric(), CI_Lower = numeric(), CI_Upper = numeric(),
  p_roc_vs_mod1 = numeric(), p_roc_vs_mod5 = numeric(), p_roc_vs_mod6 = numeric(),
  stringsAsFactors = FALSE
)

# --- AUC calculation function -------------------------------------------------
calculate_auc <- function(data, model_formula, model_name, ethnicity_label,
                          roc_mod1 = NULL, roc_mod5 = NULL, roc_mod6 = NULL) {
  tryCatch({
    reg_formula <- as.formula(paste("as.factor(at_risk_75c)", "~", model_formula))
    model_glm <- glm(reg_formula, data = data, family = "binomial", na.action = na.exclude)
    data$predicted_prob <- predict(model_glm, type = "response")

    mod_roc   <- roc(as.factor(data$at_risk_75c) ~ data$predicted_prob, quiet = TRUE)
    ci_result <- ci.auc(mod_roc, conf.level = 0.95)

    p_roc_vs_mod1 <- if (!is.null(roc_mod1)) roc.test(mod_roc, roc_mod1, method = "delong")$p.value else NA
    p_roc_vs_mod5 <- if (!is.null(roc_mod5)) roc.test(mod_roc, roc_mod5, method = "delong")$p.value else NA
    p_roc_vs_mod6 <- if (!is.null(roc_mod6)) roc.test(mod_roc, roc_mod6, method = "delong")$p.value else NA

    data.frame(
      Model = model_name, Ethnicity = ethnicity_label,
      AUC = as.numeric(mod_roc$auc), CI_Lower = ci_result[1], CI_Upper = ci_result[3],
      p_roc_vs_mod1 = p_roc_vs_mod1, p_roc_vs_mod5 = p_roc_vs_mod5, p_roc_vs_mod6 = p_roc_vs_mod6,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat("Error in model", model_name, "for ethnicity", ethnicity_label, ":", e$message, "\n")
    data.frame(
      Model = model_name, Ethnicity = ethnicity_label,
      AUC = NA, CI_Lower = NA, CI_Upper = NA,
      p_roc_vs_mod1 = NA, p_roc_vs_mod5 = NA, p_roc_vs_mod6 = NA,
      stringsAsFactors = FALSE
    )
  })
}

# --- Fit reference ROC objects (for DeLong comparisons) -----------------------

# Overall reference ROCs
roc_mod1_overall <- roc(as.factor(mrs_all$at_risk_75c) ~ predict(glm(as.formula(paste("as.factor(at_risk_75c)", "~", mod_1_RHS_equation_overall)), data = mrs_all, family = "binomial", na.action = na.exclude), type = "response"), quiet = TRUE)
roc_mod5_overall <- roc(as.factor(mrs_all$at_risk_75c) ~ predict(glm(as.formula(paste("as.factor(at_risk_75c)", "~", mod_5_RHS_equation_overall)), data = mrs_all, family = "binomial", na.action = na.exclude), type = "response"), quiet = TRUE)
roc_mod6_overall <- roc(as.factor(mrs_all$at_risk_75c) ~ predict(glm(as.formula(paste("as.factor(at_risk_75c)", "~", mod_6_RHS_equation_overall)), data = mrs_all, family = "binomial", na.action = na.exclude), type = "response"), quiet = TRUE)

# Ethnic-specific reference ROCs
roc_mod1_ethnicities <- list()
roc_mod5_ethnicities <- list()
roc_mod6_ethnicities <- list()

for (ethnicity in desired_order[!desired_order %in% "Overall"]) {
  subset_data <- mrs_all %>% filter(Ethnicity == ethnicity)
  roc_mod1_ethnicities[[ethnicity]] <- roc(as.factor(subset_data$at_risk_75c) ~ predict(glm(as.formula(paste("as.factor(at_risk_75c)", "~", mod_1_RHS_equation_base)), data = subset_data, family = "binomial", na.action = na.exclude), type = "response"), quiet = TRUE)
  roc_mod5_ethnicities[[ethnicity]] <- roc(as.factor(subset_data$at_risk_75c) ~ predict(glm(as.formula(paste("as.factor(at_risk_75c)", "~", mod_5_RHS_equation_base)), data = subset_data, family = "binomial", na.action = na.exclude), type = "response"), quiet = TRUE)
  roc_mod6_ethnicities[[ethnicity]] <- roc(as.factor(subset_data$at_risk_75c) ~ predict(glm(as.formula(paste("as.factor(at_risk_75c)", "~", mod_6_RHS_equation_base)), data = subset_data, family = "binomial", na.action = na.exclude), type = "response"), quiet = TRUE)
}

# --- Compute AUCs per ethnicity -----------------------------------------------

# Standardise ethnicity labels if needed
mrs_all$Ethnicity[mrs_all$Ethnicity == "Chinese"] <- "C"
mrs_all$Ethnicity[mrs_all$Ethnicity == "Malay"]   <- "M"
mrs_all$Ethnicity[mrs_all$Ethnicity == "Indian"]  <- "I"

for (ethnicity in desired_order[!desired_order %in% "Overall"]) {
  subset_data <- mrs_all %>% filter(Ethnicity == ethnicity)

  results <- rbind(results, calculate_auc(subset_data, mod_1_RHS_equation_base, "Model_1", ethnicity, roc_mod1_ethnicities[[ethnicity]], roc_mod5_ethnicities[[ethnicity]], roc_mod6_ethnicities[[ethnicity]]))
  results <- rbind(results, calculate_auc(subset_data, mod_2_RHS_equation_base, "Model_2", ethnicity, roc_mod1_ethnicities[[ethnicity]], roc_mod5_ethnicities[[ethnicity]], roc_mod6_ethnicities[[ethnicity]]))
  results <- rbind(results, calculate_auc(subset_data, mod_3_RHS_equation_base, "Model_3", ethnicity, roc_mod1_ethnicities[[ethnicity]], roc_mod5_ethnicities[[ethnicity]], roc_mod6_ethnicities[[ethnicity]]))
  results <- rbind(results, calculate_auc(subset_data, mod_4_RHS_equation_base, "Model_4", ethnicity, roc_mod1_ethnicities[[ethnicity]], roc_mod5_ethnicities[[ethnicity]], roc_mod6_ethnicities[[ethnicity]]))
  results <- rbind(results, calculate_auc(subset_data, mod_5_RHS_equation_base, "Model_5", ethnicity, roc_mod1_ethnicities[[ethnicity]], roc_mod5_ethnicities[[ethnicity]], roc_mod6_ethnicities[[ethnicity]]))
  results <- rbind(results, calculate_auc(subset_data, mod_6_RHS_equation_base, "Model_6", ethnicity, roc_mod1_ethnicities[[ethnicity]], roc_mod5_ethnicities[[ethnicity]], roc_mod6_ethnicities[[ethnicity]]))
}

# --- Compute AUCs overall (all ethnicities pooled, with Ethnicity covariate) --

results <- rbind(results, calculate_auc(mrs_all, mod_1_RHS_equation_overall, "Model_1", "Overall", roc_mod1_overall, roc_mod5_overall, roc_mod6_overall))
results <- rbind(results, calculate_auc(mrs_all, mod_2_RHS_equation_overall, "Model_2", "Overall", roc_mod1_overall, roc_mod5_overall, roc_mod6_overall))
results <- rbind(results, calculate_auc(mrs_all, mod_3_RHS_equation_overall, "Model_3", "Overall", roc_mod1_overall, roc_mod5_overall, roc_mod6_overall))
results <- rbind(results, calculate_auc(mrs_all, mod_4_RHS_equation_overall, "Model_4", "Overall", roc_mod1_overall, roc_mod5_overall, roc_mod6_overall))
results <- rbind(results, calculate_auc(mrs_all, mod_5_RHS_equation_overall, "Model_5", "Overall", roc_mod1_overall, roc_mod5_overall, roc_mod6_overall))
results <- rbind(results, calculate_auc(mrs_all, mod_6_RHS_equation_overall, "Model_6", "Overall", roc_mod1_overall, roc_mod5_overall, roc_mod6_overall))

results$Ethnicity <- factor(results$Ethnicity, levels = rev(desired_order))

# --- Reshape for downstream use -----------------------------------------------
results_long <- results %>%
  tidyr::pivot_longer(
    cols = starts_with("p_roc_vs_"),
    names_to = "Comparison", values_to = "P_Value",
    names_prefix = "p_roc_vs_"
  )

bonferroni_alpha <- 0.05 / length(unique(results_long$Comparison))

# --- Plot: All 6 models -------------------------------------------------------
ggplot(results, aes(x = AUC, y = Ethnicity, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge2(reverse = TRUE), width = 0.7) +
  geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper),
                position = position_dodge2(width = 2, reverse = TRUE),
                width = 0.7, color = "black") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_fill_manual(values = model_colors) +
  labs(title = "Model Performance by Ethnicity and Overall",
       x = "AUC (95% CI)", y = "Ethnicity Group", fill = "Model") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10), legend.position = "bottom")

ggsave("plot_auc_bymod_ethnicity.pdf", width = 10, height = 8)

print(results)
saveRDS(results, "summary_roc_by_ethnicity.rds")

# --- Plot: Models 1, 5, 6 only ------------------------------------------------
models_to_plot <- c("Model_1", "Model_5", "Model_6")

results_156 <- results %>%
  filter(Model %in% models_to_plot) %>%
  mutate(Ethnicity = factor(Ethnicity, levels = rev(desired_order)),
         Model = factor(Model, levels = models_to_plot))

model_colors_156 <- model_colors[models_to_plot]

p_156 <- ggplot(results_156, aes(x = AUC, y = Ethnicity, fill = Model)) +
  geom_col(position = position_dodge2(reverse = TRUE), width = 0.7) +
  geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper),
                position = position_dodge2(width = 2, reverse = TRUE),
                width = 0.7, color = "black") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_fill_manual(values = model_colors_156) +
  labs(title = "Model Performance (Models 1, 5, 6) by Ethnicity and Overall",
       x = "AUC (95% CI)", y = "Ethnicity Group", fill = "Model") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10), legend.position = "bottom")

print(p_156)

# --- Plot: Models 2, 3, 4, 5, 6 only ------------------------------------------
models_to_plot <- c("Model_2", "Model_3", "Model_4", "Model_5", "Model_6")

results_23456 <- results %>%
  filter(Model %in% models_to_plot) %>%
  mutate(Ethnicity = factor(Ethnicity, levels = rev(desired_order)),
         Model = factor(Model, levels = models_to_plot))

model_colors_23456 <- model_colors[models_to_plot]

p_23456 <- ggplot(results_23456, aes(x = AUC, y = Ethnicity, fill = Model)) +
  geom_col(position = position_dodge2(reverse = TRUE), width = 0.7) +
  geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper),
                position = position_dodge2(width = 2, reverse = TRUE),
                width = 0.7, color = "black") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_fill_manual(values = model_colors_23456) +
  labs(title = "Model Performance (Models 2-6) by Ethnicity and Overall",
       x = "AUC (95% CI)", y = "Ethnicity Group", fill = "Model") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10), legend.position = "bottom")

print(p_23456)
ggsave("plot_auc_models_2_3_4_5_6_bymod_ethnicity.pdf", plot = p_23456, width = 10, height = 8)


# ==============================================================================
# STEP 7: MRS-quartile association with cIMT risk (logistic regression)
# ==============================================================================
#
# For each model, fit a logistic regression of risk status on MRS quartiles
# (Q1 = reference) + covariates. Extract odds ratios and 95% CIs for Q2-Q4.
# Note: exp(beta) from logistic regression yields odds ratios.
# ==============================================================================

library(writexl)

mrs_combined_test <- readRDS("mrs_combined_test.rds")
all_cov <- colnames(mrs_combined_test)

# ------------------------------------------------------------------------------
# 7a. Model 1: at_risk_75c ~ MRS_age_gender_quartiles + Age + Gender + Ethnicity
# ------------------------------------------------------------------------------

risk_by_mrsquartile <- mrs_combined_test %>%
  group_by(MRS_age_gender_quartiles) %>%
  summarise(count_at_risk = sum(at_risk_75c == "Y"),
            count_non_risk = sum(at_risk_75c == "N"))

mod_1_RHS_equation <- 'as.factor(MRS_age_gender_quartiles) + Age + as.factor(Gender) + as.factor(Ethnicity)'
reg_formula <- as.formula(paste("as.factor(at_risk_75c)", "~", mod_1_RHS_equation))
model_glm <- glm(reg_formula, data = mrs_combined_test, family = "binomial", na.action = na.exclude)
coef_matrix <- summary(model_glm)$coefficients[2:4, ]

regout_mod1 <- data.frame(
  MRS_age_gender_quartiles = paste0("Q", 2:4),
  beta = coef_matrix[, 1], se = coef_matrix[, 2],
  t = coef_matrix[, 3], p = coef_matrix[, 4]
) %>%
  mutate(
    relative_risk = exp(beta),
    ci_lower = exp(beta - 1.96 * se), ci_upper = exp(beta + 1.96 * se),
    ci_95 = sprintf("%.2f (%.2f-%.2f)", relative_risk, ci_lower, ci_upper)
  ) %>%
  rename(`Relative Risk (95% CI)` = ci_95, `P-value` = p) %>%
  left_join(risk_by_mrsquartile, by = "MRS_age_gender_quartiles") %>%
  relocate(count_at_risk, count_non_risk, .after = MRS_age_gender_quartiles)

rownames(regout_mod1) <- c("mod_1_Q2", "mod_1_Q3", "mod_1_Q4")

ggplot(regout_mod1, aes(x = relative_risk, y = factor(MRS_age_gender_quartiles, levels = c("Q4", "Q3", "Q2")))) +
  geom_bar(stat = "identity", fill = "#FFA500", width = 0.5) +
  labs(x = "Relative Risk", y = "MRS Quartile", title = "Relative Risk by MRS Quartile (Model 1)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(), panel.spacing.y = unit(0.02, "lines")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) + scale_y_discrete()

# ------------------------------------------------------------------------------
# 7b. Model 2: at_risk_75c ~ MRS_fullcov_quartiles + Age + Gender + Ethnicity
# ------------------------------------------------------------------------------

risk_by_mrsquartile <- mrs_combined_test %>%
  group_by(MRS_fullcov_quartiles) %>%
  summarise(count_at_risk = sum(at_risk_75c == "Y"),
            count_non_risk = sum(at_risk_75c == "N"))

mod_2_RHS_equation <- 'as.factor(MRS_fullcov_quartiles) + Age + as.factor(Gender) + as.factor(Ethnicity)'
reg_formula <- as.formula(paste("as.factor(at_risk_75c)", "~", mod_2_RHS_equation))
model_glm <- glm(reg_formula, data = mrs_combined_test, family = "binomial", na.action = na.exclude)
coef_matrix <- summary(model_glm)$coefficients[2:4, ]

regout_mod2 <- data.frame(
  MRS_fullcov_quartiles = paste0("Q", 2:4),
  beta = coef_matrix[, 1], se = coef_matrix[, 2],
  t = coef_matrix[, 3], p = coef_matrix[, 4]
) %>%
  mutate(
    relative_risk = exp(beta),
    ci_lower = exp(beta - 1.96 * se), ci_upper = exp(beta + 1.96 * se),
    ci_95 = sprintf("%.2f (%.2f-%.2f)", relative_risk, ci_lower, ci_upper)
  ) %>%
  rename(`Relative Risk (95% CI)` = ci_95, `P-value` = p) %>%
  left_join(risk_by_mrsquartile, by = "MRS_fullcov_quartiles") %>%
  relocate(count_at_risk, count_non_risk, .after = MRS_fullcov_quartiles)

rownames(regout_mod2) <- c("mod_2_Q2", "mod_2_Q3", "mod_2_Q4")

ggplot(regout_mod2, aes(x = relative_risk, y = factor(MRS_fullcov_quartiles, levels = c("Q4", "Q3", "Q2")))) +
  geom_bar(stat = "identity", fill = "#FF0000", width = 0.5) +
  labs(x = "Relative Risk", y = "MRS Quartile", title = "Relative Risk by MRS Quartile (Model 2)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(), panel.spacing.y = unit(0.02, "lines")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) + scale_y_discrete()

# ------------------------------------------------------------------------------
# 7c. Model 3: + Smoking + WBC
# ------------------------------------------------------------------------------

risk_by_mrsquartile <- mrs_combined_test %>%
  group_by(MRS_fullcov_quartiles) %>%
  summarise(count_at_risk = sum(at_risk_75c == "Y"),
            count_non_risk = sum(at_risk_75c == "N"))

cov_wbc <- c(all_cov[38:43])
mod_3_RHS_equation <- 'as.factor(MRS_fullcov_quartiles) + Age + as.factor(Gender) + as.factor(Ethnicity) + as.factor(smoking)'
for (element_loop in cov_wbc) {
  mod_3_RHS_equation <- paste(mod_3_RHS_equation, '+', element_loop, sep = "")
}

reg_formula <- as.formula(paste("as.factor(at_risk_75c)", "~", mod_3_RHS_equation))
model_glm <- glm(reg_formula, data = mrs_combined_test, family = "binomial", na.action = na.exclude)
coef_matrix <- summary(model_glm)$coefficients[2:4, ]

regout_mod3 <- data.frame(
  MRS_fullcov_quartiles = paste0("Q", 2:4),
  beta = coef_matrix[, 1], se = coef_matrix[, 2],
  t = coef_matrix[, 3], p = coef_matrix[, 4]
) %>%
  mutate(
    relative_risk = exp(beta),
    ci_lower = exp(beta - 1.96 * se), ci_upper = exp(beta + 1.96 * se),
    ci_95 = sprintf("%.2f (%.2f-%.2f)", relative_risk, ci_lower, ci_upper)
  ) %>%
  rename(`Relative Risk (95% CI)` = ci_95, `P-value` = p) %>%
  left_join(risk_by_mrsquartile, by = "MRS_fullcov_quartiles") %>%
  relocate(count_at_risk, count_non_risk, .after = MRS_fullcov_quartiles)

rownames(regout_mod3) <- c("mod_3_Q2", "mod_3_Q3", "mod_3_Q4")

ggplot(regout_mod3, aes(x = relative_risk, y = factor(MRS_fullcov_quartiles, levels = c("Q4", "Q3", "Q2")))) +
  geom_bar(stat = "identity", fill = "#C80000", width = 0.5) +
  labs(x = "Relative Risk", y = "MRS Quartile", title = "Relative Risk by MRS Quartile (Model 3)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(), panel.spacing.y = unit(0.02, "lines")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) + scale_y_discrete()

# ------------------------------------------------------------------------------
# 7d. Model 4: + Smoking + WBC + Control probe PCs
# ------------------------------------------------------------------------------

risk_by_mrsquartile <- mrs_combined_test %>%
  group_by(MRS_fullcov_quartiles) %>%
  summarise(count_at_risk = sum(at_risk_75c == "Y"),
            count_non_risk = sum(at_risk_75c == "N"))

cov_wbc_pc <- c(all_cov[8:43])
mod_4_RHS_equation <- 'as.factor(MRS_fullcov_quartiles) + Age + as.factor(Gender) + as.factor(Ethnicity) + as.factor(smoking)'
for (element_loop in cov_wbc_pc) {
  mod_4_RHS_equation <- paste(mod_4_RHS_equation, '+', element_loop, sep = "")
}

reg_formula <- as.formula(paste("as.factor(at_risk_75c)", "~", mod_4_RHS_equation))
model_glm <- glm(reg_formula, data = mrs_combined_test, family = "binomial", na.action = na.exclude)
coef_matrix <- summary(model_glm)$coefficients[2:4, ]

regout_mod4 <- data.frame(
  MRS_fullcov_quartiles = paste0("Q", 2:4),
  beta = coef_matrix[, 1], se = coef_matrix[, 2],
  t = coef_matrix[, 3], p = coef_matrix[, 4]
) %>%
  mutate(
    relative_risk = exp(beta),
    ci_lower = exp(beta - 1.96 * se), ci_upper = exp(beta + 1.96 * se),
    ci_95 = sprintf("%.2f (%.2f-%.2f)", relative_risk, ci_lower, ci_upper)
  ) %>%
  rename(`Relative Risk (95% CI)` = ci_95, `P-value` = p) %>%
  left_join(risk_by_mrsquartile, by = "MRS_fullcov_quartiles") %>%
  relocate(count_at_risk, count_non_risk, .after = MRS_fullcov_quartiles)

rownames(regout_mod4) <- c("mod_4_Q2", "mod_4_Q3", "mod_4_Q4")

ggplot(regout_mod4, aes(x = relative_risk, y = factor(MRS_fullcov_quartiles, levels = c("Q4", "Q3", "Q2")))) +
  geom_bar(stat = "identity", fill = "#7B0000", width = 0.5) +
  labs(x = "Relative Risk", y = "MRS Quartile", title = "Relative Risk by MRS Quartile (Model 4)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(), panel.spacing.y = unit(0.02, "lines")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) + scale_y_discrete()

# ------------------------------------------------------------------------------
# 7e. Combine and export all quartile results
# ------------------------------------------------------------------------------

colnames(regout_mod1)[1] <- "MRS_quartile"
colnames(regout_mod2)[1] <- "MRS_quartile"
colnames(regout_mod3)[1] <- "MRS_quartile"
colnames(regout_mod4)[1] <- "MRS_quartile"

combined_regout_MRSquartile <- rbind(regout_mod1, regout_mod2, regout_mod3, regout_mod4)
print(combined_regout_MRSquartile)
write_xlsx(combined_regout_MRSquartile, "combined_regout_MRSquartile.xlsx")

# ------------------------------------------------------------------------------
# 7f. Combined horizontal plot (all models x quartiles)
# ------------------------------------------------------------------------------

regout_mod1$model <- "Model 1 (Basic)"
regout_mod2$model <- "Model 2 (+Covariates)"
regout_mod3$model <- "Model 3 (+WBC)"
regout_mod4$model <- "Model 4 (+PCs)"
combined_plot_data <- rbind(regout_mod1, regout_mod2, regout_mod3, regout_mod4)

model_colors_risk <- c("Model 1 (Basic)"       = "#FCA500",
                       "Model 2 (+Covariates)"  = "#F90102",
                       "Model 3 (+WBC)"         = "#C80000",
                       "Model 4 (+PCs)"         = "#7B0000")

combined_plot_data$MRS_quartile <- factor(combined_plot_data$MRS_quartile,
                                          levels = c("Q4", "Q3", "Q2"), ordered = TRUE)

reversed_quartile_plot <- ggplot(combined_plot_data,
                                 aes(x = relative_risk,
                                     y = factor(model, levels = rev(names(model_colors_risk))),
                                     fill = model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, aes(group = MRS_quartile)) +
  scale_fill_manual(values = model_colors_risk, name = "Model") +
  labs(x = "Relative Risk", y = NULL,
       title = "Relative Risk by Model and MRS Quartile (Grouped)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.position = "right", legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), panel.grid.major.y = element_blank()) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))

print(reversed_quartile_plot)

sessionInfo()
