################################################################################
#
#  Epigenome-Wide Association Study (EWAS) Regression Pipeline
#
#  Description:
#  This script performs an epigenome-wide association study testing the
#  association between DNA methylation (beta values) at each CpG site and
#  log-transformed cIMT, adjusting for covariates. The analysis is run
#  separately for each ethnic group.
#
#  Model:
#    log(cIMT) ~ CpG_beta + Age + Gender + Smoking + WBC (6) + 30 ctrl PCs
#
#  This script shows the Chinese ("C") analysis. To run for Malay or Indian,
#  change the ethnicity subset from phe_Ethnicity$'C' to $'M' or $'I',
#  update variable names accordingly (e.g. phe_CHN -> phe_MLY), and adjust
#  the WBC/PC column indices if needed.
#
#  Inputs:
#    - ctrlprobes.RData :
#        Contains `ctrl.all` [samples x control probes], control probe
#        intensities from the preprocessing pipeline.
#    - housemanEstimates_QN_detPrelax_betas.RData :
#        Contains `constrainedCoefs` [samples x cell types], estimated
#        white blood cell proportions (Houseman method).
#    - beta_QN_rmGSwap_rmDup_HELIOS.RData :
#        Contains `betaRmDup` [CpG probes x samples], quantile-normalised
#        beta values after gender-swap and duplicate removal.
#    - cIMT_samplesheet_wsmoking.rds :
#        Phenotype data frame [samples x variables] with columns:
#          SentrixID, FREG5_Ethnic_Group, FREG8_Age, FREG7_Gender,
#          CIMT_distal_mean, smoking
#
#  Output:
#    - PC_cp_CHN.rds :
#        Ethnic-specific control probe PCA scores [samples x PCs].
#    - regout_cIMTmean_CHN_AgeGenderWBCsmokingAdj.RData :
#        Contains `regout_CHN`, a matrix [CpG probes x 4] with the
#        methylation coefficient summary statistics:
#          Column 1: Estimate (regression coefficient)
#          Column 2: Std. Error
#          Column 3: t value
#          Column 4: Pr(>|t|) (p-value)
#
#  Dependencies:
#    install.packages(c("dplyr", "forcats"))
#
################################################################################


# ==============================================================================
# 1. Load libraries and data
# ==============================================================================

library(dplyr)
library(forcats)

# Control probe intensities [samples x control probes]
load("ctrlprobes.RData")

# WBC estimates [samples x cell types]
load("housemanEstimates_QN_detPrelax_betas.RData")

# Methylation betas [CpG probes x samples]
load("beta_QN_rmGSwap_rmDup_HELIOS.RData")
# betaRmDup <- betaRmDup[c(1:10),] # run test on subset


# ==============================================================================
# 2. Prepare phenotype data
# ==============================================================================

phe <- readRDS("cIMT_samplesheet_wsmoking.rds")
dim(phe) #1363; 17

# Create log-transformed outcome
phe$log_CIMT_distal_mean <- log(as.numeric(phe$CIMT_distal_mean))

# Rename columns and recode smoking into 3 categories:
#   Original codes: 1 = never, 2 = passive/minimal -> non-smoker
#                   3 = former                      -> ex-smoker
#                   4 = occasional, 5 = daily       -> smoker
phe <- phe %>%
  rename(
    Age = FREG8_Age,
    Gender = FREG7_Gender
  ) %>%
  mutate(smoking = case_when(
    smoking == 1 ~ 1,
    smoking == 2 ~ 1,
    smoking == 3 ~ 2,
    smoking == 4 ~ 3,
    smoking == 5 ~ 3,
    is.na(smoking) ~ NA_real_
  )) %>%
  mutate(smoking = factor(smoking,
                          levels = c(1, 2, 3),
                          labels = c('non-smoker', 'ex-smoker', 'smoker')))


# ==============================================================================
# 3. Subset to target ethnic group
# ==============================================================================

# Split phenotype by ethnicity and select Chinese
phe_Ethnicity <- split(phe, phe$FREG5_Ethnic_Group)
phe_CHN <- phe_Ethnicity$'C'   # change to $'M' for Malay or $'I' for Indian
dim(phe_CHN)


# ==============================================================================
# 4. Compute ethnic-specific control probe PCs
# ==============================================================================

ctrl.all_CHN <- ctrl.all[rownames(ctrl.all) %in% phe_CHN[, 4], ]
dim(ctrl.all_CHN)

pca_cp_CHN <- prcomp(na.omit(ctrl.all_CHN))
PC_cp_CHN <- predict(pca_cp_CHN)
saveRDS(PC_cp_CHN, "PC_cp_CHN.rds")
#PC_cp_CHN <- readRDS("PC_cp_CHN.rds")


# ==============================================================================
# 5. Subset methylation betas to ethnic group
# ==============================================================================

dim(betaRmDup) #837772 markers

# Remove CpG sites that are entirely NA across all samples
betaRmDup <- betaRmDup[rowSums(is.na(betaRmDup)) != ncol(betaRmDup), ]
dim(betaRmDup)

# Select ethnic-specific betas
betaRmDup_CHN <- betaRmDup[, colnames(betaRmDup) %in% phe_CHN[, 4]]
dim(betaRmDup_CHN)


# ==============================================================================
# 6. Build master covariate table
# ==============================================================================

# Format control probe PCs for merging
PC_cp_CHN <- as.data.frame(PC_cp_CHN)
PC_cp_CHN$SentrixID <- row.names(PC_cp_CHN)

# Format WBC proportions for merging
constrainedCoefs <- as.data.frame(constrainedCoefs)
constrainedCoefs$SentrixID <- row.names(constrainedCoefs)

# Merge phenotype, control probe PCs, and WBC proportions
phe_CHN_all_variables <- merge(phe_CHN, PC_cp_CHN, by = "SentrixID")
phe_CHN_all_variables <- merge(phe_CHN_all_variables, constrainedCoefs, by = "SentrixID")
dim(phe_CHN_all_variables)


# ==============================================================================
# 7. Align beta matrix columns to covariate table rows
# ==============================================================================

beta_order <- phe_CHN_all_variables$SentrixID
betaRmDup_CHN <- betaRmDup_CHN[, match(beta_order, colnames(betaRmDup_CHN))]
colnames(betaRmDup_CHN) == phe_CHN_all_variables$SentrixID


# ==============================================================================
# 8. Construct regression formula
# ==============================================================================
#
# Additional covariates (control probe PCs and WBC proportions) are identified
# by column index in the merged covariate table.
# NOTE: Column indices may differ across ethnicities due to merging order.
# Inspect `Header_List` and adjust indices accordingly.
# ==============================================================================

Header_List <- c()
Header_List <- colnames(phe_CHN_all_variables)
Independent_Variable_List <- c(Header_List[21:50], Header_List[225:230])
print(Independent_Variable_List)

# Build the right-hand side of the regression formula
# Loop through covariates one by one and add to the equation
RHS_equation <- '+Age+as.factor(Gender)+as.factor(smoking)'
for (element_loop in Independent_Variable_List) {
  RHS_equation <- paste(RHS_equation, '+', element_loop, sep = "")
}
print(RHS_equation)


# ==============================================================================
# 9. Run EWAS: probe-by-probe linear regression
# ==============================================================================
#
# For each CpG site, fit:
#   log(cIMT) ~ beta_at_CpG + Age + Gender + Smoking + control_probe_PCs + WBC
# and extract the coefficient summary for the methylation term (row 2).
#
# Output: regout_CHN [CpG probes x 4] with columns:
#   Estimate, Std. Error, t value, Pr(>|t|)
# ==============================================================================

#betaRmDup_CHN <- as.matrix(betaRmDup_CHN)
regout_CHN <- matrix(0, nrow(betaRmDup_CHN), 4)

for (i in 1:nrow(betaRmDup_CHN)) {
  reg_formula <- ''
  reg_formula <- as.formula(as.character(paste('log_CIMT_distal_mean~', 'betaRmDup_CHN[i,]', RHS_equation, sep = "")))
  regout_CHN[i, ] <- summary(lm(reg_formula, data = phe_CHN_all_variables, na.action = na.exclude))$coefficient[2, ]
}

rownames(regout_CHN) <- rownames(betaRmDup_CHN)
dim(regout_CHN)
print(reg_formula)
save(regout_CHN, file = "regout_cIMTmean_CHN_AgeGenderWBCsmokingAdj.RData")
