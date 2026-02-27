################################################################################
#
#  Matched-Background Permutation for Pathway / Enrichment Analysis
#
#  Description:
#  This script generates null-matched background CpG sets for use in
#  permutation-based enrichment or pathway analysis. For each CpG of
#  interest (e.g. EWAS sentinel hits), it identifies background CpGs
#  that are matched on mean methylation level and variability (SD) but
#  located at independent genomic loci (>5 kb away).
#
#  This approach controls for the known relationship between methylation
#  level/variability and the probability of being associated with a trait,
#  avoiding inflation in enrichment tests.
#
#  Pipeline steps:
#    1. Load EWAS-relevant CpGs and manifest with hg38 coordinates
#    2. Compute pooled mean and SD of methylation across ethnic groups
#    3. Build a filtered background table (excluding CpGs of interest
#       and any CpG within 5 kb of a CpG of interest)
#    4. For each CpG of interest, find matched background CpGs using a
#       sliding threshold on mean and SD differences, then sample
#       N permuted matches without replacement
#    5. Clean and save permuted background CpG sets
#
#  Matching criteria (per CpG of interest):
#    i.   |mean_background - mean_hit| < threshold_mean
#    ii.  |sd_background   - sd_hit|   < threshold_sd
#    iii. Distance > 5 kb (i.e. not at the same genomic locus)
#    Thresholds start at mean < 0.025 and SD < 0.0025, increasing
#    incrementally until the required number of matches is reached.
#
#  Inputs:
#    - EPIC.hg38.manifest.tsv.gz :
#        Illumina EPIC manifest (hg38 coordinates). Key columns used:
#          Probe_ID  -> CpG identifier
#          CpG_chrm  -> chromosome
#          CpG_beg   -> start position (0-based)
#
#    - asian_relevant_cpgs.rds :
#        Character vector of CpG IDs of interest (e.g. EWAS hits)
#
#    - backg_cpgs.rds :
#        Character vector of all background CpG IDs (autosomal probes
#        that passed QC)
#
#    - beta_QN_rmGSwap_rmDup_HELIOS.RData :
#        Contains `betaRmDup`, a matrix [CpG probes x samples] of
#        quantile-normalised beta values
#
#    - cIMT_samplesheet_wsmoking.rds :
#        Phenotype file with columns: SentrixID, FREG5_Ethnic_Group
#
#  Outputs:
#    - msd.rds :
#        Per-CpG, per-ethnicity summary statistics (n, mean, variance)
#    - msd_heliosMA_EPIC.rds :
#        Pooled mean and SD per CpG with genomic positions
#    - backg_<N>_permutations_asianrel.RData :
#        Contains `matches` [hits x N_permutations] matrix of matched
#        background CpG IDs, and `matchParameters` [hits x 3] matrix
#        recording the thresholds used for each hit
#    - backg_cpgs_perm_<N>_cleaned.rds :
#        Cleaned version of `matches` with NA rows removed
#
#  Dependencies:
#    install.packages(c("dplyr", "tibble", "tidyr", "stringr"))
#
################################################################################


# ==============================================================================
# 0. User-configurable parameters
# ==============================================================================

# -- File paths --
manifest_file   <- "EPIC.hg38.manifest.tsv.gz"
hits_file       <- "asian_relevant_cpgs.rds"
background_file <- "backg_cpgs.rds"
beta_file       <- "beta_QN_rmGSwap_rmDup_HELIOS.RData"
phenotype_file  <- "cIMT_samplesheet_wsmoking.rds"

# -- Matching parameters --
n_permutations  <- 20       # number of matched background CpGs per hit
min_match_pool  <- 20       # minimum pool size required before sampling
locus_distance  <- 5000     # minimum distance (bp) between hit and match


# ==============================================================================
# 1. Load data
# ==============================================================================

library(dplyr)
library(tibble)
library(tidyr)
library(stringr)

# Manifest: CpG -> chromosome + position (hg38)
manifest <- read.delim(manifest_file)
manifest <- manifest %>%
  dplyr::rename(CpG = Probe_ID, Chr_hg38 = CpG_chrm, Start_hg38_0based = CpG_beg) %>%
  dplyr::select(CpG, Chr_hg38, Start_hg38_0based)

# CpGs of interest and background set
asian_relevant_cpgs <- readRDS(hits_file)
backg_cpgs          <- readRDS(background_file)


# ==============================================================================
# 2. Compute pooled mean and SD of methylation across ethnic groups
# ==============================================================================
#
# Pooled statistics account for differing sample sizes per ethnic group:
#   pooled_mean = sum(n_k * mean_k) / sum(n_k)
#   pooled_sd   = sqrt( sum((n_k - 1) * var_k) / (sum(n_k) - K) )
# ==============================================================================

# Load methylation betas
load(beta_file)  # loads: betaRmDup [CpGs x samples]
betaRmDup <- data.frame(betaRmDup, check.names = FALSE)

# Read phenotype and exclude ethnic group "O" (not analysed in EWAS)
cIMT_pheno  <- readRDS(phenotype_file)
samplesheet <- cIMT_pheno[!cIMT_pheno$FREG5_Ethnic_Group == "O", ]
sampleid    <- samplesheet$SentrixID

# Subset methylation betas to background CpGs and EWAS-analysed samples
helios_meth_betas <- betaRmDup %>%
  dplyr::filter(rownames(betaRmDup) %in% backg_cpgs) %>%
  dplyr::select(dplyr::all_of(sampleid))

# Per-CpG, per-ethnicity summary statistics (n, mean, variance)
msd <- helios_meth_betas %>%
  as.data.frame() %>%
  rownames_to_column("CpG") %>%
  pivot_longer(cols = -CpG, names_to = "SentrixID", values_to = "value") %>%
  left_join(samplesheet[, c("SentrixID", "FREG5_Ethnic_Group")], by = "SentrixID") %>%
  group_by(CpG, FREG5_Ethnic_Group) %>%
  summarise(
    n    = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    var  = var(value, na.rm = TRUE)
  ) %>%
  ungroup()

# Calculate pooled mean and SD, then append genomic positions
msd_pooled <- msd %>%
  group_by(CpG) %>%
  summarise(
    meanMeth = sum(n * mean) / sum(n),
    sdMeth   = sqrt(sum((n - 1) * var) / (sum(n) - n()))
  ) %>%
  left_join(manifest[, c("CpG", "Chr_hg38", "Start_hg38_0based")], by = "CpG") %>%
  dplyr::rename(pos = Start_hg38_0based) %>%
  filter(!is.na(pos))

rownames(msd_pooled) <- msd_pooled$CpG

saveRDS(msd, "msd.rds")
saveRDS(msd_pooled, "msd_heliosMA_EPIC.rds")


# ==============================================================================
# 3. Build filtered background table
# ==============================================================================
#
# Remove from the background pool:
#   (a) The CpGs of interest themselves
#   (b) Any CpG within `locus_distance` bp of a CpG of interest
#       (to ensure genomic independence of matched controls)
# ==============================================================================

asian_relevant_cpgs <- readRDS(hits_file)
msd <- readRDS("msd_heliosMA_EPIC.rds")

# Subset pooled stats for the CpGs of interest
hits      <- data.frame(msd[as.character(asian_relevant_cpgs), ])
hits$CG   <- hits$CpG

# Identify all CpGs within `locus_distance` of any hit
for (c in 1:length(hits$CG)) {

  print(c)
  cg          <- hits$CG[c]
  cg.coord    <- msd[cg, ]
  cis.cg.coord <- msd[msd$Chr_hg38 == cg.coord$Chr_hg38 &
                       abs(msd$pos - cg.coord$pos) < locus_distance, ]

  if (c == 1) { cis.cgs <- cis.cg.coord$CpG }
  if (c > 1)  { cis.cgs <- c(cis.cgs, cis.cg.coord$CpG) }
}

cis.cgs <- unique(cis.cgs)

# Background = all CpGs minus hits and their cis neighbours
backg <- data.frame(msd[!(rownames(msd) %in% hits$CG |
                           rownames(msd) %in% cis.cgs), ])
rownames(backg) <- backg$CpG


# ==============================================================================
# 4. Generate matched-background permutations
# ==============================================================================
#
# For each CpG of interest, find background CpGs matched on mean and SD
# using a sliding threshold grid:
#   - SD thresholds:   0.0025, 0.005, ..., 0.025
#   - Mean thresholds: 0.025,  0.050, ..., 0.250
# Start with the tightest combination; use the first pair that yields
# >= `min_match_pool` eligible background CpGs. Then sample
# `n_permutations` unique matches without replacement.
#
# Output matrices:
#   matches        : [n_hits x n_permutations] matched background CpG IDs
#   matchParameters: [n_hits x 3] mean_threshold, sd_threshold, pool_size
# ==============================================================================

# Expand all threshold combinations
parameters <- expand.grid(
  seq(0.0025, 0.025, 0.0025),   # SD thresholds
  seq(0.025, 0.25, 0.025),      # Mean thresholds
  stringsAsFactors = FALSE
)
colnames(parameters) <- c("SD", "Mean")

# Initialise output matrices
matches <- matrix(nrow = nrow(hits), ncol = n_permutations)
colnames(matches) <- as.character(seq(from = 1, to = n_permutations, 1))
rownames(matches) <- hits$CG

matchParameters <- matrix(nrow = nrow(hits), ncol = 3)
colnames(matchParameters) <- c("mean_threshold", "sd_threshold", "n_matches")
rownames(matchParameters) <- hits$CG

for (h in 1:nrow(hits)) {

  cpg.match <- NA
  len <- NA
  hit <- hits[h, ]
  print(h)

  # -- 4a. Calculate pool size across all threshold combinations --
  for (para in 1:nrow(parameters)) {
    n <- length(sample(rownames(backg[
      abs(backg$meanMeth - hit$meanMeth) < parameters[para, 2] &
      abs(backg$sdMeth   - hit$sdMeth)   < parameters[para, 1], ])))
    if (para == 1) { sel <- c(parameters[para, 2], parameters[para, 1], n) }
    if (para > 1)  { sel <- rbind(sel, c(parameters[para, 2], parameters[para, 1], n)) }
  }
  colnames(sel) <- c("mean", "sd", "n_matches")

  # -- 4b. Select tightest thresholds yielding >= min_match_pool matches --
  sel <- sel[which(sel[, 3] > min_match_pool), , drop = FALSE]
  if (nrow(sel) > 0) {
    meanthres <- sel[1, 1]
    sdthres   <- sel[1, 2]
    nthresh   <- sel[1, 3]
  } else {
    message("No rows meet the condition for hit ", h)
    meanthres <- sdthres <- nthresh <- NA
  }

  # -- 4c. Sample matched CpGs without replacement --
  for (p in 1:n_permutations) {
    set.seed(p)
    if (exists("cgm")) { rm(cgm) }
    cgms <- sample(rownames(backg[
      abs(backg$meanMeth - hit$meanMeth) < meanthres &
      abs(backg$sdMeth   - hit$sdMeth)   < sdthres, ]))
    cgm <- setdiff(cgms, cpg.match)[1]
    if (!is.na(cgm)) {
      cpg.match <- c(cpg.match, cgm)
    } else {
      cpg.match <- c(cpg.match, (as.character(rownames(hit))))  # fallback
      print("no matching marker found")
    }
  }

  matches[h, ]         <- na.omit(cpg.match)
  matchParameters[h, ] <- c(meanthres, sdthres, nthresh)
}

save(matches, matchParameters,
     file = paste0("backg_", n_permutations, "_permutations_asianrel.RData"))


# ==============================================================================
# 5. Clean and save permuted CpG sets
# ==============================================================================
#
# Remove rows where any permutation column contains an NA-like string
# (arises when the fallback was used or matching failed).
# ==============================================================================

load(paste0("backg_", n_permutations, "_permutations_asianrel.RData"))

backg_cpgs_df_cleaned <- matches %>%
  data.frame() %>%
  # Convert "NA..." strings to real NA
  mutate(across(everything(),
                ~ ifelse(str_starts(as.character(.), fixed("NA")), NA, .))) %>%
  # Remove rows with any NA values
  filter(if_all(everything(), ~ !is.na(.)))

saveRDS(backg_cpgs_df_cleaned,
        paste0("backg_cpgs_perm_", n_permutations, "_cleaned.rds"))

sessionInfo()
