################################################################################
#
#  DNA Methylation Array Preprocessing Pipeline
#  Platform: Illumina EPIC (850K) array
#
#  Description:
#  This script performs quality control, filtering, and normalisation of
#  Illumina EPIC array methylation data starting from raw IDAT files.
#
#  Pipeline steps:
#    1. Read raw IDAT files and apply Illumina background correction
#    2. Extract signal intensities for Type I and Type II probes
#    3. Extract control probe intensities and summarise via PCA
#    4. Compute detection p-values
#    5. Filter probes and samples by detection p-value and call rate
#    6. Restrict to autosomal CpG and CpH probes
#    7. Quantile-normalise signal intensities
#    8. Calculate beta values
#
#  Inputs:
#    - IDAT files (.idat) in the working directory, one pair per sample
#    - Manifest file "manifestb2.csv" (Illumina EPIC manifest CSV with a
#      7-row header; columns must include: Name, Infinium_Design_Type,
#      Color_Channel, CHR, MAPINFO)
#
#  Outputs:
#    - intensities.RData         : raw probe signal intensities (probes x samples)
#    - ctrlprobes.RData          : control probe matrix and PCA scores
#    - detectionPvalue.RData     : detection p-value matrix (probes x samples)
#    - callRates_detP001.RData   : per-sample and per-marker call rates
#    - beta_raw_detP001.RData    : raw beta values before normalisation
#    - beta_QN_detP001_marker095.RData : final quantile-normalised beta values
#
#  Key data structures:
#    - TypeII.{Green,Red}.All    : matrices [Type II probes x samples],
#                                  unmethylated (Green) and methylated (Red)
#    - TypeI.{Green,Red}.{M,U}.All : matrices [Type I probes x samples],
#                                  methylated (M) and unmethylated (U) split
#                                  by colour channel (Green or Red)
#    - ctrl.all                  : matrix [samples x control probes], intensities
#                                  for QC control probes (BS conversion, staining,
#                                  extension, hybridisation, etc.)
#    - dp.all                    : matrix [probes x samples], detection p-values
#    - beta                      : matrix [autosomal probes x samples], beta
#                                  values ranging from 0 (unmethylated) to 1
#                                  (fully methylated)
#
#  Dependencies:
#    install.packages("BiocManager")
#    BiocManager::install(c("minfi", "IlluminaHumanMethylationEPICmanifest",
#                           "S4Vectors", "limma"))
#
################################################################################


# ==============================================================================
# 0. Load required packages
# ==============================================================================

require(minfi)
require(IlluminaHumanMethylationEPICmanifest)
require(S4Vectors)
require(limma)


# ==============================================================================
# 1. Read IDAT files, background-correct, and extract probe intensities
# ==============================================================================

# Identify unique sample basenames from IDAT filenames in the working directory.
# Each sample has two IDAT files (*_Grn.idat and *_Red.idat) sharing a common
# 19-character prefix (Sentrix barcode + position).
filenames <- unique(substr(dir(), 1, 19))
dim(filenames) <- c(length(filenames), 1)

for (i in 1:nrow(filenames)) {
  message("Processing sample ", i, " of ", nrow(filenames), ": ", filenames[i, ])

  # Read a single sample into an RGChannelSet and apply Illumina background
  # subtraction (NOOB-style negative control correction).
  RGset <- read.metharray(file.path(paste0("./", filenames[i, ])), verbose = TRUE)
  RGset <- bgcorrect.illumina(RGset)

  # --------------------------------------------------------------------------
  # 1a. Type II probe intensities
  #     Type II probes use a single bead: Green channel = unmethylated signal,
  #     Red channel = methylated signal.
  # --------------------------------------------------------------------------
  TypeII.Name  <- getProbeInfo(RGset, type = "II")$Name
  TypeII.Green <- as.matrix(getGreen(RGset)[getProbeInfo(RGset, type = "II")$AddressA, ])
  TypeII.Red   <- as.matrix(getRed(RGset)[getProbeInfo(RGset, type = "II")$AddressA, ])
  rownames(TypeII.Green) <- rownames(TypeII.Red) <- TypeII.Name
  colnames(TypeII.Green) <- colnames(TypeII.Red) <- sampleNames(RGset)

  # --------------------------------------------------------------------------
  # 1b. Type I probe intensities
  #     Type I probes use two beads (methylated = AddressB, unmethylated =
  #     AddressA), measured in the same colour channel.
  # --------------------------------------------------------------------------

  # -- Green-channel Type I probes --
  TypeI.Green.Name <- getProbeInfo(RGset, type = "I-Green")$Name
  TypeI.Green.M    <- as.matrix(getGreen(RGset)[getProbeInfo(RGset, type = "I-Green")$AddressB, ])
  TypeI.Green.U    <- as.matrix(getGreen(RGset)[getProbeInfo(RGset, type = "I-Green")$AddressA, ])
  rownames(TypeI.Green.M) <- rownames(TypeI.Green.U) <- TypeI.Green.Name
  colnames(TypeI.Green.M) <- colnames(TypeI.Green.U) <- sampleNames(RGset)

  # -- Red-channel Type I probes --
  TypeI.Red.Name <- getProbeInfo(RGset, type = "I-Red")$Name
  TypeI.Red.M    <- as.matrix(getRed(RGset)[getProbeInfo(RGset, type = "I-Red")$AddressB, ])
  TypeI.Red.U    <- as.matrix(getRed(RGset)[getProbeInfo(RGset, type = "I-Red")$AddressA, ])
  rownames(TypeI.Red.M) <- rownames(TypeI.Red.U) <- TypeI.Red.Name
  colnames(TypeI.Red.M) <- colnames(TypeI.Red.U) <- sampleNames(RGset)

  # --------------------------------------------------------------------------
  # 1c. Control probe intensities
  #     These are used for QC assessment. Control probe indices below are
  #     specific to the EPIC manifest; verify against your manifest version.
  # --------------------------------------------------------------------------
  control <- getProbeInfo(RGset, type = "Control")

  # Helper: extract control probe intensities into a named matrix
  extract_ctrl <- function(rows, channel = c("Green", "Red")) {
    channel <- match.arg(channel)
    probe_names <- control[rows, ]$ExtendedType
    addresses   <- control[rows, ]$Address
    mat <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(probe_names),
                  dimnames = list(probe_names, sampleNames(RGset)))
    if (channel == "Green") {
      mat[probe_names, ] <- as.matrix(getGreen(RGset)[addresses, ])
    } else {
      mat[probe_names, ] <- as.matrix(getRed(RGset)[addresses, ])
    }
    return(mat)
  }

  # Bisulfite conversion I (converted probes only; excludes BS Conversion I-U)
  BSCI.Green <- extract_ctrl(16:17, "Green")
  BSCI.Red   <- extract_ctrl(18:20, "Red")

  # Bisulfite conversion II
  BSCII.Red  <- extract_ctrl(26:29, "Red")

  # Staining
  stain.Red   <- extract_ctrl(3, "Red")
  stain.Green <- extract_ctrl(5, "Green")

  # Extension (single-base)
  extensionA.Red   <- extract_ctrl(9,  "Red")
  extensionT.Red   <- extract_ctrl(8,  "Red")
  extensionC.Green <- extract_ctrl(10, "Green")
  extensionG.Green <- extract_ctrl(7,  "Green")

  # Hybridisation (High / Medium / Low)
  hybridH.Green <- extract_ctrl(13, "Green")
  hybridM.Green <- extract_ctrl(12, "Green")
  hybridL.Green <- extract_ctrl(11, "Green")

  # Target removal
  target.Green <- extract_ctrl(14:15, "Green")

  # Specificity I
  specI.Green <- extract_ctrl(30:32, "Green")
  specI.Red   <- extract_ctrl(36:38, "Red")

  # Specificity II
  specII.Red <- extract_ctrl(42:44, "Red")

  # Non-polymorphic
  np.Red   <- extract_ctrl(45:46, "Red")
  np.Green <- extract_ctrl(47:48, "Green")

  # Normalisation control probes (NORM_A, NORM_T, NORM_C, NORM_G)
  normC.Green.Name <- control[control[, 2] == "NORM_C", 4]
  normC.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(normC.Green.Name),
                        dimnames = list(normC.Green.Name, sampleNames(RGset)))
  normC.Green[normC.Green.Name, ] <- as.matrix(getGreen(RGset)[control[control[, 2] == "NORM_C", 1], ])

  normG.Green.Name <- control[control[, 2] == "NORM_G", 4]
  normG.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(normG.Green.Name),
                        dimnames = list(normG.Green.Name, sampleNames(RGset)))
  normG.Green[normG.Green.Name, ] <- as.matrix(getGreen(RGset)[control[control[, 2] == "NORM_G", 1], ])

  normA.Red.Name <- control[control[, 2] == "NORM_A", 4]
  normA.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(normA.Red.Name),
                      dimnames = list(normA.Red.Name, sampleNames(RGset)))
  normA.Red[normA.Red.Name, ] <- as.matrix(getRed(RGset)[control[control[, 2] == "NORM_A", 1], ])

  normT.Red.Name <- control[control[, 2] == "NORM_T", 4]
  normT.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(normT.Red.Name),
                      dimnames = list(normT.Red.Name, sampleNames(RGset)))
  normT.Red[normT.Red.Name, ] <- as.matrix(getRed(RGset)[control[control[, 2] == "NORM_T", 1], ])

  # Combine all control probe intensities into a single matrix
  # Result: ctrl [control probes x 1], transposed later to [samples x control probes]
  ctrl <- rbind(
    as.matrix(BSCI.Green), as.matrix(BSCI.Red), as.matrix(BSCII.Red),
    stain.Red, stain.Green,
    extensionA.Red, extensionT.Red, extensionC.Green, extensionG.Green,
    hybridH.Green, hybridM.Green, hybridL.Green,
    as.matrix(target.Green),
    as.matrix(specI.Green), as.matrix(specI.Red), as.matrix(specII.Red),
    np.Red[1, , drop = FALSE], np.Red[2, , drop = FALSE],
    np.Green[1, , drop = FALSE], np.Green[2, , drop = FALSE],
    as.matrix(normC.Green), as.matrix(normG.Green),
    as.matrix(normA.Red), as.matrix(normT.Red)
  )

  # --------------------------------------------------------------------------
  # 1d. Detection p-values (m+u method)
  #     Low p-values indicate signal reliably above background.
  # --------------------------------------------------------------------------
  dp <- detectionP(RGset, type = "m+u")

  # --------------------------------------------------------------------------
  # 1e. Accumulate across samples
  # --------------------------------------------------------------------------
  if (exists("TypeII.Red.All")) {
    TypeII.Red.All    <- cbind(TypeII.Red.All, TypeII.Red)
    TypeII.Green.All  <- cbind(TypeII.Green.All, TypeII.Green)
    TypeI.Red.M.All   <- cbind(TypeI.Red.M.All, TypeI.Red.M)
    TypeI.Red.U.All   <- cbind(TypeI.Red.U.All, TypeI.Red.U)
    TypeI.Green.M.All <- cbind(TypeI.Green.M.All, TypeI.Green.M)
    TypeI.Green.U.All <- cbind(TypeI.Green.U.All, TypeI.Green.U)
    ctrl.all          <- rbind(ctrl.all, t(ctrl))
    dp.all            <- cbind(dp.all, dp)
  } else {
    TypeII.Red.All    <- TypeII.Red
    TypeII.Green.All  <- TypeII.Green
    TypeI.Red.M.All   <- TypeI.Red.M
    TypeI.Red.U.All   <- TypeI.Red.U
    TypeI.Green.M.All <- TypeI.Green.M
    TypeI.Green.U.All <- TypeI.Green.U
    ctrl.all          <- t(ctrl)
    dp.all            <- dp
  }
}


# ==============================================================================
# 2. PCA on control probe intensities (for QC / batch assessment)
# ==============================================================================

# ctrl.all: [samples x control probes]. PCA yields scores that can be used
# as covariates or to identify outlier/batch effects.
pca <- prcomp(na.omit(ctrl.all))
ctrlprobes.scores <- pca$x
colnames(ctrlprobes.scores) <- paste0(colnames(ctrlprobes.scores), "_cp")


# ==============================================================================
# 3. Save intermediate data
# ==============================================================================

save(TypeII.Red.All, TypeII.Green.All,
     TypeI.Red.M.All, TypeI.Red.U.All,
     TypeI.Green.M.All, TypeI.Green.U.All,
     file = "./intensities.RData")

save(ctrl.all, ctrlprobes.scores, file = "./ctrlprobes.RData")
save(dp.all, file = "./detectionPvalue.RData")


# ==============================================================================
# 4. Define autosomal probe set from manifest
# ==============================================================================

# The manifest CSV has a 7-line header (skip = 7). Key columns:
#   Name                 : probe ID (e.g. "cg00000029", "ch.1.1234")
#   Infinium_Design_Type : "I" or "II"
#   Color_Channel        : "Grn" or "Red" (Type I only)
#   CHR                  : chromosome (1-22, X, Y)
#   MAPINFO              : genomic coordinate
anno <- read.csv("manifestb2.csv", as.is = TRUE, skip = 7)
anno <- anno[, c("Infinium_Design_Type", "Color_Channel", "CHR", "MAPINFO", "Name")]

# Select autosomal CpG ("cg") and CpH ("ch.") probes, excluding sex chromosomes
cgs  <- anno[substr(anno$Name, 1, 2) == "cg"  & !(anno$CHR %in% c("X", "Y")), ]
cas  <- anno[substr(anno$Name, 1, 3) == "ch." & !(anno$CHR %in% c("X", "Y")), ]
auto <- c(cgs$Name, cas$Name)


# ==============================================================================
# 5. Apply detection p-value filter
# ==============================================================================

# Probes with detection p >= threshold are set to NA (unreliable signal).
det_pval_threshold <- 0.01

load("intensities.RData")
load("detectionPvalue.RData")

# Type II probes
d <- dp.all[rownames(TypeII.Green.All), colnames(TypeII.Green.All)]
TypeII.Green.All.d <- ifelse(d < det_pval_threshold, TypeII.Green.All, NA)
TypeII.Red.All.d   <- ifelse(d < det_pval_threshold, TypeII.Red.All, NA)

# Type I Green probes
d <- dp.all[rownames(TypeI.Green.M.All), colnames(TypeI.Green.M.All)]
TypeI.Green.M.All.d <- ifelse(d < det_pval_threshold, TypeI.Green.M.All, NA)
TypeI.Green.U.All.d <- ifelse(d < det_pval_threshold, TypeI.Green.U.All, NA)

# Type I Red probes
d <- dp.all[rownames(TypeI.Red.M.All), colnames(TypeI.Red.M.All)]
TypeI.Red.M.All.d <- ifelse(d < det_pval_threshold, TypeI.Red.M.All, NA)
TypeI.Red.U.All.d <- ifelse(d < det_pval_threshold, TypeI.Red.U.All, NA)

rm(dp.all, d)


# ==============================================================================
# 6. Compute raw beta values and call rates (before normalisation)
# ==============================================================================

# Subset to autosomal probes
samples  <- colnames(TypeI.Red.M.All)

markers_II  <- intersect(rownames(TypeII.Green.All.d), auto)
markers_IG  <- intersect(rownames(TypeI.Green.M.All.d), auto)
markers_IR  <- intersect(rownames(TypeI.Red.M.All.d), auto)

# Beta = Methylated / (Methylated + Unmethylated + offset)
# The offset of 100 stabilises beta values when total intensity is low.
TypeII.betas       <- TypeII.Green.All.d[markers_II, samples] /
                      (TypeII.Red.All.d[markers_II, samples] +
                       TypeII.Green.All.d[markers_II, samples] + 100)
TypeI.Green.betas  <- TypeI.Green.M.All.d[markers_IG, samples] /
                      (TypeI.Green.M.All.d[markers_IG, samples] +
                       TypeI.Green.U.All.d[markers_IG, samples] + 100)
TypeI.Red.betas    <- TypeI.Red.M.All.d[markers_IR, samples] /
                      (TypeI.Red.M.All.d[markers_IR, samples] +
                       TypeI.Red.U.All.d[markers_IR, samples] + 100)

beta_raw <- as.matrix(rbind(TypeII.betas, TypeI.Green.betas, TypeI.Red.betas))

# Call rates: proportion of non-NA values
sample.call <- colSums(!is.na(beta_raw)) / nrow(beta_raw)
marker.call <- rowSums(!is.na(beta_raw)) / ncol(beta_raw)

save(sample.call, marker.call, file = "callRates_detP001.RData")
save(beta_raw, file = "beta_raw_detP001.RData")

rm(TypeII.betas, TypeI.Green.betas, TypeI.Red.betas, beta_raw)


# ==============================================================================
# 7. Call-rate filtering
# ==============================================================================

# Remove samples with < 98% call rate and markers with < 95% call rate
load("callRates_detP001.RData")

sample_callrate_threshold <- 0.98
marker_callrate_threshold <- 0.95

samples <- names(sample.call[sample.call > sample_callrate_threshold])

# Subset each probe type to passing autosomal markers and passing samples
markers_II <- intersect(rownames(TypeII.Green.All.d), auto)
markers_II <- intersect(markers_II, names(marker.call[marker.call >= marker_callrate_threshold]))
TypeII.Green <- TypeII.Green.All.d[markers_II, samples]
TypeII.Red   <- TypeII.Red.All.d[markers_II, samples]

markers_IG <- intersect(rownames(TypeI.Green.M.All.d), auto)
markers_IG <- intersect(markers_IG, names(marker.call[marker.call >= marker_callrate_threshold]))
TypeI.Green.M <- TypeI.Green.M.All.d[markers_IG, samples]
TypeI.Green.U <- TypeI.Green.U.All.d[markers_IG, samples]

markers_IR <- intersect(rownames(TypeI.Red.M.All.d), auto)
markers_IR <- intersect(markers_IR, names(marker.call[marker.call >= marker_callrate_threshold]))
TypeI.Red.M <- TypeI.Red.M.All.d[markers_IR, samples]
TypeI.Red.U <- TypeI.Red.U.All.d[markers_IR, samples]

rm(TypeII.Green.All.d, TypeII.Red.All.d,
   TypeI.Green.M.All.d, TypeI.Green.U.All.d,
   TypeI.Red.M.All.d, TypeI.Red.U.All.d)


# ==============================================================================
# 8. Quantile normalisation and final beta calculation
# ==============================================================================

# Quantile normalisation is applied separately within each probe-type /
# channel group to bring sample intensity distributions into alignment.
TypeII.Green  <- normalizeQuantiles(TypeII.Green)
TypeII.Red    <- normalizeQuantiles(TypeII.Red)
TypeI.Green.M <- normalizeQuantiles(TypeI.Green.M)
TypeI.Green.U <- normalizeQuantiles(TypeI.Green.U)
TypeI.Red.M   <- normalizeQuantiles(TypeI.Red.M)
TypeI.Red.U   <- normalizeQuantiles(TypeI.Red.U)

# Compute beta values from normalised intensities
TypeII.betas      <- TypeII.Green  / (TypeII.Red + TypeII.Green + 100)
TypeI.Green.betas <- TypeI.Green.M / (TypeI.Green.M + TypeI.Green.U + 100)
TypeI.Red.betas   <- TypeI.Red.M   / (TypeI.Red.M + TypeI.Red.U + 100)

# Final beta matrix: [autosomal probes x samples]
beta <- as.matrix(rbind(TypeII.betas, TypeI.Green.betas, TypeI.Red.betas))

rm(TypeII.Green, TypeII.Red, TypeI.Green.M, TypeI.Green.U,
   TypeI.Red.M, TypeI.Red.U, TypeII.betas, TypeI.Green.betas, TypeI.Red.betas)

save(beta, file = "beta_QN_detP001_marker095.RData")

message("Pipeline complete. Final beta matrix: ",
        nrow(beta), " probes x ", ncol(beta), " samples.")

sessionInfo()
