################################################################################
#
#  Trans-Ancestry Meta-Analysis of EWAS Results
#  Software: METAL (https://genome.sph.umich.edu/wiki/METAL)
#
#  Description:
#  This script performs a sample-size-weighted p-value-based meta-analysis
#  (Stouffer's method) across ethnic-specific EWAS summary statistics using
#  METAL. It combines results from multiple ancestry groups and tests for
#  heterogeneity of effects across groups.
#
#  Method:
#    Stouffer's method combines z-scores (derived from p-values and effect
#    directions) weighted by sample size. This is preferred over inverse-
#    variance weighting when effect sizes may not be directly comparable
#    across studies (e.g. different ancestry groups or platforms).
#
#  Input files:
#    Tab-separated text files, one per ancestry group. Each file must contain
#    (at minimum) the following columns:
#      CpG   : CpG probe identifier (e.g. "cg00000029")
#      beta  : regression coefficient (effect size; used for direction only)
#      p     : p-value from the EWAS regression
#      N     : sample size (used as weight in Stouffer's method)
#
#  Output file:
#    A tab-separated results table (.tbl) with meta-analysis statistics
#    for each CpG, including:
#      - Combined p-value and z-score (Stouffer-weighted)
#      - Direction of effect across input studies (e.g. "+++" or "+-+")
#      - Heterogeneity test (Cochran's Q) p-value and I-squared
#
#  Usage:
#    metal meta_analysis.metal
#
#  Notes:
#    - Input files should NOT be genomic-inflation corrected prior to
#      meta-analysis if post-hoc correction is preferred.
#    - All input files must share the same column names and format.
#    - Update file paths below to match your directory structure.
#
################################################################################


# ==============================================================================
# 1. Global settings
# ==============================================================================

# Column delimiter in all input files
SEPARATOR   TAB

# Column name mapping:
#   MARKERLABEL  -> column containing the CpG probe ID
#   EFFECT       -> column containing the effect estimate (for direction)
#   PVALUE       -> column containing the association p-value
#   WEIGHTLABEL  -> column containing the sample size (used as weight)
MARKERLABEL CpG
EFFECT      beta
PVALUE      p
WEIGHTLABEL N

# Use sample-size-weighted Stouffer's method (as opposed to STDERR for
# inverse-variance-weighted fixed-effects meta-analysis)
SCHEME      SAMPLESIZE


# ==============================================================================
# 2. Process input files (one per ancestry group)
# ==============================================================================

# Each PROCESS command reads one ancestry-specific EWAS summary statistics file.
# Add or remove PROCESS lines to match the number of ancestry groups.

# -- Chinese --
PROCESS HELIOS_ewas_results_CHN.txt

# -- Malay --
PROCESS HELIOS_ewas_results_MLY.txt

# -- Indian --
PROCESS HELIOS_ewas_results_IND.txt


# ==============================================================================
# 3. Output and analysis
# ==============================================================================

# Output file prefix and suffix (METAL appends the suffix to the prefix)
OUTFILE MA_cIMTmean_ASIAN_HELIOS_CHNMLYIND_stouffer .tbl

# Perform heterogeneity analysis (Cochran's Q test) to assess consistency
# of effects across ancestry groups
ANALYZE HETEROGENEITY

QUIT
