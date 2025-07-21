######################
### INITIALIZATION ###
######################

# Load libraries
library(MethylCallR, lib.loc = "/home/FCAM/arork/software/rlibs/v4.3.3/")
library(IlluminaHumanMethylationEPICv2manifest, lib.loc = "/home/FCAM/arork/software/rlibs/v4.3.3/")
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38, lib.loc = "/home/FCAM/arork/software/rlibs/v4.3.3/")
library(gmqn, lib.loc = "/home/FCAM/arork/software/rlibs/v4.3.3/")

# Set seed and working directory #
set.seed(123)
setwd("/core/cbc/core_projects/2025/jbrown2025/07_differential_methylation/")

#########################################
### DIFFERENTIAL METHYLATION ANALYSES ###
#########################################

#############
# Read Data #
#############

# Read normalized MethylCallR data #
normalized <- readRDS("../05_normalization/normalized_data.rds")

# Extract the familial group with one mother and two infants #
bvals <- normalized$beta[,1:3]
pd <- normalized$pd[1:3,]

# Create a design matrix #
design <- MeCall.MakeModelMatrix(pd = pd, covariate = NULL, interest = "Status", PCA.result = NULL,
  n.pc = 2, Cell.comp = NULL, celltypes = NULL)

# Conduct a differential methylation analysis at the probe level. Note: Is is atypical to use beta #
# values over M values in differential methylation analyses, but this appars to be the recommendation. #
dm.probes <- MeCall.DMP(meth = bvals, pd = pd, interest = "Status", design = design, cutoff.P = 1,
  multi.P = 'BH', arraytype = "EPICv2")

# Conduct a differential methylation analysis at the region level #
dm.regions <- MeCall.DMR(meth = bvals, pd = pd, interest = "Status", arraytype = "EPICv2")

#########################
### CONCLUDE ANALYSIS ###
#########################

# Write data to a RDS files #
saveRDS(dm.probes, "dm_probe_data.rds")
saveRDS(dm.regions, "dm_region_data.rds")

# Print any warnings and sessionInfo #
writeLines(capture.output(warnings()), "../logs/jbrown_s07_R_warnings.txt")
writeLines(capture.output(sessionInfo()), "../logs/jbrown_s07_R_sessionInfo.txt")

# Clear the environment for the next project #
rm(list = ls())

