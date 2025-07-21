######################
### INITIALIZATION ###
######################

# Load libraries #
library(MethylCallR, lib.loc = "/home/FCAM/arork/software/rlibs/v4.3.3/")
library(IlluminaHumanMethylationEPICv2manifest, lib.loc = "/home/FCAM/arork/software/rlibs/v4.3.3/")
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38, lib.loc = "/home/FCAM/arork/software/rlibs/v4.3.3/")
library(gmqn, lib.loc = "/home/FCAM/arork/software/rlibs/v4.3.3/")

# Set seed and working directory #
set.seed(123)
setwd("/core/cbc/core_projects/2025/jbrown2025/05_normalization/")

##########################
### DATA NORMALIZATION ###
##########################

# Read filtered idat data into a variable #
filtered.data <- readRDS("../03_filtering/filtered_data.rds")

# Normalize filtered methylation levels w/ Noob + BMIQ #
normalized.data <- MeCall.Norm(data = filtered.data, method = c("Noob", "BMIQ"))

#########################
### CONCLUDE ANALYSIS ###
#########################

# Write data to a RDS files #
saveRDS(normalized.data, "normalized_data.rds")

# Print any warnings and sessionInfo #
writeLines(capture.output(warnings()), "../logs/jbrown_s05_R_warnings.txt")
writeLines(capture.output(sessionInfo()), "../logs/jbrown_s05_R_sessionInfo.txt")

# Clear the environment for the next project #
rm(list = ls())

