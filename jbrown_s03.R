######################
### INITIALIZATION ###
######################

# Load libraries #
library(MethylCallR, lib.loc = "/home/FCAM/arork/software/rlibs/v4.3.3/")
library(IlluminaHumanMethylationEPICv2manifest, lib.loc = "/home/FCAM/arork/software/rlibs/v4.3.3/")

# Set seed and working directory #
set.seed(123)
setwd("/core/cbc/core_projects/2025/jbrown2025/03_filtering/")

######################
### DATA FILTERING ###
######################

# Read idat files into a variable #
unfiltered.data <- MeCall.ReadIdat(idat.dir = "../raw_data/",
  pd.file = "../supplement/methycallr_samples.csv",
  platform = "EPICv2")

# Filter out any low quality probes. Note: setting "detP.probe.cutoff" too close to 0.05 throws an #
# error. It is unclear why. We will simply approximate to keep as many quality probes as possible. #
filtered.data <- MeCall.Filtering(data = unfiltered.data,
  badSample.minfi = FALSE, # Do not remove samples
  detP.probe.cutoff = 0.0498, # Remove probes w/ detection p-values >= ~0.05
  detP.sample.cutoff = 1, # Do not remove samples
  NoCG = FALSE) # There is one non-CG sub-telomeric probe in our list

#########################
### CONCLUDE ANALYSIS ###
#########################

# Write data to RDS files #
saveRDS(unfiltered.data, "unfiltered_data.rds")
saveRDS(filtered.data, "filtered_data.rds")

# Print any warnings and sessionInfo #
writeLines(capture.output(warnings()), "../logs/jbrown_s03_R_warnings.txt")
writeLines(capture.output(sessionInfo()), "../logs/jbrown_s03_R_sessionInfo.txt")

# Clear the environment for the next project #
rm(list = ls())

