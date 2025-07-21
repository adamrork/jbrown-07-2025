######################
### INITIALIZATION ###
######################

# Load libraries
library(MethylCallR, lib.loc = "/home/FCAM/arork/software/rlibs/v4.3.3/")
library(sesame)
library(dplyr)

# Set seed and working directory #
set.seed(123)
setwd("/core/cbc/core_projects/2025/jbrown2025/01_write_rds/")

# Read in samples #
samples <- read.csv("../supplement/methycallr_samples.csv")

########################
### CREATE RDS FILES ###
########################

#########################
# Gather Data w/ SeSAMe #
#########################

# Obtain idat file paths #
idats <- searchIDATprefixes("../raw_data/")

# Interpolate beta values, pvalue with SeSAMe #
bvals <- openSesame(idats, func = getBetas) %>% data.frame(check.names = FALSE)
pvals <- openSesame(idats, func = pOOBAH, return.pval = TRUE) %>% data.frame(check.names = FALSE)
mvals <- BetaValueToMValue(bvals)

# Create a list for easier downstream handling #
sesame.data <- list(bvals, mvals, pvals)
names(sesame.data) <- c("bvals", "mvals", "pvals")

##############################
# Gather Data w/ MethylCallR #
##############################

# Interpolate beta values, pvalue with MethylCallR #
methylcallr.data <- MeCall.ReadIdat(idat.dir = "../raw_data/",
  pd.file = "../supplement/methycallr_samples.csv",
  platform = "EPICv2")

methylcallr.data <- list(methylcallr.data$beta, methylcallr.data$M, methylcallr.data$detP)
methylcallr.data <- lapply(methylcallr.data, data.frame, check.names = FALSE)
names(methylcallr.data) <- c("bvals", "mvals", "pvals")

#############################
# Clean and Sort Dataframes #
#############################

# Rename SeSAMe dataframe column names and sort columns to be in alphanumeric order by name #
for ( i in 1:3 ) {
  colnames(sesame.data[[i]]) <- samples[,2][match(colnames(sesame.data[[i]]), samples[,1])]
  sesame.data[[i]] <- sesame.data[[i]][, order(colnames(sesame.data[[i]]))]

  methylcallr.data[[i]] <- methylcallr.data[[i]][, order(colnames(methylcallr.data[[i]]))]
}

###################
# Write RDS files #
###################

# Write RDS files #
saveRDS(sesame.data, "unfiltered_sesame_data.rds")
saveRDS(methylcallr.data, "unfiltered_methylcallr_data.rds")

#########################
### CONCLUDE ANALYSIS ###
#########################

# Print any warnings and sessionInfo #
writeLines(capture.output(warnings()), "../logs/jbrown_s01_R_warnings.txt")
writeLines(capture.output(sessionInfo()), "../logs/jbrown_s01_R_sessionInfo.txt")

# Clear the environment for the next project #
rm(list = ls())

