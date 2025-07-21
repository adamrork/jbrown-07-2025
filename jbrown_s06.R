######################
### INITIALIZATION ###
######################

# Load libraries
library(tidyverse)

# Set seed and working directory #
set.seed(123)
setwd("~/Desktop/general/projects/clients/2025/003_judy_brown/rds_files/normalized_data/")

########################################
### NORMALIZATION QUALITY ASSESSMENT ###
########################################

#############
# Read Data #
#############

# Read sub-telomeric probes list #
probes <- as.vector(read.table("../../data/subtelomeric_probes.txt", header = T)[,1])

# Read normalized MethylCallR data #
normalized <- readRDS("normalized_data.rds")

all.norm <- list(normalized$beta, normalized$M)
names(all.norm) <- c("bvals", "mvals")

# Gather all the subtelomeric probes #
sbtl.norm <- list()

sbtl.norm[["bvals"]] <- all.norm[["bvals"]][which(rownames(all.norm[["bvals"]]) %>% gsub("_.*", "", .) %in% probes),]
sbtl.norm[["mvals"]] <- all.norm[["mvals"]][which(rownames(all.norm[["bvals"]]) %>% gsub("_.*", "", .) %in% probes),]

# Convert data from wide to long format #
all.norm <- lapply(all.norm, data.frame, check.names = FALSE) %>%
  lapply(pivot_longer, cols = 1:7) %>% lapply(data.frame)

sbtl.norm <- lapply(sbtl.norm, data.frame, check.names = FALSE) %>%
  lapply(pivot_longer, cols = 1:7) %>% lapply(data.frame)

# Rename columns #
colnames(all.norm[["bvals"]]) <- c("Samples", "BVals")
colnames(sbtl.norm[["bvals"]]) <- c("Samples", "BVals")

colnames(all.norm[["mvals"]]) <- c("Samples", "MVals")
colnames(sbtl.norm[["mvals"]]) <- c("Samples", "MVals")

##########################
### Data Vizualization ###
##########################

# Generate density plots and a fairly granular histogram of the raw beta-values from all probes #

# Create some simple plot titles and filename suffixes #
all.norm.plot.title <- paste0("All Normalized Probes")
all.norm.bval.fname <- paste0(all.norm.plot.title, " - Beta Value Density.pdf")
all.norm.mval.fname <- paste0(all.norm.plot.title, " - M Value Density.pdf")

sbtl.norm.plot.title <- paste0("All Normalized Subtelomeric Probes")
sbtl.norm.bval.fname <- paste0(sbtl.norm.plot.title, " - Beta Value Density.pdf")
sbtl.norm.mval.fname <- paste0(sbtl.norm.plot.title, " - M Value Density.pdf")

###############
# Beta Values #
###############

# Create density plots of the beta values #
all.norm.bval.dens <- ggplot(all.norm[["bvals"]], aes(x = BVals, fill = Samples, color = Samples)) +
  geom_density(alpha = 0.05) + xlab("Beta Value") + ylab("Density") + ggtitle(all.norm.plot.title) +
  theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5))

sbtl.norm.bval.dens <- ggplot(sbtl.norm[["bvals"]], aes(x = BVals, fill = Samples, color = Samples)) +
  geom_density(alpha = 0.05) + xlab("Beta Value") + ylab("Density") + ggtitle(sbtl.norm.plot.title) +
  theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5))

############
# M Values #
############

# Create density plots of the M values #
all.norm.mval.dens <- ggplot(all.norm[["mvals"]], aes(x = MVals, fill = Samples, color = Samples)) +
  geom_density(alpha = 0.05) + xlab("M Value") + ylab("Density") + ggtitle(all.norm.plot.title) +
  theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5))

sbtl.norm.mval.dens <- ggplot(sbtl.norm[["mvals"]], aes(x = MVals, fill = Samples, color = Samples)) +
  geom_density(alpha = 0.05) + xlab("M Value") + ylab("Density") + ggtitle(sbtl.norm.plot.title) +
  theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5))

# Save each density plot to a PDF #
ggsave(all.norm.bval.fname, plot = all.norm.bval.dens, device = "pdf", path = "../../figures/normalized_qa/",
       width = 8, height = 6, units = "in", dpi = 400)

ggsave(sbtl.norm.bval.fname, plot = sbtl.norm.bval.dens, device = "pdf", path = "../../figures/normalized_qa/",
       width = 8, height = 6, units = "in", dpi = 400)

ggsave(all.norm.mval.fname, plot = all.norm.mval.dens, device = "pdf", path = "../../figures/normalized_qa/",
       width = 8, height = 6, units = "in", dpi = 400)

ggsave(sbtl.norm.mval.fname, plot = sbtl.norm.mval.dens, device = "pdf", path = "../../figures/normalized_qa/",
       width = 8, height = 6, units = "in", dpi = 400)

#########################
### CONCLUDE ANALYSIS ###
#########################

# Print any warnings and sessionInfo #
writeLines(capture.output(warnings()), "../../logs/jbrown_s06_R_warnings.txt")
writeLines(capture.output(sessionInfo()), "../../logs/jbrown_s06_R_sessionInfo.txt")

# Clear the environment for the next project #
rm(list = ls())

#
# NOTE: This script was ran locally and interactively
#

