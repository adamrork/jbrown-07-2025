######################
### INITIALIZATION ###
######################

# Load libraries
library(tidyverse)

# Set seed and working directory #
set.seed(123)
setwd("~/Desktop/general/projects/clients/2025/003_judy_brown/rds_files/filtered_data/")

####################################
### FILTERING QUALITY ASSESSMENT ###
####################################

#############
# Read Data #
#############

# Read sub-telomeric probes list #
probes <- as.vector(read.table("../../data/subtelomeric_probes.txt", header = T)[,1])

# Read filtered MethylCallR data #
filtered <- readRDS("filtered_data.rds")

all.filt <- list(filtered$beta, filtered$M)
names(all.filt) <- c("bvals", "mvals")

# Gather all the subtelomeric probes #
sbtl.filt <- list()

sbtl.filt[["bvals"]] <- all.filt[["bvals"]][which(rownames(all.filt[["bvals"]]) %>% gsub("_.*", "", .) %in% probes),]
sbtl.filt[["mvals"]] <- all.filt[["mvals"]][which(rownames(all.filt[["bvals"]]) %>% gsub("_.*", "", .) %in% probes),]

# Gather some statistics on the number of filtered subtelomeric probes identified #
dim(sbtl.filt[["bvals"]])[1] # 765 probes

# Convert data from wide to long format #
all.filt <- lapply(all.filt, data.frame, check.names = FALSE) %>%
  lapply(pivot_longer, cols = 1:7) %>% lapply(data.frame)
  
sbtl.filt <- lapply(sbtl.filt, data.frame, check.names = FALSE) %>%
  lapply(pivot_longer, cols = 1:7) %>% lapply(data.frame)
  
# Rename columns #
colnames(all.filt[["bvals"]]) <- c("Samples", "BVals")
colnames(sbtl.filt[["bvals"]]) <- c("Samples", "BVals")

colnames(all.filt[["mvals"]]) <- c("Samples", "MVals")
colnames(sbtl.filt[["mvals"]]) <- c("Samples", "MVals")

##########################
### Data Vizualization ###
##########################

# Generate density plots and a fairly granular histogram of the raw beta-values from all probes #

# Create some simple plot titles and filename suffixes #
all.filt.plot.title <- paste0("All Filtered Probes")
all.filt.bval.fname <- paste0(all.filt.plot.title, " - Beta Value Density.pdf")
all.filt.mval.fname <- paste0(all.filt.plot.title, " - M Value Density.pdf")

sbtl.filt.plot.title <- paste0("All Filtered Subtelomeric Probes")
sbtl.filt.bval.fname <- paste0(sbtl.filt.plot.title, " - Beta Value Density.pdf")
sbtl.filt.mval.fname <- paste0(sbtl.filt.plot.title, " - M Value Density.pdf")

###############
# Beta Values #
###############
  
# Create density plots of the beta values #
all.filt.bval.dens <- ggplot(all.filt[["bvals"]], aes(x = BVals, fill = Samples, color = Samples)) +
  geom_density(alpha = 0.05) + xlab("Beta Value") + ylab("Density") + ggtitle(all.filt.plot.title) +
  theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5))

sbtl.filt.bval.dens <- ggplot(sbtl.filt[["bvals"]], aes(x = BVals, fill = Samples, color = Samples)) +
  geom_density(alpha = 0.05) + xlab("Beta Value") + ylab("Density") + ggtitle(sbtl.filt.plot.title) +
  theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5))
  
############
# M Values #
############
  
# Create density plots of the M values #
all.filt.mval.dens <- ggplot(all.filt[["mvals"]], aes(x = MVals, fill = Samples, color = Samples)) +
  geom_density(alpha = 0.05) + xlab("M Value") + ylab("Density") + ggtitle(all.filt.plot.title) +
  theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5))
  
sbtl.filt.mval.dens <- ggplot(sbtl.filt[["mvals"]], aes(x = MVals, fill = Samples, color = Samples)) +
  geom_density(alpha = 0.05) + xlab("M Value") + ylab("Density") + ggtitle(sbtl.filt.plot.title) +
  theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5))

# Save each density plot to a PDF #
ggsave(all.filt.bval.fname, plot = all.filt.bval.dens, device = "pdf", path = "../../figures/filtered_qa/",
       width = 8, height = 6, units = "in", dpi = 400)
  
ggsave(sbtl.filt.bval.fname, plot = sbtl.filt.bval.dens, device = "pdf", path = "../../figures/filtered_qa/",
        width = 8, height = 6, units = "in", dpi = 400)
  
ggsave(all.filt.mval.fname, plot = all.filt.mval.dens, device = "pdf", path = "../../figures/filtered_qa/",
       width = 8, height = 6, units = "in", dpi = 400)
  
ggsave(sbtl.filt.mval.fname, plot = sbtl.filt.mval.dens, device = "pdf", path = "../../figures/filtered_qa/",
        width = 8, height = 6, units = "in", dpi = 400)

#########################
### CONCLUDE ANALYSIS ###
#########################

# Print any warnings and sessionInfo #
writeLines(capture.output(warnings()), "../../logs/jbrown_s04_R_warnings.txt")
writeLines(capture.output(sessionInfo()), "../../logs/jbrown_s04_R_sessionInfo.txt")

# Clear the environment for the next project #
rm(list = ls())

#
# NOTE: This script was ran locally and interactively
#

