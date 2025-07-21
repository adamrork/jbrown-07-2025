######################
### INITIALIZATION ###
######################

# Load libraries
library(tidyverse)

# Set seed and working directory #
set.seed(123)
setwd("~/Desktop/general/projects/clients/2025/003_judy_brown/rds_files/unfiltered_data/")

##################################
### INITIAL QUALITY ASSESSMENT ###
##################################

#############
# Read Data #
#############

# Read sub-telomeric probes list #
probes <- as.vector(read.table("../../data/subtelomeric_probes.txt", header = T)[,1])

# Read SeSAMe and MethylCallR data #
sme <- readRDS("unfiltered_sesame_data.rds")
mcr <- readRDS("unfiltered_methylcallr_data.rds")

# Gather all the best-perfoming probes #
best.sme <- best.mcr <- list()

# Remove probes not having interpolated beta values in all samples #
best.sme[["bvals"]] <- sme[["bvals"]][complete.cases(sme[["bvals"]]),]
best.mcr[["bvals"]] <- mcr[["bvals"]][complete.cases(mcr[["bvals"]]),]

# Identify probes which have detection p-values < 0.05 in all samples #
sig.sme.probes <- rownames(sme[["pvals"]])[apply(sme[["pvals"]] < 0.05, 1, all)]
sig.mcr.probes <- rownames(mcr[["pvals"]])[apply(mcr[["pvals"]] < 0.05, 1, all)]

# Keep only those probes which have detection p-values < 0.05 in all samples #
best.sme[["bvals"]] <- best.sme[["bvals"]][which(rownames(best.sme[["bvals"]]) %in% sig.sme.probes),]
best.mcr[["bvals"]] <- best.mcr[["bvals"]][which(rownames(best.mcr[["bvals"]]) %in% sig.mcr.probes),]

for ( i in 2:3 ) {
  best.sme[[i]] <- sme[[i]][which(rownames(sme[[i]]) %in% rownames(best.sme[["bvals"]])),]
  best.mcr[[i]] <- mcr[[i]][which(rownames(mcr[[i]]) %in% rownames(best.mcr[["bvals"]])),]
}

names(best.sme) <- names(best.mcr) <- c("bvals", "mvals", "pvals")

# Gather all the subtelomeric probes #
sbtl.sme <- sbtl.mcr <- list()

for ( i in 1:3 ) {
  sbtl.sme[[i]] <- sme[[i]][which(rownames(sme[[i]]) %>% gsub("_.*", "", .) %in% probes),]
  sbtl.mcr[[i]] <- mcr[[i]][which(rownames(mcr[[i]]) %>% gsub("_.*", "", .) %in% probes),]
}

names(sbtl.sme) <- names(sbtl.mcr) <- c("bvals", "mvals", "pvals")

# Gather all the best-perfoming subtelomeric probes #
best.sbtl.sme <- best.sbtl.mcr <- list()

for ( i in 1:3 ) {
  best.sbtl.sme[[i]] <- best.sme[[i]][which(rownames(best.sme[[i]]) %>% gsub("_.*", "", .) %in% probes),]
  best.sbtl.mcr[[i]] <- best.mcr[[i]][which(rownames(best.mcr[[i]]) %>% gsub("_.*", "", .) %in% probes),]
}

names(best.sbtl.sme) <- names(best.sbtl.mcr) <- c("bvals", "mvals", "pvals")

# Create a list of lists for easier handling #
sme.lists <- list(sme, best.sme, sbtl.sme, best.sbtl.sme)
mcr.lists <- list(mcr, best.mcr, sbtl.mcr, best.sbtl.mcr)

names(sme.lists) <- c("All", "Best", "Subtelomeric", "Best Subtelomeric")
names(mcr.lists) <- c("All", "Best", "Subtelomeric", "Best Subtelomeric")

# Convert data from wide to long format #
for ( i in 1:4 ) {
  sme.lists[[i]] <- lapply(sme.lists[[i]], data.frame, check.names = FALSE) %>%
    lapply(pivot_longer, cols = 1:7) %>% lapply(data.frame)
  
  mcr.lists[[i]] <- lapply(mcr.lists[[i]], data.frame, check.names = FALSE) %>%
    lapply(pivot_longer, cols = 1:7) %>% lapply(data.frame)
  
  # Rename columns #
  colnames(sme.lists[[i]][["bvals"]]) <- colnames(mcr.lists[[i]][["bvals"]]) <- c("Samples", "BVals")
  colnames(sme.lists[[i]][["mvals"]]) <- colnames(mcr.lists[[i]][["mvals"]]) <- c("Samples", "MVals")
  colnames(sme.lists[[i]][["pvals"]]) <- colnames(mcr.lists[[i]][["pvals"]]) <- c("Samples", "PVals")
}

# Gather some statistics on the number of subtelomeric probes identified #
dim(sbtl.sme$bvals)[1] # 1266 probes
dim(sbtl.mcr$bvals)[1] # 1266 probes

##########################
### Data Vizualization ###
##########################

# Generate density plots and a fairly granular histogram of the raw beta-values from all probes #
for (lst in names(sme.lists)) {
  
  # Create some simple plot titles and filename suffixes #
  sme.plot.title <- paste0("SeSAMe - ", lst, " Probes")
  sme.bval.fname <- paste0(sme.plot.title, " - Beta Value Density.pdf")
  sme.mval.fname <- paste0(sme.plot.title, " - M Value Density.pdf")
  sme.pval.fname <- paste0(sme.plot.title, " - P Value Histogram.pdf")
  
  mcr.plot.title <- paste0("MethylCallR - ", lst, " Probes")
  mcr.bval.fname <- paste0(mcr.plot.title, " - Beta Value Density.pdf")
  mcr.mval.fname <- paste0(mcr.plot.title, " - M Value Density.pdf")
  mcr.pval.fname <- paste0(mcr.plot.title, " - P Value Histogram.pdf")
  
  ###############
  # Beta Values #
  ###############
  
  # Create density plots of the beta values #
  sme.bval.dens <- ggplot(sme.lists[[lst]][["bvals"]], aes(x = BVals, fill = Samples, color = Samples)) +
    geom_density(alpha = 0.05) + xlab("Beta Value") + ylab("Density") + ggtitle(sme.plot.title) +
    theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  mcr.bval.dens <- ggplot(mcr.lists[[lst]][["bvals"]], aes(x = BVals, fill = Samples, color = Samples)) +
    geom_density(alpha = 0.05) + xlab("Beta Value") + ylab("Density") + ggtitle(mcr.plot.title) +
    theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  ############
  # M Values #
  ############
  
  # Create density plots of the M values #
  sme.mval.dens <- ggplot(sme.lists[[lst]][["mvals"]], aes(x = MVals, fill = Samples, color = Samples)) +
    geom_density(alpha = 0.05) + xlab("M Value") + ylab("Density") + ggtitle(sme.plot.title) +
    theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  mcr.mval.dens <- ggplot(mcr.lists[[lst]][["mvals"]], aes(x = MVals, fill = Samples, color = Samples)) +
    geom_density(alpha = 0.05) + xlab("M Value") + ylab("Density") + ggtitle(mcr.plot.title) +
    theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  ############
  # p-values #
  ############
  
  # Create density plots of the detection p-values #
  sme.pval.dens <- ggplot(sme.lists[[lst]][["pvals"]], aes(x = PVals, fill = Samples)) +
    geom_histogram(position = "identity", alpha = 0.5, bins = 100) +
    xlab("p-value") + ylab("Density") + ggtitle(sme.plot.title) +
    theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  mcr.pval.dens <- ggplot(mcr.lists[[lst]][["pvals"]], aes(x = PVals, fill = Samples)) +
    geom_histogram(position = "identity", alpha = 0.5, bins = 100) +
    xlab("p-value") + ylab("Density") + ggtitle(mcr.plot.title) +
    theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  # Save each density plot to a PDF #
  ggsave(sme.bval.fname, plot = sme.bval.dens, device = "pdf", path = "../../figures/initial_qa/",
         width = 8, height = 6, units = "in", dpi = 400)
  
  ggsave(mcr.bval.fname, plot = mcr.bval.dens, device = "pdf", path = "../../figures/initial_qa/",
         width = 8, height = 6, units = "in", dpi = 400)
  
  ggsave(sme.mval.fname, plot = sme.mval.dens, device = "pdf", path = "../../figures/initial_qa/",
         width = 8, height = 6, units = "in", dpi = 400)
  
  ggsave(mcr.mval.fname, plot = mcr.mval.dens, device = "pdf", path = "../../figures/initial_qa/",
         width = 8, height = 6, units = "in", dpi = 400)
  
  ggsave(sme.pval.fname, plot = sme.pval.dens, device = "pdf", path = "../../figures/initial_qa/",
         width = 8, height = 6, units = "in", dpi = 400)
  
  ggsave(mcr.pval.fname, plot = mcr.pval.dens, device = "pdf", path = "../../figures/initial_qa/",
         width = 8, height = 6, units = "in", dpi = 400)
}

#########################
### CONCLUDE ANALYSIS ###
#########################

# Print any warnings and sessionInfo #
writeLines(capture.output(warnings()), "../../logs/jbrown_s02_R_warnings.txt")
writeLines(capture.output(sessionInfo()), "../../logs/jbrown_s02_R_sessionInfo.txt")

# Clear the environment for the next project #
rm(list = ls())

#
# NOTE: This script was ran locally and interactively
#

