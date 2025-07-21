######################
### INITIALIZATION ###
######################

# Load libraries
library(tidyverse)

# Set seed and working directory #
set.seed(123)
setwd("~/Desktop/general/projects/clients/2025/003_judy_brown/rds_files/normalized_data/")

########################################
### COMPARATIVE METHYLATION ANALYSES ###
########################################

#############
# Read Data #
#############

# Read samples #
samples <- read.csv("../../data/methylcallr_samples.csv", header = T)

# Read sub-telomeric probes list #
probes <- as.vector(read.table("../../data/subtelomeric_probes.txt", header = T)[,1])

# Read normalized MethylCallR data #
normalized <- readRDS("normalized_data.rds")

# Read differential methylation resutls #
dmp.results <- readRDS("../differential_methylation/dm_probe_data.rds")@Main.info
dmr.results <- readRDS("../differential_methylation/dm_region_data.rds")

# It's not clear that the p-values in the DMR results are adjusted, so do that here to be safe #
dmr.results$adjusted.p <- p.adjust(dmr.results$p.value, method = "BH")
dmr.results[which(dmr.results$seqnames == "chr19"),]

# Set variables containing sample names for each familial group #
group.a.ids <- c("11787M", "11787IP1", "71663IPII")
group.b.ids <- c("36149M", "36149IP")
group.c.ids <- c("50543M", "50543IPI")

# Split data into three subsets, by familial group #
bvals.a <- normalized$beta[,group.a.ids]
bvals.b <- normalized$beta[,group.b.ids]
bvals.c <- normalized$beta[,group.c.ids]

mvals.a <- normalized$M[,group.a.ids]
mvals.b <- normalized$M[,group.b.ids]
mvals.c <- normalized$M[,group.c.ids]

# Combine data into lists, one per familial group #
group.a <- list(bvals.a, mvals.a)
group.b <- list(bvals.b, mvals.b)
group.c <- list(bvals.c, mvals.c)

# Add names to the list elements #
names(group.a) <- names(group.b) <- names(group.c) <- c("bvals", "mvals")

# Group A - probe-wise delta beta and delta M values of Mother vs Infants #
group.a.deltaB <- group.a$bvals[,1] - rowMeans(group.a$bvals[,2:3]) %>% data.frame()
group.a.deltaM <- group.a$mvals[,1] - rowMeans(group.a$mvals[,2:3]) %>% data.frame()

# Group B - probe-wise delta beta and delta M values of Mother vs Infant #
group.b.deltaB <- group.b$bvals[,1] - group.b$bvals[,2] %>% data.frame()
group.b.deltaM <- group.b$mvals[,1] - group.b$mvals[,2] %>% data.frame()

# Group C - probe-wise delta beta and delta M values of Mother vs Infant #
group.c.deltaB <- group.c$bvals[,1] - group.c$bvals[,2] %>% data.frame()
group.c.deltaM <- group.c$mvals[,1] - group.c$mvals[,2] %>% data.frame()

# Gather results into lists, one per familial group #
group.a.comparisons <- list(group.a.deltaB, group.a.deltaM)
group.b.comparisons <- list(group.b.deltaB, group.b.deltaM)
group.c.comparisons <- list(group.c.deltaB, group.c.deltaM)

# Add names to the list elements #
names(group.a.comparisons[[1]]) <- names(group.b.comparisons[[1]]) <- names(group.c.comparisons[[1]]) <- "deltaB"
names(group.a.comparisons[[2]]) <- names(group.b.comparisons[[2]]) <- names(group.c.comparisons[[2]]) <- "deltaM"

# Rank each set of results (largest = 1, smallest = nrows) #
for ( i in 1:2 ) {
  group.a.comparisons[[i]]$rank <- rank(-group.a.comparisons[[i]][,1], ties.method = "first")
  group.a.comparisons[[i]]<- group.a.comparisons[[i]][order(group.a.comparisons[[i]]$rank),]
  rownames(group.a.comparisons[[i]]) <- rownames(group.a.comparisons[[i]]) %>% gsub("_.*", "", .)

  group.b.comparisons[[i]]$rank <- rank(-group.b.comparisons[[i]][,1], ties.method = "first")
  group.b.comparisons[[i]]<- group.b.comparisons[[i]][order(group.b.comparisons[[i]]$rank),]
  rownames(group.b.comparisons[[i]]) <- rownames(group.b.comparisons[[i]]) %>% gsub("_.*", "", .)

  group.c.comparisons[[i]]$rank <- rank(-group.c.comparisons[[i]][,1], ties.method = "first")
  group.c.comparisons[[i]]<- group.c.comparisons[[i]][order(group.c.comparisons[[i]]$rank),]
  rownames(group.c.comparisons[[i]]) <- rownames(group.c.comparisons[[i]]) %>% gsub("_.*", "", .)
}

# Subset the data to include only subtelomeric probes #
group.a.str.comparisons <- group.b.str.comparisons <- group.c.str.comparisons <- list()

# Susbet the probes and fix the ranks to be subtelomeric subset-specific. Note: Probes should already be #
# ordered and probably don't need to be re-ranked, but we'll do so here in case it is useful downstream. #
for ( i in 1:2 ) {
  group.a.str.comparisons[[i]] <- group.a.comparisons[[i]][which(rownames(group.a.comparisons[[i]]) %in% probes),]
  group.a.str.comparisons[[i]]$st.rank <- rank(-group.a.str.comparisons[[i]][,1], ties.method = "first")
  group.a.str.comparisons[[i]] <- group.a.str.comparisons[[i]][order(group.a.str.comparisons[[i]]$st.rank),]

  group.b.str.comparisons[[i]] <- group.b.comparisons[[i]][which(rownames(group.b.comparisons[[i]]) %in% probes),]
  group.b.str.comparisons[[i]]$st.rank <- rank(-group.b.str.comparisons[[i]][,1], ties.method = "first")
  group.b.str.comparisons[[i]] <- group.b.str.comparisons[[i]][order(group.b.str.comparisons[[i]]$st.rank),]

  group.c.str.comparisons[[i]] <- group.c.comparisons[[i]][which(rownames(group.c.comparisons[[i]]) %in% probes),]
  group.c.str.comparisons[[i]]$st.rank <- rank(-group.c.str.comparisons[[i]][,1], ties.method = "first")
  group.c.str.comparisons[[i]] <- group.c.str.comparisons[[i]][order(group.c.str.comparisons[[i]]$st.rank),]

}

# Add names to the list elements #
names(group.a.str.comparisons) <- names(group.b.str.comparisons) <- names(group.c.str.comparisons) <- c("deltaB", "deltaM")

# For Group A, gather significantlly differentially methylated probes #
sig.dmp.results <- dmp.results[which(dmp.results$adjusted.p < 0.05),]
sig.dmr.results <- dmr.results[which(dmr.results$adjusted.p < 0.05),]

# Add ranks to the differentially methylated probe results
sig.dmp.results$rank <- rank(-sig.dmp.results$deltabeta, ties.method = "first")
sig.dmp.results <- sig.dmp.results[order(sig.dmp.results$rank),]

# Extract subtelomeric probes from the DM probe list #
sig.st.dmp.results <- sig.dmp.results[which(sig.dmp.results$CpGid %in% probes),]

# Transform the differentially methylated probe results to fit the ranked list format #
group.a.dmp.results <- data.frame(sig.st.dmp.results$deltabeta)
colnames(group.a.dmp.results) <- "deltaB"
rownames(group.a.dmp.results) <- rownames(sig.st.dmp.results) %>% gsub("_.*", "", .)

group.a.dmp.results$st.rank <- rank(-group.a.dmp.results[,1], ties.method = "first")
group.a.dmp.results <- group.a.dmp.results[order(group.a.dmp.results$st.rank),]

##########################
### DATA VISUALIZATION ###
##########################

# Before we do any data visualization, add some factors to the data to keep plotting functions clean #

# A = Δβ > 0.1
# B = 0.1 >= Δβ > 0.0
# C = 0.0 < Δβ <= -0.1
# D = Δβ < -0.1

group.a.str.comparisons[[1]]$Status <- ifelse(group.a.str.comparisons[[1]][,1] > 0,
  ifelse(group.a.str.comparisons[[1]][,1] > 0.1, "A", "B"),
  ifelse(group.a.str.comparisons[[1]][,1] < -0.1, "D", "C"))

group.b.str.comparisons[[1]]$Status <- ifelse(group.b.str.comparisons[[1]][,1] > 0,
  ifelse(group.b.str.comparisons[[1]][,1] > 0.1, "A", "B"),
  ifelse(group.b.str.comparisons[[1]][,1] < -0.1, "D", "C"))

group.c.str.comparisons[[1]]$Status <- ifelse(group.c.str.comparisons[[1]][,1] > 0,
  ifelse(group.c.str.comparisons[[1]][,1] > 0.1, "A", "B"),
  ifelse(group.c.str.comparisons[[1]][,1] < -0.1, "D", "C"))

group.a.dmp.results$Status <- ifelse(group.a.dmp.results[,1] > 0,
  ifelse(group.a.dmp.results[,1] > 0.1, "A", "B"),
  ifelse(group.a.dmp.results[,1] < -0.1, "D", "C"))

# Gather the number of "differentially methylated" subtelomeric probes #
sum(group.a.str.comparisons[[1]]$Status == "A" | group.a.str.comparisons[[1]]$Status == "D")
sum(group.b.str.comparisons[[1]]$Status == "A" | group.b.str.comparisons[[1]]$Status == "D")
sum(group.c.str.comparisons[[1]]$Status == "A" | group.c.str.comparisons[[1]]$Status == "D")

# Hypermethylated in X (X vs. Y) #
sum(group.a.str.comparisons[[1]]$Status == "A")
sum(group.b.str.comparisons[[1]]$Status == "A")
sum(group.c.str.comparisons[[1]]$Status == "A")

# Hypomethylated in X (X vs. Y) #
sum(group.a.str.comparisons[[1]]$Status == "D")
sum(group.b.str.comparisons[[1]]$Status == "D")
sum(group.c.str.comparisons[[1]]$Status == "D")

# Create a list of annotations #
hypr.lab <- "Top 5 Hypermethylated Subtelomeric Probes"
hypo.lab <- "Top 5 Hypomethylated Subtelomeric Probes"

# a.top and a.btm have been checked against the DMP results to ensure a 1:1 match #
a.top <- c("cg09959355 : Δβ = 0.814", "cg02365900 : Δβ = 0.807", "cg07876904 : Δβ = 0.801",
           "cg22503684 : Δβ = 0.784", "cg16032694 : Δβ = 0.774")
a.btm <- c("cg00169401 : Δβ = -0.779", "cg00804354 : Δβ = -0.694", "cg17855943 : Δβ = -0.669",
           "cg23825768 : Δβ = -0.669", "cg02345255 : Δβ = -0.647")
b.top <- c("cg09959355 : Δβ = 0.654", "cg12284090 : Δβ = 0.588", "cg07013148 : Δβ = 0.426",
           "cg01129115 : Δβ = 0.406", "cg17709815 : Δβ = 0.400")
b.btm <- c("cg21767790 : Δβ = -0.454", "cg17053016 : Δβ = -0.436", "cg24565620 : Δβ = -0.360",
           "cg17001765 : Δβ = -0.287", "cg09642964 : Δβ = -0.258")
c.top <- c("cg23652859 : Δβ = 0.843", "cg09959355 : Δβ = 0.806", "cg02152128 : Δβ = 0.786",
           "cg09477755 : Δβ = 0.763", "cg21042248 : Δβ = 0.732")
c.btm <- c("cg03611028 : Δβ = -0.853", "cg26617156 : Δβ = -0.851", "cg17053016 : Δβ = -0.833",
           "cg25722686 : Δβ = -0.823", "cg10120070 : Δβ = -0.819")

# Plot ranked list of probes by delta beta values. Note: the method of adding #
# annotations employed here may be a tad cumbersome, but it is at least explicit #
group.a.dm.hist <- ggplot(group.a.dmp.results, aes(x = st.rank, y = deltaB, color = Status, fill = Status)) +
  geom_col() +
  scale_color_manual(labels = c("Δβ > 0.1", "Δβ < -0.1"),
                     values = c("#FF6363", "#65B3DB")) +
  scale_fill_manual(labels = c("Δβ > 0.1", "Δβ < -0.1"),
                    values = c("#FF6363", "#65B3DB")) +
  xlab("Rank") + ylab ("Δβ") + ggtitle("DM Subtelomeric Probes - 11787M vs. 11787IP1 & 71663IPII") +
  annotate(geom = "text", x = 90, y = 0.72, size = 4, color = "#CC1B1B", label = hypr.lab) +
  annotate(geom = "text", x = 90, y = 0.65, size = 4, color = "#CC1B1B", label = a.top[1]) +
  annotate(geom = "text", x = 90, y = 0.58, size = 4, color = "#CC1B1B", label = a.top[2]) +
  annotate(geom = "text", x = 90, y = 0.51, size = 4, color = "#CC1B1B", label = a.top[3]) +
  annotate(geom = "text", x = 90, y = 0.44, size = 4, color = "#CC1B1B", label = a.top[4]) +
  annotate(geom = "text", x = 90, y = 0.37, size = 4, color = "#CC1B1B", label = a.top[5]) +
  annotate(geom = "text", x = 40, y = -0.48, size = 4, color = "#316FCC", label = hypo.lab) +
  annotate(geom = "text", x = 40, y = -0.55, size = 4, color = "#316FCC", label = a.btm[1]) +
  annotate(geom = "text", x = 40, y = -0.62, size = 4, color = "#316FCC", label = a.btm[2]) +
  annotate(geom = "text", x = 40, y = -0.69, size = 4, color = "#316FCC", label = a.btm[3]) +
  annotate(geom = "text", x = 40, y = -0.76, size = 4, color = "#316FCC", label = a.btm[4]) +
  annotate(geom = "text", x = 40, y = -0.83, size = 4, color = "#316FCC", label = a.btm[5]) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 16),
                     axis.title = element_text(size = 14), axis.text = element_text(size = 10),
                     legend.title = element_text(size = 12), legend.text = element_text(size = 10))

group.a.hist <- ggplot(group.a.str.comparisons[[1]], aes(x = st.rank, y = deltaB, color = Status, fill = Status)) +
  geom_col() +
  scale_color_manual(labels = c("Δβ > 0.1", "0.1 ≥ Δβ > 0.0", "0.0 < Δβ ≤ -0.1", "Δβ < -0.1"),
                     values = c("#FF6363", "#FB9A99", "#A6CEE3", "#65B3DB")) +
  scale_fill_manual(labels = c("Δβ > 0.1", "0.1 ≥ Δβ > 0.0", "0.0 < Δβ ≤ -0.1", "Δβ < -0.1"),
                    values = c("#FF6363", "#FB9A99", "#A6CEE3", "#65B3DB")) +
  xlab("Rank") + ylab ("Δβ") + ggtitle("Subtelomeric Probes - 11787M vs. 11787IP1 & 71663IPII") + labs(tag = "5") +
  annotate(geom = "text", x = 525, y = 0.72, size = 4, color = "#CC1B1B", label = hypr.lab) +
  annotate(geom = "text", x = 525, y = 0.65, size = 4, color = "#CC1B1B", label = a.top[1]) +
  annotate(geom = "text", x = 525, y = 0.58, size = 4, color = "#CC1B1B", label = a.top[2]) +
  annotate(geom = "text", x = 525, y = 0.51, size = 4, color = "#CC1B1B", label = a.top[3]) +
  annotate(geom = "text", x = 525, y = 0.44, size = 4, color = "#CC1B1B", label = a.top[4]) +
  annotate(geom = "text", x = 525, y = 0.37, size = 4, color = "#CC1B1B", label = a.top[5]) +
  annotate(geom = "text", x = 225, y = -0.48, size = 4, color = "#316FCC", label = hypo.lab) +
  annotate(geom = "text", x = 225, y = -0.55, size = 4, color = "#316FCC", label = a.btm[1]) +
  annotate(geom = "text", x = 225, y = -0.62, size = 4, color = "#316FCC", label = a.btm[2]) +
  annotate(geom = "text", x = 225, y = -0.69, size = 4, color = "#316FCC", label = a.btm[3]) +
  annotate(geom = "text", x = 225, y = -0.76, size = 4, color = "#316FCC", label = a.btm[4]) +
  annotate(geom = "text", x = 225, y = -0.83, size = 4, color = "#316FCC", label = a.btm[5]) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 16),
                     axis.title = element_text(size = 14), axis.text = element_text(size = 10),
                     legend.title = element_text(size = 12), legend.text = element_text(size = 10),
                     plot.tag = element_text(size = 25), plot.tag.position = "topleft")

group.b.hist <- ggplot(group.b.str.comparisons[[1]], aes(x = st.rank, y = deltaB, color = Status, fill = Status)) +
  geom_col() +
  scale_color_manual(labels = c("Δβ > 0.1", "0.1 ≥ Δβ > 0.0", "0.0 < Δβ ≤ -0.1", "Δβ < -0.1"),
                     values = c("#FF6363", "#FB9A99", "#A6CEE3", "#65B3DB")) +
  scale_fill_manual(labels = c("Δβ > 0.1", "0.1 ≥ Δβ > 0.0", "0.0 < Δβ ≤ -0.1", "Δβ < -0.1"),
                    values = c("#FF6363", "#FB9A99", "#A6CEE3", "#65B3DB")) +
  xlab("Rank") + ylab ("Δβ") + ggtitle("Subtelomeric Probes - 36149M vs. 36149IP") + labs(tag = "6") +
  annotate(geom = "text", x = 525, y = 0.60, size = 4, color = "#CC1B1B", label = hypr.lab) +
  annotate(geom = "text", x = 525, y = 0.55, size = 4, color = "#CC1B1B", label = b.top[1]) +
  annotate(geom = "text", x = 525, y = 0.50, size = 4, color = "#CC1B1B", label = b.top[2]) +
  annotate(geom = "text", x = 525, y = 0.45, size = 4, color = "#CC1B1B", label = b.top[3]) +
  annotate(geom = "text", x = 525, y = 0.40, size = 4, color = "#CC1B1B", label = b.top[4]) +
  annotate(geom = "text", x = 525, y = 0.35, size = 4, color = "#CC1B1B", label = b.top[5]) +
  annotate(geom = "text", x = 225, y = -0.20, size = 4, color = "#316FCC", label = hypo.lab) +
  annotate(geom = "text", x = 225, y = -0.25, size = 4, color = "#316FCC", label = b.btm[1]) +
  annotate(geom = "text", x = 225, y = -0.30, size = 4, color = "#316FCC", label = b.btm[2]) +
  annotate(geom = "text", x = 225, y = -0.35, size = 4, color = "#316FCC", label = b.btm[3]) +
  annotate(geom = "text", x = 225, y = -0.40, size = 4, color = "#316FCC", label = b.btm[4]) +
  annotate(geom = "text", x = 225, y = -0.45, size = 4, color = "#316FCC", label = b.btm[5]) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 16),
                     axis.title = element_text(size = 14), axis.text = element_text(size = 10),
                     legend.title = element_text(size = 12), legend.text = element_text(size = 10),
                     plot.tag = element_text(size = 25), plot.tag.position = "topleft")

group.c.hist <- ggplot(group.c.str.comparisons[[1]], aes(x = st.rank, y = deltaB, color = Status, fill = Status)) +
  geom_col() +
  scale_color_manual(labels = c("Δβ > 0.1", "0.1 ≥ Δβ > 0.0", "0.0 < Δβ ≤ -0.1", "Δβ < -0.1"),
                     values = c("#FF6363", "#FB9A99", "#A6CEE3", "#65B3DB")) +
  scale_fill_manual(labels = c("Δβ > 0.1", "0.1 ≥ Δβ > 0.0", "0.0 < Δβ ≤ -0.1", "Δβ < -0.1"),
                    values = c("#FF6363", "#FB9A99", "#A6CEE3", "#65B3DB")) +
  xlab("Rank") + ylab ("Δβ") + ggtitle("Subtelomeric Probes - 50543M vs. 50543IPI") + labs(tag = "7") +
  annotate(geom = "text", x = 525, y = 0.72, size = 4, color = "#CC1B1B", label = hypr.lab) +
  annotate(geom = "text", x = 525, y = 0.65, size = 4, color = "#CC1B1B", label = c.top[1]) +
  annotate(geom = "text", x = 525, y = 0.58, size = 4, color = "#CC1B1B", label = c.top[2]) +
  annotate(geom = "text", x = 525, y = 0.51, size = 4, color = "#CC1B1B", label = c.top[3]) +
  annotate(geom = "text", x = 525, y = 0.44, size = 4, color = "#CC1B1B", label = c.top[4]) +
  annotate(geom = "text", x = 525, y = 0.37, size = 4, color = "#CC1B1B", label = c.top[5]) +
  annotate(geom = "text", x = 225, y = -0.48, size = 4, color = "#316FCC", label = hypo.lab) +
  annotate(geom = "text", x = 225, y = -0.55, size = 4, color = "#316FCC", label = c.btm[1]) +
  annotate(geom = "text", x = 225, y = -0.62, size = 4, color = "#316FCC", label = c.btm[2]) +
  annotate(geom = "text", x = 225, y = -0.69, size = 4, color = "#316FCC", label = c.btm[3]) +
  annotate(geom = "text", x = 225, y = -0.76, size = 4, color = "#316FCC", label = c.btm[4]) +
  annotate(geom = "text", x = 225, y = -0.83, size = 4, color = "#316FCC", label = c.btm[5]) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 16),
                     axis.title = element_text(size = 14), axis.text = element_text(size = 10),
                     legend.title = element_text(size = 12), legend.text = element_text(size = 10),
                     plot.tag = element_text(size = 25), plot.tag.position = "topleft")

ggsave("DM Subtelomeric Probes - 11787M vs. 11787IP1 & 71663IPII.png", path = "../../figures/differential_methylation/",
       plot = group.a.dm.hist, device = "png", width = 11, height = 5, units = "in", dpi = 500)
ggsave("Subtelomeric Probes - 11787M vs. 11787IP1 & 71663IPII.png", path = "../../figures/differential_methylation/",
       plot = group.a.hist, device = "png", width = 11, height = 5, units = "in", dpi = 500)
ggsave("Subtelomeric Probes - 36149M vs. 36149IP.png", path = "../../figures/differential_methylation/",
       plot = group.b.hist, device = "png", width = 11, height = 5, units = "in", dpi = 500)
ggsave("Subtelomeric Probes - 50543M vs. 50543IPI.png", path = "../../figures/differential_methylation/",
       plot = group.c.hist, device = "png", width = 11, height = 5, units = "in", dpi = 500)

#########################
### CONCLUDE ANALYSIS ###
#########################

# Write several relevant tables of results to files #
write.table(dmp.results, "../../tsv_files/differential_methylation/Group A - Differentially Methylated Probes.tsv")
write.table(dmr.results, "../../tsv_files/differential_methylation/Group A - Differentially Methylated Regions.tsv")
write.table(group.a.str.comparisons[[1]][,1:3], "../../tsv_files/differential_methylation/Group A - Ranked delta betas.tsv")
write.table(group.b.str.comparisons[[1]][,1:3], "../../tsv_files/differential_methylation/Group B - Ranked delta betas.tsv")
write.table(group.c.str.comparisons[[1]][,1:3], "../../tsv_files/differential_methylation/Group C - Ranked delta betas.tsv")

# Print any warnings and sessionInfo #
writeLines(capture.output(warnings()), "../../logs/jbrown_s08_R_warnings.txt")
writeLines(capture.output(sessionInfo()), "../../logs/jbrown_s08_R_sessionInfo.txt")

# Clear the environment for the next project #
rm(list = ls())

#
# NOTE: This script was ran locally and interactively
#

