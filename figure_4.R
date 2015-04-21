
# script to reproduce the analysis and figures from Bulgarelli et al., 2015
#
# if you use any of the following code, please cite:
#
# Davide Bulgarelli, Ruben Garrido-Oter, Philipp C. Münch, Aaron Weiman,
# Johannes Dröge, Yao Pan, Alice C. McHardy, Paul Schulze-Lefert
# Structure and Function of the Bacterial Root Microbiota in Wild
# and Domesticated Barley. Cell Host and Microbe, 2015
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

### root-enriched barley and A. thaliana microbiota comparison

library(ape)

# cleanup
rm(list = ls())

# source functions
source("tree.plot.func.R")

# compile c function and load
system("R CMD SHLIB reorder_phylo.c")
dyn.load("reorder_phylo.so")

# input files
barley.design.file <- "./data/mapping.txt"
barley.count.matrix.file <- "./data/count_matrix.txt"
at.design.file <- "./data/at_mapping.txt"
at.count.matrix.file <- "./data/at_count_matrix.txt"
at_rotus_ids.file <- "./data/at_golm_rotus_NCBI_IDs.txt"
at_rotus_names.file <- "./data/at_golm_rotus_NCBI_names.txt"
barley_rotus_ids.file <- "./data/barley_rotus_NCBI_IDs.txt"
barley_rotus_names.file <- "./data/barley_rotus_NCBI_names.txt"
#
# load data
barley.design <- read.table(file=barley.design.file, sep="\t", header=T)
barley.count_matrix <- read.table(file=barley.count.matrix.file, sep="\t")
at.design <- read.table(file=at.design.file, sep="\t", header=T)
at.count_matrix <- read.table(file=at.count.matrix.file, sep="\t")
barley_rotus_ids <- read.table(file=barley_rotus_ids.file, sep="\t")
barley_rotus_names <- read.table(file=barley_rotus_names.file, sep="\t")
at_rotus_ids <- read.table(file=at_rotus_ids.file, sep="\t")
at_rotus_names <- read.table(file=at_rotus_names.file, sep="\t")

### BARLEY

## get OTU IDs
barley_rotus <- as.character(barley_rotus_names[, 1])

# normalize count matrix
barley_abundance_matrix <- apply(barley.count_matrix, 2, function(xx) xx / sum(xx))

# subset root samples and root OTUs
barley_root_samples <- as.character(barley.design$SampleID[barley.design$Compartment=="Root"])
idx <- colnames(barley_abundance_matrix) %in% barley_root_samples
barley_rotus_abundance_matrix <- barley_abundance_matrix[, idx]
idx <- rownames(barley_abundance_matrix) %in% barley_rotus
barley_rotus_abundance_matrix <- barley_rotus_abundance_matrix[idx, ]

# perc. of rOTUs w.r.t. to the whole root community
barley_rotus_perc <- mean(colSums(barley_rotus_abundance_matrix))

# calculate mean OTU abundances
barley_mean_rotus_abundances <- apply(barley_rotus_abundance_matrix, 1, mean)

barley_rotus <- data.frame(ncbi_name=as.character(barley_rotus_names[, 2]),
                           ncbi_id=barley_rotus_ids[, 2],
                           abundance=barley_mean_rotus_abundances)

# aggregate mean abundances by taxon
idx <- sort(unique(barley_rotus$ncbi_id), index.return=T)$ix
a <- aggregate(.~barley_rotus$ncbi_id, barley_rotus, FUN=sum)$abundance
barley_rotus <- data.frame(ncbi_id=unique(barley_rotus$ncbi_id)[idx],
                           ncbi_name=unique(barley_rotus$ncbi_name)[idx],
                           abundance=a)

# add all Rhizobiales abundances
barley_rhizobiales <- barley_rotus$abundance[barley_rotus$ncbi_name=="Hyphomicrobiaceae"] +
                      barley_rotus$abundance[barley_rotus$ncbi_name=="Rhizobiaceae"]
barley_rotus <- barley_rotus[!barley_rotus$ncbi_name %in% c("Hyphomicrobiaceae", "Rhizobiaceae"), ]
barley_rotus <- rbind(barley_rotus, data.frame(ncbi_id=356,
                                               ncbi_name="Rhizobiales",
                                               abundance=barley_rhizobiales))

# ARABIDOPSIS

# get OTU IDs
at_rotus <- as.character(at_rotus_names[, 1])

# normalize count matrix
at_abundance_matrix <- apply(at.count_matrix, 2, function(xx) xx / sum(xx))

# subset root samples and root OTUs
at_root_samples <- as.character(at.design$Sample_ID[at.design$Compartment=="Root"])
at_golm_samples <- as.character(at.design$Sample_ID[at.design$Site=="Golm"])
idx <- colnames(at_abundance_matrix) %in% intersect(at_root_samples, at_golm_samples)
at_rotus_abundance_matrix <- at_abundance_matrix[, idx]
at_rotus_abundance_matrix <- at_rotus_abundance_matrix[rownames(at_abundance_matrix) %in% at_rotus, ]

# perc. of rOTUs w.r.t. to the whole root community
at_rotus_perc <- mean(colSums(at_rotus_abundance_matrix))

# calculate mean OTU abundance
at_mean_rotus_abundances <- apply(at_rotus_abundance_matrix, 1, mean)
idx <- match(at_rotus_ids[, 1], names(at_mean_rotus_abundances))
at_mean_rotus_abundances <- at_mean_rotus_abundances[idx]

at_rotus <- data.frame(ncbi_name=as.character(at_rotus_names[, 2]),
                           ncbi_id=at_rotus_ids[, 2],
                           abundance=at_mean_rotus_abundances)

# add up all Burkholderiales abundances
burk_names <-  c("Burkholderiales", "Burkholderiaceae",
                 "uncultured Aquabacterium sp.", "Methylibium sp. DR10")
at_burkholderiales <- sum(at_rotus$abundance[at_rotus$ncbi_name %in% burk_names])
at_rotus <- at_rotus[!at_rotus$ncbi_name %in% c(burk_names, "uncultured Sorangiineae bacterium"), ]

# aggregate mean abundances by taxon
idx <- sort(unique(at_rotus$ncbi_id), index.return=T)$ix
at_rotus <- data.frame(ncbi_id=unique(at_rotus$ncbi_id)[idx],
                           ncbi_name=unique(at_rotus$ncbi_name)[idx],
                           abundance=aggregate(.~at_rotus$ncbi_id, at_rotus, FUN=sum)$abundance)

at_rotus <- rbind(at_rotus, data.frame(ncbi_id=80840,
                                       ncbi_name="Burkholderiales",
                                       abundance=at_burkholderiales))

# get all unique NCBI ids and names
all_ids <- unique(union(barley_rotus$ncbi_id, at_rotus$ncbi_id))
all_names <- unique(union(barley_rotus$ncbi_name, at_rotus$ncbi_name))

# tree plotting (Fig. 4)

# colors
c_blue <- rgb(0, 0, 1, .6)
c_red <- rgb(1, 0, 0, .6)
c_green <- rgb(0, 1, 0, .6)

# tree files
ids.tree.file <- "./data/NCBI_IDs.tree"
names.tree.file <- "./data/NCBI_names.tree"

# read newick trees 
tree_ids <- read.newick(file=ids.tree.file)
tree_names <- read.newick(file=names.tree.file)

# get node labels from the tree with numeric ids
node_labels <- data.frame(ncbi_id=c(tree_ids$tip.label, tree_ids$node.label))

# generate cex matrix 
barley_rotus <- barley_rotus[, -2]
at_rotus <- at_rotus[, -2]
cex <- merge(node_labels, at_rotus, by="ncbi_id", all.x=T)
cex <- merge(cex, barley_rotus, by="ncbi_id", all.x=T)
cex[is.na(cex)] <- 0
# cex <- cex[-1, ]
rownames(cex) <- cex$ncbi_id
cex <- cex[, 2:3]

# reorder according to the tree node labels
cex <- cex[match(node_labels$ncbi_id, rownames(cex)), ]

# log_transformed cex
barley_cex <- log2(cex[, 2] * barley_rotus_perc * 500 + 1)
at_cex <- log2(cex[, 1] * at_rotus_perc * 500 + 1)

plot_btree(floating_pies=F, tree=tree_names, cex_labels=2, type="cladogram",
           abundances=at_cex, abundances2=barley_cex,
           cex_abundances=.7, cex_abundances2=.7,
           colors=c_red, colors2=c_blue, show.node.label=T,
           offset_labels=9, output_file_name="./figures/Fig_4.pdf",
           points_pch=16, margins=c(1, 1, 1, 1), edge.width=1)

