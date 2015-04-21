
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

### comparison of 16S rRNA amplicon v. shotgun metagenome abundances

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
ampl.count.matrix.file <- "./data/count_matrix.txt"
ampl.ncbi.ids.file <- "./data/barley_ncbi_ids_order_or_above.txt"
ampl.ncbi.names.file <- "./data/barley_ncbi_names_order_or_above.txt"
metagenome.ncbi.ids.file <- "./data/metagenome_ncbi_ids_order_or_above.txt"
metagenome.ncbi.names.file <- "./data/metagenome_ncbi_names_order_or_above.txt"
metagenome.abundance.matrix.file <- "./data/abundances_coverage_per_sample_order.tsv"

# load data
barley.design <- read.table(file=barley.design.file, sep="\t", header=T)
ampl.count.matrix <- read.table(file=ampl.count.matrix.file, sep="\t")
ampl.ncbi.ids <- read.table(file=ampl.ncbi.ids.file, sep="\t")
ampl.ncbi.names <- read.table(file=ampl.ncbi.names.file, sep="\t")
metagenome.ncbi.ids <- read.table(file=metagenome.ncbi.ids.file, sep="\t")
metagenome.ncbi.names <- read.table(file=metagenome.ncbi.names.file, sep="\t")
metagenome.abundance.matrix <- read.table(file=metagenome.abundance.matrix.file, sep="\t")[-1, ]

# 16S ampl DATA

# subset root samples
barley.root.samples <- as.character(barley.design$SampleID[barley.design$Compartment=="Rhizosphere"])
ampl.count.matrix <- ampl.count.matrix[, colnames(ampl.count.matrix) %in% barley.root.samples]

# normalize count matrix
ampl.abundance.matrix <- apply(ampl.count.matrix, 2, function(x) x/sum(x))

# calculate mean OTU abundances
ampl.mean.abundances <- apply(ampl.abundance.matrix, 1, mean)

ampl.abundances <- data.frame(ncbi.name=as.character(ampl.ncbi.names[, 2]),
                                  ncbi.id=ampl.ncbi.ids[, 2],
                                  abundance=ampl.mean.abundances)

# aggregate abundances up to order level
idx <- sort(unique(ampl.ncbi.ids[, 2]), index.return=T)$ix
a <- aggregate(.~ampl.abundances$ncbi.id, ampl.abundances, FUN=sum)$abundance
ampl.abundances <- data.frame(ncbi.id=unique(ampl.abundances$ncbi.id)[idx],
                                  ncbi.name=unique(ampl.abundances$ncbi.name)[idx],
                                  abundance=a)


# METAGENOME DATA

# normalize
metagenome.abundance.matrix <- apply(metagenome.abundance.matrix, 2, function(x) x/sum(x))

# calculate mean taxon abundances
metagenome.mean.abundances <- apply(metagenome.abundance.matrix, 1, mean)

metagenome.abundances <- data.frame(ncbi.id=metagenome.ncbi.ids[, 1],
                                    ncbi.name=metagenome.ncbi.names[, 1],
                                    abundance=metagenome.mean.abundances)

# tree plotting (Fig. 5)

# colors
c_blue <- rgb(0, 0, 1, .6)
c_red <- rgb(1, 0, 0, .6)
c_green <- rgb(0, 1, 0, .6)
c_yellow <- rgb(1, 1, 0, .6)

# tree files
ids.tree.file <- "./data/NCBI_IDs_order.tree"
names.tree.file <- "./data/NCBI_names_order.tree"

# read newick trees 
tree_ids <- read.newick(file=ids.tree.file)
tree_names <- read.newick(file=names.tree.file)

# get node labels from the tree with numeric ids
node_labels <- data.frame(ncbi.id=c(tree_ids$tip.label, tree_ids$node.label))

# generate cex matrix

ampl.abundances <- ampl.abundances[, -2]
metagenome.abundances <- metagenome.abundances[, -2]
cex <- merge(node_labels, ampl.abundances, by="ncbi.id", all.x=T)
cex <- merge(cex, metagenome.abundances, by="ncbi.id", all.x=T)
cex[is.na(cex)] <- 0
rownames(cex) <- cex$ncbi.id
cex <- cex[, 2:3]

# rescale ampl cex (some taxon are missing from the reference tree)
cex[, 1] <- cex[, 1] / sum(cex[, 1])

# reorder according to the tree node labels
cex <- cex[match(node_labels$ncbi.id, rownames(cex)), ] * 10

# log_transformed cex
ampl.cex <- log2(cex[, 1] * 1e6 + 1)
metagenome.cex <- log2(cex[, 2] * 1e6 + 1)

ampl.colors <- rep(c_green, length(ampl.cex))
ampl.colors[ampl.cex > metagenome.cex] <- c_red

metagenome.colors <- rep(c_green, length(metagenome.cex))
metagenome.colors[metagenome.cex > ampl.cex] <- c_blue

plot_btree(floating_pies=F, tree=tree_names,
           abundances=ampl.cex, abundances2=metagenome.cex,
           cex_abundances=.08, cex_abundances2=.08,
           colors=ampl.colors, colors2=metagenome.colors,
           cex_labels=.4, tip.color="transparent", type="fan",
           offset_labels=6.5, output_file_name="./figures/Fig_5.pdf",
           points_pch=16, margins=c(7, 7, 7, 7), edge.width=1)

