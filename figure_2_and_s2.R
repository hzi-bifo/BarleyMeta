
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

### ananlysis of beta-diversity (CPCoA)

# cleanup
rm(list=ls())

# load functions
source("cpcoa.plot.func.R")

# input files

# QIIME mapping
design.file <- "./data/mapping.txt"
# OTU taxonomy
taxonomy.file <- "./data/taxonomy.txt"
# OTU count table
count.matrix.file <- "./data/count_matrix.txt"
# UniFrac distance matrices
wunifrac.file <- "./data/weighted_unifrac.txt"
uwunifrac.file <- "./data/unweighted_unifrac.txt"


# load the design table, the taxonomy information,
# the count table and distance matrices from QIIME

d <- read.table(file=design.file, sep="\t", header=T)
d <- data.frame(Compartment=d$Compartment,
                Experiment=d$Experiment,
                Description=d$Description,
                row.names=d$SampleID)
taxonomy <- read.table(file=taxonomy.file, fill=T, header=T, sep="\t")
count_matrix <- t(read.table(file=count.matrix.file, sep="\t"))
wunifrac <- read.table(wunifrac.file)
uwunifrac <- read.table(uwunifrac.file)

# apply transformations to OTU count table

matrix_log_transformed <- log_transform(count_matrix, threshold=2)
count_matrix <- count_matrix[, match(colnames(matrix_log_transformed), colnames(count_matrix))]
matrix_w_transformed <- wisconsin_transform(count_matrix)

# generate subset without soil samples for genotype-constrained PCoA

t_nosoil <- matrix_log_transformed[which(!grepl("Bulk", rownames(count_matrix))), ]
d_nosoil <- d[which(!grepl("Soil", d$Compartment)), ]
wunifrac_nosoil <- wunifrac[which(!grepl("Bulk", rownames(count_matrix))), ]
uwunifrac_nosoil <- uwunifrac[which(!grepl("Bulk", rownames(count_matrix))), ]

###
### Figs. 2A, 2B, S2A
###

# CPCoA + compartment and genotype using bray-curtis distance

sqrt_transform <- T

# bray-curtis

capscale.gen_bc_nosoil <- capscale(t_nosoil ~ Description + Condition(Experiment + Compartment),
                                   data=d_nosoil, add=F, sqrt.dist=sqrt_transform, dist="bray")

capscale.comp <- capscale(matrix_log_transformed ~ Compartment + Condition(Experiment),
                          data=d, add=F, sqrt.dist=sqrt_transform, dist="bray")

# ANOVA-like permutation analysis
# anova.cca is a vegan wrapper for CCA, which uses the function permutest

perm_anova.gen_bc_nosoil <- anova.cca(capscale.gen_bc_nosoil)
print(perm_anova.gen_bc_nosoil)

perm_anova.comp <- anova.cca(capscale.comp)
print(perm_anova.comp)

# generate variability tables and calculate confidence intervals for the variance

var_tbl.gen_bc_nosoil <- variability_table(capscale.gen_bc_nosoil)
cap_var.gen_bc_nosoil <- cap_var_props(capscale.gen_bc_nosoil)
ci.gen_bc_nosoil <-  cca_ci(capscale.gen_bc_nosoil) 

var_tbl.comp <- variability_table(capscale.comp)
cap_var.comp <- cap_var_props(capscale.comp)
ci.comp <-  cca_ci(capscale.comp) 

# extract the weighted average (sample) scores

wa.gen_bc_nosoil <- capscale.gen_bc_nosoil$CCA$wa
wa.comp <- capscale.comp$CCA$wa

# extract OTU scores from CCA object

otu.scores.gen_bc_nosoil <- capscale.gen_bc_nosoil$CCA$v[, 1:2]
otu.scores.comp <- capscale.comp$CCA$v[, 1:2]

# extract centroids of constrained factor

centroids.gen_bc_nosoil <- capscale.gen_bc_nosoil$CCA$centroids[, 1:2]
centroids.comp <- capscale.comp$CCA$centroids[, 1:2]

# plot CPCoA (sample scores)

plot_cap(p=wa.gen_bc_nosoil, d=d, col_var=c("spontaneum", "landrace", "modern"),
         pch_var=c("Root", "Rhizosphere"),
         col_comp=c("Description", "Description", "Description"),
         pch_comp=c("Compartment", "Compartment"), shapes=c(16, 17),
         colors=c("darkgreen", "red", "blue"), file_name="./figures/Fig_2B.pdf",
         svg=F, constraint="Bray-Curtis - constrained by Genotype",
         ci=ci.gen_bc_nosoil, var_tbl=var_tbl.gen_bc_nosoil,
         cap_var=cap_var.gen_bc_nosoil, perm_anova=perm_anova.gen_bc_nosoil,
         cex.main=0.8, cex=1.3)

plot_cap(p=wa.comp, d=d, col_var=c("Soil", "Root", "Rhizosphere"),
         pch_var=c("Soil", "Root", "Rhizosphere"),
         col_comp=c("Compartment", "Compartment", "Compartment"),
         pch_comp=c("Compartment", "Compartment", "Compartment"), shapes=c(0, 1, 2),
         colors=c("black", "black", "black"), file_name="./figures/Fig_2A.pdf",
         svg=F, constraint="Bray-Curtis - constrained by Compartment",
         ci=ci.comp, var_tbl=var_tbl.comp,
         cap_var=cap_var.comp, perm_anova=perm_anova.comp,
         cex.main=0.8, cex=1.3)


# plot OTU scores and sample biplots

plot_biplots(p=otu.scores.gen_bc_nosoil,
             centroids=centroids.gen_bc_nosoil,
             pch=16, col=c("red", "blue", "darkgreen"),
             file_name="./figures/Fig_S2A_right.pdf", svg=F,
             constraint="Bray-Curtis - cosntrained by Genotype (log scale)",
             ci=ci.gen_bc_nosoil, var_tbl=var_tbl.gen_bc_nosoil,
             cap_var=cap_var.gen_bc_nosoil, perm_anova=perm_anova.gen_bc_nosoil,
             abundances=matrix_log_transformed, otu_subset=NULL,
             taxonomy=taxonomy, cex.main=0.8, cex=0.7)

plot_biplots(p=otu.scores.comp,
             centroids=centroids.comp,
             pch=16, col=c("red", "blue", "darkgreen"),
             file_name="./figures/FigS_2A_left.pdf", svg=F,
             constraint="Bray-Curtis - cosntrained by Genotype (log scale)",
             ci=ci.comp, var_tbl=var_tbl.comp,
             cap_var=cap_var.comp, perm_anova=perm_anova.comp,
             abundances=matrix_log_transformed, otu_subset=NULL,
             taxonomy=taxonomy, cex.main=0.8, cex=0.7)

###
### Figs. S2B
###

# CPCoA + compartment and genotype using bray-curtis distance

sqrt_transform <- T

# unweighted unifrac

capscale.gen_bc_nosoil <- capscale(uwunifrac_nosoil ~ Description + Condition(Experiment + Compartment),
                                     data=d_nosoil, add=T, sqrt.dist=sqrt_transform)

capscale.comp <- capscale(uwunifrac ~ Compartment + Condition(Experiment), 
                          data=d, add=T, sqrt.dist=sqrt_transform)

# ANOVA-like permutation analysis
# anova.cca is a vegan wrapper for CCA, which uses the function permutest

perm_anova.gen_bc_nosoil <- anova.cca(capscale.gen_bc_nosoil)
print(perm_anova.gen_bc_nosoil)

perm_anova.comp <- anova.cca(capscale.comp)
print(perm_anova.comp)

# generate variability tables and calculate confidence intervals for the variance

var_tbl.gen_bc_nosoil <- variability_table(capscale.gen_bc_nosoil)
cap_var.gen_bc_nosoil <- cap_var_props(capscale.gen_bc_nosoil)
ci.gen_bc_nosoil <-  cca_ci(capscale.gen_bc_nosoil) 

var_tbl.comp <- variability_table(capscale.comp)
cap_var.comp <- cap_var_props(capscale.comp)
ci.comp <-  cca_ci(capscale.comp) 

# extract the weighted average (sample) scores

wa.gen_bc_nosoil <- capscale.gen_bc_nosoil$CCA$wa
wa.comp <- capscale.comp$CCA$wa

# extract OTU scores from CCA object

otu.scores.gen_bc_nosoil <- capscale.gen_bc_nosoil$CCA$v[, 1:2]
otu.scores.comp <- capscale.comp$CCA$v[, 1:2]

# extract centroids of constrained factor

centroids.gen_bc_nosoil <- capscale.gen_bc_nosoil$CCA$centroids[, 1:2]
centroids.comp <- capscale.comp$CCA$centroids[, 1:2]

# plot CPCoA (sample scores)

plot_cap(p=wa.gen_bc_nosoil, d=d, col_var=c("spontaneum", "landrace", "modern"),
         pch_var=c("Root", "Rhizosphere"),
         col_comp=c("Description", "Description", "Description"),
         pch_comp=c("Compartment", "Compartment"), shapes=c(16, 17),
         colors=c("darkgreen", "red", "blue"), file_name="./figures/Fig_S2B_right.pdf",
         svg=F, constraint="Bray-Curtis - constrained by Genotype",
         ci=ci.gen_bc_nosoil, var_tbl=var_tbl.gen_bc_nosoil,
         cap_var=cap_var.gen_bc_nosoil, perm_anova=perm_anova.gen_bc_nosoil,
         cex.main=0.8, cex=1.3)

plot_cap(p=wa.comp, d=d, col_var=c("Soil", "Root", "Rhizosphere"),
         pch_var=c("Soil", "Root", "Rhizosphere"),
         col_comp=c("Compartment", "Compartment", "Compartment"),
         pch_comp=c("Compartment", "Compartment", "Compartment"), shapes=c(0, 1, 2),
         colors=c("black", "black", "black"), file_name="./figures/Fig_S2B_left.pdf",
         svg=F, constraint="Bray-Curtis - constrained by Compartment",
         ci=ci.comp, var_tbl=var_tbl.comp,
         cap_var=cap_var.comp, perm_anova=perm_anova.comp,
         cex.main=0.8, cex=1.3)

###
### Figs. S2C
###

# CPCoA + compartment and genotype using bray-curtis distance

sqrt_transform <- T

# unweighted unifrac

# weighted unifrac

capscale.gen_bc_nosoil <- capscale(wunifrac_nosoil ~ Description + Condition(Experiment + Compartment),
                                   data=d_nosoil, add=T, sqrt.dist=sqrt_transform)

capscale.comp <- capscale(wunifrac ~ Compartment + Condition(Experiment), data=d,
                          add=T, sqrt.dist=sqrt_transform)

# ANOVA-like permutation analysis
# anova.cca is a vegan wrapper for CCA, which uses the function permutest

perm_anova.gen_bc_nosoil <- anova.cca(capscale.gen_bc_nosoil)
print(perm_anova.gen_bc_nosoil)

perm_anova.comp <- anova.cca(capscale.comp)
print(perm_anova.comp)

# generate variability tables and calculate confidence intervals for the variance

var_tbl.gen_bc_nosoil <- variability_table(capscale.gen_bc_nosoil)
cap_var.gen_bc_nosoil <- cap_var_props(capscale.gen_bc_nosoil)
ci.gen_bc_nosoil <-  cca_ci(capscale.gen_bc_nosoil) 

var_tbl.comp <- variability_table(capscale.comp)
cap_var.comp <- cap_var_props(capscale.comp)
ci.comp <-  cca_ci(capscale.comp) 

# extract the weighted average (sample) scores

wa.gen_bc_nosoil <- capscale.gen_bc_nosoil$CCA$wa
wa.comp <- capscale.comp$CCA$wa

# extract OTU scores from CCA object

otu.scores.gen_bc_nosoil <- capscale.gen_bc_nosoil$CCA$v[, 1:2]
otu.scores.comp <- capscale.comp$CCA$v[, 1:2]

# extract centroids of constrained factor

centroids.gen_bc_nosoil <- capscale.gen_bc_nosoil$CCA$centroids[, 1:2]
centroids.comp <- capscale.comp$CCA$centroids[, 1:2]

# plot CPCoA (sample scores)

plot_cap(p=wa.gen_bc_nosoil, d=d, col_var=c("spontaneum", "landrace", "modern"),
         pch_var=c("Root", "Rhizosphere"),
         col_comp=c("Description", "Description", "Description"),
         pch_comp=c("Compartment", "Compartment"), shapes=c(16, 17),
         colors=c("darkgreen", "red", "blue"), file_name="./figures/Fig_S2C_right.pdf",
         svg=F, constraint="Bray-Curtis - constrained by Genotype",
         ci=ci.gen_bc_nosoil, var_tbl=var_tbl.gen_bc_nosoil,
         cap_var=cap_var.gen_bc_nosoil, perm_anova=perm_anova.gen_bc_nosoil,
         cex.main=0.8, cex=1.3)

plot_cap(p=wa.comp, d=d, col_var=c("Soil", "Root", "Rhizosphere"),
         pch_var=c("Soil", "Root", "Rhizosphere"),
         col_comp=c("Compartment", "Compartment", "Compartment"),
         pch_comp=c("Compartment", "Compartment", "Compartment"), shapes=c(0, 1, 2),
         colors=c("black", "black", "black"), file_name="./figures/Fig_S2C_left.pdf",
         svg=F, constraint="Bray-Curtis - constrained by Compartment",
         ci=ci.comp, var_tbl=var_tbl.comp,
         cap_var=cap_var.comp, perm_anova=perm_anova.comp,
         cex.main=0.8, cex=1.3)

