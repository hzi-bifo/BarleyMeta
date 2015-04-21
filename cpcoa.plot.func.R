
library(vegan)
library(calibrate)
library(Biostrings)
library(ape)

variability_table <- function(cca){

        chi <- c(cca$tot.chi,
                       cca$CCA$tot.chi, cca$CA$tot.chi)
        variability_table <- cbind(chi, chi/chi[1])
        colnames(variability_table) <- c("inertia", "proportion")
        rownames(variability_table) <- c("total", "constrained", "unconstrained")
        return(variability_table)

}

cap_var_props <- function(cca){

        eig_tot <- sum(cca$CCA$eig)
        var_propdf <- cca$CCA$eig/eig_tot
        return(var_propdf)
}

pca_var_props <- function(cca){

        eig_tot <- sum(cca$CA$eig)
        var_propdf <- cca$CA$eig/eig_tot
        return(var_propdf)
}

cca_ci <- function(cca, permutations=5000){

        var_tbl <- variability_table(cca)
        p <- permutest(cca, permutations=permutations)
        ci <- quantile(p$F.perm, c(.05,.95))*p$chi[1]/var_tbl["total", "inertia"]
        return(ci)

}

plot_cap <- function(p, d, col_var, pch_var, col_comp, pch_comp,
                     shapes, colors, file_name, svg=F, constraint, 
                     ci, var_tbl, cap_var, perm_anova, cex.main=0.8, cex=0.8) {
	
	col <- rep("black", dim(p)[1])
	pch <- rep(19, dim(p)[1])
	for (i in 1:length(col_var)){
		index <- grep(col_var[i], d[rownames(p), col_comp[i]]) 
		col[index] <- colors[i]
	}
	for (i in 1:length(pch_var)){
		index <- grep(pch_var[i], d[rownames(p), pch_comp[i]]) 
		pch[index] <- shapes[i]
	}
		xlab <- paste("Constrained PCoA 1 (", format(cap_var[1]*100, digits=4), " %)", sep="")
		ylab <- paste("Constrained PCoA 1 (", format(cap_var[2]*100, digits=4), " %)", sep="")
		main <- paste(constraint, ": [", format(var_tbl["constrained", "proportion"]*100, digits=2),
				"% of variance; P < ", format(perm_anova[1,4], digits=2),
                "; 95% CI = ", format(ci[1]*100, digits=2), 
				"%, ", format(ci[2]*100, digits=2), "%]", sep="")	
		if(svg) svg(file=file_name) else pdf(file=file_name)
		plot(p, col=col, pch=pch, xlab=xlab, ylab=ylab, main=main, cex.main=cex.main, cex=cex)
		abline(v=0, h=0, lty="dotted")
		dev.off()
	
}

plot_pcoa <- function(p, d, variables, col_var, pch_var, col_comp, pch_comp,
                      shapes, colors, file_name, svg=F, pcoa_var, cex.main=0.8, cex=0.8) {
	
	col <- rep("black", dim(p)[1])
	pch <- rep(19, dim(p)[1])
	for (i in 1:length(col_var)){
		index <- grep(col_var[i], d[rownames(p), col_comp[i]]) 
		col[index] <- colors[i]
	}
	for (i in 1:length(pch_var)){
		index <- grep(pch_var[i], d[rownames(p), pch_comp[i]]) 
		pch[index] <- shapes[i]
	}
		xlab <- paste("PCo 1 (", format(pcoa_var[1]*100, digits=4), " %)", sep="")
		ylab <- paste("PCo 1 (", format(pcoa_var[2]*100, digits=4), " %)", sep="")
		main <- "Unconstrained PCoA"	
		if(svg) svg(file=file_name) else pdf(file=file_name)
		plot(p, col=col, pch=pch, xlab=xlab, ylab=ylab, main=main, cex.main=cex.main, cex=cex)
		abline(v=0, h=0, lty="dotted")	
		dev.off()
	
}

plot_pcoa_spp <- function(p, centroids=NULL, pch, file_name, svg=F, pcoa_var,
                          abundances=NULL, otu_subset=NULL, taxonomy=NULL, cex.main=0.8, cex=0.8){

	col <- rep("black", dim(p)[1])

	if(!is.null(taxonomy)){
		col[which(taxonomy[rownames(p),]$Phylum=="Proteobacteria")] <- "green4"
		col[which(taxonomy[rownames(p),]$Phylum=="Bacteroidetes")] <- "blue1"
		col[which(taxonomy[rownames(p),]$Phylum=="Actinobacteria")] <- "red1"
	}

	if(!is.null(abundances)){
		cex <- colMeans(abundances[, rownames(p)] * 0.5)
	}


		xlab <- paste("PCo 1 (", format(pcoa_var[1]*100, digits=4), " %)", sep="")
		ylab <- paste("PCo 1 (", format(pcoa_var[2]*100, digits=4), " %)", sep="")
		main <- "Unconstrained PCoA"	
		if(svg) svg(file=file_name) else pdf(file=file_name)



		plot(p, col=col, pch=pch, xlab=xlab, ylab=ylab, main=main, cex.main=cex.main, cex=cex)
		abline(v=0, h=0, lty="dotted")	

		if(!is.null(otu_subset)){
			idx <- which(rownames(p) %in% otu_subset)
			points(p[idx, 1:2], pch=16, col=col[idx], cex=cex[idx])
		}


		dev.off()
	
}

plot_biplots <- function(p, centroids, pch, col, file_name, svg=F, constraint, ci,
                         var_tbl, cap_var, perm_anova, abundances=NULL, otu_subset=NULL,
                         taxonomy=NULL, cex.main=0.8, cex=10) {
	
	colors <- rep("black", dim(p)[1]) 

	xlab <- paste("Constrained PCoA 1 (", format(cap_var[1]*100, digits=4), " %)", sep="")
	ylab <- paste("Constrained PCoA 1 (", format(cap_var[2]*100, digits=4), " %)", sep="")
	main <- paste(constraint, ": [", format(var_tbl["constrained", "proportion"]*100, digits=2),
			"% of variance; P < ", format(perm_anova[1,4], digits=2), "; 95% CI = ", format(ci[1]*100, digits=2), 
			"%, ", format(ci[2]*100, digits=2), "%]", sep="")	
	if(svg) svg(file=file_name) else pdf(file=file_name)

	if(!is.null(taxonomy)){
		colors[which(taxonomy[rownames(p),]$Phylum=="Proteobacteria")] <- "green4"
		colors[which(taxonomy[rownames(p),]$Phylum=="Bacteroidetes")] <- "blue1"
		colors[which(taxonomy[rownames(p),]$Phylum=="Actinobacteria")] <- "red1"
	}

	if(!is.null(abundances)){
		cex <- colMeans(abundances[, rownames(p)] * 0.5)
	}

	plot(rbind(centroids, p), col="white", pch=19, xlab=xlab, ylab=ylab, main=main, cex.main=cex.main)
	points(p, col=colors, cex=cex, pch=pch)
	arrows(x0=c(0,0,0,0), y0=c(0,0,0,0), x1=centroids[,1], y1=centroids[,2], col=col, pch=17, length=0.1, cex=1.5, lwd=2)
	abline(v=0, h=0, lty="dotted")

	if(!is.null(otu_subset)){
		idx <- which(rownames(p) %in% otu_subset)
		points(p[idx, 1:2], pch=16, col=colors[idx], cex=cex[idx])
	} 
	dev.off()

}

log_transform <- function(matrix, scale=1000, threshold=2) {

        transformed <- matrix / rowSums(matrix, na=T) * scale
        transformed <- transformed[,apply(transformed, 2, max) > threshold]
        transformed <- log2(transformed + 1)

        return(transformed)

}

wisconsin_transform <- function(m) {

	m_transform <- apply(m, 2, function(xx) xx / sum(xx))
	m_transform <- t(apply(m_transform, 1, function(xx) xx / sum(xx))) * 1000
	
	return(m_transform)

}

cpcoa_analysis <- function(count_matrix, d, taxonomy, label){

    # apply transformations to count matrix
    
    matrix_log_transformed <- log_transform(count_matrix, threshold=2)
    count_matrix <- count_matrix[, match(colnames(matrix_log_transformed),
                                         colnames(count_matrix))]
    matrix_w_transformed <- wisconsin_transform(count_matrix)
    
    # generate subset without soil samples for genotype-constrained PCoA
    
    #t_nosoil <- matrix_w_transformed[which(!grepl("Bulk", rownames(count_matrix))), ]
    t_nosoil <- matrix_log_transformed[which(!grepl("Bulk", rownames(count_matrix))), ]
    d_nosoil <- d[which(!grepl("Soil", d$Compartment)), ]
    
    # cpcoa by compartment, experiment + compartment and genotype using bray-curtis distance
    
    sqrt_transform <- F
    
    capscale.gen_bc_nosoil <- capscale(t_nosoil ~ Description + Condition(Experiment + Compartment),
                                       data=d_nosoil, add=F, sqrt.dist=sqrt_transform, dist="bray")
    
    # anova-like permutation analysis -anova.cca is a vegan
    # wrapper for cca which uses the function permutest-
    
    perm_anova.gen_bc_nosoil <- anova.cca(capscale.gen_bc_nosoil,
                                          strata=c("Experiment", "Compartment"))
    print(perm_anova.gen_bc_nosoil)
    
    # generate variability tables and calculate confidence intervals for the variance
    
    var_tbl.gen_bc_nosoil <- variability_table(capscale.gen_bc_nosoil)
    cap_var.gen_bc_nosoil <- cap_var_props(capscale.gen_bc_nosoil)
    ci.gen_bc_nosoil <-  cca_ci(capscale.gen_bc_nosoil) 
    
    # extract the weighted average (sample) scores
    
    wa.gen_bc_nosoil <- capscale.gen_bc_nosoil$CCA$wa
    
    # extract OTU scores from CCA object
    
    otu.scores.gen_bc_nosoil <- capscale.gen_bc_nosoil$CCA$v[, 1:2]
    
    # extract centroids of constrained factor
    
    centroids.gen_bc_nosoil <- capscale.gen_bc_nosoil$CCA$centroids[, 1:2]
    
    # plot CPCoA
    
    plot_cap(p=wa.gen_bc_nosoil, d=d, col_var=c("spontaneum", "landrace", "modern"),
             pch_var=c("Root", "Rhizosphere"),
             col_comp=c("Description", "Description", "Description"),
             pch_comp=c("Compartment", "Compartment"), shapes=c(16, 17),
             colors=c("darkgreen", "red", "blue"),
             file_name=paste("./figures/gen_bc_nosoil_", label, ".pdf", sep=""),
             svg=F, constraint=label,
             ci=ci.gen_bc_nosoil, var_tbl=var_tbl.gen_bc_nosoil,
             cap_var=cap_var.gen_bc_nosoil, perm_anova=perm_anova.gen_bc_nosoil,
             cex.main=0.8, cex=1.3)
    
    # plot OTU scores and sample biplots
    
    plot_biplots(p=otu.scores.gen_bc_nosoil,
                 centroids=centroids.gen_bc_nosoil,
                 pch=16, col=c("red", "blue", "darkgreen"),
                 file_name=paste("./figures/gen_bc_nosoil_otus_log_", label, ".pdf", sep=""), svg=F,
                 constraint=label,
                 ci=ci.gen_bc_nosoil, var_tbl=var_tbl.gen_bc_nosoil,
                 cap_var=cap_var.gen_bc_nosoil, perm_anova=perm_anova.gen_bc_nosoil,
                 abundances=matrix_log_transformed, otu_subset=NULL,
                 taxonomy=taxonomy, cex.main=0.8, cex=0.7)
    
    plot_biplots(p=otu.scores.gen_bc_nosoil,
                 centroids=centroids.gen_bc_nosoil,
                 pch=16, col=c("red", "blue", "darkgreen"),
                 file_name=paste("./figures/gen_bc_nosoil_otus_", label, ".pdf", sep=""), svg=F,
                 constraint=label,
                 ci=ci.gen_bc_nosoil, var_tbl=var_tbl.gen_bc_nosoil,
                 cap_var=cap_var.gen_bc_nosoil, perm_anova=perm_anova.gen_bc_nosoil,
                 abundances=count_matrix * 0.08, otu_subset=NULL,
                 taxonomy=taxonomy, cex.main=0.8, cex=0.7)

    return(paste(label, ": [", format(var_tbl.gen_bc_nosoil["constrained", "proportion"]*100, digits=5),
                 "% of variance; P < ", format(perm_anova.gen_bc_nosoil[1,5], digits=5), "; 95% CI = ",
                 format(ci.gen_bc_nosoil[1]*100, digits=5),
                 "%, ", format(ci.gen_bc_nosoil[2]*100, digits=5), "%]", sep=""))

}

