
# given a clustering function and a distance matrix
# gives you a tree with bootstra support values
tree_bt <- function(func, data, main, type, direction="downwards", ...) {
    tree <- func((t(data)))
    tree$node.label <- boot.phylo(tree, t(data), FUN=func, B=100)
    plot.phylo(tree, show.node.label=T, type=type, direction=direction, ...)
    title(main=main)
}

# function to read a reference tree from a newick file
# and collapse internal nodes with only one child
read.ref.tree <- function(file) {
    tree <- read.newick(file=file)
    tree <- collapse.singles(tree)
    tree$edge.length[] <- 1
    return(tree)
}

# function to plot pie charts on the nodes of the tree
floating.pie <- function (xpos, ypos, x, edges = 200, radius = .1, 
                          col = NULL, cex.p=0.1, startpos = 0, ...) {
    u <- par("usr")
    user.asp <- diff(u[3:4]) / diff(u[1:2])
    p <- par("pin")
    inches.asp <- p[2] / p[1]
    asp <- user.asp/inches.asp
    if (!is.numeric(x) || any(is.na(x) | x < 0)) {
    	stop("floating.pie: x values must be non-negative")
    }
    x <- c(0, cumsum(x) / sum(x))
    dx <- diff(x)
    nx <- length(dx)
    if (is.null(col)) 
        col <- rainbow(nx)
    else if (length(col) < nx) 
        col <- rep(col, nx)
    bc <- 2 * pi * (x[1:nx] + dx/2) + startpos
    for (i in 1:nx) {
        n <- max(2, floor(edges * dx[i]))
        t2p <- 2 * pi * seq(x[i], x[i + 1], length = n) + startpos
        xc <- c(cos(t2p) * radius + xpos, xpos)
        yc <- c(sin(t2p) * radius * asp + ypos, ypos)
        polygon(xc, yc, col = col[i], cex=0.01, border="transparent", ...)
    }
}

# reordering of the tree (for plotting)
reorder.phylo <- function (x, order = "cladewise", index.only = FALSE, ...){
    ORDER <- c("cladewise", "postorder", "pruningwise")
    io <- pmatch(order, ORDER)
    if (is.na(io)) 
        stop("ambiguous order")
    order <- ORDER[io]
    nb.edge <- dim(x$edge)[1]
    if (!is.null(attr(x, "order"))) 
        if (attr(x, "order") == order) 
            if (index.only) 
                return(1:nb.edge)
            else return(x)
    nb.node <- x$Nnode
    if (nb.node == 1) 
        return(x)
    nb.tip <- length(x$tip.label)
    if (io == 3) {
        x <- reorder(x)
        neworder <- .C("neworder_pruningwise", as.integer(nb.tip), as.integer(nb.node),
                       as.integer(x$edge[, 1]), as.integer(x$edge[,2]),
                       as.integer(nb.edge), integer(nb.edge))[[6]]
    }
    else {
        neworder <- .C("neworder_phylo", as.integer(nb.tip), as.integer(x$edge[, 1]),
                       as.integer(x$edge[, 2]), as.integer(nb.edge), integer(nb.edge),
                       io, DUP = FALSE, NAOK = TRUE)[[5]]
    }
    if (index.only) 
        return(neworder)
    x$edge <- x$edge[neworder, ]
    if (!is.null(x$edge.length)) 
        x$edge.length <- x$edge.length[neworder]
    if (sum(is.na(x$edge.length)))
      x$edge.length[] <- 1
    attr(x, "order") <- order
    x
}

# function to plot the tree with node charts
# uses plot.phylo to obtain the node coordinates and then redraws the tree
# cex values need to be adjusted
plot_btree <- function(tree, abundances, abundances2=NULL, pie_sizes, colors, colors2,
                       output_file_name=NULL, svg=F,
                       tip.color="black", show.node.label=F,
                       cex_abundances=1, cex_abundances2=1, points_pch=16, cex_labels=1,
                       offset_labels=5.3, floating_pies=F, 
                       colors_pie=c("blue", "green", "red"), type="fan",
                       edge=T, color_labels="black", color_alg_lines="dark grey",
                       margins=c(7, 7, 7, 7), edge.color="black", ...) {

    if(svg) {
        svg(file=output_file_name) 
    } else {
        pdf(file=output_file_name)
    }
    
    par(mar=margins)
    plot.phylo2(tree, type=type, tip.color=tip.color,
                use.edge.length=T, label.offset=0,  show.node.label=show.node.label,
                cex=.5, edge.color=edge.color, ...)
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    XX <- lastPP$xx
    YY <- lastPP$yy
    xx.tips <- lastPP$xx[1:lastPP$Ntip]
    yy.tips <- lastPP$yy[1:lastPP$Ntip]
    angle <- atan2(yy.tips, xx.tips)
    xx.labels <- offset_labels * cos(angle)
    yy.labels <- offset_labels * sin(angle)
    if(type == "fan") {
        par(xpd=NA)
        sapply(1:length(xx.labels), function(i) lines(c(xx.tips[i], xx.labels[i]),
                                                      c(yy.tips[i], yy.labels[i]),
                                                      lty=3, col=color_alg_lines))
        s <- xx.labels < 0
        angle <- angle * 180 / pi
        angle[s] <- angle[s] + 180
        adj <- rep(0, length(tree$tip.label))
        adj[which(xx.labels<0)] <- 1 
        sapply(1:length(tree$tip.label), 
               function(i) text(x=xx.labels[i], y=yy.labels[i],
                                labels=tree$tip.label[i], srt=angle[i],
                                offset=1, cex=cex_labels, adj=adj[i],
                                col=color_labels))
    }
    if(floating_pies) {
    	sapply(1:length(pie_sizes), 
                function(i) {
                    if(rowSums(abundances)[i]!=0)
                        floating.pie(col=colors_pie, xpos=XX[i],
                                     ypos=YY[i], x=abundances[i,], 
                                     radius=pie_sizes[i] * cex_abundances)
                            })
    } else {
    	points(x=XX, y=YY, col=colors, cex=abundances * cex_abundances, pch=points_pch)
        points(x=XX, y=YY, col=colors2, cex=abundances2 * cex_abundances2, pch=points_pch)
    }
    dev.off()
}

# function to read a Newick string with node labels & (possible) singles
# written by Liam J. Revell 2013
read.newick<-function(file="",text){
	# check to see if reading from file
	if(file!="") text<-scan(file,sep="\n",what="character")
	if(length(text)>1){
		tree<-lapply(text,newick)
		class(tree)<-"multiPhylo"
	} else tree<-newick(text)
	return(tree)
}

# main Newick string function
# written by Liam J. Revell 2013
newick<-function(text){
	text<-unlist(strsplit(text, NULL))
	tip.label<-vector(mode="character")
	node.label<-vector(mode="character") 
	edge<-matrix(c(1,NA),1,2) 
	edge.length<-vector()
	currnode<-1
	Nnode<-currnode
	i<-j<-k<-1
	while(text[i]!=";"){
		if(text[i]=="("){
			if(j>nrow(edge)) edge<-rbind(edge,c(NA,NA))
			edge[j,1]<-currnode
			i<-i+1
			# is the next element a label?
			if(is.na(match(text[i],c("(",")",",",":",";")))){
				temp<-getLabel(text,i)
				tip.label[k]<-temp$label
				i<-temp$end
				edge[j,2]<--k
				k<-k+1
				# is there a branch length?
				if(text[i]==":"){
					temp<-getEdgeLength(text,i)
					edge.length[j]<-temp$edge.length
					i<-temp$end
				}	
			} else if(text[i]=="("){
				Nnode<-Nnode+1 # creating a new internal node
				currnode<-Nnode
				edge[j,2]<-currnode # move to new internal node
			}
			j<-j+1
		} else if(text[i]==")"){
			i<-i+1
			# is the next element a label?
			if(is.na(match(text[i],c("(",")",",",":",";")))){
				temp<-getLabel(text,i)
				node.label[currnode]<-temp$label
				i<-temp$end
			}
			# is there a branch length?
			if(text[i]==":"){
				temp<-getEdgeLength(text,i)
				if(currnode>1){ 
					ii<-match(currnode,edge[,2])
					edge.length[ii]<-temp$edge.length
				} else root.edge<-temp$edge.length
				i<-temp$end
			}	
			if(currnode>1) currnode<-edge[match(currnode,edge[,2]),1] # move down the tree
		} else if(text[i]==","){
			if(j>nrow(edge)) edge<-rbind(edge,c(NA,NA))
			edge[j,1]<-currnode
			i<-i+1
			# is the next element a label?
			if(is.na(match(text[i],c("(",")",",",":",";")))){
				temp<-getLabel(text,i)
				tip.label[k]<-temp$label
				i<-temp$end
				edge[j,2]<--k
				k<-k+1
				# is there a branch length?
				if(text[i]==":"){
					temp<-getEdgeLength(text,i)
					edge.length[j]<-temp$edge.length
					i<-temp$end
				}
			} else if(text[i]=="("){
				Nnode<-Nnode+1 # creating a new internal node
				currnode<-Nnode
				edge[j,2]<-currnode # move to internal node
			}
			j<-j+1
		}
	}
	Ntip<-k-1
	edge[edge>0]<-edge[edge>0]+Ntip
	edge[edge<0]<--edge[edge<0]
	edge.length[is.na(edge.length)]<-0
	node.label[is.na(node.label)]<-""
	if(length(node.label)==0) node.label<-NULL
	# assemble into "phylo" object
	tree<-list(edge=edge,Nnode=as.integer(Nnode),tip.label=tip.label,edge.length=edge.length,node.label=node.label)
	class(tree)<-"phylo"
	return(tree)
}

# function gets label
# written by Liam J. Revell 2011-2013
getLabel<-function(text,start,stop.char=c(",",":",")",";")){
	i<-0
	label<-vector()
	while(is.na(match(text[i+start],stop.char))){
		label[i+1]<-text[i+start]
		i<-i+1
	}
	return(list(label=paste(label,collapse=""),end=i+start))
}

# function gets branch length
# written by Liam J. Revell 2011-2013
getEdgeLength<-function(text,start){
	i<-start+1; m<-1
	temp<-vector()
	while(is.na(match(text[i],c(",",")",";")))){
		temp[m]<-text[i]
		i<-i+1
		m<-m+1
	}
	return(list(edge.length=as.numeric(paste(temp,collapse="")),end=i))
}

# adapted verson of plot.phylo from ape package
plot.phylo2 <- function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
    show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
    edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
    adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
    label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
    direction = "rightwards", lab4ut = "horizontal", tip.color = "black", 
    plot = TRUE, rotate.tree = 0, open.angle = 0, ...)  {
    Ntip <- length(x$tip.label)
    if (Ntip < 2) {
        warning("found less than 2 tips in the tree")
        return(NULL)
    }
#    if (any(tabulate(x$edge[, 1]) == 1)) 
#        stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles()")
    .nodeHeight <- function(Ntip, Nnode, edge, Nedge, yy) .C("node_height", 
        as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
            1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(yy), 
        DUP = FALSE, PACKAGE = "ape")[[6]]
    .nodeDepth <- function(Ntip, Nnode, edge, Nedge) .C("node_depth", 
        as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
            1]), as.integer(edge[, 2]), as.integer(Nedge), double(Ntip + 
            Nnode), DUP = FALSE, PACKAGE = "ape")[[6]]
    .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
        edge.length) .C("node_depth_edgelength", as.integer(Ntip), 
        as.integer(Nnode), as.integer(edge[, 1]), as.integer(edge[, 
            2]), as.integer(Nedge), as.double(edge.length), double(Ntip + 
            Nnode), DUP = FALSE, PACKAGE = "ape")[[7]]
    Nedge <- dim(x$edge)[1]
    Nnode <- x$Nnode
    ROOT <- Ntip + 1
    type <- match.arg(type, c("phylogram", "cladogram", "fan", 
        "unrooted", "radial"))
    direction <- match.arg(direction, c("rightwards", "leftwards", 
        "upwards", "downwards"))
    if (is.null(x$edge.length)) 
        use.edge.length <- FALSE
    if (type %in% c("unrooted", "radial") || !use.edge.length || 
        is.null(x$root.edge) || !x$root.edge) 
        root.edge <- FALSE
    if (type == "fan" && root.edge) {
        warning("drawing root edge with type = 'fan' is not yet supported")
        root.edge <- FALSE
    }
    phyloORclado <- type %in% c("phylogram", "cladogram")
    horizontal <- direction %in% c("rightwards", "leftwards")
    xe <- x$edge
    if (phyloORclado) {
        phyOrder <- attr(x, "order")
        if (is.null(phyOrder) || phyOrder != "cladewise") {
            x <- reorder(x)
            if (!identical(x$edge, xe)) {
                ereorder <- match(x$edge[, 2], xe[, 2])
                if (length(edge.color) > 1) {
                  edge.color <- rep(edge.color, length.out = Nedge)
                  edge.color <- edge.color[ereorder]
                }
                if (length(edge.width) > 1) {
                  edge.width <- rep(edge.width, length.out = Nedge)
                  edge.width <- edge.width[ereorder]
                }
                if (length(edge.lty) > 1) {
                  edge.lty <- rep(edge.lty, length.out = Nedge)
                  edge.lty <- edge.lty[ereorder]
                }
            }
        }
        yy <- numeric(Ntip + Nnode)
        TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
        yy[TIPS] <- 1:Ntip
    }
    z <- reorder(x, order = "pruningwise")
    if (phyloORclado) {
        if (is.null(node.pos)) {
            node.pos <- 1
            if (type == "cladogram" && !use.edge.length) 
                node.pos <- 2
        }
        if (node.pos == 1) 
            yy <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
        else {
            ans <- .C("node_height_clado", as.integer(Ntip), 
                as.integer(Nnode), as.integer(z$edge[, 1]), as.integer(z$edge[, 
                  2]), as.integer(Nedge), double(Ntip + Nnode), 
                as.double(yy), DUP = FALSE, PACKAGE = "ape")
            xx <- ans[[6]] - 1
            yy <- ans[[7]]
        }
        if (!use.edge.length) {
            if (node.pos != 2) 
                xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge) - 
                  1
            xx <- max(xx) - xx
        }
        else {
            xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
                z$edge.length)
        }
    }
    else {
        twopi <- 2 * pi
        rotate.tree <- twopi * rotate.tree/360
        switch(type, fan = {
            TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
            xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360, 
                length.out = Ntip)
            theta <- double(Ntip)
            theta[TIPS] <- xx
            theta <- c(theta, numeric(Nnode))
            theta <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, 
                theta)
            if (use.edge.length) {
                r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, 
                  Nedge, z$edge.length)
            } else {
                r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
                r <- 1/r
            }
            theta <- theta + rotate.tree
            xx <- r * cos(theta)
            yy <- r * sin(theta)
        }, unrooted = {
            nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
            XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, 
                z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip, 
                Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
            xx <- XY$M[, 1] - min(XY$M[, 1])
            yy <- XY$M[, 2] - min(XY$M[, 2])
        }, radial = {
            X <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
            X[X == 1] <- 0
            X <- 1 - X/Ntip
            yy <- c((1:Ntip) * twopi/Ntip, rep(0, Nnode))
            Y <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
            xx <- X * cos(Y + rotate.tree)
            yy <- X * sin(Y + rotate.tree)
        })
    }
    if (phyloORclado) {
        if (!horizontal) {
            tmp <- yy
            yy <- xx
            xx <- tmp - min(tmp) + 1
        }
        if (root.edge) {
            if (direction == "rightwards") 
                xx <- xx + x$root.edge
            if (direction == "upwards") 
                yy <- yy + x$root.edge
        }
    }
    if (no.margin) 
        par(mai = rep(0, 4))
    if (is.null(x.lim)) {
        if (phyloORclado) {
            if (horizontal) {
                x.lim <- c(0, NA)
                pin1 <- par("pin")[1]
                strWi <- strwidth(x$tip.label, "inches")
                xx.tips <- xx[1:Ntip] * 1.04
                alp <- try(uniroot(function(a) max(a * xx.tips + 
                  strWi) - pin1, c(0, 1e+06))$root, silent = TRUE)
                if (is.character(alp)) 
                  tmp <- max(xx.tips) * 1.5
                else {
                  tmp <- if (show.tip.label) 
                    max(xx.tips + strWi/alp)
                  else max(xx.tips)
                }
                x.lim[2] <- tmp
            }
            else x.lim <- c(1, Ntip)
        }
        else switch(type, fan = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                x.lim <- c(min(xx) - offset, max(xx) + offset)
            } else x.lim <- c(min(xx), max(xx))
        }, unrooted = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                x.lim <- c(0 - offset, max(xx) + offset)
            } else x.lim <- c(0, max(xx))
        }, radial = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.03 * cex)
                x.lim <- c(-1 - offset, 1 + offset)
            } else x.lim <- c(-1, 1)
        })
    }
    else if (length(x.lim) == 1) {
        x.lim <- c(0, x.lim)
        if (phyloORclado && !horizontal) 
            x.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
            x.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
                cex)
        if (type == "radial") 
            x.lim[1] <- if (show.tip.label) 
                -1 - max(nchar(x$tip.label) * 0.03 * cex)
            else -1
    }
    if (phyloORclado && direction == "leftwards") 
        xx <- x.lim[2] - xx
    if (is.null(y.lim)) {
        if (phyloORclado) {
            if (horizontal) 
                y.lim <- c(1, Ntip)
            else {
                y.lim <- c(0, NA)
                pin2 <- par("pin")[2]
                strWi <- strwidth(x$tip.label, "inches")
                yy.tips <- yy[1:Ntip] * 1.04
                alp <- try(uniroot(function(a) max(a * yy.tips + 
                  strWi) - pin2, c(0, 1e+06))$root, silent = TRUE)
                if (is.character(alp)) 
                  tmp <- max(yy.tips) * 1.5
                else {
                  tmp <- if (show.tip.label) 
                    max(yy.tips + strWi/alp)
                  else max(yy.tips)
                }
                y.lim[2] <- tmp
            }
        }
        else switch(type, fan = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                y.lim <- c(min(yy) - offset, max(yy) + offset)
            } else y.lim <- c(min(yy), max(yy))
        }, unrooted = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                y.lim <- c(0 - offset, max(yy) + offset)
            } else y.lim <- c(0, max(yy))
        }, radial = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.03 * cex)
                y.lim <- c(-1 - offset, 1 + offset)
            } else y.lim <- c(-1, 1)
        })
    }
    else if (length(y.lim) == 1) {
        y.lim <- c(0, y.lim)
        if (phyloORclado && horizontal) 
            y.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
            y.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
                cex)
        if (type == "radial") 
            y.lim[1] <- if (show.tip.label) 
                -1 - max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
            else -1
    }
    if (phyloORclado && direction == "downwards") 
        yy <- max(yy) - yy
    if (phyloORclado && root.edge) {
        if (direction == "leftwards") 
            x.lim[2] <- x.lim[2] + x$root.edge
        if (direction == "downwards") 
            y.lim[2] <- y.lim[2] + x$root.edge
    }
    asp <- if (type %in% c("fan", "radial", "unrooted")) 
        1
    else NA
    plot(0, type = "n", xlim = x.lim, ylim = y.lim, ann = FALSE, 
        axes = FALSE, asp = asp, ...)
    if (plot) {
        if (is.null(adj)) 
            adj <- if (phyloORclado && direction == "leftwards") 
                1
            else 0
        if (phyloORclado && show.tip.label) {
            MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
            loy <- 0
            if (direction == "rightwards") {
                lox <- label.offset + MAXSTRING * 1.05 * adj
            }
            if (direction == "leftwards") {
                lox <- -label.offset - MAXSTRING * 1.05 * (1 - 
                  adj)
            }
            if (!horizontal) {
                psr <- par("usr")
                MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
                  psr[1])
                loy <- label.offset + MAXSTRING * 1.05 * adj
                lox <- 0
                srt <- 90 + srt
                if (direction == "downwards") {
                  loy <- -loy
                  srt <- 180 + srt
                }
            }
        }
        if (type == "phylogram") {
            phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal, 
                edge.color, edge.width, edge.lty)
        }
        else {
            if (type == "fan") {
                ereorder <- match(z$edge[, 2], x$edge[, 2])
                if (length(edge.color) > 1) {
                  edge.color <- rep(edge.color, length.out = Nedge)
                  edge.color <- edge.color[ereorder]
                }
                if (length(edge.width) > 1) {
                  edge.width <- rep(edge.width, length.out = Nedge)
                  edge.width <- edge.width[ereorder]
                }
                if (length(edge.lty) > 1) {
                  edge.lty <- rep(edge.lty, length.out = Nedge)
                  edge.lty <- edge.lty[ereorder]
                }
                circular.plot(z$edge, Ntip, Nnode, xx, yy, theta, 
                  r, edge.color, edge.width, edge.lty)
            }
            else cladogram.plot(x$edge, xx, yy, edge.color, edge.width, 
                edge.lty)
        }
        if (root.edge) 
            switch(direction, rightwards = segments(0, yy[ROOT], 
                x$root.edge, yy[ROOT]), leftwards = segments(xx[ROOT], 
                yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT]), 
                upwards = segments(xx[ROOT], 0, xx[ROOT], x$root.edge), 
                downwards = segments(xx[ROOT], yy[ROOT], xx[ROOT], 
                  yy[ROOT] + x$root.edge))
        if (show.tip.label) {
            if (is.expression(x$tip.label)) 
                underscore <- TRUE
            if (!underscore) 
                x$tip.label <- gsub("_", " ", x$tip.label)
            if (phyloORclado) 
                text(xx[1:Ntip] + lox, yy[1:Ntip] + loy, x$tip.label, 
                  adj = adj, font = font, srt = srt, cex = cex, 
                  col = tip.color)
            if (type == "unrooted") {
                if (lab4ut == "horizontal") {
                  y.adj <- x.adj <- numeric(Ntip)
                  sel <- abs(XY$axe) > 0.75 * pi
                  x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
                    1.05
                  sel <- abs(XY$axe) > pi/4 & abs(XY$axe) < 0.75 * 
                    pi
                  x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
                    (2 * abs(XY$axe)[sel]/pi - 0.5)
                  sel <- XY$axe > pi/4 & XY$axe < 0.75 * pi
                  y.adj[sel] <- strheight(x$tip.label)[sel]/2
                  sel <- XY$axe < -pi/4 & XY$axe > -0.75 * pi
                  y.adj[sel] <- -strheight(x$tip.label)[sel] * 
                    0.75
                  text(xx[1:Ntip] + x.adj * cex, yy[1:Ntip] + 
                    y.adj * cex, x$tip.label, adj = c(adj, 0), 
                    font = font, srt = srt, cex = cex, col = tip.color)
                }
                else {
                  adj <- abs(XY$axe) > pi/2
                  srt <- 180 * XY$axe/pi
                  srt[adj] <- srt[adj] - 180
                  adj <- as.numeric(adj)
                  xx.tips <- xx[1:Ntip]
                  yy.tips <- yy[1:Ntip]
                  if (label.offset) {
                    xx.tips <- xx.tips + label.offset * cos(XY$axe)
                    yy.tips <- yy.tips + label.offset * sin(XY$axe)
                  }
                  font <- rep(font, length.out = Ntip)
                  tip.color <- rep(tip.color, length.out = Ntip)
                  cex <- rep(cex, length.out = Ntip)
                  for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                    cex = cex[i], x$tip.label[i], adj = adj[i], 
                    font = font[i], srt = srt[i], col = tip.color[i])
                }
            }
            if (type %in% c("fan", "radial")) {
                xx.tips <- xx[1:Ntip]
                yy.tips <- yy[1:Ntip]
                angle <- atan2(yy.tips, xx.tips)
                if (label.offset) {
                  xx.tips <- xx.tips + label.offset * cos(angle)
                  yy.tips <- yy.tips + label.offset * sin(angle)
                }
                s <- xx.tips < 0
                angle <- angle * 180/pi
                angle[s] <- angle[s] + 180
                adj <- as.numeric(s)
                font <- rep(font, length.out = Ntip)
                tip.color <- rep(tip.color, length.out = Ntip)
                cex <- rep(cex, length.out = Ntip)
                for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                  x$tip.label[i], font = font[i], cex = cex[i], 
                  srt = angle[i], adj = adj[i], col = tip.color[i])
            }
        }
        if (show.node.label) 
            text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)], 
                x$node.label, adj = adj, font = font, srt = srt, 
                cex = cex)
    }
    L <- list(type = type, use.edge.length = use.edge.length, 
        node.pos = node.pos, show.tip.label = show.tip.label, 
        show.node.label = show.node.label, font = font, cex = cex, 
        adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
        x.lim = x.lim, y.lim = y.lim, direction = direction, 
        tip.color = tip.color, Ntip = Ntip, Nnode = Nnode)
    assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)), 
        envir = .PlotPhyloEnv)
    invisible(L)


}

