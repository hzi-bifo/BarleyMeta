#' This function create the boxplot graphic Figure 6A (top)
#' Top-ranking protein families under positive selection with significantly 
#' increased Dn/Ds statistic. The distribution at the top shows the density 
#' function over all protein families smoothed with a Gaussian kernel function. 
#' The green bar indicates the average ∼Dn/Ds over all the samples, the blue 
#' bar the average ∼Dn/Ds for all TIGFRAMS annotated with the term “patho” 
#' and/or “secretion.” The boxplot shows the distribution of the ∼Dn/Ds across 
#' all samples for the top 50 ranked TIGRFAM families under positive selection,
#' with families sorted by their median ∼Dn/Ds in descending order. TIGRFAMs 
#' annotated with “repeat” or with a mean repetitive value of more than 50% 
#' were discarded
#' 
#' @param input_path : data folder where dnds values for samples A-F are stored
#' @keywords dnds, figure
#' @export
#' @examples 
#' create_figure6A2()

create_figure_6A2 <- function(data_folder = "data/figure6/FigureA/protein_familiy_data.csv",
                           rep = 0.5,
                           use_gap_corrected=T,
                           seperate=F){
  require(ggplot2)
  require(reshape2)
  dat <- read.csv2(data_folder, sep=";", header=T)
  # filter rep
  dat <- dat[which(as.numeric(as.matrix(dat$meanRep)) < rep),]
  # make df
  dat$median.dnds <- as.numeric(as.matrix(dat$median.dnds))
  if(use_gap_corrected){
    big <- data.frame(dat[,42:47])
    big <- cbind(dat[,1], big)
    colnames(big) <- c("Name", "A", "B", "C", "D", "E", "F")
  } else {
    big <- data.frame(dat[,1:7])
  }
  
  big_m <- melt(big, id.vars='Name')
  big_m$value <- as.numeric(as.matrix(big_m$value))/2

  the_mean <- mean(big_m[is.finite(big_m$value),]$value, na.rm=T)
  cat(paste("Mean value = ", the_mean,"\n"))
  # plot with median data
  m <- ggplot(dat, aes(x = median.dnds)) + geom_density(position="identity") + xlim(0,1.8)
  
  # plot with all data
  
  if(seperate){
    m <- ggplot(big_m, aes(x = value, fill=variable)) 
    m <- m + geom_density(position="identity", alpha=0.5) + theme_bw() + xlim(0,1.8)
  } else {
    m <- ggplot(big_m, aes(x = value))
    m <- m + geom_density(position="identity", fill="grey80") + theme_bw() + xlim(0,1.8)
  }
   
  m <- m + geom_vline(xintercept = the_mean, color="darkgreen", linetype=5)
  return(m)
}
