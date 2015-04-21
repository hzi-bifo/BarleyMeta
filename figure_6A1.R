#' script to reproduce the analysis and figures from Bulgarelli et al., 2015
#'
#' if you use any of the following code, please cite:
#'
#' Davide Bulgarelli, Ruben Garrido-Oter, Philipp C. Münch, Aaron Weiman,
#' Johannes Dröge, Yao Pan, Alice C. McHardy, Paul Schulze-Lefert
#' Structure and Function of the Bacterial Root Microbiota in Wild
#' and Domesticated Barley. Cell Host and Microbe, 2015
#'
# originally by Philipp C. Münch

#' This function create the boxplot graphic Figure 6A.
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
#' create_figureA()
#' 
create_figure_6A1 <- function(input_path = "data/figure6/FigureA/figure_6a_data.csv",
                           annotation_path = "data/figure6/FigureA/figure_6a_data.csv"){
  #### load dependencies ####
  require(ggplot2)
  require(reshape2)
  require(matrixStats)
  
  #### Import preprocessed dnds data ####
  dnds_data <- read.csv2(input_path, sep=";", header=T)
  dnds_data <- dnds_data[,-1]
  # process annotation files
  annotation_files <- list.files(annotation_path)
  annotation_names <- sub("^([^.]*).*", "\\1", annotation_files)
#  i <- 1
  # iterate over each annotation file
#  for (annotation_file in annotation_files){
#    dnds_data[[annotation_names[i]]] <- FALSE
#    annotation_data <- read.table(paste(annotation_path,annotation_file,sep=""), sep="", header=F)
#    hits <- match(as.character(as.matrix(dnds_data$Name)),
#          as.character(as.matrix(annotation_data$V1)))
#    dnds_data[[annotation_names[i]]][hits] <- TRUE
#    i <- i + 1
#  }
 
  # filter out repeats
  dnds_data <- dnds_data[-grep("repeat",dnds_data$annotation),]

  # create dataframe for plotting
  df <- melt(dnds_data)

  p <- ggplot(df, aes(x= annotation, y=value/2)) + geom_boxplot(fill="grey80")+ coord_flip() + theme_bw()
  p <- p + xlim(as.character(as.matrix(dnds_data$annotation)))

#  p <- p + xlim(0,1.8) # + ylab("dN/dS value") + xlab(" ")
  p <- p + theme(axis.ticks=element_blank(), panel.border=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank()) # for 2 plots 2
  p <- p + geom_hline(yintercept = 0.759840522, color="darkgreen", linetype=5) 
  

  return(p)
}