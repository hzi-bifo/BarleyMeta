#' This function create the sequence graphic for Figure 6.
#' Sequence clusters of residues under positive selection in selected protein families. 
#' Top: dots indicate ∼Dn/Ds for a given position in the protein sequence, and 
#' their color corresponds to the proportion of gaps in the multiple sequence 
#' alignment (MSA). Gray-shaded areas indicate significant clusters of residues 
#' under positive selection. Gray-shaded horizontal lines indicate repetitive 
#' elements. Bottom: Jensen-Shannon divergence as a function of the positions 
#' in the MSA.
#' 
#' @param data_folder path to protein specific folder where the 
#' data for each protein familiy is stored. Defaults to example_data/TIGR00450/.
#' @param only_legend outputs only the figure legend Default FALSE
#' @param window_size size for the sliding window to detect residues under positive
#' selection DEFAULT 10
#' @param gap_threshold positions with a gap threshold higher than this variable
#' will be removed from sliding window statistics
#' @keywords dnds, figure
#' @export
#' @examples 
#' create_figure_6B()
create_figure_6B <- function(data_folder = "data/figure6/FigureB/TIGR01573/",
                          only_legend = FALSE,
                          window_size = 10,
                          gap_threshold = 0.6){
  #### load dependencies ####
  require(SpatioTemporal)
  require(ggplot2)
  require(gridExtra)
  require(zoo)
  source("g_legend.R")
  source("create_cluster_matrix.R")
  source("sliding_window.R")
  
  #### get dataset information ####
  data_files <- basename(list.files(data_folder))
  sample_name <-sub("^([^.]*).*", "\\1", data_files[1] )

  #### import data #### 
  gap_data <- try(read.table(paste(data_folder, data_files[grep(".gapdata", data_files)],sep="")))
  if (inherits(gap_data, 'try-error')){  cat("No gap information provided. Please ad .gapdata \n")} 
  dnds_data <- try(read.table(paste(data_folder, data_files[grep(".dndsdata", data_files)],sep=""),header=T))
  if (inherits(dnds_data, 'try-error')){  cat("No dnds information provided. Please ad .dndsdata \n")} 
  js_data <- try(read.table(paste(data_folder, data_files[grep(".jsdata", data_files)] ,sep="")))
  if (inherits(js_data, 'try-error')){  cat("No js information provided. Please ad .jsdata \n")} 
  repeat_data <- try(read.table(paste(data_folder, data_files[grep(".repeatdata", data_files)],sep="")))
  if (inherits(repeat_data, 'try-error')){  cat("No repeat information provided. Please ad .repeatdata \n")} 
  tm_data <- try(read.table(paste(data_folder, data_files[grep(".tmdata", data_files)],sep=""), sep=";"))
  if (inherits(tm_data, 'try-error')){  cat("No tm information provided. Please ad .tmdata \n")} 

  #### check data ####
  # check for length
  if (nrow(gap_data) != nrow(dnds_data)){ stop("Number of rows in dnds data is not equal number of rows of gap data")}  
  
  #### prepare data frame ####
  # process the dnds data
  dN <- dnds_data$tree_nonsyn + dnds_data$tree_nonsyn 
  dS <- dnds_data$leaf_syn + dnds_data$tree_syn
  dNdS <- rep(0, length(dN))
  dNdS[which(dS > 0)] <- dN[which(dS > 0)]/dS[which(dS > 0)]

  # process the dnds data with sliding window
  pvalue_window <- sliding_window(dN, dS, gap_data =  gap_data, 
                                  w_size = window_size,
                                  g_threshold = gap_threshold)
  
  # generate a matrix with start and end positions for plotting of significant windows
  epitopes <- create_cluster_matrix(pvalues = pvalue_window, w_size = window_size, 
                                    significance_level = 0.05)
  
  # create data frame with all informations
  df <- data.frame(position=dnds_data$pos,
                   selection=dNdS,
                   selection_smooth=smooth.spline(dnds_data$pos, dNdS, spar=0.35)$y,
                   gap=gap_data$V1,
                   gap_smooth=smooth.spline(dnds_data$pos, gap_data$V1, spar=0.35)$y,
                   js=js_data$V1, 
                   js_smooth=smooth.spline(dnds_data$pos, js_data$V1, spar=0.35)$y,
                   stringsAsFactors=FALSE) 
  
  # process repeat informations if supplied 
  # generating a boolean matrix with TRUE if there is a repeat element on this position
  if (!inherits(repeat_data, 'try-error')){ # if annotation is provided
    df$is_repeat <- rep(FALSE, length(dNdS))
    if (!is.null(repeat_data)){ # if there are repetitive elements
      for (segment in 1:nrow(repeat_data)){
        try(df$is_repeat[repeat_data$V1[segment]:repeat_data$V2[segment]] <- TRUE,TRUE)
      }
    }
  } 
  
  # process membrane informations if supplied 
  if (!inherits(tm_data, 'try-error')){ # if annotation is provided
    df$is_inside_membrane <- rep(FALSE, length(dNdS))
    df$is_outside_membrane <- rep(FALSE, length(dNdS))
    df$is_tm_helix <- rep(FALSE, length(dNdS))
    if (!is.null(tm_data)){ # if there are repetitive elements
      for (segment in 1:nrow(tm_data)){
        #TODO: add tm infomation to data frame 
      }
    }
  } 
  
  #### genereate plot #### 
  # create sequence plot figure with dnds and gap informations
  fig_a <- ggplot(df, aes(position, selection))
  # highlight significant aereas
  fig_a <- fig_a + geom_rect(data=epitopes, aes(NULL, NULL, xmin = start, xmax = end), 
                             ymin = -Inf, ymax = Inf,  fill="grey80")
  
  fig_a <- fig_a + geom_ribbon(aes(x=position, ymax=selection_smooth,ymin=1), fill="grey60")
  fig_a <- fig_a + geom_point(aes(colour=gap), size = 1.7, alpha=3/4)
  fig_a <- fig_a + scale_colour_gradient(low = "green", high="red")
  fig_a <- fig_a + theme_bw() + xlab(" ") + ylab("dN/dS ratio") 
  fig_a <- fig_a + geom_line(aes(y=selection_smooth))
  fig_a <- fig_a + labs(title=paste(sample_name,"; g_threshold=", gap_threshold,"; window_size=", window_size,sep=""))
  fig_a <- fig_a + geom_hline(aes(yintercept=1))
 
  # plot repeat informations
  if (!inherits(repeat_data, 'try-error')){
    fig_a <- fig_a + geom_segment(data=repeat_data,aes(x = V1, y = -1.2, xend = V2, yend = -1.2), color = "grey")
  }
 
  # create sequence plot figure with JS informations
  fig_b <- ggplot(df, aes(position, js_smooth)) 
  fig_b <- fig_b + geom_line(aes(y=js_smooth))
  fig_b <- fig_b + theme_bw() + xlab("Sequence Position (nt)") + ylab("JS score") 
  fig_b <- fig_b + scale_y_continuous(breaks=c(0,1),limits=c(0,1)) 

  # plot membrane informations
  if (!inherits(tm_data, 'try-error')){
    fig_b <- fig_b + geom_rect(data=tm_data,aes(NULL, NULL, xmin = V2, xmax = V3, fill=V1), ymin = -Inf, ymax = Inf,alpha=0.2)
  }
  
  if (only_legend){
    fig_b_legend <- g_legend(fig_b)  
    fig_a_legend <- g_legend(fig_a)
    fig <- grid.arrange(fig_b_legend, fig_a_legend, nrow=2, heights=1:3)
    
  } else {
    fig_b <- fig_b + theme(legend.position="none")
    fig_a <- fig_a + theme(legend.position="none")    
    fig <- grid.arrange(fig_a, fig_b, nrow=2, heights=c(3/4,1/4), ncol=1)    
  } 
  return(fig_a)
}