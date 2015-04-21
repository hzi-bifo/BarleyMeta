#' This function creates start and end coordinates of clusters froma list of pvalues
#' @param pvalues: Some p-values
#' @param significance_level: p-values below this value are considered as significant
#' @param w_size: size of the sliding window
#' @keywords dnds, figure
#' @export
#' @examples 
#' create_cluster_matrix()
#' 
create_cluster_matrix <- function(pvalues,
                          significance_level=0.05,
                          w_size=10){
  iteration <- 1 
  i <- 1
  clusters <- NULL
  
  while (!is.na (which(pvalues[i:length(pvalues)] <= significance_level)[1])) { # check if there are some significant elements
    if (i==1){
      i <- which(pvalues[i:length(pvalues)] <= significance_level)[1]  # first sig elem
      clusters$start[iteration] <- i
      i <- which(pvalues[i:length(pvalues)] > significance_level)[1] + clusters$start[iteration] - 1 # last sig. elem.
      clusters$end[iteration] <- i 
      iteration <- iteration + 1
    } else {
      i <- which(pvalues[i:length(pvalues)] <= significance_level)[1] + clusters$end[iteration-1] 
      clusters$start[iteration] <- i
      i <- which(pvalues[i:length(pvalues)] > significance_level)[1] + clusters$start[iteration] - 1
      clusters$end[iteration] <- i 
      iteration <- iteration + 1
    }
  }
  if(is.null(clusters)){ # add some zero values to prevent plotting issues
    clusters$start[1] <- 0
    clusters$end[1] <- 0
    clusters <- as.data.frame(clusters)
  } else {
    clusters <- as.data.frame(clusters)
    clusters$start <- clusters$start - w_size/2
    clusters$end <- clusters$end + w_size/2
  }
  return(clusters)
}