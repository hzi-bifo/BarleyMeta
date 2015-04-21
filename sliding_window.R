#' This function creates a p-value matrix indicating significant positions based
#' on dN/dS statistic
#' @param dN: dN values per position in the alignment
#' @param dS: dS values per position in the alignment
#' @param gap_data: gap proportion per position in the alignment
#' @param w_size: size of the sliding window DEFAULT 10 
#' @param g_threshold: positions with a gap proportion higher than this threshold
#' will be removed. DEFAULT 75
#' @keywords dnds, figure
#' @export
#' @examples 
#' sliding_window()
sliding_window <- function(dN, dS, gap_data, w_size=10, g_threshold=0.2){
  require(zoo)
  dNdS_window <- rep(0,nrow(dnds_data))
  if (w_size>0){
    dN_window <- rollsum(dN, w_size, fill = list(NA, NULL, NA))
    dS_window <- rollsum(dS, w_size, fill = list(NA, NULL, NA))
    dNdS_window[which(dS_window > 0)] <- dN_window[which(dS_window > 0)] / 
      dS_window[which(dS_window > 0)]
    dNdS_gap <- rollmean(gap_data$V1, w_size, fill=list(NA,NULL,NA))
    dNdS_gap[which(is.na(dNdS_gap))] <- 1
  }

  # iterate over window and find significant windows with FDR correction
  dN_all <- sum(dN, na.rm=T)
  dS_all <- sum(dS, na.rm=T)
  window <- rep(1, nrow(dnds_data))
  for (i in 1:length(dNdS_window)){
    if (dNdS_gap[i] < g_threshold){ 
      if (!is.na(dN_window[i]) && !is.na(dS_window[i])){
        test_matrix <-  matrix(c(dN_window[i], dN_all,
                                 dS_window[i], dS_all),nrow=2)
        pval <- fisher.test(test_matrix,alternative=c("greater"))[["p.value"]]
        window[i] <- p.adjust(pval, method="fdr",n = length(window))
      } else {
        window[i] <- 1
      }
    } else { # g_threshold
      window[i] <- 1
    }
  }
  return(window)
}