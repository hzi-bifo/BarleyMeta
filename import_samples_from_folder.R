#' This function import all samples from a folder and create a data frame with
#' all information needed for the boxplot
#' 
#' @param data_folder : data folder where dnds values for samples A-F are stored
#' @keywords dnds, figure
#' @export
#' @examples 
#' import_samples_from_folder()
#' 
import_samples_from_folder <- function(data_folder = "example_data/FigureA/"){
  
  files <-  list.files(data_folder)
  sample_names <-sub("^([^.]*).*", "\\1", files )
  
  # create matrix of protein names
  merged_names <- ""
  for (file in files){
    cat(file)
    file_data <- read.table(paste(data_folder,file,sep=""), sep="", header=F)
    merged_names <- c(merged_names, as.character(as.matrix(file_data$V1)))
  }
  names <- unique(merged_names)
  
  # create matrix
  samples <- matrix(ncol=length(files) + 1, nrow=length(names))
  samples[,1] <- names
  colnames(samples) <- c("names",sample_names)
  
  # fill out the matrix with sample information
  i <- 2
  for (file in files){
    file_data <- read.table(paste(data_folder,file,sep=""), sep="", header=F)
    file_data$dnds <- file_data$V2 / file_data$V3
    samples[,i] <- file_data$dnds[match(samples[,1], as.character(as.matrix(file_data$V1)))]
    i <- i + 1 
  }
  samples <- samples[which(samples[,1] !=""),] #remove empty row
  return(data.frame(samples))
}