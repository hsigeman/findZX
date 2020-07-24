library(doBy)
library(data.table)

# Reads in a bed-file with coverage values for each sample
# Removes outliers and makes a histogram for each individual
# over the mean genome coverage in 5kbp windows.

################################################################################
################################## FUNCTIONS ###################################
################################################################################

remove_outliers <- function(column) {
  
  outliers <- boxplot(column, plot = FALSE)$out
  column <- column[!(column %in% outliers)]
  
  return(column)
}

normalize <- function(column) {
  
  column_median <- median(column)
  column_normalized <- column/column_median
  
  return(column_normalized)
}

################################################################################
################################################################################
################################################################################

args <- commandArgs(trailingOnly = TRUE)

filename = args[1]
outfilename = args[2]
synteny = args[3]
sample_names = args[4:length(args)]

# read gencov
cov = read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)

if (synteny == "with-synteny") {
  cov <- cov[-1:-13]
} else {
  cov <- cov[-1:-3]
}

colnames(cov) <- c(sample_names)

# remove outliers for each sample, save resulting array in a list
# samples should be able to have different lengths
cov_no_outliers <- apply(cov, 2, remove_outliers)

# normalize on median of each sample
cov_norm <- lapply(cov_no_outliers, normalize)

# plot histogram for each sample
pdf(outfilename, width = 14)

par(mfrow=c(2,(length(cov_norm)/2)), mar=c(5,4,4,2), oma=c(0,1,0,0))

args <- list(xlab="gen.cov. normalized on median")
titles <- colnames(cov)
mapply(hist, cov_norm, breaks = 100, main=titles, MoreArgs=args)

dev.off()
