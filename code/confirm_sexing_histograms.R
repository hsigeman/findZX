library(doBy)
library(data.table)

# Reads in a bed-file with coverage values for each sample
# Calculates the mean coverage on chr Z and 4 for each individual
# and plots the coverage and the ratio chr z/4

# remove outliers
remove_outliers <- function(data_table) {
  
  outliers <- boxplot(data_table$S)$out
  data_table <- data_table[-which(data_table$S %in% outliers),]
  
  return(data_table)
}



args <- commandArgs(trailingOnly = TRUE)


filename = args[1]
outfilename = args[2]
synteny = args[3]
sample_names <- args[4:length(args)]

# read gencov
cov = read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)

if (synteny == "with-synteny") {
  cov <- cov[-1:-10][-1:-3]
}

outliers <- boxplot(cov$V14)$out
cov$V14 <- cov[-which(cov$V14 %in% outliers),]

# normalize

par(mfrow=c(2,5), mar=c(1,1,1,1), oma=c(0,0,0,0), xpd=TRUE)

hist(cov$V14, breaks = 50)
hist(cov$V15, breaks = 50)
hist(cov$V16, breaks = 50)
hist(cov$V17, breaks = 50)
hist(cov$V18, breaks = 50)
hist(cov$V19, breaks = 50)
hist(cov$V20, breaks = 50)
hist(cov$V21, breaks = 50)
hist(cov$V22, breaks = 50)
hist(cov$V23, breaks = 50)
