library(doBy)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

filename = args[1]
outfilename = args[2]

het = read.table(filename,header=TRUE,fill=TRUE,stringsAsFactor=FALSE)

sample_names <- colnames(het)[-1:-3]

het <- transform(het, range=round(end/5000))

het_mean <- summaryBy(sample_names ~ chr + range, data=het, keep.names=TRUE)

het_mean <- het_mean[-1:-2]

het_list <- list(het_mean[1:6])


pdf(outfilename, width = 14)

par(mfrow=c(2,(length(het_mean)/2)), mar=c(5,4,4,2), oma=c(0,1,0,0))

args <- list(xlab="mean heterozygosity 5kbp windows")
titles <- colnames(het_mean)

mapply(hist, het_mean, breaks = 100, main=titles, xlab="mean heterozygosity 5kbp windows")

dev.off()
