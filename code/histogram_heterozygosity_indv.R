library(doBy)
library(data.table)

# Reads in a file which for each sample and snp have a 1 if the individual is
# heterozygot and 1 if homozygot and NA if the genotype is missing.
# Calculates the mean heterozygosity per snp in 5kbp windows and shows 
# it in a historgram for each individual.

args <- commandArgs(trailingOnly = TRUE)

filename = args[1]
outfilename = args[2]
synteny = args[3]
sample_names = args[4:length(args)]

het <- read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)

if (synteny == "with-synteny") {
  het <- het[-1:-7][-4:-6]
}

sample_names <- gsub("-", "_", sample_names)

colnames(het) <- c("chr", "start", "end", sample_names)

Thet <- transform(het, range=round(end/5000))

f <- as.formula(paste(paste(sample_names, collapse = "+"), "~", "chr + range"))
het_mean <- summaryBy(f, data=Thet, keep.names=TRUE)

het_mean <- het_mean[-1:-2]


pdf(outfilename, width = 14)

par(mfrow=c(2,(length(het_mean)/2)), mar=c(5,4,4,2), oma=c(0,1,0,0))

args <- list(xlab="mean heterozygosity 5kbp windows")
titles <- colnames(het_mean)

mapply(hist, het_mean, breaks = 100, main=titles, xlab="mean heterozygosity 5kbp windows")

dev.off()
