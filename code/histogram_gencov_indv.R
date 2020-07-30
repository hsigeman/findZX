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
  column_normalized <- log(column/column_median)
  
  return(column_normalized)
}

new_name <- function(n,sufix) { 
  
  n <- paste(n,"", sep = sufix) 
  return(n) 
}

################################################################################
################################################################################
################################################################################

args <- commandArgs(trailingOnly = TRUE)

file_gencov = args[1]
file_snp = args[2]
outfile_gencov = args[3]
outfile_snp = args[4]
outfile_scatt = args[5]
synteny = args[6]
sample_names = args[7:length(args)]

# read gencov
cov = read.table(file_gencov,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
het <- read.table(file_snp,header=TRUE,fill=TRUE,stringsAsFactor=FALSE)


if (synteny == "with-synteny") {
  cov <- cov[-1:-10]
  het <- het[-1:-7][-4:-6]
}

colnames(cov) <- c("chr","start","end",sample_names)
colnames(het) <- c("chr","start","end",sample_names)


Thet <- transform(het, range=floor(end/5000))

f <- as.formula(paste(paste(sample_names, collapse = "+"), "~", "chr + range"))
het_mean <- summaryBy(f, data=Thet, keep.names=TRUE, na.rm = TRUE)

cov$range <- cov$start/5000

cov_het <- merge(cov, het_mean, by = c("chr","range"))
cov_het <- cov_het[-3:-4]

#remove outliers and normalize before plotting

pdf(outfile_scatt, width = 14)
par(mfrow=c(2,((length(cov_het)-2)/4)), mar=c(5,4,4,2), oma=c(0,1,0,0))

nr_samples <- length(sample_names)

for (i in 1:nr_samples) {
  x = i + 2
  y = i + 2 + nr_samples
  plot(log(cov_het[,x]), log(cov_het[,y]), xlab = "log of gen.cov. [10^x]", pch = 20,
       ylab = "log of mean heterozygosity", main = paste(sample_names[i], ", 5kbp windows", sep = ""))
}

dev.off()


# remove columns: chr and range (start, stop)
cov <- cov[-length(cov)][-1:-3]
het_mean <- het_mean[-1:-2]

# remove outliers for each sample, save resulting array in a list
# samples should be able to have different lengths
cov_no_outliers <- apply(cov, 2, remove_outliers)

# normalize on median of each sample and take the logarithm
cov_norm <- lapply(cov_no_outliers, normalize)



# plot histogram for each sample
pdf(outfile_gencov, width = 14)
par(mfrow=c(2,(length(cov_norm)/2)), mar=c(5,4,4,2), oma=c(0,1,0,0))

args <- list(xlab="log of gen.cov. normalized on median [10^x]",
             ylab="Frequency of 5kbp windows")
titles <- colnames(cov)
mapply(hist, cov_norm, breaks = 100, main=titles, MoreArgs=args)

dev.off()



pdf(outfile_snp, width = 14)
par(mfrow=c(2,(length(het_mean)/2)), mar=c(5,4,4,2), oma=c(0,1,0,0))

args <- list(xlab="log of mean heterozygosity [10^x]",
             ylab="Frequency of 5kbp windows")
titles <- colnames(het_mean)

mapply(hist, log(het_mean), breaks = 100, main=titles, MoreArgs=args)

dev.off()
