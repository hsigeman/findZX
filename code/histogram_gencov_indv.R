library(doBy)
library(data.table)
library(ggplot2)
library(cowplot)
require(scales)

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
  column_normalized <- as.data.frame(column_normalized)
  
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
outfile = args[3]
synteny = args[4]
sample_names = args[5:length(args)]

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

nr_samples <- length(sample_names)

# remove columns: chr and range (start, stop)
cov <- cov[-length(cov)][-1:-3]
het_mean <- het_mean[-1:-2]

# remove outliers for each sample, save resulting array in a list
# samples should be able to have different lengths
cov_no_outliers <- apply(cov, 2, remove_outliers)

# normalize on median of each sample and take the logarithm
cov_norm <- lapply(cov_no_outliers, normalize)

################################################################################
################################### PLOTTING ###################################
################################################################################

plist <- list()

for (i in 1:nr_samples) {

###################### HETEROZYGOSITY VS GENOME COVERAGE #######################
  
  x = i + 2
  y = i + 2 + nr_samples
  
  outliers <- boxplot(cov_het[,x], plot = FALSE)$out
  no_outliers <- cov_het[!(cov_het[,x] %in% outliers),]
  
  p <- ggplot(data = no_outliers, aes(x = no_outliers[,x], y = no_outliers[,y])) + geom_bin2d(bins = 10)
  p <- p + labs(x = "genome coverage", y = "heterozygosity") + theme_bw()
  
############################### GENOME COVERAGE ################################
  
  df <- cov_norm[[i]]
  hg <- ggplot(df, aes(x = column_normalized)) + geom_histogram(bins = 100)
  hg <- hg + labs(x="genome coverage, normalized on median", y="Frequency")
  hg <- hg + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))
  hg <- hg + annotation_logticks(sides="b") + theme_bw()
  
################################ HETEROZYGOSITY ################################
  
  hh <- ggplot(het_mean, aes(x = het_mean[,1])) + geom_histogram(bins = 100)
  hh <- hh + labs(x="heterozygosity", y="Frequency")
  hh <- hh + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))
  hh <- hh + annotation_logticks(sides="b") + theme_bw()
  
  
  pg <- plot_grid("", "", "", p,hg,hh,ncol = 3, rel_heights = c(1,20),
                  labels = c(sample_names[i], "", "", "A", "B", "C"))
  plist[[i]] <- pg
  
}


pdf(outfile, width = 20)

for (i in 1:nr_samples) {
  
  print(plist[[i]])
  
}

dev.off()

