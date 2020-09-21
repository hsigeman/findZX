library(doBy)
library(data.table)
library(ggplot2)
library(cowplot)
require(scales)

# Reads in a bed-file with coverage values and a file with heterozygosity for 
# each site for each sample, values for each sample are in separate columns.
# Removes outliers and makes two histograms and one scatterplot for each individual:
# histograms for genome coverage and heterozygosity and a scatterplot of genome
# coverage vs heterozygosity.
# Heterozygosity is given for each variable site. 1 = heterozygot, 0 = homozygot,
# na = missing data.

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

length_chr <- function(data_table) {
  # Returns the length of each chr/scaffold as a dataframe.
  
  max_per_chr <- setDT(data_table)[, .SD[which.max(end)], by=chr]
  
  max_per_chr <- as.data.frame(max_per_chr[,1:2])
  
  return(max_per_chr)
}

################################################################################
################################################################################
################################################################################

set.seed(999)

args <- commandArgs(trailingOnly = TRUE)

file_gencov = args[1]
file_snp = args[2]
outfile = args[3]
synteny = args[4]
sample_names = args[5:length(args)]

################################################################################
################################# READ FILES ###################################
################################################################################

cov <- read.table(file_gencov,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
het <- read.table(file_snp,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)


if (synteny == "with-synteny") {
  cov <- cov[-1:-10]
  het <- het[-1:-10]
}

colnames(cov)[1:3] <- c("chr","start","end")
colnames(het)[1:3] <- c("chr","start","end")

################################################################################
################################ CALCULATIONS ##################################
################################################################################

col_names <- colnames(cov)[4:length(cov)]
colnames(cov) <- c("chr","start","end",sample_names)

Thet <- transform(het, range=floor(end/5000))
len_chr <- length_chr(Thet)
colnames(len_chr) <- c("chr","length")

f <- as.formula(paste(paste(col_names, collapse = "+"), "~", "chr + range"))
het_mean <- summaryBy(f, data=Thet, keep.names=TRUE, na.rm = TRUE)

colnames(het_mean) <- c("chr","range",sample_names)

Tcov <- transform(cov, range=floor(start/5000))

cov_het <- merge(Tcov, het_mean, by = c("chr","range"))
cov_het <- cov_het[-3:-4]
cov_het <- merge(cov_het, len_chr, by = "chr")

nr_samples <- length(sample_names)

# remove columns: chr and range (start, stop)
cov <- cov[-1:-3]
het_mean <- het_mean[,-1:-2]
het_mean <- as.data.frame(het_mean)

# remove outliers for each sample, save resulting array in a list
# samples should be able to have different lengths
cov_no_outliers <- apply(cov, 2, remove_outliers)

# normalize on median of each sample
#cov_norm <- lapply(cov_no_outliers, normalize)
cov_norm <- cov_no_outliers

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
  
  p <- ggplot(data = no_outliers, aes(x = no_outliers[,x], y = no_outliers[,y])) + 
    geom_bin2d() + 
    labs(x = "genome coverage", y = "heterozygosity") + 
    theme_bw() +
    scale_fill_gradient(low="white",high="darkblue",trans="log10") +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels = trans_format("log10", math_format(10^.x))) + 
    annotation_logticks(sides="b")
  
  gl <- ggplot(data = no_outliers, aes(x = no_outliers[,x], y = length)) + 
    geom_bin2d() + 
    labs(x = "genome coverage", y = "scaffold length") + 
    theme_bw() +
    scale_fill_gradient(low="white",high="darkblue",trans="log10")
  
  hl <- ggplot(data = no_outliers, aes(x = no_outliers[,y], y = length)) + 
    geom_bin2d() + 
    labs(x = "heterozygosity", y = "scaffold length") + 
    theme_bw() +
    scale_fill_gradient(low="white",high="darkblue",trans="log10")
  
############################### GENOME COVERAGE ################################
  
  df <- as.data.frame(cov_norm[[i]])
  hg <- ggplot(df, aes(x = cov_norm[[i]])) + 
     geom_histogram() + 
     labs(x="genome coverage", y="Frequency") + 
     scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                      labels = trans_format("log10", math_format(10^.x))) + 
     annotation_logticks(sides="b") + 
     theme_bw()
  
################################ HETEROZYGOSITY ################################
  
  hh <- ggplot(het_mean, aes(x = het_mean[,i])) + 
     geom_histogram() + 
     labs(x="heterozygosity", y="Frequency") + 
     theme_bw()
  
  
  pg <- plot_grid("","","","","",p,gl,hl,hg,hh,ncol = 5, rel_heights = c(1,20),
                  labels = c(sample_names[i],"","","","","A","B","C","D","E"))
  plist[[i]] <- pg
  
}


pdf(outfile, width = 30)

for (i in 1:nr_samples) {
  
  print(plist[[i]])
  
}

dev.off()

