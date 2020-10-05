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
chr_file = args[5]
sample_names = args[6:length(args)]

################################################################################
################################# READ FILES ###################################
################################################################################

cov <- read.table(file_gencov,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
het <- read.table(file_snp,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)


if (synteny == "with-synteny") {
  cov <- cov[-c(1:7,11:13)]
  het <- het[-c(1:7,11:13)]
}

colnames(cov)[1:3] <- c("chr","start","end")
colnames(het)[1:3] <- c("chr","start","end")

################################################################################
################################ CALCULATIONS ##################################
################################################################################

col_names <- colnames(cov)[4:length(cov)]

Thet <- transform(het, range=floor(end/5000))

len_chr <- length_chr(cov)
colnames(len_chr) <- c("chr","length")

f <- as.formula(paste(paste(col_names, collapse = "+"), "~", "chr + range"))
Thet <- summaryBy(f, data=Thet, keep.names=TRUE, na.rm = TRUE)
colnames(Thet) <- c("chr","range",sample_names)

Tcov <- transform(cov, range=floor(start/5000))
Tcov <- summaryBy(f, data=Tcov, keep.names=TRUE, na.rm = TRUE)
colnames(Tcov) <- c("chr","range",sample_names)

cov_het <- merge(Tcov, Thet, by = c("chr","range"))
cov_het <- merge(cov_het, len_chr, by = "chr")

nr_samples <- length(sample_names)
cov_het <- as.data.frame(cov_het)

for (i in 1:nr_samples) {
  
  outliers <- boxplot(cov_het[,(i+2)], plot = FALSE)$stats[5]
  cov_het[,(i+2)][cov_het[,(i+2)] > outliers] = NA
  
}

if (file.exists(chr_file)) {
  
  chromosome <- read.csv(chr_file, header = FALSE, sep = ",")
  chromosome <- as.factor(chromosome)
  
  cov_het <- cov_het[cov_het$chr %in% chromosome, ]
  
}

################################################################################
################################### PLOTTING ###################################
################################################################################

plist <- list()

for (i in 1:nr_samples) {

###################### HETEROZYGOSITY VS GENOME COVERAGE #######################
  
  x = i + 2
  y = i + 2 + nr_samples
  
  gh <- ggplot(data = cov_het, aes(x = cov_het[,x], y = cov_het[,y])) + 
    geom_bin2d(na.rm = TRUE) + 
    labs(x = "genome coverage", y = "heterozygosity", title = sample_names[i]) + 
    theme_bw() +
    scale_fill_gradient(low="white",high="darkblue",trans="log10")
  
  gl <- ggplot(data = cov_het, aes(x = cov_het[,x], y = length)) + 
    geom_bin2d(na.rm = TRUE) + 
    labs(x = "genome coverage", y = "scaffold length", title = sample_names[i]) + 
    theme_bw() +
    scale_fill_gradient(low="white",high="darkblue",trans="log10")
  
  hl <- ggplot(data = cov_het, aes(x = cov_het[,y], y = length)) + 
    geom_bin2d(na.rm = TRUE) + 
    labs(x = "heterozygosity", y = "scaffold length", title = sample_names[i]) + 
    theme_bw() +
    scale_fill_gradient(low="white",high="darkblue",trans="log10")
  
############################### GENOME COVERAGE ################################
  
  g <- ggplot(cov_het, aes(x = cov_het[,x])) + 
     geom_histogram(bins = 50, na.rm = TRUE) + 
     labs(x="genome coverage", y="Frequency", title = sample_names[i]) + 
     theme_bw()
  
  # Find the x value for the bin with the highest count
  hist_stats_hg <- ggplot_build(g)$data[[1]]
  max_bin <- which(hist_stats_hg$density == max(hist_stats_hg$density))
  
  # Calcualte half of the x value for the highest bin
  max_x <- hist_stats_hg$x[max_bin]
  halfMax_x <- max_x / 2
  
  g <- g + geom_vline(xintercept=halfMax_x, color = "blue") + 
    geom_vline(xintercept=max_x, color = "red")
  
################################ HETEROZYGOSITY ################################
  
  h <- ggplot(cov_het, aes(x = cov_het[,y])) + 
     geom_histogram(bins = 50, na.rm = TRUE) + 
     labs(x="heterozygosity", y="Frequency", title = sample_names[i]) + 
     theme_bw()
  
################################# CHROMOSOMES ##################################
  
  if (dim(len_chr)[1] >= 50) {
    nr_chr <- 50
  } else {
    nr_chr <- dim(len_chr)[1]
  }
  
  top_chr <- as.factor(len_chr[with(len_chr, order(-length)),][1:nr_chr,1])
  
  cov_het_subset <- cov_het[cov_het$chr %in% top_chr,]
  cov_het_subset$chr <- factor(cov_het_subset$chr, levels = top_chr)
  
  c <- ggplot(cov_het_subset, aes(x=range, y=cov_het_subset[,x])) + 
    geom_smooth(na.rm = TRUE) + 
    facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
    labs(y = "genome coverage", x = "position [5kbp window]", title = sample_names[i]) +
    theme(axis.text.x = element_blank())
  
##################################### ALL ######################################
  
  pg <- plot_grid(gh, gl, hl, g, h, ncol = 5)
  
  p <- plot_grid(pg, c, ncol = 1)

  plist[[i]] <- p
  
}


pdf(outfile, width = 30)

for (i in 1:nr_samples) {
  
  print(plist[[i]])
  
}

dev.off()

