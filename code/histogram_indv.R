library(doBy)
library(data.table)
library(ggplot2)
library(cowplot)
require(scales)

# Reads in a bed-file with coverage values and a file with heterozygosity for 
# each site for each sample, values for each sample are in separate columns.
# Removes outliers and makes two histograms and three heatmaps for each individual:
# histograms for genome coverage and heterozygosity and heatmaps of genome
# coverage vs heterozygosity, genome coverage vs length, heterozygosity vs length.
#
# Heterozygosity is given for each variable site. 1 = heterozygot, 0 = homozygot,
# na = missing data.

# Outputs a coverage file where the windows with a coverage in the range of N1
# is set to 1 and the other windows set to 0

set.seed(999)

args <- commandArgs(trailingOnly = TRUE)

file_gencov = args[1]
file_snp = args[2]
outPdf = args[3]
outCov = args[4]
synteny = args[5]
chr_file = args[6]
sample_names = args[7:length(args)]

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

nr_samples <- length(sample_names)
col_names <- colnames(cov)[4:length(cov)]
f <- as.formula( paste( paste( col_names, 
                               collapse = "+"), 
                        "~", "chr + range"))

max_per_chr <- setDT(cov)[, .SD[which.max(end)], by=chr]
len_chr <- as.data.frame( max_per_chr[,1:2])
colnames(len_chr) <- c("chr","length")

Tcov <- transform(cov, 
                  range=floor(start/5000))
Tcov <- summaryBy(f, 
                  data=Tcov, 
                  keep.names=TRUE, 
                  na.rm = TRUE)
colnames(Tcov) <- c("chr","range",sample_names)

Thet <- transform(het, 
                  range=floor(end/5000))
Thet <- summaryBy(f, 
                  data=Thet, 
                  keep.names=TRUE, 
                  na.rm = TRUE)
colnames(Thet) <- c("chr","range",sample_names)

cov_het <- merge(Tcov, 
                 Thet, 
                 by = c("chr","range"))
cov_het <- merge(cov_het, 
                 len_chr, 
                 by = "chr")
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

simple_cov <- cov_het[,1:(2+nr_samples)]

for (i in 1:nr_samples) {

############## HEATMAP HETEROZYGOSITY VS GENOME COVERAGE VS LENGTH #############
  
  x = i + 2
  y = i + 2 + nr_samples
  
  gh <- ggplot(data = cov_het, 
               aes(x = cov_het[,x], 
                   y = cov_het[,y])) + 
    geom_bin2d(na.rm = TRUE) + 
    labs(x = "genome coverage", 
         y = "heterozygosity", 
         title = sample_names[i]) + 
    theme_bw() +
    scale_fill_gradient(low="white",
                        high="darkblue",
                        trans="log10")
  
  gl <- ggplot(data = cov_het, 
               aes(x = cov_het[,x], 
                   y = length)) + 
    geom_bin2d(na.rm = TRUE) + 
    labs(x = "genome coverage", 
         y = "scaffold length", 
         title = sample_names[i]) + 
    theme_bw() +
    scale_fill_gradient(low="white",
                        high="darkblue",
                        trans="log10")
  
  hl <- ggplot(data = cov_het, 
               aes(x = cov_het[,y], 
                   y = length)) + 
    geom_bin2d(na.rm = TRUE) + 
    labs(x = "heterozygosity", 
         y = "scaffold length", 
         title = sample_names[i]) + 
    theme_bw() +
    scale_fill_gradient(low="white",
                        high="darkblue",
                        trans="log10")
  
############## HISTOGRAM GENOME COVERAGE WITH LINES FOR N1 AND N2 ##############
  
  g <- ggplot(cov_het, 
              aes(x = cov_het[,x])) + 
     geom_histogram(bins = 50, 
                    na.rm = TRUE) + 
     labs(x="genome coverage", 
          y="Frequency", 
          title = sample_names[i]) + 
     theme_bw()
  
  # Find the x value for the bin with the highest count
  hist_stats_hg <- ggplot_build(g)$data[[1]]
  max_bin <- which(hist_stats_hg$density == max(hist_stats_hg$density))
  
  # Calcualte half of the x value for the highest bin
  max_x <- hist_stats_hg$x[max_bin]
  halfMax_x <- max_x / 2
  
  g <- g + geom_vline(xintercept = halfMax_x, 
                      color = "blue") + 
    geom_vline(xintercept = max_x, 
               color = "red")
  
  # All windows with a coverage in the N1 range are set to 1, all other 0
  window_range <- 0.25
  
  simple_cov[,x][simple_cov[,x] < halfMax_x*(1-window_range)] <- 0
  simple_cov[,x][simple_cov[,x] > halfMax_x*(1+window_range)] <- 0
  simple_cov[,x][simple_cov[,x] > 0] <- 1
  
########################### HISTOGRAM HETEROZYGOSITY ###########################
  
  h <- ggplot(cov_het, 
              aes(x = cov_het[,y])) + 
     geom_histogram(bins = 50, 
                    na.rm = TRUE) + 
     labs(x="heterozygosity", 
          y="Frequency", 
          title = sample_names[i]) + 
     theme_bw()
  
################################# CHROMOSOMES ##################################
  
  # Pick out the 50 largest chromosomes/scaffolds or all if there are not more 
  # than 50 chromosomes
  if (dim(len_chr)[1] >= 50) {
    nr_chr <- 50
  } else {
    nr_chr <- dim(len_chr)[1]
  }
  
  top_chr <- as.factor(len_chr[with(len_chr, order(-length)),][1:nr_chr,1])
  
  cov_het_subset <- cov_het[cov_het$chr %in% top_chr,]
  cov_het_subset$chr <- factor(cov_het_subset$chr, levels = top_chr)
  
  c <- ggplot(cov_het_subset, 
              aes(x=range, 
                  y=cov_het_subset[,x])) + 
    geom_smooth(na.rm = TRUE) + 
    facet_grid(. ~ chr, 
               scales = "free_x", 
               space = "free_x") +
    labs(y = "genome coverage", 
         x = "position [5kbp window]", 
         title = sample_names[i]) +
    theme(axis.text.x = element_blank())
  
##################################### ALL ######################################
  
  pg <- plot_grid(gh, gl, hl, g, h, ncol = 5)
  
  p <- plot_grid(pg, c, ncol = 1)

  plist[[i]] <- p
  
}


pdf(outPdf, width = 30)

for (i in 1:nr_samples) {
  
  print(plist[[i]])
  
}
dev.off()

# Write coverage table where windows around N1 is 1 and all other windows are 0
write.table(simple_cov, outCov, 
            quote=FALSE, 
            sep="\t", 
            row.names = F, col.names = F, 
            na = "NA")
