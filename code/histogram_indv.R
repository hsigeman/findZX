library(doBy)
library(data.table)
library(ggplot2)
library(cowplot)
require(scales)
library(ggExtra)

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
file_read_len = args[3]

outPdf = args[4]
synteny = args[5]
chr_file = args[6]

sample_names = args[7:length(args)]

data <- ggdraw() + draw_label(paste0("Data points from tables: \n", file_gencov, " \n ",
   file_snp))


################################################################################
################################# READ FILES ###################################
################################################################################

cov <- read.table(file_gencov, header=FALSE, fill=TRUE, stringsAsFactor=FALSE)
het <- read.table(file_snp, header=FALSE, fill=TRUE, stringsAsFactor=FALSE)
read_len <- read.csv(file_read_len, header = FALSE)

colnames(cov)[1:3] <- c("chr","start","end")
colnames(het)[1:3] <- c("chr","start","end")

nr_samples <- length(sample_names)

for (i in 1:nr_samples) {
  sample_read_len <- read_len[read_len == sample_names[i],][1,2]
  
  cov[,(i+3)] <- cov[,(i+3)]*sample_read_len/5000
}

################################################################################
################################ CALCULATIONS ##################################
################################################################################

col_names <- colnames(cov)[4:length(cov)]

max_per_chr <- setDT(cov)[, .SD[which.max(end)], by=chr]
len_chr <- as.data.frame( max_per_chr[,1:2])
colnames(len_chr) <- c("chr","length")

cov_het <- merge(cov, 
                 het, 
                 by = c("chr", "start", "end"))
cov_het <- merge(cov_het, 
                 len_chr, 
                 by = "chr")
cov_het <- as.data.frame(cov_het)
cov <- as.data.frame(cov)

for (i in 1:nr_samples) {
  
  outliers <- boxplot(cov_het[,(i+3)], plot = FALSE)$stats[5]
  cov_het[,(i+3)][cov_het[,(i+3)] > outliers] = NA
  
  outliers <- boxplot(cov[,(i+3)], plot = FALSE)$stats[5]
  cov[,(i+3)][cov[,(i+3)] > outliers] = NA
  
}

if (file.exists(chr_file)) {
  
#  chromosome <- read.csv(chr_file, header = FALSE, sep = ",")
#  chromosome <- as.factor(chromosome)
  chromosome <- read.delim(chr_file, header = FALSE)
  chromosome$V1 <- trimws(chromosome$V1, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  chromosome <- as.factor(chromosome$V1)
  cov_het <- cov_het[cov_het$chr %in% chromosome, ]
  
}

################################################################################
################################### PLOTTING ###################################
################################################################################

text_size_colour = list(theme_bw(base_family="Courier", base_size = 12) + 
    theme(axis.text.x= element_text(colour="black", size=12)) +
    theme(axis.text.y= element_text(colour="black", size=12)) +
    theme(axis.title.x = element_text(colour="black",size=14)) + 
    theme(axis.title.y = element_text(colour="black",size=14)) + 
    theme(axis.ticks = element_line(colour = "black", size = 0.5)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    theme(plot.title=element_text(family="Courier", size=14, colour="black", hjust = 0.5)) +
    theme(panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                colour = "lightgrey"), 
          panel.grid.minor = element_line(size = 0.15, linetype = 'solid',
                                colour = "lightgrey")))


plist <- list()

for (i in 1:nr_samples) {

############## HEATMAP HETEROZYGOSITY VS GENOME COVERAGE VS LENGTH #############
  
  x = i + 3
  y = i + 3 + nr_samples
  
  gh <- ggplot(data = cov_het, 
               aes(x = cov_het[,x], 
                   y = cov_het[,y])) + 
    geom_bin2d(na.rm = TRUE) + 
    labs(x = "genome coverage", 
         y = "heterozygosity") + 
    text_size_colour +
    scale_fill_gradient(low="white",
                        high="black",
                        trans="log10") + 
    theme(legend.position = "bottom")
  
  gh <- ggMarginal(gh, type="histogram", color="black", fill="blue")

  gl <- ggplot(data = cov_het, 
               aes(x = cov_het[,x], 
                   y = length/1000000)) + 
    geom_bin2d(na.rm = TRUE) + 
    labs(x = "genome coverage", 
         y = "scaffold length (Mb)") + 
    text_size_colour +
    scale_fill_gradient(low="white",
                        high="black",
                        trans="log10") + 
    theme(legend.position = "bottom")
  
  gl <- ggMarginal(gl, type="histogram", color="black", fill="blue")
  hl <- ggplot(data = cov_het, 
               aes(x = cov_het[,y], 
                   y = length/1000000)) + 
    geom_bin2d(na.rm = TRUE) + 
    labs(x = "heterozygosity", 
         y = "scaffold length (Mb)") + 
    text_size_colour +
    scale_fill_gradient(low="white",
                        high="black",
                        trans="log10")  + 
    theme(legend.position = "bottom")

  hl <- ggMarginal(hl, type="histogram", color="black", fill="blue") 

############## HISTOGRAM GENOME COVERAGE WITH LINES FOR N1 AND N2 ##############
  
  g <- ggplot(cov_het, 
              aes(x = cov_het[,x])) + 
     geom_histogram(bins = 50, 
                    na.rm = TRUE, color="black", size=0.1, fill="blue") + 
     labs(x="genome coverage", 
          y="Frequency") + 
     text_size_colour 
  
  # Find the x value for the bin with the highest count
  hist_stats_hg <- ggplot_build(g)$data[[1]]
  max_bin <- which(hist_stats_hg$density == max(hist_stats_hg$density))
  
  # Calcualte half of the x value for the highest bin
  max_x <- hist_stats_hg$x[max_bin]
  halfMax_x <- max_x / 2
  
  #g <- g + geom_vline(xintercept = halfMax_x, 
  #                    color = "blue") + 
  #  geom_vline(xintercept = max_x, 
  #             color = "red")
  
########################### HISTOGRAM HETEROZYGOSITY ###########################
  
  h <- ggplot(cov_het, 
              aes(x = cov_het[,y])) + 
     geom_histogram(bins = 50, 
                    na.rm = TRUE, color="black", size=0.1, fill="blue") + 
     labs(x="heterozygosity", 
          y="Frequency") + 
     text_size_colour
  
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
              aes(x=end/1000000, 
                  y=cov_het_subset[,x])) + 
    geom_bin2d(na.rm = TRUE) + 
    scale_fill_gradient(low="white",
                        high="black",
                        trans="log10") + 
  #  geom_point(alpha=0.1, size=0.2) +
    geom_smooth(na.rm = TRUE, span = 0.3) + 
    facet_grid(. ~ chr, 
               scales = "free_x", 
               space = "free_x") +
    labs(y = "genome coverage", 
         x = "Chromosome/scaffold position (Mb)") +
    scale_x_continuous(limits = c(0, NA)) +
    theme(axis.text.x = element_blank()) +
    text_size_colour
  
    title <- ggdraw() +
    draw_label(paste0("\"Confirm sexing\" plots based on 5kb window values \t \t \t Sample: ", sample_names[i]),
    fontface = 'bold',
    x = 0,
    hjust = 0
    ) +
    theme(
    plot.margin = margin(0, 0, 0, 7))
 
##################################### ALL ######################################
  
  pg <- plot_grid(g,h, gh, gl, hl, 
                  ncol = 5,
                  labels = c("A","B","C","D","E"))
  
  p <- plot_grid(pg, c, ncol = 1, labels = c(" ","F"))
  
  p2 <- plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1))
  
  plist[[i]] <- p2
  
}


pdf(outPdf, width = 20, height = 10)

for (i in 1:nr_samples) {
  
  print(plist[[i]]) 
}
print(data)
dev.off()
