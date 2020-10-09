library(doBy)
library(data.table)
library(plot3D)
library(ggplot2)
require(scales)
library(cowplot)
library(viridis)

# Reads in files ending with '.chr.out' produced by calculate_gencov_windows.R 
# and calculate_heterozygosityDiff_windows.R
# Makes scatterplots for genome coverage vs difference in heterozygosity and
# chromosomes/scaffolds length. Two values are plotted against each other and
# the third is shown through coloring of the points. 
# Makes a 3D scatterplot of all three variables and colors after length.

set.seed(999)
source("code/functions.R")

args <- commandArgs(trailingOnly = TRUE)

file1 = args[1]
file2 = args[2]
file3 = args[3]
filesnp = args[4]
scatter2D_out = args[5]
scatter3D_out = args[6]
chr_file = args[7]

################################################################################
################################# READ FILES ###################################
################################################################################

cov_00_table <- read.table(file1, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
cov_02_table <- read.table(file2, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
cov_04_table <- read.table(file3, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
snp_table <- read.table(filesnp, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)

if (file.exists(chr_file)) {
  
  chromosome <- read.csv(chr_file, header = FALSE, sep = ",")
  chromosome <- as.factor(chromosome)
  
  cov_00_table <- cov_00_table[cov_00_table$chr %in% chromosome, ]
  cov_02_table <- cov_02_table[cov_02_table$chr %in% chromosome, ]
  cov_04_table <- cov_04_table[cov_04_table$chr %in% chromosome, ]
  snp_table <- snp_table[snp_table$chr %in% chromosome, ]
  
}

cov_00 <- gen_data_4plotting(cov_00_table, c("chr", "length", "diff"))
cov.select.00 <- cov_00$df
max.cov.00 <- cov_00$max
min.cov.00 <- cov_00$min
median.cov.00 <- cov_00$median

cov_02 <- gen_data_4plotting(cov_02_table, c("chr", "length", "diff"))
cov.select.02 <- cov_02$df
max.cov.02 <- cov_02$max
min.cov.02 <- cov_02$min
median.cov.02 <- cov_02$median

cov_04 <- gen_data_4plotting(cov_04_table, c("chr", "length", "diff"))
cov.select.04 <- cov_04$df
max.cov.04 <- cov_04$max
min.cov.04 <- cov_04$min
median.cov.04 <- cov_04$median

snp <- gen_data_4plotting(snp_table, c("chr", "length", "diff"))
snp.select <- snp$df
max.snp <- snp$max
min.snp <- snp$min
median.snp <- snp$median

################################################################################
################################# MERGE DATA ###################################
################################################################################

cov.select <- merge(cov.select.00, cov.select.02, by = c("factor"))[c(1,3,5)]
colnames(cov.select) <- c("factor", "cov00", "cov02")

cov.select <- merge(cov.select, cov.select.04, by = c("factor"))[c(1,2,3,5)]
colnames(cov.select) <- c("factor", "cov00", "cov02", "cov04")

cov.select <- merge(cov.select, snp.select, by = c("factor"))
colnames(cov.select) <- c("Chromosome", "cov00", "cov02", "cov04", "length", "hetDiff")

################################################################################
############################# SCATTER PLOT LENGTH ##############################
################################################################################

cov.select <- cov.select[order(cov.select$length),]

l0 <- ggplot(cov.select, aes(x = cov00, y = hetDiff, size=length)) + 
  labs(title = file1, x = "", y = "difference in heterozygosity") + 
  theme_bw() + 
  geom_point(aes(color = length)) +
  theme(legend.position="none") +
  scale_color_gradient(low="lightgrey", high="red")

l2 <- ggplot(cov.select, aes(x = cov02, y = hetDiff, size=length)) + 
  labs(title = file2, x = "difference in normalized genome coverage", y = "") + 
  theme_bw() + 
  geom_point(aes(color = length)) + 
  theme(legend.position="none") +
  scale_color_gradient(low="lightgrey", high="red")

l4 <- ggplot(cov.select, aes(x = cov04, y = hetDiff, size=length)) + 
  labs(title = file3, x = "", color = "scaffold length", y = "") + 
  theme_bw() + 
  geom_point(aes(color = length)) + 
  scale_color_gradient(low="lightgrey", high="red")

legend <- get_legend(l4)

l4 <- l4 + theme(legend.position="none")

l <- plot_grid(l0, l2, l4, legend, ncol = 4, rel_widths = c(3,3,3,1))

################################################################################
############################ SCATTER PLOT COVERAGE #############################
################################################################################

cov.select <- cov.select[order(cov.select$cov00),]

c0 <- ggplot(cov.select, aes(y = length, x = hetDiff)) + 
  labs(title = "", y = "scaffold length", x = "", color = "coverage") + 
  theme_bw() + 
  geom_point(aes(color = cov00)) + 
  scale_color_viridis()
  
legend <- get_legend(c0)
  
c0 <- c0 + theme(legend.position="none")

cov.select <- cov.select[order(cov.select$cov02),]

c2 <- ggplot(cov.select, aes(y = length, x = hetDiff)) + 
  labs(title = "", y = "", x = "difference in heterozygosity") + 
  theme_bw() + 
  geom_point(aes(color = cov02)) +
  scale_color_viridis() +
  theme(legend.position="none")

cov.select <- cov.select[order(cov.select$cov04),]

c4 <- ggplot(cov.select, aes(y = length, x = hetDiff)) + 
  labs(title = "", x = "", y = "") + 
  theme_bw() + 
  geom_point(aes(color = cov04)) + 
  scale_color_viridis() + 
  theme(legend.position="none")

c <- plot_grid(c0, c2, c4, legend, ncol = 4, rel_widths = c(3,3,3,1))

################################################################################
############################ SCATTER PLOT HETDIFF ##############################
################################################################################

cov.select <- cov.select[order(cov.select$hetDiff),]
hetDiff_low <- cov.select[cov.select$hetDiff < 0,]

mid <- 0

hl0 <- ggplot(hetDiff_low, aes(y = length, x = cov00)) + 
  labs(title = "", y = "scaffold length", x = "") + 
  theme_bw() + 
  geom_point(aes(color = hetDiff)) + 
  scale_color_gradient2(midpoint = mid, mid = "lightgrey",low="blue", high="red") + 
  theme(legend.position="none")

hl2 <- ggplot(hetDiff_low, aes(y = length, x = cov02)) + 
  labs(title = "", y = "", x = "difference in normalized genome coverage") + 
  theme_bw() + 
  geom_point(aes(color = hetDiff)) + 
  scale_color_gradient2(midpoint = mid, mid = "lightgrey",low="blue", high="red") + 
  theme(legend.position="none")

hl4 <- ggplot(hetDiff_low, aes(y = length, x = cov04)) + 
  labs(title = "", x = "", color = "heterozygosity", y = "") + 
  theme_bw() + 
  geom_point(aes(color = hetDiff)) + 
  scale_color_gradient2(midpoint = mid, mid = "lightgrey",low="blue", high="red")

legend <- get_legend(hl4)

hl4 <- hl4 + theme(legend.position="none")

hl <- plot_grid(hl0, hl2, hl4, legend, ncol = 4, rel_widths = c(3,3,3,1))

################################################################################

hetDiff_high <- cov.select[cov.select$hetDiff >= 0,]

hh0 <- ggplot(hetDiff_high, aes(y = length, x = cov00)) + 
  labs(title = "", y = "scaffold length", x = "") + 
  theme_bw() + 
  geom_point(aes(color = hetDiff)) + 
  scale_color_gradient2(midpoint = mid, mid = "lightgrey",low="blue", high="red") + 
  theme(legend.position="none")

hh2 <- ggplot(hetDiff_high, aes(y = length, x = cov02)) + 
  labs(title = "", y = "", x = "difference in normalized genome coverage") + 
  theme_bw() + 
  geom_point(aes(color = hetDiff)) + 
  scale_color_gradient2(midpoint = mid, mid = "lightgrey",low="blue", high="red") + 
  theme(legend.position="none")

hh4 <- ggplot(hetDiff_high, aes(y = length, x = cov04)) + 
  labs(title = "", x = "", color = "heterozygosity", y = "") + 
  theme_bw() + 
  geom_point(aes(color = hetDiff)) + 
  scale_color_gradient2(midpoint = mid, mid = "lightgrey",low="blue", high="red")

legend <- get_legend(hh4)

hh4 <- hh4 + theme(legend.position="none")

hh <- plot_grid(hh0, hh2, hh4, legend, ncol = 4, rel_widths = c(3,3,3,1))

################################################################################

pg <- plot_grid(l,c,hl,hh, ncol = 1, labels = 'AUTO')

ggsave(scatter2D_out, plot = pg, device = pdf(), width = 14, height = 14)

################################################################################
############################### SCATTER PLOT 3D ################################
################################################################################

pdf(file=scatter3D_out, width = 14, height = 8)

par(mfrow=c(1,3), mar=c(2,1,2,0), oma=c(0,0,0,0), xpd=TRUE)

scatter3D(cov.select$length, cov.select$cov00, cov.select$hetDiff,
          colvar = cov.select$length, pch = 19, xlab = "Scaffold length", 
          ylab = "difference in normalized genome coverage", main = "nm = 0",
          zlab = "difference in heterozygosity", bty = "b2",
          colkey = FALSE, col = viridis(length(cov.select$length), direction = -1),
          phi = 20, theta = 60)

scatter3D(cov.select$length, cov.select$cov02, cov.select$hetDiff,
          colvar = cov.select$length, pch = 19, xlab = "Scaffold length", 
          ylab = "difference in normalized genome coverage", main = "nm = 2",
          zlab = "difference in heterozygosity", bty = "b2",
          colkey = FALSE, col = viridis(length(cov.select$length), direction = -1),
          phi = 20, theta = 60)

scatter3D(cov.select$length, cov.select$cov04, cov.select$hetDiff,
          colvar = cov.select$length, pch = 19, xlab = "Scaffold length", 
          ylab = "difference in normalized genome coverage", main = "nm = 4",
          zlab = "difference in heterozygosity", bty = "b2",
          colkey = FALSE, col = viridis(length(cov.select$length), direction = -1),
          phi = 20, theta = 60)

dev.off()

