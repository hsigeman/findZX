library(doBy)
library(data.table)
library(plot3D)
library(ggplot2)
require(scales)
library(cowplot)
library(viridisLite)
library(ggpubr)

# Reads in files ending with '.chr.out' produced by calculate_gencov_windows.R 
# and calculate_heterozygosityDiff_windows.R
# Makes scatterplots for genome coverage vs difference in heterozygosity and
# chromosomes/scaffolds length. Two values are plotted against each other and
# the third is shown through coloring of the points. 
# Makes a 3D scatterplot of all three variables and colors after length.

set.seed(999)
source("workflow/scripts/functions.R")

args <- commandArgs(trailingOnly = TRUE)

file1 = args[1]
file2 = args[2]
file3 = args[3]
filesnp = args[4]
scatter2D_out = args[5]
#scatter3D_out = args[6]
chr_file = args[6]
ED1 = args[7]
ED2 = args[8]
ED3 = args[9]

################################################################################
################################# READ FILES ###################################
################################################################################

ED1 = gsub("\\.", "-", ED1)
ED2 = gsub("\\.", "-", ED2)
ED3 = gsub("\\.", "-", ED3)
ED1 = gsub("$", " mismatches", ED1)
ED2 = gsub("$", " mismatches", ED2)

ED1 = gsub("0-0", "0", ED1)
ED1 = gsub("0-", "<= ", ED1)
ED2 = gsub("0-", "<= ", ED2)


scatter2D_out_base = gsub("\\.pdf", "", scatter2D_out)

cov_00_table <- read.table(file1, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
cov_02_table <- read.table(file2, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
cov_04_table <- read.table(file3, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
snp_table <- read.table(filesnp, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)

if (file.exists(chr_file)) {
  
#  chromosome <- read.csv(chr_file, header = FALSE, sep = ",")
#  chromosome <- as.factor(chromosome)
  chromosome <- read.delim(chr_file, header = FALSE)
  chromosome$V1 <- trimws(chromosome$V1, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  chromosome <- chromosome$V1 
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

x_axis <- expression(Delta*~genome~coverage)
y_axis <- expression(Delta*~heterozygosity)

scaleFUN <- function(x) sprintf("%.2f", x)

text_size_colour = list(theme_bw(base_family="Courier", base_size = 12) + 
    theme(axis.text.x= element_text(colour="black", size=12)) +
    theme(axis.text.y= element_text(colour="black", size=12)) +
    theme(axis.title.x = element_text(colour="black",size=15)) + 
    theme(axis.title.y = element_text(colour="black",size=15)) + 
    theme(plot.title=element_text(family="Courier", size=15, colour="black", hjust = 0.5)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    theme(plot.margin= margin(1, 1, 1, 1, "mm")) +
    theme(panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                colour = "lightgrey"), 
          panel.grid.minor = element_line(size = 0.15, linetype = 'solid',
                                colour = "lightgrey")))

theme_title <- function(...) {
  theme_gray(base_family = "Courier") + 
    theme(plot.title = element_text(face = "bold"))
}
title_theme <- calc_element("plot.title", theme_title())

theme_description <- function(...) {
  theme_gray(base_family = "Courier") + 
    theme(plot.title = element_text(face = "plain",size = 12, hjust = 0))
}
theme_description <- calc_element("plot.title", theme_description())



cov.select <- cov.select[order(abs(cov.select$length)),]



l0 <- ggplot(cov.select, aes(x = cov00, y = hetDiff, size=length/1000000)) + 
  labs(title = sprintf("%s", ED1), x = "", y = y_axis) + 
  geom_hline(yintercept = round(median.snp, digits=2), linetype="dashed", size=0.2) +
  geom_vline(xintercept = median.cov.00, linetype="dashed", size=0.2) +
  geom_point(aes(color = round(length/1000000, digits = 2)), alpha = 0.5) +
  text_size_colour +
  theme(legend.position="none") +
  scale_color_gradient(low="gray90", high="blue") #+
 # scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
 # scale_x_continuous(labels = scales::number_format(accuracy = 0.01))

l2 <- ggplot(cov.select, aes(x = cov02, y = hetDiff, size=length/1000000)) + 
  #labs(title = sprintf("%s mismatches", ED2), x = "difference in normalized genome coverage", y = "") + 
  labs(title = sprintf("%s", ED2), x = " ", y = "") + 
  geom_hline(yintercept = median.snp, linetype="dashed", size=0.2) +
  geom_vline(xintercept = median.cov.02, linetype="dashed", size=0.2) +
  geom_point(aes(color = round(length/1000000, digits = 2)), alpha = 0.5) + 
  text_size_colour +
  theme(legend.position="none") +
  scale_color_gradient(low="gray90", high="blue") #+
 # scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
 # scale_x_continuous(labels = scales::number_format(accuracy = 0.01))

l4 <- ggplot(cov.select, aes(x = cov04, y = hetDiff, size=length/1000000)) + 
  labs(title = sprintf("%s", ED3), x = "", size = "scaffold length (Mb)", color = "scaffold length (Mb)", y = "") + 
  geom_hline(yintercept = median.snp, linetype="dashed", size=0.2) +
  geom_vline(xintercept = median.cov.04, linetype="dashed", size=0.2) + 
  geom_point(aes(color = round(length/1000000, digits = 2)), alpha = 0.5) + 
  text_size_colour +
  scale_color_gradient(low="gray90", high="blue") #+
 # scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
 # scale_x_continuous(labels = scales::number_format(accuracy = 0.01))

legend_l <- get_legend(l4)

l4 <- l4 + theme(legend.position="none")

#l <- plot_grid(l0, l2, l4, legend_l, ncol = 4, rel_widths = c(3,3,3,2))


################################################################################
############################ SCATTER PLOT HETDIFF ##############################
################################################################################


cov.select <- cov.select[order(abs(cov.select$hetDiff)),]

mid <- 0

hl0 <- ggplot(cov.select, aes(y = length/1000000, x = cov00)) + 
  labs(title = "", y = "scaffold length (Mb)", x = "") + 
  geom_vline(xintercept = median.cov.00, linetype="dashed", size=0.2) + 
  geom_point(aes(color = hetDiff)) + 
  scale_color_gradient2(midpoint = mid, mid = "lightgrey",low="blue", high="red") + 
  text_size_colour +
  theme(legend.position="none") #+
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  #scale_x_continuous(labels = scales::number_format(accuracy = 0.01))

hl2 <- ggplot(cov.select, aes(y = length/1000000, x = cov02)) + 
  labs(title = "", y = "", x = x_axis) + 
  geom_vline(xintercept = median.cov.02, linetype="dashed", size=0.2) + 
  geom_point(aes(color = hetDiff)) + 
  scale_color_gradient2(midpoint = mid, mid = "lightgrey",low="blue", high="red") + 
  text_size_colour +
  theme(legend.position="none") #+
#  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
 # scale_x_continuous(labels = scales::number_format(accuracy = 0.01))

hl4 <- ggplot(cov.select, aes(y = length/1000000, x = cov04)) + 
  labs(title = "", x = "", color = "heterozygosity", y = "") + 
  geom_vline(xintercept = median.cov.04, linetype="dashed", size=0.2) + 
  geom_point(aes(color = hetDiff)) + 
  scale_color_gradient2(midpoint = mid, mid = "lightgrey",low="blue", high="red") + 
  text_size_colour #+
 # scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
 # scale_x_continuous(labels = scales::number_format(accuracy = 0.01))

legend_hl <- get_legend(hl4)

hl4 <- hl4 + theme(legend.position="none")

#hl <- ggarrange(hl0, hl2, hl4, ncol = 3, align="hv")
#hl <- plot_grid(hl, legend_hl, ncol = 2, rel_widths = c(9,2))

pg <- ggarrange(l0, l2, l4, hl0, hl2, hl4, ncol = 3, nrow= 2, align="hv",
  labels =c("A", "B", "C", "D","E", "F"))

pg_legend <- plot_grid(legend_l, legend_hl, ncol = 1)

pg <- plot_grid(pg, pg_legend, ncol = 2, rel_widths = c(1, 0.25))


################################################################################

#pg <- plot_grid(pg,c,hl, ncol = 1, labels = 'AUTO')

#ggsave(scatter2D_out, plot = pg, device = pdf(), width = 14, height = 14)

#pg <- plot_grid(l,hl, ncol = 1, labels = 'AUTO')
title <- ggdraw() + draw_label("Sex differences (heterogametic-homogametic) per chromosome/scaffold", 
    fontfamily = title_theme$family,
    fontface = title_theme$face,
    size = title_theme$size
  )


pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1.5))

data <- ggdraw() + draw_label(paste0("Data points from tables: \n", file1, " \n ",
   file2, " \n ", file3, " \n ", filesnp))


  pdf(file=scatter2D_out, width = 14, height = 8)
  print(pg)
  print(data)
  dev.off()
  

#ggsave(scatter2D_out, plot = pg, device = pdf(), width = 14, height = 8)

#ggsave(sprintf("%s.png", scatter2D_out_base), plot = pg, device = png(), width = 14, height = 8, dpi = 900)


################################################################################
############################### SCATTER PLOT 3D ################################
################################################################################

#pdf(file=scatter3D_out, width = 14, height = 8)

#par(mfrow=c(1,3), mar=c(2,1,2,0), oma=c(0,0,0,0), xpd=TRUE)

#scatter3D(cov.select$length, cov.select$cov00, cov.select$hetDiff,
#          colvar = cov.select$length, pch = 19, xlab = "Scaffold length", 
#          ylab = "difference in normalized genome coverage", main = sprintf("%s mismatches", ED1),
#          zlab = "difference in heterozygosity", bty = "b2",
#          colkey = FALSE, col = viridis(length(cov.select$length), direction = -1),
#          phi = 20, theta = 60)

#scatter3D(cov.select$length, cov.select$cov02, cov.select$hetDiff,
#          colvar = cov.select$length, pch = 19, xlab = "Scaffold length", 
#          ylab = "difference in normalized genome coverage", main = sprintf("%s mismatches", ED2),
#          zlab = "difference in heterozygosity", bty = "b2",
#          colkey = FALSE, col = viridis(length(cov.select$length), direction = -1),
#          phi = 20, theta = 60)

#scatter3D(cov.select$length, cov.select$cov04, cov.select$hetDiff,
#          colvar = cov.select$length, pch = 19, xlab = "Scaffold length", 
#          ylab = "difference in normalized genome coverage", main = sprintf("%s mismatches", ED3),
#          zlab = "difference in heterozygosity", bty = "b2",
#          colkey = FALSE, col = viridis(length(cov.select$length), direction = -1),
#          phi = 20, theta = 60)
#dev.off()
