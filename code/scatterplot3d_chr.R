library(doBy)
library(data.table)
library(plot3D)
library(ggplot2)
require(scales)
library(cowplot)

# reads in chr.out files and plots cov vs het and colors chr dots dependen on
# chr length
# 3d scatterplot cov vs het vs chr

source("code/functions.R")

set.seed(999)

args <- commandArgs(trailingOnly = TRUE)

file00 = args[1]
file02 = args[2]
file04 = args[3]
filesnp = args[4]
scatter2D_out = args[5]
scatter3D_out = args[6]

################################################################################
################################# READ FILES ###################################
################################################################################

cov_00_table <- read.table(file00, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
cov_00 <- gen_data_4plotting(cov_00_table, c("chr", "length", "ratio"))
cov.select.00 <- cov_00$df
max.cov.00 <- cov_00$max
min.cov.00 <- cov_00$min
median.cov.00 <- cov_00$median

cov_02_table <- read.table(file02, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
cov_02 <- gen_data_4plotting(cov_02_table, c("chr", "length", "ratio"))
cov.select.02 <- cov_02$df
max.cov.02 <- cov_02$max
min.cov.02 <- cov_02$min
median.cov.02 <- cov_02$median

cov_04_table <- read.table(file04, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
cov_04 <- gen_data_4plotting(cov_04_table, c("chr", "length", "ratio"))
cov.select.04 <- cov_04$df
max.cov.04 <- cov_04$max
min.cov.04 <- cov_04$min
median.cov.04 <- cov_04$median

snp_table <- read.table(filesnp, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
snp <- gen_data_4plotting(snp_table, c("chr", "length", "diff"))
snp.select <- snp$df
max.snp <- snp$max
min.snp <- snp$min
median.snp <- snp$median

################################################################################
################################# MERGE DATA ###################################
################################################################################

cov.select <- merge(cov.select.00, cov.select.02, by = c("factor", "x"))
cov.select <- merge(cov.select, cov.select.04, by = c("factor", "x"))
cov.select <- merge(cov.select, snp.select, by = c("factor"))

cov.select <- cov.select[,-6]
colnames(cov.select) <- c("scaffold", "length", "cov00", "cov02", "cov04", "hetDiff")

cov.select <- cov.select[order(cov.select$length),]

################################################################################
############################### SCATTER PLOT 2D ################################
################################################################################

c0 <- ggplot(cov.select, aes(x = cov00, y = hetDiff)) + 
  labs(title = "nm = 0", x = "", 
       y = "mean difference in heterozygosity") + theme_bw() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) + 
  annotation_logticks(sides="b") + theme_bw() + geom_point() + 
  geom_point(aes(color = length)) + scale_color_gradient(low = "blue", high = "red") + 
  theme(legend.position="none")

c2 <- ggplot(cov.select, aes(x = cov02, y = hetDiff)) + 
  labs(title = "nm = 2", x = "mean genome coverage", 
       y = "") + theme_bw() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) + 
  annotation_logticks(sides="b") + theme_bw() + geom_point() + 
  geom_point(aes(color = length)) + scale_color_gradient(low = "blue", high = "red") + 
  theme(legend.position="none")

c4 <- ggplot(cov.select, aes(x = cov04, y = hetDiff)) + 
  labs(title = "nm = 4", x = "", 
       y = "") + theme_bw() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) + 
  annotation_logticks(sides="b") + theme_bw() + geom_point() + 
  geom_point(aes(color = length)) + scale_color_gradient(low = "blue", high = "red")
legend <- get_legend(c4)
c4 <- c4 + theme(legend.position="none")

pg <- plot_grid(c0, c2, c4, legend, ncol = 4, rel_widths = c(3,3,3,1), 
                labels = c("A", "B", "C", ""))

ggsave(scatter2D_out, plot = pg, device = pdf(), width = 14, height = 9)

################################################################################
############################### SCATTER PLOT 3D ################################
################################################################################

pdf(file=scatter3D_out, width = 14, height = 5)

par(mfrow=c(1,3), mar=c(2,1,2,0), oma=c(0,0,0,0), xpd=TRUE)

scatter3D(cov.select$length, log(cov.select$cov00), cov.select$hetDiff,
          colvar = cov.select$length, pch = 19, xlab = "Scaffold length", 
          ylab = "log of normalized genome coverage", main = "nm = 0",
          zlab = "difference in heterozygosity", bty = "b2",
          colkey = FALSE)

scatter3D(cov.select$length, log(cov.select$cov02), cov.select$hetDiff,
          colvar = cov.select$length, pch = 19, xlab = "Scaffold length", 
          ylab = "log of normalized genome coverage", main = "nm = 2",
          zlab = "difference in heterozygosity", bty = "b2",
          colkey = FALSE)

scatter3D(cov.select$length, log(cov.select$cov04), cov.select$hetDiff,
          colvar = cov.select$length, pch = 19, xlab = "Scaffold length", 
          ylab = "log of normalized genome coverage", main = "nm = 4",
          zlab = "difference in heterozygosity", bty = "b2",
          colkey = FALSE)

dev.off()

