library(doBy)
library(data.table)
library(scatterplot3d)

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


cov_00 <- gen_data_4plotting(file00, c("chr", "length", "ratio"))
cov.select.00 <- cov_00$df
max.cov.00 <- cov_00$max
max.cov.00
min.cov.00 <- cov_00$min
min.cov.00
median.cov.00 <- cov_00$median
median.cov.00

cov_02 <- gen_data_4plotting(file02, c("chr", "length", "ratio"))
cov.select.02 <- cov_02$df
max.cov.02 <- cov_02$max
max.cov.02
min.cov.02 <- cov_02$min
min.cov.02
median.cov.02 <- cov_02$median
median.cov.02

cov_04 <- gen_data_4plotting(file04, c("chr", "length", "ratio"))
cov.select.04 <- cov_04$df
max.cov.04 <- cov_04$max
max.cov.04
min.cov.04 <- cov_04$min
min.cov.04
median.cov.04 <- cov_04$median
median.cov.04

snp <- gen_data_4plotting(filesnp, c("chr", "length", "diff"))
snp.select <- snp$df
max.snp <- snp$max
max.snp
min.snp <- snp$min
min.snp
median.snp <- snp$median
median.snp


cov.select <- merge(cov.select.00, cov.select.02, by = c("factor", "x"))
cov.select <- merge(cov.select, cov.select.04, by = c("factor", "x"))
cov.select <- merge(cov.select, snp.select, by = c("factor"))

cov.select <- cov.select[-6]
colnames(cov.select) <- c("factor", "x", "cov00", "cov02", "cov04", "hetDiff")

cov.select <- cov.select[order(cov.select$x),]


c <- colorRampPalette(colors = c('blue','green','yellow', 'red'))(length(cov.select$x))



pdf(file=scatter2D_out, width = 18, height = 9)
par(mfrow=c(1,3), mar=c(4,4,4,1), oma=c(1,1,0,0), xpd=TRUE)

plot(cov.select$cov00, cov.select$hetDiff, col = c, pch = 20, 
     main = "Statistics for each scaffold, nm = 0", xlab = "Mean genome coverage",
     ylab = "Mean difference in heterozygosity")
lgd_ = rep(NA, 11)
lgd_[c(1,11)] = c(min(cov.select$x),max(cov.select$x))
legend("bottomright",
       legend = lgd_,
       fill = colorRampPalette(colors = c('blue','green','yellow', 'red'))(11),
       border = NA,
       y.intersp = 0.8,
       title = "Scaffold length")

plot(cov.select$cov02, cov.select$hetDiff, col = c, pch = 20, 
     main = "Statistics for each scaffold, nm = 2", xlab = "Mean genome coverage",
     ylab = "Mean difference in heterozygosity")

plot(cov.select$cov04, cov.select$hetDiff, col = c, pch = 20, 
     main = "Statistics for each scaffold, nm = 4", xlab = "Mean genome coverage",
     ylab = "Mean difference in heterozygosity")

dev.off()



pdf(file=scatter3D_out, width = 18, height = 9)
par(mfrow=c(1,3), mar=c(4,4,4,1), oma=c(1,1,0,0), xpd=TRUE)
# 3d plot het/cov/len
scatterplot3d(cov.select$x, cov.select$cov00, cov.select$hetDiff,  color = c, 
              pch = 20, angle = 150, main = "Statistics for each scaffold, nm = 0", 
              xlab = "Scaffold length", ylab = "Mean genome coverage", 
              zlab = "Mean difference in heterozygosity")
legend("topright",
       legend = lgd_,
       fill = colorRampPalette(colors = c('blue','green','yellow', 'red'))(11),
       border = NA,
       y.intersp = 0.8,
       title = "Scaffold length")

scatterplot3d(cov.select$x, cov.select$cov02, cov.select$hetDiff, color = c, 
              pch = 20, angle = 150, main = "Statistics for each scaffold, nm = 2", 
              xlab = "Scaffold length", ylab = "Mean genome coverage", 
              zlab = "Mean difference in heterozygosity")

scatterplot3d(cov.select$x, cov.select$cov04, cov.select$hetDiff, color = c, 
              pch = 20, angle = 150, main = "Statistics for each scaffold, nm = 4", 
              xlab = "Scaffold length", ylab = "Mean genome coverage", 
              zlab = "Mean difference in heterozygosity")

dev.off()

