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



cov.select <- merge(cov.select.00, cov.select.02, by = c("factor", "x"))
cov.select <- merge(cov.select, cov.select.04, by = c("factor", "x"))
cov.select <- merge(cov.select, snp.select, by = c("factor"))

cov.select <- cov.select[-6]
colnames(cov.select) <- c("factor", "x", "cov00", "cov02", "cov04", "hetDiff")

cov.select <- cov.select[order(cov.select$x),]


cr <- colorRamp(c('blue','green','red'), space = "rgb")


pdf(file=scatter2D_out, width = 18, height = 9)
par(mfrow=c(1,3), mar=c(4,4,4,1), oma=c(1,1,0,0), xpd=TRUE)

plot(log(cov.select$cov00), cov.select$hetDiff, col = rgb(cr(cov.select$x / max(cov.select$x))/255), pch = 20, 
     main = "Statistics for each scaffold, nm = 0", xlab = "Mean genome coverage [log-scale 10^x]",
     ylab = "Mean difference in heterozygosity")
lgd_ = rep(NA, 11)
lgd_[c(1,11)] = c(min(cov.select$x),max(cov.select$x))
legend("bottomright", legend = lgd_, fill = colorRampPalette(colors = c('blue','green','red'))(11),
       border = NA, y.intersp = 0.8, title = "Scaffold length")

plot(log(cov.select$cov02), cov.select$hetDiff, col = rgb(cr(cov.select$x / max(cov.select$x))/255), pch = 20, 
     main = "Statistics for each scaffold, nm = 2", xlab = "Mean genome coverage [log-scale 10^x]",
     ylab = "Mean difference in heterozygosity")

plot(log(cov.select$cov04), cov.select$hetDiff, col = rgb(cr(cov.select$x / max(cov.select$x))/255), pch = 20, 
     main = "Statistics for each scaffold, nm = 4", xlab = "Mean genome coverage [log-scale 10^x]",
     ylab = "Mean difference in heterozygosity")

dev.off()



pdf(file=scatter3D_out, width = 18, height = 9)
par(mfrow=c(1,3), mar=c(4,4,4,1), oma=c(1,1,0,0), xpd=TRUE)
# 3d plot het/cov/len
scatterplot3d(log(cov.select$x), log(cov.select$cov00), cov.select$hetDiff,  color = rgb(cr(cov.select$x / max(cov.select$x))/255), 
              pch = 20, main = "Statistics for each scaffold, nm = 0", angle = 150,
              xlab = "Scaffold length [log-scale 10^x bp]", ylab = "Mean normalized genome coverage [log-scale 10^x]",
              zlab = "Mean difference in heterozygosity")

legend("topright", legend = lgd_, fill = colorRampPalette(colors = c('blue','green','red'))(11),
       border = NA, y.intersp = 0.8, title = "Scaffold length", bg = "white", cex = 1.5)

scatterplot3d(log(cov.select$x), log(cov.select$cov02), cov.select$hetDiff, color = rgb(cr(cov.select$x / max(cov.select$x))/255), 
              pch = 20, main = "Statistics for each scaffold, nm = 2", angle = 150,
              xlab = "Scaffold length [log-scale 10^x bp]", ylab = "Mean normalized genome coverage [log-scale 10^x]", 
              zlab = "Mean difference in heterozygosity")

scatterplot3d(log(cov.select$x), log(cov.select$cov04), cov.select$hetDiff, color = rgb(cr(cov.select$x / max(cov.select$x))/255), 
              pch = 20, main = "Statistics for each scaffold, nm = 4", angle = 150,
              xlab = "Scaffold length [log-scale 10^x bp]", ylab = "Mean normalized genome coverage [log-scale 10^x]", 
              zlab = "Mean difference in heterozygosity")

dev.off()

