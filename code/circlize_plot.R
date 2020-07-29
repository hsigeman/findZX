library(circlize)
library(doBy)
library(data.table)

# Reads genome coverage files and snp file and makes a circlized plot
# and a 2D scatterplot between gencov and difference in heterozygosity

source("code/functions.R")

########## CALCULATION ##########

set.seed(999)

args <- commandArgs(trailingOnly = TRUE)

file00 = args[1]
file02 = args[2]
file04 = args[3]
filesnp = args[4]

circlize_out = args[5]
scatter_out = args[6]
chr_file = args[7]


cov_00_table <- read.table(file00, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)

if (dim(cov_00_table)[1] == 0) {
  print("Warning: No chromsome/scaffold over 1Mbp! Check output based on chromosome instead.")
  
  pdf(file=scatter_out, width = 9, height = 9)
  dev.off()
  
  pdf(file=circlize_out, width = 9, height = 9)
  dev.off()
} else {
  cov_00 <- gen_data_4plotting(cov_00_table, c("chr", "range", "ratio"))
  cov.select.00 <- cov_00$df
  max.cov.00 <- cov_00$max
  min.cov.00 <- cov_00$min
  median.cov.00 <- cov_00$median

  cov_02_table <- read.table(file02, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
  cov_02 <- gen_data_4plotting(cov_02_table, c("chr", "range", "ratio"))
  cov.select.02 <- cov_02$df
  max.cov.02 <- cov_02$max
  min.cov.02 <- cov_02$min
  median.cov.02 <- cov_02$median

  cov_04_table <- read.table(file04, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
  cov_04 <- gen_data_4plotting(cov_04_table, c("chr", "range", "ratio"))
cov.select.04 <- cov_04$df
max.cov.04 <- cov_04$max
min.cov.04 <- cov_04$min
median.cov.04 <- cov_04$median

snp_table <- read.table(filesnp, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
snp <- gen_data_4plotting(snp_table, c("chr", "range", "diff"))
snp.select <- snp$df
max.snp <- snp$max
min.snp <- snp$min
median.snp <- snp$median

########## SCATTER PLOTTING ##########

cov.select <- merge(cov.select.00, cov.select.02, by = c("factor", "x"))
cov.select <- merge(cov.select, cov.select.04, by = c("factor", "x"))
cov.select <- merge(cov.select, snp.select, by = c("factor", "x"))

colnames(cov.select) <- c("factor", "x", "cov00", "cov02", "cov04", "hetDiff")


if (file.exists(chr_file)) {
  chromosome <- read.csv(chr_file, header = FALSE, sep = ",")
} else {
  max_per_chr <- setDT(cov.select.00)[, .SD[which.max(x)], by=factor]
  max_per_chr <- as.data.frame(max_per_chr[,1:2])
  chromosome <- as.data.frame(matrix((as.character(max_per_chr[order(-max_per_chr$x),][,1])), nrow = 1))
}
  
    
cov.select$factor <- ordered(cov.select$factor, chromosome)
cov.select.00$factor <- ordered(cov.select.00$factor, levels = chromosome)
cov.select.02$factor <- ordered(cov.select.02$factor, levels = chromosome)
cov.select.04$factor <- ordered(cov.select.04$factor, levels = chromosome)
snp.select$factor <- ordered(snp.select$factor, levels = chromosome)


# Color each chromosome with a unique color
c <- colorRampPalette(colors = c('red','yellow','green','blue'))(length(unique(cov.select$factor)))
palette(c)

pdf(file=scatter_out, width = 9, height = 9)

par(mfrow=c(3,1), mar=c(4,4,4,1), oma=c(0,10,0,0), xpd=TRUE)

plot(cov.select$cov00, cov.select$hetDiff, col = cov.select$factor, xlim = c(min.cov.00, max.cov.00),
     xlab = "Normalized genome coverage, nm = 0", ylab = "Difference in proportion of heterozygosity",
     pch = 20)
abline(h = median.snp)
abline(v = median.cov.00)
plot(cov.select$cov02, cov.select$hetDiff, col = cov.select$factor, xlim = c(min.cov.00, max.cov.00),
     xlab = "Normalized genome coverage, nm = 2", ylab = "Difference in proportion of heterozygosity",
     pch = 20)
abline(h = median.snp)
abline(v = median.cov.02)
par(xpd=NA)
legend("left", legend=levels(unique(cov.select$factor)),fill=1:length(cov.select$factor),inset=c(-0.2,0.6))
par(xpd=TRUE)
plot(cov.select$cov04, cov.select$hetDiff, col = cov.select$factor, xlim = c(min.cov.00, max.cov.00),
     xlab = "Normalized genome coverage, nm = 4", ylab = "Difference in proportion of heterozygosity",
     pch = 20)
abline(h = median.snp)
abline(v = median.cov.04)

dev.off()


########## CIRCLIZE PLOTTING ##########

factor.nr <- as.numeric(length(unique(cov.select.04$factor)))

circos.clear()
pdf(file=circlize_out, width = 9, height = 9)
par(mfrow=c(1,1), mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7) 
circos.par(cell.padding = c(0, 0, 0, 0), "track.height" = 0.15, gap.after = c(rep(1, factor.nr-1), 10))
circos.initialize(factors = cov.select.04$factor, x = cov.select.04$x)


circos.trackPlotRegion(factors = snp.select$factor, ylim = c(min.snp,max.snp), x = snp.select$x, y = snp.select$y, panel.fun = function(x, y) {
  
  grey = c("#FFFFFF", "#CCCCCC", "#999999")
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim) + uy(9, "mm"), cex = 1.1, sector.index) 
  circos.lines(x, y, col = "blue", area = TRUE, baseline = median.snp)
  circos.yaxis(side = "left", sector.index = levels(cov.select.00$factor)[1])
})

circos.trackPlotRegion(factors = cov.select.00$factor, ylim = c(min.cov.00,max.cov.00), x = cov.select.00$x, y = cov.select.00$y, panel.fun = function(x, y) {
  
  grey = c("#FFFFFF", "#CCCCCC", "#999999")
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.lines(x, y, col = "red", area = TRUE, baseline = median.cov.00)
  circos.yaxis(side = "left", sector.index = levels(cov.select.00$factor)[1])
})

circos.trackPlotRegion(factors = cov.select.02$factor, ylim = c(min.cov.02,max.cov.02), x = cov.select.02$x, y = cov.select.02$y, panel.fun = function(x, y) {
  
  grey = c("#FFFFFF", "#CCCCCC", "#999999")
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.lines(x, y, col = "red", area = TRUE, baseline = median.cov.02)
  circos.yaxis(side = "left", sector.index = levels(cov.select.00$factor)[1])
})

circos.trackPlotRegion(factors = cov.select.04$factor, ylim = c(min.cov.04,max.cov.04), x = cov.select.04$x, y = cov.select.04$y, panel.fun = function(x, y) {
  
  grey = c("#FFFFFF", "#CCCCCC", "#999999")
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.lines(x, y, col = "red", area = TRUE, baseline = median.cov.04)
  circos.yaxis(side = "left", sector.index = levels(cov.select.00$factor)[1])
  circos.axis("bottom", direction = "inside", labels.facing = "reverse.clockwise")
})

dev.off()

}
