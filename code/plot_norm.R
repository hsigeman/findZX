library(circlize)
library(doBy)
library(data.table)

set.seed(999)

args <- commandArgs(trailingOnly = TRUE)

species= args[1]
file00= args[2]
file02= args[3]
file04= args[4]
filesnp= args[5]

cov_00 =read.table(file00, header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov_00$V1 <- as.factor(cov_00$V1)

cov_00$V1 <- ordered(cov_00$V1,
                  levels = c("1", "1A", "1B", "2", "3", "4", "4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
                             "15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "Z"))

myvars <- c("V1", "V2", "V3")
setDF(cov_00)
cov.select.00 <- cov_00[myvars]
cov.select.00 <- na.omit(cov.select.00)
cov.select.00 <- plyr::rename(cov.select.00, c("V1"="factor", "V2"="x", "V3"="y"))

max.cov.00 <- max(cov.select.00$y)
max.cov.00
min.cov.00 <- min(cov.select.00$y)
min.cov.00

cov_chr_sum.00 <- summaryBy(y ~ factor, data=cov.select.00, keep.names=TRUE)
cov_chr_sum.00$type <- "cov.00"



cov_02 =read.table(file02, header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov_02$V1 <- as.factor(cov_02$V1)

cov_02$V1 <- ordered(cov_02$V1,
                  levels = c("1", "1A", "1B", "2", "3", "4", "4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
                             "15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "Z"))

myvars <- c("V1", "V2", "V3")
setDF(cov_02)
cov.select.02 <- cov_02[myvars]
cov.select.02 <- na.omit(cov.select.02)
cov.select.02 <- plyr::rename(cov.select.02, c("V1"="factor", "V2"="x", "V3"="y"))

max.cov.02 <- max(cov.select.02$y)
max.cov.02
min.cov.02 <- min(cov.select.02$y)
min.cov.02

cov_chr_sum.02 <- summaryBy(y ~ factor, data=cov.select.02, keep.names=TRUE)
cov_chr_sum.02$type <- "cov.02"



cov_04 =read.table(file04, header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov_04$V1 <- as.factor(cov_04$V1)

cov_04$V1 <- ordered(cov_04$V1,
                  levels = c("1", "1A", "1B", "2", "3", "4", "4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
                             "15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "Z"))

myvars <- c("V1", "V2", "V3")
setDF(cov_04)
cov.select.04 <- cov_04[myvars]
cov.select.04 <- na.omit(cov.select.04)
cov.select.04 <- plyr::rename(cov.select.04, c("V1"="factor", "V2"="x", "V3"="y"))

max.cov.04 <- max(cov.select.04$y)
max.cov.04
min.cov.04 <- min(cov.select.04$y)
min.cov.04

cov_chr_sum.04 <- summaryBy(y ~ factor, data=cov.select.04, keep.names=TRUE)
cov_chr_sum.04$type <- "cov.04"




snp =read.table(filesnp, header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
snp$V1 <- as.factor(snp$V1)

snp$V1 <- ordered(snp$V1,
                  levels = c("1", "1A", "1B", "2", "3", "4", "4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
                             "15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "Z"))

myvars <- c("V1", "V2", "V3")
setDF(snp)
snp.select <- snp[myvars]
snp.select <- na.omit(snp.select)
snp.select <- plyr::rename(snp.select, c("V1"="factor", "V2"="x", "V3"="y"))

max.snp <- max(snp.select$y)
max.snp
min.snp <- min(snp.select$y)
min.snp

snp_chr_sum <- summaryBy(y ~ factor, data=snp.select, keep.names=TRUE)
snp_chr_sum$type <- "snv"



######### PLOTTING

factor.nr <- as.numeric(length(unique(cov.select.04$factor)))

circos.clear()
pdf(file=args[6], width = 9, height = 9)
par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7) 
circos.par(cell.padding = c(0, 0, 0, 0), "track.height" = 0.15, gap.after = c(rep(1, factor.nr-1), 10))
circos.initialize(factors = cov.select.04$factor, x = cov.select.04$x)


circos.trackPlotRegion(factors = snp.select$factor, ylim = c(min.snp,max.snp), x = snp.select$x, y = snp.select$y, panel.fun = function(x, y) {
  
  grey = c("#FFFFFF", "#CCCCCC", "#999999")
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim) + uy(9, "mm"), cex = 1.1, sector.index) 
  circos.lines(x, y, col = "blue", area = TRUE, baseline = 1)
  circos.yaxis(side = "left", sector.index = 1)
  # circos.axis("top")
  
})
circos.trackPlotRegion(factors = cov.select.00$factor, ylim = c(min.cov.00,max.cov.00), x = cov.select.00$x, y = cov.select.00$y, panel.fun = function(x, y) {
  
  grey = c("#FFFFFF", "#CCCCCC", "#999999")
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.lines(x, y, col = "red", area = TRUE, baseline = 1)
  circos.yaxis(side = "left", sector.index = 1)
  # circos.axis("bottom", direction = "inside", labels.facing = "reverse.clockwise")
})

circos.trackPlotRegion(factors = cov.select.02$factor, ylim = c(min.cov.02,max.cov.02), x = cov.select.02$x, y = cov.select.02$y, panel.fun = function(x, y) {
  
  grey = c("#FFFFFF", "#CCCCCC", "#999999")
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.lines(x, y, col = "red", area = TRUE, baseline = 1)
  circos.yaxis(side = "left", sector.index = 1)
  # circos.axis("bottom", direction = "inside", labels.facing = "reverse.clockwise")
})

circos.trackPlotRegion(factors = cov.select.04$factor, ylim = c(min.cov.04,max.cov.04), x = cov.select.04$x, y = cov.select.04$y, panel.fun = function(x, y) {
  
  grey = c("#FFFFFF", "#CCCCCC", "#999999")
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.lines(x, y, col = "red", area = TRUE, baseline = 1)
  circos.yaxis(side = "left", sector.index = 1)
  circos.axis("bottom", direction = "inside", labels.facing = "reverse.clockwise")
})
dev.off()

