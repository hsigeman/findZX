library(circlize)
library(doBy)
library(data.table)

# Reads genome coverage files and snp file and makes a circlized plot
# Hard coded, if you change something with the input file, use with caution
# Zebra finch chromosomes

########## FUNCTIONS ##########

gen_data_4plotting <- function(filename, y_column) {
	data_table =read.table(filename, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
	data_table$chr <- as.factor(data_table$chr)

	#data_table$chr <- ordered(data_table$chr,
	#			levels = c("1", "1A", "1B", "2", "3", "4", "4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
	#				"15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "Z"))

	myvars <- c("chr", "range", y_column)
	setDF(data_table)
	data_table.select <- data_table[myvars]
	data_table.select <- na.omit(data_table.select)
	colnames(data_table.select) <- c("factor", "x", "y")

	max.data_table <- max(data_table.select$y)
	min.data_table <- min(data_table.select$y)
	median.data_table <- median(data_table.select$y)

	return(list(df = data_table.select, max = max.data_table, min = min.data_table, median = median.data_table))
}


########## CALCULATION ##########

set.seed(999)

args <- commandArgs(trailingOnly = TRUE)

file00 = args[1]
file02 = args[2]
file04 = args[3]
filesnp = args[4]


cov_00 <- gen_data_4plotting(file00, "ratio")
cov.select.00 <- cov_00$df
max.cov.00 <- cov_00$max
max.cov.00
min.cov.00 <- cov_00$min
min.cov.00
median.cov.00 <- cov_00$median
median.cov.00

cov_02 <- gen_data_4plotting(file02, "ratio")
cov.select.02 <- cov_02$df
max.cov.02 <- cov_02$max
max.cov.02
min.cov.02 <- cov_02$min
min.cov.02
median.cov.02 <- cov_02$median
median.cov.02

cov_04 <- gen_data_4plotting(file04, "ratio")
cov.select.04 <- cov_04$df
max.cov.04 <- cov_04$max
max.cov.04
min.cov.04 <- cov_04$min
min.cov.04
median.cov.04 <- cov_04$median
median.cov.04

snp <- gen_data_4plotting(filesnp, "diff")
snp.select <- snp$df
max.snp <- snp$max
max.snp
min.snp <- snp$min
min.snp
median.snp <- snp$median
median.snp

########## PLOTTING ##########

factor.nr <- as.numeric(length(unique(cov.select.04$factor)))

circos.clear()
pdf(file=args[5], width = 9, height = 9)
par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7) 
circos.par(cell.padding = c(0, 0, 0, 0), "track.height" = 0.15, gap.after = c(rep(1, factor.nr-1), 10))
circos.initialize(factors = cov.select.04$factor, x = cov.select.04$x)


circos.trackPlotRegion(factors = snp.select$factor, ylim = c(min.snp,max.snp), x = snp.select$x, y = snp.select$y, panel.fun = function(x, y) {
  
  grey = c("#FFFFFF", "#CCCCCC", "#999999")
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim) + uy(9, "mm"), cex = 1.1, sector.index) 
  circos.lines(x, y, col = "blue", area = TRUE, baseline = median.snp)
  circos.yaxis(side = "left", sector.index = 1)
})

circos.trackPlotRegion(factors = cov.select.00$factor, ylim = c(min.cov.00,max.cov.00), x = cov.select.00$x, y = cov.select.00$y, panel.fun = function(x, y) {
  
  grey = c("#FFFFFF", "#CCCCCC", "#999999")
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.lines(x, y, col = "red", area = TRUE, baseline = median.cov.00)
  circos.yaxis(side = "left", sector.index = 1)
})

circos.trackPlotRegion(factors = cov.select.02$factor, ylim = c(min.cov.02,max.cov.02), x = cov.select.02$x, y = cov.select.02$y, panel.fun = function(x, y) {
  
  grey = c("#FFFFFF", "#CCCCCC", "#999999")
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.lines(x, y, col = "red", area = TRUE, baseline = median.cov.02)
  circos.yaxis(side = "left", sector.index = 1)
})

circos.trackPlotRegion(factors = cov.select.04$factor, ylim = c(min.cov.04,max.cov.04), x = cov.select.04$x, y = cov.select.04$y, panel.fun = function(x, y) {
  
  grey = c("#FFFFFF", "#CCCCCC", "#999999")
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.lines(x, y, col = "red", area = TRUE, baseline = median.cov.04)
  circos.yaxis(side = "left", sector.index = 1)
  circos.axis("bottom", direction = "inside", labels.facing = "reverse.clockwise")
})

dev.off()

