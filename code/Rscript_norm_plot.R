#install.packages("circlize")
library(circlize)
library(doBy)
library(data.table)

set.seed(999)

args <- commandArgs(trailingOnly = TRUE)

species=args[1]
female=args[2]
male=args[3]


cov=read.table(args[4],header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov <- plyr::rename(cov, c("V1"="contig", "V2"="contig_start","V3"="contig_end", "V4"="chr", "V7"="F_cov", "V8"="M_cov", "V5"="start", "V6"="end"))

# Remove chromosomes smaller than 1 Mb #
max_per_chr <- setDT(cov)[, .SD[which.max(end)], by=chr]
chr_over_1mb <- subset(max_per_chr$chr, max_per_chr$end>1000000)
cov <- cov[ cov$chr %in% chr_over_1mb, ]

cov <- subset(cov, cov$chr!="Un")
random <- unique(cov$chr[grep("random", cov$chr)])
cov <- cov[ ! cov$chr %in% random, ]
random <- unique(cov$chr[grep("LG", cov$chr)])
cov <- cov[ ! cov$chr %in% random, ]
random <- unique(cov$chr[grep("MT", cov$chr)])
cov <- cov[ ! cov$chr %in% random, ]

cov$chr <- ordered(cov$chr,
                   levels = c("1", "1A", "1B", "2", "3", "4", "4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
                              "15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "Z"))

cov$cov_ratio <- cov$F_cov/cov$M_cov
cov$cov_ratio_scaled <- cov$cov_ratio / median(cov$cov_ratio, na.rm = TRUE)
cov <- subset(cov, cov$cov_ratio != "Inf")
cov <- subset(cov, cov$cov_ratio_scaled != "Inf")

# Remove outliers
outliers <- boxplot(cov_ratio_scaled ~ chr, data=cov)$out
cov <- cov[-which(cov$cov_ratio_scaled %in% outliers),]

cov_1Mb_ranges <- transform(cov, range=round(end/1000000))
cov_1Mb_ranges_sum <- summaryBy(cov_ratio_scaled ~ chr + range, data=cov_1Mb_ranges, keep.names=TRUE)
cov_100kb_ranges <- transform(cov, range=round(end/100000))
cov_100kb_ranges_sum <- summaryBy(cov_ratio_scaled ~ chr + range, data=cov_100kb_ranges, keep.names=TRUE)
cov <- cov_1Mb_ranges_sum
cov$chr <- as.factor(cov$chr)

myvars <- c("chr", "range", "cov_ratio_scaled")
setDF(cov)
cov.select <- cov[myvars]
cov.select <- na.omit(cov.select)
cov.select <- plyr::rename(cov.select, c("chr"="factor", "range"="x", "cov_ratio_scaled"="y"))

max.cov.00 <- max(cov.select$y)
max.cov.00
min.cov.00 <- min(cov.select$y)
min.cov.00

cov.select.00 <- cov.select

cov_chr_sum.00 <- summaryBy(y ~ factor, data=cov.select.00, keep.names=TRUE)
cov_chr_sum.00$type <- "cov.00"

outfile <- sprintf(paste("results/",species,"/",species,".1Mbp_cov_scaled.nm.00.txt",sep=""))
write.table(cov_1Mb_ranges_sum, outfile, quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")

outfile <- sprintf(paste("results/",species,"/",species,".100kbp_cov_scaled.nm.00.txt",sep=""))
write.table(cov_100kb_ranges_sum, outfile, quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")


#####################

cov=read.table(args[5],header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov <- plyr::rename(cov, c("V1"="contig", "V2"="contig_start","V3"="contig_end", "V4"="chr", "V7"="F_cov", "V8"="M_cov", "V5"="start", "V6"="end"))

# Remove chromosomes smaller than 1 Mb #
max_per_chr <- setDT(cov)[, .SD[which.max(end)], by=chr]
chr_over_1mb <- subset(max_per_chr$chr, max_per_chr$end>1000000)
cov <- cov[ cov$chr %in% chr_over_1mb, ]

cov <- subset(cov, cov$chr!="Un")
random <- unique(cov$chr[grep("random", cov$chr)])
cov <- cov[ ! cov$chr %in% random, ]
random <- unique(cov$chr[grep("LG", cov$chr)])
cov <- cov[ ! cov$chr %in% random, ]
random <- unique(cov$chr[grep("MT", cov$chr)])
cov <- cov[ ! cov$chr %in% random, ]

cov$chr <- ordered(cov$chr,
                   levels = c("1", "1A", "1B", "2", "3", "4", "4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
                              "15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "Z"))

cov$cov_ratio <- cov$F_cov/cov$M_cov
cov$cov_ratio_scaled <- cov$cov_ratio / median(cov$cov_ratio, na.rm = TRUE)
cov <- subset(cov, cov$cov_ratio != "Inf")
cov <- subset(cov, cov$cov_ratio_scaled != "Inf")

# Remove outliers
outliers <- boxplot(cov_ratio_scaled ~ chr, data=cov)$out
cov <- cov[-which(cov$cov_ratio_scaled %in% outliers),]

cov_1Mb_ranges <- transform(cov, range=round(end/1000000))
cov_1Mb_ranges_sum <- summaryBy(cov_ratio_scaled ~ chr + range, data=cov_1Mb_ranges, keep.names=TRUE)
cov_100kb_ranges <- transform(cov, range=round(end/100000))
cov_100kb_ranges_sum <- summaryBy(cov_ratio_scaled ~ chr + range, data=cov_100kb_ranges, keep.names=TRUE)
cov <- cov_1Mb_ranges_sum
cov$chr <- as.factor(cov$chr)

myvars <- c("chr", "range", "cov_ratio_scaled")
setDF(cov)
cov.select <- cov[myvars]
cov.select <- na.omit(cov.select)
cov.select <- plyr::rename(cov.select, c("chr"="factor", "range"="x", "cov_ratio_scaled"="y"))

max.cov.02 <- max(cov.select$y)
max.cov.02
min.cov.02 <- min(cov.select$y)
min.cov.02
cov.select.02 <- cov.select

cov_chr_sum.02 <- summaryBy(y ~ factor, data=cov.select.02, keep.names=TRUE)
cov_chr_sum.02$type <- "cov.02"

outfile <- sprintf(paste("results/",species,"/",species,".1Mbp_cov_scaled.nm.02.txt",sep=""))
write.table(cov_1Mb_ranges_sum, outfile, quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")

outfile <- sprintf(paste("results/",species,"/",species,".100kbp_cov_scaled.nm.02.txt",sep=""))
write.table(cov_100kb_ranges_sum, outfile, quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")


#####################

cov=read.table(args[6],header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov <- plyr::rename(cov, c("V1"="contig", "V2"="contig_start","V3"="contig_end", "V4"="chr", "V7"="F_cov", "V8"="M_cov", "V5"="start", "V6"="end"))

# Remove chromosomes smaller than 1 Mb #
max_per_chr <- setDT(cov)[, .SD[which.max(end)], by=chr]
chr_over_1mb <- subset(max_per_chr$chr, max_per_chr$end>1000000)
cov <- cov[ cov$chr %in% chr_over_1mb, ]

cov <- subset(cov, cov$chr!="Un")
random <- unique(cov$chr[grep("random", cov$chr)])
cov <- cov[ ! cov$chr %in% random, ]
random <- unique(cov$chr[grep("LG", cov$chr)])
cov <- cov[ ! cov$chr %in% random, ]
random <- unique(cov$chr[grep("MT", cov$chr)])
cov <- cov[ ! cov$chr %in% random, ]

cov$chr <- ordered(cov$chr,
                   levels = c("1", "1A", "1B", "2", "3", "4", "4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
                              "15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "Z"))

cov$cov_ratio <- cov$F_cov/cov$M_cov
cov$cov_ratio_scaled <- cov$cov_ratio / median(cov$cov_ratio, na.rm = TRUE)
cov <- subset(cov, cov$cov_ratio != "Inf")
cov <- subset(cov, cov$cov_ratio_scaled != "Inf")

# Remove outliers
outliers <- boxplot(cov_ratio_scaled ~ chr, data=cov)$out
cov <- cov[-which(cov$cov_ratio_scaled %in% outliers),]

cov_1Mb_ranges <- transform(cov, range=round(end/1000000))
cov_1Mb_ranges_sum <- summaryBy(cov_ratio_scaled ~ chr + range, data=cov_1Mb_ranges, keep.names=TRUE)
cov_100kb_ranges <- transform(cov, range=round(end/100000))
cov_100kb_ranges_sum <- summaryBy(cov_ratio_scaled ~ chr + range, data=cov_100kb_ranges, keep.names=TRUE)
cov <- cov_1Mb_ranges_sum
cov$chr <- as.factor(cov$chr)

myvars <- c("chr", "range", "cov_ratio_scaled")
setDF(cov)
cov.select <- cov[myvars]
cov.select <- na.omit(cov.select)
cov.select <- plyr::rename(cov.select, c("chr"="factor", "range"="x", "cov_ratio_scaled"="y"))

max.cov.04 <- max(cov.select$y)
max.cov.04
min.cov.04 <- min(cov.select$y)
min.cov.04

cov.select.04 <- cov.select

cov_chr_sum.04 <- summaryBy(y ~ factor, data=cov.select.04, keep.names=TRUE)
cov_chr_sum.04$type <- "cov.04"

outfile <- sprintf(paste("results/",species,"/",species,".1Mbp_cov_scaled.nm.04.txt",sep=""))
write.table(cov_1Mb_ranges_sum, outfile, quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")

outfile <- sprintf(paste("results/",species,"/",species,".100kbp_cov_scaled.nm.04.txt",sep=""))
write.table(cov_100kb_ranges_sum, outfile, quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")


###############


###### SNP 

snp=read.table(args[7],header=FALSE,fill=TRUE,stringsAsFactor=FALSE)


max_per_chr <- setDT(snp)[, .SD[which.max(V6)], by=V4]
chr_over_1mb <- subset(max_per_chr$V4, max_per_chr$V6>1000000)
snp <- snp[ snp$V4 %in% chr_over_1mb, ]
snp <- subset(snp, snp$V4!="Un")
random <- unique(snp$V4[grep("random", snp$V4)])
snp <- snp[ ! snp$V4 %in% random, ]
random <- unique(snp$V4[grep("LG", snp$V4)])
snp <- snp[ ! snp$V4 %in% random, ]
random <- unique(snp$V4[grep("MT", snp$V4)])
snp <- snp[ ! snp$V4 %in% random, ]

snp$V4 <- ordered(snp$V4,
                  levels = c("1", "1A", "1B", "2", "3", "4", "4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
                             "15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "Z"))

snp$count <- 1
snp <- subset(snp, snp$V1=="S")

### Count private female and male alleles per chromosome and Mb
snp_1Mb_ranges <- transform(snp, range=round(V5/1000000))
snp_1Mb_ranges_count <- aggregate(count ~ V4 + V3 + range, data = snp_1Mb_ranges[which(snp_1Mb_ranges$V1=="S"),], sum)
### Transform to wide format
snp_1Mb_ranges_count_wide <- dcast(snp_1Mb_ranges_count, V4 + range ~ V3, value.var="count")
snp_1Mb_ranges_count_wide[is.na(snp_1Mb_ranges_count_wide)] <- 0
names(snp_1Mb_ranges_count_wide)[names(snp_1Mb_ranges_count_wide)==female] <- "female"
names(snp_1Mb_ranges_count_wide)[names(snp_1Mb_ranges_count_wide)==male] <- "male"
snp_1Mb_ranges_count_wide$snp_diff <- snp_1Mb_ranges_count_wide$female - snp_1Mb_ranges_count_wide$male
snp_1Mb_ranges_count_wide$snp_ratio <- snp_1Mb_ranges_count_wide$female/snp_1Mb_ranges_count_wide$male
names(snp_1Mb_ranges_count_wide)[names(snp_1Mb_ranges_count_wide)=="V4"] <- "chr"
head(snp_1Mb_ranges_count_wide)
snp_1Mb_ranges_count_wide <- subset(snp_1Mb_ranges_count_wide, snp_1Mb_ranges_count_wide$snp_ratio != "Inf")
snp_1Mb_ranges_count_wide$snp_ratio_scaled <- snp_1Mb_ranges_count_wide$snp_ratio / median(snp_1Mb_ranges_count_wide$snp_ratio, na.rm = TRUE)
snp_1Mb_ranges_count_wide$snp_diff_scaled <- snp_1Mb_ranges_count_wide$snp_diff / median(snp_1Mb_ranges_count_wide$snp_diff, na.rm = TRUE)

#myvars <- c("chr", "range", "snp_ratio_scaled")
myvars <- c("chr", "range", "snp_diff")
snp.select <- snp_1Mb_ranges_count_wide[myvars]
snp.select <- na.omit(snp.select)
snp.select <- plyr::rename(snp.select, c("chr"="factor", "range"="x", "snp_diff"="y"))

head(snp.select)


max.snp <- max(snp.select$y)
max.snp
min.snp <- min(snp.select$y)
min.snp

snp_chr_sum <- summaryBy(y ~ factor, data=snp.select, keep.names=TRUE)
snp_chr_sum$type <- "snv"

outfile <- sprintf(paste("results/",species,"/",species,".1Mbp_snp_diff.txt",sep=""))
write.table(snp_1Mb_ranges_count_wide, outfile, quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")

snp_100kbp_ranges_count <- transform(snp, range=round(V5/100000))
snp_100kbp_ranges_count <- aggregate(count ~ V4 + V3 + range, data = snp_100kbp_ranges_count[which(snp_100kbp_ranges_count$V1=="S"),], sum)
snp_100kbp_ranges_count_wide <- dcast(snp_100kbp_ranges_count, V4 + range ~ V3, value.var="count")
snp_100kbp_ranges_count_wide[is.na(snp_100kbp_ranges_count_wide)] <- 0
names(snp_100kbp_ranges_count_wide)[names(snp_100kbp_ranges_count_wide)==female] <- "female"
names(snp_100kbp_ranges_count_wide)[names(snp_100kbp_ranges_count_wide)==male] <- "male"
snp_100kbp_ranges_count_wide$snp_diff <- snp_100kbp_ranges_count_wide$female - snp_100kbp_ranges_count_wide$male
snp_100kbp_ranges_count_wide$snp_ratio <- snp_100kbp_ranges_count_wide$female/snp_100kbp_ranges_count_wide$male
names(snp_100kbp_ranges_count_wide)[names(snp_100kbp_ranges_count_wide)=="V4"] <- "chr"
head(snp_100kbp_ranges_count_wide)
snp_100kbp_ranges_count_wide <- subset(snp_100kbp_ranges_count_wide, snp_100kbp_ranges_count_wide$snp_ratio != "Inf")
snp_100kbp_ranges_count_wide$snp_ratio_scaled <- snp_100kbp_ranges_count_wide$snp_ratio / median(snp_100kbp_ranges_count_wide$snp_ratio, na.rm = TRUE)
snp_100kbp_ranges_count_wide$snp_diff_scaled <- snp_100kbp_ranges_count_wide$snp_diff / median(snp_100kbp_ranges_count_wide$snp_diff, na.rm = TRUE)

outfile <- sprintf(paste("results/",species,"/",species,".100kbp_snp_diff.txt",sep=""))
write.table(snp_100kbp_ranges_count_wide, outfile, quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")

######### PLOTTING

factor.nr <- as.numeric(length(unique(cov.select$factor)))

circos.clear()
pdf(file=args[8], width = 9, height = 9)
par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7) 
circos.par(cell.padding = c(0, 0, 0, 0), "track.height" = 0.15, gap.after = c(rep(1, factor.nr-1), 10))
circos.initialize(factors = cov.select$factor, x = cov.select$x)


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


circos.clear()
pdf(file=paste("results/",species,"/",species,".pdf",sep=""), width = 9, height = 9)
par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7) 
circos.par(cell.padding = c(0, 0, 0, 0), "track.height" = 0.15, gap.after = c(rep(1, factor.nr-1), 10))
circos.initialize(factors = cov.select$factor, x = cov.select$x)


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

