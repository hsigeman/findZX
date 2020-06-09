#install.packages("circlize")
library(circlize)
library(doBy)
library(data.table)

set.seed(999)

args <- commandArgs(trailingOnly = TRUE)

species=args[1]
female=args[2]
male=args[3]


Rscript cal_cov_ratio_ranges.R species 00 args[4]
Rscript cal_cov_ratio_ranges.R species 02 args[5]
Rscript cal_cov_ratio_ranges.R species 04 args[6]


cov_00 =read.table(paste("results/",species,"/",species,".1Mbp_cov_scaled.nm.00.txt",sep=""),
                   header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov_00$V1 <- as.factor(cov_00$V1)

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



cov_02 =read.table(paste("results/",species,"/",species,".1Mbp_cov_scaled.nm.02.txt",sep=""),
                   header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov_02$V1 <- as.factor(cov_02$V1)

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



cov_04 =read.table(paste("results/",species,"/",species,".1Mbp_cov_scaled.nm.00.txt",sep=""),
                   header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov_04$V1 <- as.factor(cov_04$V1)

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