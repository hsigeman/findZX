library(doBy)
library(data.table)

# Takes a bed file with the difference in heterozygosity for sites. 
# Calculates the mean in 1Mbp and 100kbp windows for each sex.

set.seed(999)

source("code/functions.R")
args <- commandArgs(trailingOnly = TRUE)

filename = args[1]

snp = read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)

snp <- plyr::rename(snp, c("V1"="chr", "V2"="start",
                           "V3"="end", "V4"="diff"))

snp_mean_chr <- mean_diff_chr(snp)

write.table(snp_mean_chr, args[4], quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")


snp <- remove_chr_less_than_1mb(snp)


snp_1Mb_ranges_mean <- mean_diff_win(snp, 1000000)

write.table(snp_1Mb_ranges_mean, args[2], quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")


snp_100kb_ranges_mean <- mean_diff_win(snp, 100000)

write.table(snp_100kb_ranges_mean, args[3], quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")
