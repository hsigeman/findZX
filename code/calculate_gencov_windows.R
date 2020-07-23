library(doBy)
library(data.table)

# Calculates the mean genome of a statistica in 1Mbp and 100kbp windows from 5kbp windows
# Can calculate the mean of genome coverage or allele diversity (for 5kbp windows)


set.seed(999)

args <- commandArgs(trailingOnly = TRUE)

filename = args[1]


source("code/functions.R")

cov = read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov <- plyr::rename(cov, c("V1"="chr", "V2"="start", "V3"="end", "V4"="heterogametic", 
                           "V5"="homogametic"))

cov <- calculate_ratio(cov)
cov <- remove_outliers(cov)


mean_whole_chr <- mean_chr(cov)

len_chr <- length_chr(cov)
colnames(len_chr) <- c("chr", "length")

mean_whole_chr <- merge(mean_whole_chr, len_chr, by = "chr")


cov <- remove_chr_less_than_1mb(cov)

mean_1Mb_ranges <- mean_win(cov, 1000000)
mean_100kb_ranges <- mean_win(cov, 100000)

write.table(mean_1Mb_ranges, args[2], quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")

write.table(mean_100kb_ranges, args[3], quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")

write.table(mean_whole_chr, args[4], quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")
