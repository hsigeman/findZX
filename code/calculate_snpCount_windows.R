library(doBy)
library(data.table)

# Takes a bed file with sites. Counts the number of sites in 1Mbp and 100kbp windows
# for each sex. Calculates the ratio and difference between sexes.


set.seed(999)

source("code/functions.R")
args <- commandArgs(trailingOnly = TRUE)

filename = args[1]

snp = read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)

snp <- plyr::rename(snp, c("V1"="chr", "V2"="start",
                           "V3"="end", "V4"="sex"))

snp <- remove_chr_less_than_1mb(snp)

snp$count <- 1



snp_1Mb_ranges_count <- count_snp_win(snp, 1000000)

snp_1Mb_ranges_count_wide <- transform_wide(snp_1Mb_ranges_count)

snp_1Mb_ranges_count_wide <- calculate_ratio(snp_1Mb_ranges_count_wide)
snp_1Mb_ranges_count_wide <- calculate_diff(snp_1Mb_ranges_count_wide)

write.table(snp_1Mb_ranges_count_wide, args[2], quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")



snp_100kbp_ranges_count <- count_snp_win(snp, 100000)

snp_100kbp_ranges_count_wide <- transform_wide(snp_100kbp_ranges_count)

snp_100kbp_ranges_count_wide <- calculate_ratio(snp_100kbp_ranges_count_wide)
snp_100kbp_ranges_count_wide <- calculate_diff(snp_100kbp_ranges_count_wide)

write.table(snp_100kbp_ranges_count_wide, args[3], quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")
