library(doBy)
library(data.table)

set.seed(999)

args <- commandArgs(trailingOnly = TRUE)

species= args[1]
filename= args[2]

snp=read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)

snp <- plyr::rename(snp, c("V1"="sex", "V2"="chr",
                           "V3"="start", "V4"="end"))

snp <- remove_chr_less_than_1mb(snp)

snp$count <- 1

snp_1Mb_ranges_count <- count_snp_win(snp, 1000000)

snp_1Mb_ranges_count_wide <- transform_wide(snp_1Mb_ranges_count)

#names(snp_1Mb_ranges_count_wide)[names(snp_1Mb_ranges_count_wide)==female] <- "female"
#names(snp_1Mb_ranges_count_wide)[names(snp_1Mb_ranges_count_wide)==male] <- "male"

snp_1Mb_ranges_count_wide <- calculate_ratio(snp_1Mb_ranges_count_wide)
snp_1Mb_ranges_count_wide <- calculate_diff(snp_1Mb_ranges_count_wide)

write.table(snp_1Mb_ranges_count_wide, args[3], quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")



snp_100kbp_ranges_count <- count_snp_win(snp, 100000)

snp_100kbp_ranges_count_wide <- transform_wide(snp_100kbp_ranges_count)

#names(snp_100kbp_ranges_count_wide)[names(snp_100kbp_ranges_count_wide)==female] <- "female"
#names(snp_100kbp_ranges_count_wide)[names(snp_100kbp_ranges_count_wide)==male] <- "male"

snp_100kbp_ranges_count_wide <- calculate_ratio(snp_100kbp_ranges_count_wide)
snp_100kbp_ranges_count_wide <- calculate_diff(snp_100kbp_ranges_count_wide)

write.table(snp_100kbp_ranges_count_wide, args[4], quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")
