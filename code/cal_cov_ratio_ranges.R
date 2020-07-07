library(doBy)
library(data.table)

set.seed(999)

args <- commandArgs(trailingOnly = TRUE)

species= args[1]
filename= args[2]


source("code/functions.R")

cov=read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov <- plyr::rename(cov, c("V1"="contig", "V2"="contig_start","V3"="contig_end", 
                           "V4"="chr", "V5"="start", "V6"="end", "V7"="heterogametic_sex", 
                           "V8"="homogametic_sex"))

cov <- remove_chr_less_than_1mb(cov)

cov <- calculate_ratio(cov)

# Remove outliers
cov <- remove_outliers(cov)

cov_1Mb_ranges_sum <- count_cov_win(cov, 1000000)
cov_100kb_ranges_sum <- count_cov_win(cov, 100000)

write.table(cov_1Mb_ranges_sum, args[3], quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")

write.table(cov_100kb_ranges_sum, args[4], quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")
