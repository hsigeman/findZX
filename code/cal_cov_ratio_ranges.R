library(doBy)
library(data.table)

set.seed(999)

args <- commandArgs(trailingOnly = TRUE)

species= args[1]
file_type= args[2]
filename= args[3]


source("functions.R")

cov=read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov <- plyr::rename(cov, c("V1"="contig", "V2"="contig_start","V3"="contig_end", 
                           "V4"="chr", "V5"="start", "V6"="end", "V7"="female", 
                           "V8"="male"))

cov <- remove_chr_less_than_1mb(cov)

cov <- calculate_ratio(cov)

# Remove outliers
cov <- remove_outliers(cov)

cov_1Mb_ranges_sum <- count_cov_win(cov, 1000000)
cov_100kb_ranges_sum <- count_cov_win(cov, 100000)

outfile <- sprintf(paste("results/",species,"/",species,".1Mbp.", file_type,".txt",sep=""))
write.table(cov_1Mb_ranges_sum, outfile, quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")

outfile <- sprintf(paste("results/",species,"/",species,".100kbp.", file_type,".txt",sep=""))
write.table(cov_100kb_ranges_sum, outfile, quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")
