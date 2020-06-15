#install.packages("circlize")
#library(circlize)
#install.packages("doBy")
library(doBy)
library(data.table)

set.seed(999)

args <- commandArgs(trailingOnly = TRUE)

species= args[1]
edit_dist= args[2]
filename= args[3]


source("../snakemake-sex-chr/code/functions.R")

cov=read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov <- plyr::rename(cov, c("V1"="contig", "V2"="contig_start","V3"="contig_end", 
                           "V4"="chr", "V5"="start", "V6"="end", "V7"="female", 
                           "V8"="male"))

cov <- remove_chr_less_than_1mb(cov)


cov$chr <- ordered(cov$chr,
                   levels = c("1", "1A", "1B", "2", "3", "4", "4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
                              "15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "Z"))

cov <- calculate_ratio(cov)

# Remove outliers
cov <- remove_outliers(cov)

cov_1Mb_ranges_sum <- count_cov_win(cov, 1000000)
cov_100kb_ranges_sum <- count_cov_win(cov, 100000)

outfile <- sprintf(paste("results/",species,"/",species,".1Mbp_cov_scaled.nm.", edit_dist,".txt",sep=""))
write.table(cov_1Mb_ranges_sum, outfile, quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")

outfile <- sprintf(paste("results/",species,"/",species,".100kbp_cov_scaled.nm.", edit_dist,".txt",sep=""))
write.table(cov_100kb_ranges_sum, outfile, quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")
