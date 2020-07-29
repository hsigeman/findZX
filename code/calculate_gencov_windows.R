library(doBy)
library(data.table)

# Calculates the mean genome of a statistica in 1Mbp and 100kbp windows from 5kbp windows
# Can calculate the mean of genome coverage or allele diversity (for 5kbp windows)


set.seed(999)

args <- commandArgs(trailingOnly = TRUE)

filename = args[1]
chr_file = args[2]

source("code/functions.R")

cov = read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov <- plyr::rename(cov, c("V1"="chr", "V2"="start", "V3"="end", "V4"="heterogametic", 
                           "V5"="homogametic"))

if (file.exists(chr_file)) { 
  
  chromosomes <- read.csv(chr_file, header = FALSE, sep = ",")
  snp <- remove_chr_not_in_list(snp, chromosomes)
  
}

cov <- calculate_ratio(cov)
cov <- remove_outliers(cov)


mean_whole_chr <- mean_win(cov, ratio ~ chr)

len_chr <- length_chr(cov)
colnames(len_chr) <- c("chr", "length")

mean_whole_chr <- merge(mean_whole_chr, len_chr, by = "chr")

write.table(mean_whole_chr, args[4], quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")


cov <- remove_chr_less_than_1mb(cov)

if (dim(cov)[1] > 0) {
  
  mean_1Mb <- transform(cov, range=round(end/1000000))
  mean_1Mb <- mean_win(mean_1Mb, ratio ~ chr + range)
  
  mean_100kb <- transform(cov, range=round(end/100000))
  mean_100kb <- mean_win(mean_100kb, ratio ~ chr + range)
  
  write.table(mean_1Mb_ranges, args[2], quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")
  
  write.table(mean_100kb_ranges, args[3], quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")
  
} else {
  
  print("WARNING: No chromosomes/scaffold larger than 1Mbp")
  write.table("No chromosomes/scaffold larger than 1Mbp", args[2], quote=FALSE, row.names = F, col.names = F)
  write.table("No chromosomes/scaffold larger than 1Mbp", args[3], quote=FALSE, row.names = F, col.names = F)
  
}



