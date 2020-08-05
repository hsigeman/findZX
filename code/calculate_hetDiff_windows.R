library(doBy)
library(data.table)

# Takes as input a bed file containing the difference in heterozygosity between 
# heterogametic and homogametic individuals for each variable site. 
# Calculates the mean difference in heterozygosity in 1Mbp and 100kbp windows 
# and for each chromsome/scaffold.
# Columns in input file: chromosome/scaffold start_position end_position difference_in_heterozygosity
# Call: Rscriipt {r-script} {infile} {out 1Mbp} {out 100kbp} {out chromsome/scaffold} {chr-list}

set.seed(999)
source("code/functions.R")

args <- commandArgs(trailingOnly = TRUE)

filename = args[1]
out1Mb = args[2]
out100kb = args[3]
outChr = args[4]
chr_file = args[5]

################################################################################
################################# READ FILES ###################################
################################################################################

snp = read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)

snp <- plyr::rename(snp, c("V1"="chr", "V2"="start",
                           "V3"="end", "V4"="diff"))

if (file.exists(chr_file)) { 
  
  chromosomes <- read.csv(chr_file, header = FALSE, sep = ",")
  snp <- remove_chr_not_in_list(snp, chromosomes)
  
}

################################################################################
################################# CROMOSOMES ###################################
################################################################################

snp_mean_chr <- mean_win(snp, diff ~ chr)
len_chr <- length_chr(snp)
colnames(len_chr) <- c("chr", "length")

snp_mean_chr <- merge(snp_mean_chr, len_chr, by = "chr")

write.table(snp_mean_chr, outChr, quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")

################################################################################
################################## WINDOWS #####################################
################################################################################

snp <- remove_chr_less_than_1mb(snp)

if (dim(snp)[1] > 0) {
  
  snp_1Mb_mean <- transform(snp, range=floor(end/1000000))
  snp_1Mb_mean <- mean_win(snp_1Mb_mean, diff ~ chr + range)
  
  write.table(snp_1Mb_mean, out1Mb, quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")
  
  snp_100kbp_mean <- transform(snp, range=floor(end/100000))
  snp_100kbp_mean <- mean_win(snp_100kbp_mean, diff ~ chr + range)
  
  write.table(snp_100kb_ranges_mean, out100kb, quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")
  
} else {
  
  print("WARNING: No chromosomes/scaffold larger than 1Mbp")
  write.table("No chromosomes/scaffold larger than 1Mbp", args[2], quote=FALSE, row.names = F, col.names = F)
  write.table("No chromosomes/scaffold larger than 1Mbp", args[3], quote=FALSE, row.names = F, col.names = F)

}

