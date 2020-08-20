library(doBy)
library(data.table)

# Takes as input a bed file containing the difference in heterozygosity between 
# heterogametic and homogametic individuals for each variable site. 
# Calculates the mean difference in heterozygosity in 1Mbp and 100kbp windows 
# and for each chromosome/scaffold.
# Columns in input file: chromosome/scaffold start_position end_position difference_in_heterozygosity
# Call: Rscriipt {r-script} {infile} {out 1Mbp} {out 100kbp} {out chromosome/scaffold} {chr-list}

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
################################# CHROMOSOMES ##################################
################################################################################

snp_mean_chr <- mean_win(snp, diff ~ chr)
len_chr <- length_chr(snp)
colnames(len_chr) <- c("chr", "length")

snp_mean_chr <- merge(snp_mean_chr, len_chr, by = "chr")

write.table(snp_mean_chr, outChr, quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")

################################################################################
################################## WINDOWS #####################################
################################################################################

snp <- remove_chr_less_than_value(snp, 100000)

if (dim(snp)[1] > 0) {
  
  snp_100kbp_mean <- transform(snp, range=floor(start/100000))
  snp_100kbp_ranges_start <- summaryBy(start ~ chr + range, FUN = min, data=snp_100kbp_mean, keep.names=TRUE, na.rm = TRUE)
  snp_100kbp_ranges_end <- summaryBy(end ~ chr + range, FUN = max, data=snp_100kbp_mean, keep.names=TRUE, na.rm = TRUE)
  snp_100kbp_mean <- mean_win(snp_100kbp_mean, diff ~ chr + range)
  snp_100kbp_mean <- merge(snp_100kbp_mean,snp_100kbp_ranges_start,by=c("chr","range"))
  snp_100kbp_mean <- merge(snp_100kbp_mean,snp_100kbp_ranges_end,by=c("chr","range"))
  
  write.table(snp_100kbp_mean, out100kb, quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")
  
  
  snp <- remove_chr_less_than_value(snp, 1000000)
  
  if (dim(snp)[1] > 0) {
    
    snp_1Mb_mean <- transform(snp, range=floor(start/1000000))
    snp_1Mb_ranges_start <- summaryBy(start ~ chr + range, FUN = min, data=snp_1Mb_mean, keep.names=TRUE, na.rm = TRUE)
    snp_1Mb_ranges_end <- summaryBy(end ~ chr + range, FUN = max, data=snp_1Mb_mean, keep.names=TRUE, na.rm = TRUE)
    snp_1Mb_mean <- mean_win(snp_1Mb_mean, diff ~ chr + range)
    snp_1Mb_mean <- merge(snp_1Mb_mean,snp_1Mb_ranges_start,by=c("chr","range"))
    snp_1Mb_mean <- merge(snp_1Mb_mean,snp_1Mb_ranges_end,by=c("chr","range"))
    
    write.table(snp_1Mb_mean, out1Mb, quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")
    
  } else {
    print("WARNING: No chromosomes/scaffold larger than 1Mbp")
    write.table("No chromosomes/scaffold larger than 1Mbp", out1Mb, quote=FALSE, row.names = F, col.names = F)
  }
  
  
} else {
  
  print("WARNING: No chromosomes/scaffold larger than 100kbp")
  write.table("No chromosomes/scaffold larger than 100kbp", out1Mb, quote=FALSE, row.names = F, col.names = F)
  write.table("No chromosomes/scaffold larger than 100kbp", out100kb, quote=FALSE, row.names = F, col.names = F)

}

