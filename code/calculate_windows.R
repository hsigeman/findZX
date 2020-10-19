library(doBy)
library(data.table)

# Calculates the mean of a statistica (genome coverage or allele divergence) 
# in 1Mbp and 100kbp windows and for each chromosome/scaffold from 5kbp windows
# Columns in infile: chromosome/scaffold start_window end_window gencov_heterogametic gencov_homogametic
# Call: Rscript {r-script} {infile} {out 1Mbp} {out 100kbp} {out chromosome} {chr-list}

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

cov = read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov <- plyr::rename(cov, c("V1"="chr", "V2"="start", "V3"="end", "V4"="heterogametic", 
                           "V5"="homogametic"))

if (file.exists(chr_file)) { 
  
  chromosomes <- read.csv(chr_file, header = FALSE, sep = ",")
  cov <- remove_chr_not_in_list(cov, chromosomes)
  
}

cov <- calculate_ratio(cov)
cov <- calculate_diff(cov)
#cov <- remove_outliers(cov)

################################################################################
################################# CHROMOSOME ###################################
################################################################################

mean_whole_chr <- mean_win(cov, heterogametic + homogametic + ratio + ratio_scaled + diff + diff_scaled ~ chr)

len_chr <- length_chr(cov)
colnames(len_chr) <- c("chr", "length")

mean_whole_chr <- merge(mean_whole_chr, len_chr, by = "chr")

write.table(mean_whole_chr, outChr, quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")

################################################################################
################################## WINDOWS #####################################
################################################################################

cov <- remove_chr_less_than_value(cov, 100000)

if (dim(cov)[1] > 0) {
  
  mean_100kb <- transform(cov, range=floor(start/100000))
  mean_100kb_ranges_start <- summaryBy(start ~ chr + range, FUN = min, data=mean_100kb, keep.names=TRUE, na.rm = TRUE)
  mean_100kb_ranges_end <- summaryBy(end ~ chr + range, FUN = max, data=mean_100kb, keep.names=TRUE, na.rm = TRUE)
  mean_100kb <- mean_win(mean_100kb, heterogametic + homogametic + ratio + ratio_scaled + diff + diff_scaled ~ chr + range)
  mean_100kb <- merge(mean_100kb,mean_100kb_ranges_start,by=c("chr","range"))
  mean_100kb <- merge(mean_100kb,mean_100kb_ranges_end,by=c("chr","range"))
  
  write.table(mean_100kb, out100kb, quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")
  
  
  cov <- remove_chr_less_than_value(cov, 1000000)
  
  if (dim(cov)[1] > 0) {
    
    mean_1Mb <- transform(cov, range=floor(start/1000000))
    mean_1Mb_ranges_start <- summaryBy(start ~ chr + range, FUN = min, data=mean_1Mb, keep.names=TRUE, na.rm = TRUE)
    mean_1Mb_ranges_end <- summaryBy(end ~ chr + range, FUN = max, data=mean_1Mb, keep.names=TRUE, na.rm = TRUE)
    mean_1Mb <- mean_win(mean_1Mb, heterogametic + homogametic + ratio + ratio_scaled + diff + diff_scaled ~ chr + range)
    mean_1Mb <- merge(mean_1Mb,mean_1Mb_ranges_start,by=c("chr","range"))
    mean_1Mb <- merge(mean_1Mb,mean_1Mb_ranges_end,by=c("chr","range"))
    
    write.table(mean_1Mb, out1Mb, quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")
    
  } else {
    
    print("WARNING: No chromosomes/scaffold larger than 1Mbp")
    write.table("No chromosomes/scaffold larger than 1Mbp", out1Mb, quote=FALSE, row.names = F, col.names = F)
  }
  
} else {
  
  print("WARNING: No chromosomes/scaffold larger than 100kbp")
  write.table("No chromosomes/scaffold larger than 100kbp", out1Mb, quote=FALSE, row.names = F, col.names = F)
  write.table("No chromosomes/scaffold larger than 100kbp", out100kb, quote=FALSE, row.names = F, col.names = F)
  
}

