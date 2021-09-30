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
out = args[2]
windowsize = args[3]

windowsize <- as.integer(windowsize)

################################################################################
################################# READ FILES ###################################
################################################################################

cov = read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
cov <- plyr::rename(cov, c("V1"="chr", "V2"="start", "V3"="end", "V4"="heterogametic", 
                           "V5"="homogametic"))

#if (file.exists(chr_file)) { 
#  
#  chromosomes <- read.csv(chr_file, header = FALSE, sep = ",")
#  cov <- remove_chr_not_in_list(cov, chromosomes)
#  
#}

cov <- calculate_ratio(cov)
cov <- calculate_diff(cov)
#cov <- remove_outliers(cov)


################################################################################
################################## WINDOWS #####################################
################################################################################

cov <- remove_chr_less_than_value(cov, windowsize)

if (dim(cov)[1] > 0) {
  
  mean_bp <- transform(cov, range=floor(start/windowsize))
  mean_bp_ranges_start <- summaryBy(start ~ chr + range, FUN = min, data=mean_bp, keep.names=TRUE, na.rm = TRUE)
  mean_bp_ranges_end <- summaryBy(end ~ chr + range, FUN = max, data=mean_bp, keep.names=TRUE, na.rm = TRUE)
  mean_bp <- mean_win(mean_bp, heterogametic + homogametic + ratio + ratio_scaled + diff + diff_scaled ~ chr + range)
  mean_bp <- merge(mean_bp,mean_bp_ranges_start,by=c("chr","range"))
  mean_bp <- merge(mean_bp,mean_bp_ranges_end,by=c("chr","range"))
  
  write.table(mean_bp, out, quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")

  
} else {
  
  print("WARNING: No chromosomes/scaffold larger than specified window size")
 # write.table("No chromosomes/scaffold larger than 100kbp", out1Mb, quote=FALSE, row.names = F, col.names = F)
  write.table("No chromosomes/scaffold larger than specified window size", out, quote=FALSE, row.names = F, col.names = F)
  
}

