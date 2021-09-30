library(doBy)
library(data.table)

# Calculates the mean of a statistica (genome coverage or allele divergence) 
# in 1Mbp and 100kbp windows and for each chromosome/scaffold from 5kbp windows
# Columns in infile: chromosome/scaffold start_window end_window gencov_heterogametic gencov_homogametic
# Call: Rscript {r-script} {infile} {out 1Mbp} {out 100kbp} {out chromosome} {chr-list}

set.seed(999)
source("workflow/scripts/functions.R")

args <- commandArgs(trailingOnly = TRUE)

filename = args[1]
out = args[2]

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
################################# CHROMOSOME ###################################
################################################################################

mean_whole_chr <- mean_win(cov, heterogametic + homogametic + ratio + ratio_scaled + diff + diff_scaled ~ chr)

len_chr <- length_chr(cov)
colnames(len_chr) <- c("chr", "length")

mean_whole_chr <- merge(mean_whole_chr, len_chr, by = "chr")

write.table(mean_whole_chr, out, quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")

