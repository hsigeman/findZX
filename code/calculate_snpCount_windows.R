library(doBy)
library(data.table)

# Takes a bed file with sites. Counts the number of sites in 1Mbp and 100kbp windows for each sex. 
# Calculates the ratio and difference in counts between the sexes.
# Columns in infile: chromosome/scaffold start_position end_position sex(heterogametic/homogametic)
# Call: Rscript {r-script} {infile} {out 1Mbp} {out 100kbp} {chr-list}


set.seed(999)
source("code/functions.R")

args <- commandArgs(trailingOnly = TRUE)
filename = args[1]
out1Mb = args[2]
out100kb = args[3]
chr_file = args[4]

################################################################################
################################# READ FILES ###################################
################################################################################

snp = read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
snp <- plyr::rename(snp, c("V1"="chr", "V2"="start",
                           "V3"="end", "V4"="sex"))

if (file.exists(chr_file)) { 
  
  chromosomes <- read.csv(chr_file, header = FALSE, sep = ",")
  snp <- remove_chr_not_in_list(snp, chromosomes)

}

################################################################################
################################ CALCULATIONS ##################################
################################################################################

snp <- remove_chr_less_than_value(snp,1000000)


if (dim(snp)[1] > 0) {
  
  snp$count <- 1
  
  snp_1Mb_ranges_count <- count_snp_win(snp, 1000000)
  
  snp_1Mb_ranges_count_wide <- transform_wide(snp_1Mb_ranges_count)
  
  snp_1Mb_ranges_count_wide <- calculate_ratio(snp_1Mb_ranges_count_wide)
  snp_1Mb_ranges_count_wide <- calculate_diff(snp_1Mb_ranges_count_wide)
  
  write.table(snp_1Mb_ranges_count_wide, out1Mb, quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")
  
  
  
  snp_100kbp_ranges_count <- count_snp_win(snp, 100000)
  
  snp_100kbp_ranges_count_wide <- transform_wide(snp_100kbp_ranges_count)
  
  snp_100kbp_ranges_count_wide <- calculate_ratio(snp_100kbp_ranges_count_wide)
  snp_100kbp_ranges_count_wide <- calculate_diff(snp_100kbp_ranges_count_wide)
  
  write.table(snp_100kbp_ranges_count_wide, out100kb, quote=FALSE, sep="\t", row.names = F, col.names = T, na = "NA")
  
} else {
  
  print("WARNING: No chromosomes/scaffold larger than 1Mbp")
  write.table("No chromosomes/scaffold larger than 1Mbp", out1Mb, quote=FALSE, row.names = F, col.names = F)
  write.table("No chromosomes/scaffold larger than 1Mbp", out100kb, quote=FALSE, row.names = F, col.names = F)
  
}
