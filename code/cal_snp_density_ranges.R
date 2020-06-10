#install.packages("circlize")
#library(circlize)
library(doBy)
library(data.table)

set.seed(999)

#args <- commandArgs(trailingOnly = TRUE)

species= "LocLus" #args[1]
female= "QF-1504-LOCLUS-43_S1_L001" #args[2]
male= "QF-1504-LOCLUS-24_S3_L001" #args[3]
filename= "../../SupplementaryCode2020/code/R/LocLus.singleton.bestMatch.zf.small" #args[7]

snp=read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)

snp <- plyr::rename(snp, c("V1"="type", "V2"="allel","V3"="sample", 
                           "V4"="chr", "V5"="start", "V6"="end"))

snp <- remove_chr_less_than_1mb(snp)

snp$chr <- ordered(snp$chr,
                  levels = c("1", "1A", "1B", "2", "3", "4", "4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
                             "15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "Z"))

snp$count <- 1
snp <- subset(snp, snp$type=="S")

snp_1Mb_ranges_count <- count_snp_win(snp, 1000000)

snp_1Mb_ranges_count_wide <- transform_wide(snp_1Mb_ranges_count)

names(snp_1Mb_ranges_count_wide)[names(snp_1Mb_ranges_count_wide)==female] <- "female"
names(snp_1Mb_ranges_count_wide)[names(snp_1Mb_ranges_count_wide)==male] <- "male"

snp_1Mb_ranges_count_wide <- calculate_ratio(snp_1Mb_ranges_count_wide)
snp_1Mb_ranges_count_wide <- calculate_diff(snp_1Mb_ranges_count_wide)

outfile <- sprintf(paste("results/",species,"/",species,".1Mbp_snp_diff.txt",sep=""))
write.table(snp_1Mb_ranges_count_wide, outfile, quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")



snp_100kbp_ranges_count <- count_snp_win(snp, 100000)

snp_100kbp_ranges_count_wide <- transform_wide(snp_100kbp_ranges_count)

names(snp_100kbp_ranges_count_wide)[names(snp_100kbp_ranges_count_wide)==female] <- "female"
names(snp_100kbp_ranges_count_wide)[names(snp_100kbp_ranges_count_wide)==male] <- "male"

snp_100kbp_ranges_count_wide <- calculate_ratio(snp_100kbp_ranges_count_wide)
snp_100kbp_ranges_count_wide <- calculate_diff(snp_100kbp_ranges_count_wide)

outfile <- sprintf(paste("results/",species,"/",species,".100kbp_snp_diff.txt",sep=""))
write.table(snp_100kbp_ranges_count_wide, outfile, quote=FALSE, sep="\t", row.names = F, col.names = F, na = "NA")
