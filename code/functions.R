
remove_chr_less_than_1mb <- function(data_table) {
  # removes all chromosomes that are < 1Mbp, unkown, random, 
  # LG or MT chromosomes.
  
  max_per_chr <- setDT(data_table)[, .SD[which.max(end)], by=chr]
  chr_over_1mb <- subset(max_per_chr$chr, max_per_chr$end>1000000)
  data_table <- data_table[ data_table$chr %in% chr_over_1mb, ]
  
  data_table <- subset(data_table, data_table$chr!="Un")
  random <- unique(data_table$chr[grep("random", data_table$chr)])
  data_table <- data_table[ ! data_table$chr %in% random, ]
  random <- unique(data_table$chr[grep("LG", data_table$chr)])
  data_table <- data_table[ ! data_table$chr %in% random, ]
  random <- unique(data_table$chr[grep("MT", data_table$chr)])
  data_table <- data_table[ ! data_table$chr %in% random, ]
  
  return(data_table)
}

calculate_ratio <- function(data_table) {
  # Calculates the coverage or snp ratio between female and male normalized by the
  # median ratio. Infinite values are removed.
  
  data_table$ratio <- data_table$female/data_table$male
  data_table <- subset(data_table, data_table$ratio != "Inf")
  #data_table$ratio_scaled <- data_table$ratio / median(data_table$ratio, na.rm = TRUE)
  
  return(data_table)
}

remove_outliers <- function(cov_table) {
  
  outliers <- boxplot(ratio ~ chr, data=cov_table)$out
  cov_table <- cov_table[-which(cov_table$ratio %in% outliers),]
  
  return(cov_table)
}

count_cov_win <- function(cov_table, win_len) {
  # calculates the mean coverage ration in range windows (1Mbp, 100kbp)
  
  cov_table <- transform(cov_table, range=round(end/win_len))
  cov_table <- summaryBy(ratio ~ chr + range, data=cov_table, keep.names=TRUE)
  
  return(cov_table)
}

calculate_diff <- function(data_table) {
  # calculates the normalized difference between the number of female and male snps
  
  data_table$diff <- data_table$female - data_table$male
  data_table$diff_scaled <- data_table$diff / median(data_table$diff, na.rm = TRUE)
  
  return(data_table)
}

transform_wide <- function(snp_table) {
  # Transform to wide format and replace NA with 0
  
  snp_table <- reshape2::dcast(snp_table, chr + range ~ sample, value.var="count")
  snp_table[is.na(snp_table)] <- 0
  
  return(snp_table)
}

count_snp_win <- function(snp_table, win_len) {
  # Count private female and male alleles per chromosome and win_len
  
  snp_table <- transform(snp_table, range=round(start/win_len))
  snp_table <- aggregate(count ~ chr + sample + range, data = snp_table[which(snp_table$type=="S"),], sum)
  
  return(snp_table)
}
