
remove_chr_less_than_1mb <- function(data_table) {
  # removes all chromosomes that are < 1Mbp, unknown, random, 
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
  
  data_table$ratio <- data_table$heterogametic/data_table$homogametic
  data_table <- subset(data_table, data_table$ratio != "Inf")
  data_table$ratio_scaled <- data_table$ratio / median(data_table$ratio, na.rm = TRUE)
  
  data_table[is.na(data_table)] <- 0
  
  return(data_table)
}

remove_outliers <- function(data_table) {
  
  outliers <- boxplot(ratio ~ chr, data=data_table)$out
  data_table <- data_table[-which(data_table$ratio %in% outliers),]
  
  return(data_table)
}

mean_win <- function(data_table, win_len) {
  # calculates the mean in windows (1Mbp, 100kbp)
  
  data_table <- transform(data_table, range=round(end/win_len))
  data_table <- summaryBy(ratio ~ chr + range, data=data_table, keep.names=TRUE)
  
  return(data_table)
}

mean_diff_win <- function(data_table, win_len) {
  # calculates the mean in windows (1Mbp, 100kbp)

  data_table <- transform(data_table, range=round(end/win_len))
  data_table <- summaryBy(diff ~ chr + range, data=data_table, keep.names=TRUE)

  return(data_table)
}

calculate_diff <- function(data_table) {
  # calculates the normalized difference between the number of sites in a window
  
  data_table$diff <- data_table$heterogametic - data_table$homogametic
  data_table$diff_scaled <- data_table$diff / median(data_table$diff, na.rm = TRUE)
  
  data_table[is.na(data_table)] <- 0
  
  return(data_table)
}

transform_wide <- function(data_table) {
  # Transform to wide format and replace NA with 0
  
  data_table <- reshape2::dcast(data_table, chr + range ~ sex, value.var="count")
  data_table[is.na(data_table)] <- 0
  
  return(data_table)
}

count_snp_win <- function(snp_table, win_len) {
  # Count sites for each sex and per chromosome and window
  
  snp_table <- transform(snp_table, range=round(start/win_len))
  snp_table <- aggregate(count ~ chr + sex + range, data = snp_table, sum)
  
  return(snp_table)
}

