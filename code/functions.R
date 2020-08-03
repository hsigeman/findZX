
remove_chr_less_than_1mb <- function(data_table) {
  # removes all chromosomes that are < 1Mbp, unknown, random, 
  # LG or MT chromosomes.
  
  max_per_chr <- setDT(data_table)[, .SD[which.max(end)], by=chr]
  chr_over_1mb <- subset(max_per_chr$chr, max_per_chr$end>1000000)
  data_table <- data_table[ data_table$chr %in% chr_over_1mb, ]
  
  return(data_table)
}

remove_chr_not_in_list <- function(data_table, chr_list) {
  # only keeps chromosomes specified in a list

  data_table <- data_table[ data_table$chr %in% chr_list, ]
  
  return(data_table)
}

length_chr <- function(data_table) {
  # returns the length of each chr/scaffold
  
  max_per_chr <- setDT(data_table)[, .SD[which.max(end)], by=chr]
  
  max_per_chr <- as.data.frame(max_per_chr[,1:2])
  
  return(max_per_chr)
}

remove_outliers <- function(data_table) {
  
  outliers <- boxplot(ratio ~ chr, data=data_table, plot = FALSE)$out
  data_table <- data_table[!(data_table$ratio %in% outliers),]
  
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

calculate_diff <- function(data_table) {
  # calculates the normalized difference between the number of sites in a window
  
  data_table$diff <- data_table$heterogametic - data_table$homogametic
  data_table$diff_scaled <- data_table$diff / median(data_table$diff, na.rm = TRUE)
  
  data_table[is.na(data_table)] <- 0
  
  return(data_table)
}

mean_win <- function(data_table, formula) {
  # calculates the mean in ratio for each chr
  #formula = ratio ~ chr, diff ~ chr, diff ~ chr + range, ratio ~ chr + range
  
  data_table <- summaryBy(formula, data=data_table, keep.names=TRUE, na.rm = TRUE)
  
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
  
  snp_table <- transform(snp_table, range=floor(start/win_len))
  snp_table <- aggregate(count ~ chr + sex + range, data = snp_table, sum)
  
  return(snp_table)
}

gen_data_4plotting <- function(data_table, columns) {
  
  data_table$chr <- as.factor(data_table$chr)
  
  myvars <- columns
  setDF(data_table)
  data_table.select <- data_table[myvars]
  data_table.select <- na.omit(data_table.select)
  colnames(data_table.select) <- c("factor", "x", "y")
  
  max.data_table <- max(data_table.select$y)
  min.data_table <- min(data_table.select$y)
  median.data_table <- median(data_table.select$y)
  
  return(list(df = data_table.select, max = max.data_table, min = min.data_table, median = median.data_table))
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
