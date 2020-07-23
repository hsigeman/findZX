library(doBy)
library(data.table)

########## FUNCTIONS ##########

gen_data_4plotting <- function(filename, y_column) {
  data_table =read.table(filename, header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
  data_table$chr <- as.factor(data_table$chr)
  
  myvars <- c("chr", "range", y_column)
  setDF(data_table)
  data_table.select <- data_table[myvars]
  data_table.select <- na.omit(data_table.select)
  colnames(data_table.select) <- c("factor", "x", "y")

  max.data_table <- max(data_table.select$y)
  min.data_table <- min(data_table.select$y)
  median.data_table <- median(data_table.select$y)

  return(list(df = data_table.select, max = max.data_table, min = min.data_table, median = median.data_table))
}

color_pallete_function <- colorRampPalette(
  colors = c("red", "orange", "blue"),
  space = "Lab" # Option used when colors do not represent a quantitative scale
)

########## CALCULATION ##########

set.seed(999)

args <- commandArgs(trailingOnly = TRUE)

file00 = args[1]
file02 = args[2]
file04 = args[3]
filesnp = args[4]


cov_00 <- gen_data_4plotting(file00, "ratio")
cov.select.00 <- cov_00$df
max.cov.00 <- cov_00$max
max.cov.00
min.cov.00 <- cov_00$min
min.cov.00
median.cov.00 <- cov_00$median
median.cov.00

cov_02 <- gen_data_4plotting(file02, "ratio")
cov.select.02 <- cov_02$df
max.cov.02 <- cov_02$max
max.cov.02
min.cov.02 <- cov_02$min
min.cov.02
median.cov.02 <- cov_02$median
median.cov.02

cov_04 <- gen_data_4plotting(file04, "ratio")
cov.select.04 <- cov_04$df
max.cov.04 <- cov_04$max
max.cov.04
min.cov.04 <- cov_04$min
min.cov.04
median.cov.04 <- cov_04$median
median.cov.04

snp <- gen_data_4plotting(filesnp, "diff")
snp.select <- snp$df
max.snp <- snp$max
max.snp
min.snp <- snp$min
min.snp
median.snp <- snp$median
median.snp

cov.select <- merge(cov.select.00, cov.select.02, by = c("factor", "x"))
cov.select <- merge(cov.select, cov.select.04, by = c("factor", "x"))
cov.select <- merge(cov.select, snp.select, by = c("factor", "x"))

colnames(cov.select) <- c("factor", "x", "cov00", "cov02", "cov04", "hetDiff")

# Color each chromosome with a unique color
c <- color_pallete_function(length(unique(cov.select$factor)))


par(mfrow=c(3,1), mar=c(4,4,4,1), oma=c(0,10,0,0), xpd=TRUE)

plot(cov.select$cov00, cov.select$hetDiff, col = cov.select$factor, xlim = c(min.cov.00, max.cov.00),
     xlab = "Normalized genome coverage, nm = 0", ylab = "Difference in proportion of heterozygosity")
abline(h = median.snp)
abline(v = median.cov.00)
plot(cov.select$cov02, cov.select$hetDiff, col = cov.select$factor, xlim = c(min.cov.00, max.cov.00),
     xlab = "Normalized genome coverage, nm = 2", ylab = "Difference in proportion of heterozygosity")
abline(h = median.snp)
abline(v = median.cov.02)
par(xpd=NA)
legend("left", legend=unique(cov.select$factor),fill=1:length(cov.select$factor),inset=c(-0.2,0.6))
par(xpd=TRUE)
plot(cov.select$cov04, cov.select$hetDiff, col = cov.select$factor, xlim = c(min.cov.00, max.cov.00),
     xlab = "Normalized genome coverage, nm = 4", ylab = "Difference in proportion of heterozygosity")
abline(h = median.snp)
abline(v = median.cov.04)


# plot het/cov color after chr len
# 3d plot het/cov/len