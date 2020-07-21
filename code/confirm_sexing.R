library(doBy)
library(data.table)

# Reads in a bed-file with coverage values for each sample
# Calculates the mean coverage on chr Z and 4 for each individual
# and plots the coverage and the ratio chr z/4

args <- commandArgs(trailingOnly = TRUE)
filename = args[1]
outfilename = args[2]
synteny = args[3]

sample_names <- args[4:length(args)]

# read gencov
cov = read.table(filename,header=FALSE,fill=TRUE,stringsAsFactor=FALSE)

# look at only Z chromosome
if (synteny == "with-synteny") {
  cov <- cov[-4:-7][-7:-9]
}

colnames(cov) <- c("chr", "start", "end", sample_names)

covZ <- subset(cov, cov$chr=="Z")
cov4 <- subset(cov, cov$chr=="4")

covZ_mean <- apply(covZ[,7:length(covZ)], 2, FUN = median)
cov4_mean <- apply(cov4[,7:length(cov4)], 2, FUN = median)

cov_mean <- rbind(covZ_mean, cov4_mean)


pdf(outfilename, width = 14)
par(mfrow=c(1,2), mar=c(5,1,1,1), oma=c(0,10,0,0), xpd=TRUE)

barplot(as.matrix(cov_mean), beside = TRUE, horiz=T, las = 1, col = c("#99CCFF", 
                                                                      "#CCCCFF"),
        xlab = "Median coverage in reads", cex.names=0.7, cex.axis = 0.7)
legend("right", legend = c("chr Z", "chr 4"), fill = c("#99CCFF", "#CCCCFF"))

barplot(covZ_mean/cov4_mean, horiz=T, las = 1, 
        cex.axis = 0.7, col = "#9999CC", names.arg = FALSE, 
        xlab = "Ratio in median coverage chr Z/chr 4")

dev.off()

