library(doBy)
library(data.table)

# read gencov
cov = read.table("Data/gencov.nodup.nm.all.ZF.out",header=FALSE,fill=TRUE,stringsAsFactor=FALSE)

# look at only Z chromosome
cov <- subset(cov, cov$V8=="Z")
cov <- cov[-4:-7][-7:-9]
cov <- plyr::rename(cov, c("V1"="contig", "V2"="range_start","V3"="range_end", 
                           "V8"="chr", "V9"="start", "V10"="end"))

nr_columns <- dim(cov)
sample_name <- colnames(cov)[7:nr_columns[2]]

# for each sample column, calculate mean in 1Mbp and 100 kbp windows
cov <- transform(cov, range=round(end/1000000))
cov_mean <- summaryBy(V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22 + V23 ~ range , data=cov, keep.names=TRUE)


#cov_mean <- transform(cov_mean, range = as.factor(range))
pdf("figures/gencov.nm.all.pdf", width = 14)

# plot all samples, first half is heterogametic, second half is homogametic
par(mfrow=c(1,2), mar=c(1,2,2,1))
plot(cov_mean$range, cov_mean$V14, xaxt = 'n', type="l", ylim = c(0,2000), 
     main="Labeled heterogametic")
abline(h=median(cov_mean$V14), lty = 2)
lines(cov_mean$range, cov_mean$V15, xaxt = 'n', type="l", ylim = c(0,2000),
     col = "red")
abline(h=median(cov_mean$V15), col = "red", lty = 2)
lines(cov_mean$range, cov_mean$V16, xaxt = 'n', type="l", ylim = c(0,2000),
     col = "blue")
abline(h=median(cov_mean$V16), col = "blue", lty = 2)
lines(cov_mean$range, cov_mean$V17, xaxt = 'n', type="l", ylim = c(0,2000),
     col = "green")
abline(h=median(cov_mean$V17), col = "green", lty = 2)
lines(cov_mean$range, cov_mean$V18, xaxt = 'n', type="l", ylim = c(0,2000),
     col = "purple")
abline(h=median(cov_mean$V18), col = "purple", lty = 2)
legend("bottomleft", legend = c("SJ-2333-Aarv-249-IT_S21_L002", "SJ-2333-Aarv-252-IT_S23_L002",
                             "SJ-2333-Aarv-441-IT_S25_L002", "SJ-2333-Aarv-460-IT_S27_L002",
                             "QL-1681-19_S46_L006"), fill = c("black" ,"red","blue","green","purple"))

plot(cov_mean$range, cov_mean$V19, xaxt = 'n', type="l", ylim = c(0,2000),
     main = "Labeled homogametic")
abline(h=median(cov_mean$V19), lty = 2)
lines(cov_mean$range, cov_mean$V20, xaxt = 'n', type="l", ylim = c(0,2000),
     col = "red")
abline(h=median(cov_mean$V20), col = "red", lty = 2)
lines(cov_mean$range, cov_mean$V21, xaxt = 'n', type="l", ylim = c(0,2000),
     col = "blue")
abline(h=median(cov_mean$V21), col = "blue", lty = 2)
lines(cov_mean$range, cov_mean$V22, xaxt = 'n', type="l", ylim = c(0,2000),
     col = "green")
abline(h=median(cov_mean$V22), col = "green", lty = 2)
lines(cov_mean$range, cov_mean$V23, xaxt = 'n', type="l", ylim = c(0,2000),
     col = "purple")
abline(h=median(cov_mean$V23), col = "purple", lty = 2)
legend("bottomleft", legend = c("SJ-2333-Aarv-389-IT_S17_L002", "SJ-2333-Aarv-444-IT_S18_L002",
                             "SJ-2333-Aarv-449-IT_S19_L002", "SJ-2333-Aarv-939-IT_S20_L002",
                             "QL-1681-21_S47_L006"), fill = c("black" ,"red","blue","green","purple"))

dev.off()

