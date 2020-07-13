# read gencov

cov = read.table("Data/gencov.nodup.nm.all.ZF.out",header=FALSE,fill=TRUE,stringsAsFactor=FALSE)

# look at only Z chromosome
cov <- subset(cov, cov$V8=="Z")
cov <- cov[-4:-7][-7:-9]
cov <- plyr::rename(cov, c("V1"="contig", "V2"="contig_start","V3"="contig_end", 
                           "V8"="chr", "V9"="start", "V10"="end"))

nr_columns <- dim(cov)
sample_name <- colnames(cov)[7:nr_columns[2]]

# for each sample column, calculate mean in 1Mbp and 100 kbp windows
cov_mean <- summaryBy(V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22 + V23 ~ contig , data=cov, keep.names=TRUE)

# plot all samples, first half is heterogametic, second half is homogametic
par(mfrow=c(2,5), mar=c(1,2,2,1))
plot(cov_mean$contig, cov_mean$V14, xaxt = 'n', pch = ".", ylim = c(0,2000), 
     main="SJ-2333-Aarv-249-IT_S21_L002")
plot(cov_mean$contig, cov_mean$V15, xaxt = 'n', pch = ".", ylim = c(0,2000),
     main="SJ-2333-Aarv-252-IT_S23_L002")
plot(cov_mean$contig, cov_mean$V16, xaxt = 'n', pch = ".", ylim = c(0,2000),
     main="SJ-2333-Aarv-441-IT_S25_L002")
plot(cov_mean$contig, cov_mean$V17, xaxt = 'n', pch = ".", ylim = c(0,2000),
     main="SJ-2333-Aarv-460-IT_S27_L002")
plot(cov_mean$contig, cov_mean$V18, xaxt = 'n', pch = ".", ylim = c(0,2000),
     main="QL-1681-19_S46_L006")

plot(cov_mean$contig, cov_mean$V19, xaxt = 'n', pch = ".", ylim = c(0,2000), 
     main="SJ-2333-Aarv-389-IT_S17_L002")
plot(cov_mean$contig, cov_mean$V20, xaxt = 'n', pch = ".", ylim = c(0,2000),
     main="SJ-2333-Aarv-444-IT_S18_L002")
plot(cov_mean$contig, cov_mean$V21, xaxt = 'n', pch = ".", ylim = c(0,2000),
     main="SJ-2333-Aarv-449-IT_S19_L002")
plot(cov_mean$contig, cov_mean$V22, xaxt = 'n', pch = ".", ylim = c(0,2000),
     main="SJ-2333-Aarv-939-IT_S20_L002")
plot(cov_mean$contig, cov_mean$V23, xaxt = 'n', pch = ".", ylim = c(0,2000),
     main="QL-1681-21_S47_L006")

