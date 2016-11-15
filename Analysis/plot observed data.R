## make a plot of episizes through time
pdf(paste0(folder_out, "examine episizes.pdf"), height=8, width=12, useDingbats=F)
par(mfrow=c(2,1), mar=c(3,4,1,0.5), las=1, cex=0.9)
boxplot(t(episizes), axes=F, col="white", border="white", ylab="epidemic size")
#episizes$year <- 1986:2014
for(i in 1:29) {
  points(y=episizes[i,], x=jitter(rep(i, 19), amount=0.2), pch=19, col=transp("dodgerblue3", 0.6))  
}
box()
axis(2)
axis(1, at=seq(1:29), labels=as.character(seq(1986,2014)))

boxplot(t(episizes), axes=F, col="white", border="white", ylab="epidemic size")
#episizes$year <- 1986:2014
for(i in 1:29) {
  points(y=episizes[i,], x=jitter(rep(i, 19), amount=0.2), pch=19, col=transp("dodgerblue3", 0.6))  
}
box()
axis(2)
axis(1, at=seq(1:29), labels=as.character(seq(1986,2014)))
boxplot(t(episizes), axes=F,  add=T, col=transp("white", 0))
dev.off()

pdf(paste0(folder_out, "population sizes.pdf"), height=3, width=6, useDingbats=F)
par(mfrow=c(1,1), mar=c(4,4,0,0.5), las=1)
hist(lookup$Pop_1999, col=transp("dodgerblue3", 0.6), xlim=c(1, max(lookup$Pop_1999)), 
     breaks=30, main="", xlab="population")
dev.off()

library(fBasics)
# example distribution of epi sizes
pdf(paste0(folder_out, "example city epi size dist-plain.pdf"), height=6, width=8, useDingbats=F)
par(mfrow=c(2,1), mar=c(4,4,1,0.5), las=1)
# plain
hist(episizes[,10], main="", xlim=c(0, 10000), breaks=seq(0, 10000, 500),
     col=transp("dodgerblue", 0.7), ylab="frequency", 
     bty="n", border="dodgerblue4", xlab="epidemic size per 100,000")
# add time series
plot(episizes[,10], main="", ylim=c(0,10000),
     col=transp("dodgerblue", 0.7), pch=19, ylab="epidemic size per 100,000",
     bty="n",  xlab="year", axes=F, type="b")
axis(2, cex=0.9)
axis(1, at=seq(1:29), labels=as.character(seq(1986,2014)), cex=0.9)
dev.off()

#add means etc
hist(episizes[,10], main="", xlim=c(0, 10000), breaks=seq(0, 10000, 500),
     col=transp("dodgerblue", 0.7), ylab="frequency", 
     bty="n", border="dodgerblue4", xlab="epidemic size per 100,000")
abline(v=median(episizes[,10]), col=transp("grey45", 0.7), lty=3, lwd=2, xpd=T)
abline(v=mean(episizes[,10]), col=transp("grey45", 0.7), lty=1, lwd=2, xpd=T)
abline(v=quantile(episizes[,10], 0.25), col=transp("grey45", 0.7), lty=3, lwd=2, xpd=T)
abline(v=quantile(episizes[,10], 0.75), col=transp("grey45", 0.7), lty=3, lwd=2, xpd=T)
legend("topright", legend=c(paste0("skew: ", round(skewness(episizes[,10], method="moment"), 2)),
                            paste0("standard deviation: ", round(sd(episizes[,10]),0))),
       bty="n")
# add time series
plot(episizes[,10], main="", ylim=c(0,10000),
     col=transp("dodgerblue", 0.7), pch=19, ylab="epidemic size per 100,000",
     bty="n",  xlab="year", axes=F, type="b")
axis(2, cex=0.9)
axis(1, at=seq(1:29), labels=as.character(seq(1986,2014)), cex=0.9)
abline(h=median(episizes[,10]), col=transp("grey65", 0.7), lty=3, lwd=2)
legend("bottomright", legend=c(paste0("median crossing: ", 0.36)),
       bty="n")
dev.off()

# plot example sequence
pdf(paste0(folder_out, "example episizes.pdf"), height=8, width=12, useDingbats=F)
par(mfrow=c(2,1), mar=c(3,4,1,0.5), las=1, cex=0.9)
boxplot(t(episizes), axes=F, col="white", border="white", ylab="epidemic size")
points(episizes[,1], pch=19, col=transp("dodgerblue3", 0.6))
box()
axis(2)
axis(1, at=seq(1:29), labels=as.character(seq(1986,2014)))
legend("topleft", legend="Alsace", cex=1.7, col="white", lwd=2, bty="n")
dev.off()

# plot metric against pop size
compare.xy <- function(metric.name) {
  plot.this <- metric.name
  plot(x=data.metrics[[plot.this]][ , "population"], y=data.metrics[[plot.this]][ , "value"], main="",
       ylab="value", pch=19, col=transp("dodgerblue", 0.7),
       log="x", bty="n", xlim=c(100000, max(lookup$Pop_1999, na.rm=T)))
  legend("topleft", legend=metric.name, cex=1.7, col="white", lwd=2, bty="n")
}

pdf(paste0(folder_out, "compare pop vs metrics.pdf"), height=6, width=6, useDingbats=F)
par(mfrow=c(4,3), mar=c(3,4,1,0), las=1)
for(i in metric.names) {
  compare.xy(i)
}
dev.off()

# plot metric distribution
hist.metric <- function(metric.name) {
  plot.this <- metric.name
  hist(data.metrics[[plot.this]][ , "value"], main="",
       col=transp("dodgerblue", 0.7), ylab="", bty="n", border="dodgerblue4")
  legend("topright", legend=metric.name, cex=1.7,  bty="n")
}

pdf(paste0(folder_out, "hist of metrics.pdf"), height=6, width=6, useDingbats=F)
par(mfrow=c(4,3), mar=c(3,3,1,1), las=1)
for(i in metric.names) {
  hist.metric(i)
}
dev.off()

