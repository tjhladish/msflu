
# set paths for who is working
if(Sys.info()[["user"]] %in% c("Rosalind", "eideregg")) {
  setwd("~/Dropbox/LSHTM/Multi_season_flu_transfer/")
  # folder with outputs in
  folder_data <- "~/Dropbox/LSHTM/Multi_season_flu_transfer/"
  # folder with repo in 
  folder_repo <- "~/Documents/Influenza/msflu/Analysis/"
  # outputs folder
  folder_out <- "~/Sync/LSHTM/Collaboration/Multiseason flu/Figures/"  
} else {
  # Tom's paths
}

# # # # # # # # # # # # # # # # # # # # # # # # 
# read in location lookup
lookup <- read.csv(paste0(folder_out, "20160628 - Name_Pop_Lookup.csv"), header=T, na.strings="-")
lookup$Code <- gsub("-", "_", lookup$Code)
lookup <- lookup[lookup$Used==1,]
# read in so smallest first
pop.order <- lookup[order(lookup$Pop_1999),]
# read in epi data
source(paste0(folder_repo, "Read in epidata.R"))

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

# these are the names of the metrics and the rows we want
metric.names <- c("mean", "q0", "q25", "q50", "q75", "q100", "sd", "skew", "mc")

# # # # # # # # # # # # # # # # # # # # # # # # 
# read in json file
require(jsonlite)
data <- fromJSON(txt=paste0(folder_data, "abc_msflu.json"))
data$metrics
hist(data$metrics[grep("mean", data$metrics$name), "value"])
hist(data$metrics[grep("skew", data$metrics$name), "value"])
hist(data$metrics[grep("q0", data$metrics$name), "value"])
hist(data$metrics[grep("mc", data$metrics$name), "value"])

# add place label to metrics
data$metrics$location <- data$metrics$name
for(i in seq_along(metric.names)) {
  print(i)
  data$metrics$location <- gsub(metric.names[i], "FR", data$metrics$location)
}
# remove the letter from the metric name
data$metrics$name <- gsub("_[A-Z]", "", data$metrics$name)

# transform to a list
data.metrics <- list()
for(i in seq_along(metric.names)) {
  print(i)
  data.metrics[[metric.names[i]]] <- data$metrics[grep(metric.names[i], data$metrics$name), c("location", "value")]
  data.metrics[[metric.names[i]]]$population <- lookup$Pop_1999[match(data.metrics[[metric.names[i]]]$location, lookup$Code)]
}

# plot metric against pop size
compare.xy <- function(metric.name) {
  plot.this <- metric.name
  plot(x=data.metrics[[plot.this]][ , "population"], y=data.metrics[[plot.this]][ , "value"], main="",
       ylab="value", pch=19, col=transp("dodgerblue", 0.7),
       log="x", bty="n", xlim=c(100000, max(lookup$Pop_1999, na.rm=T)))
  legend("topleft", legend=metric.name, cex=1.7, col="white", lwd=2, bty="n")
}

pdf(paste0(folder_out, "compare pop vs metrics.pdf"), height=6, width=6, useDingbats=F)
par(mfrow=c(3,3), mar=c(3,4,1,0), las=1)
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
par(mfrow=c(3,3), mar=c(3,3,1,1), las=1)
for(i in metric.names) {
  hist.metric(i)
}
dev.off()

compare.this("q50")

# # # # # # # # # # # # # # # # # # # # # # # # 
# read in results
require("RSQLite")
drv = dbDriver("SQLite")
db = dbConnect(drv, "./msmsc_flu.sqlite")
dbListTables(db)

# # # # # # # # # # # # # # # # # # # # # # # # 
# read in all fitted metrics
basic <- dbGetQuery(db, 'select * from parameters P, metrics M, jobs J 
                          where P.serial = M.serial 
                          and P.serial = J.serial 
                          and status = \'D\'
                          and posterior > -1
                          and smcSet = 12') 
#metrics <- dbGetQuery(db, "SELECT * from metrics")
#jobs <- dbGetQuery(db, "SELECT * from jobs")
#pars <- dbGetQuery(db, "SELECT * from parameters")

# these are the rows we want
these.rows <- 1:500 #deprecated

# subset the parts of the metrics that we want (last SMC iteration)
metric.subset <- list()
for(i in seq_along(metric.names)) {
  print(i)
  metric.subset[[metric.names[i]]] <- basic[these.rows, grep(paste0(metric.names[i], "_"), colnames(basic))]
  colnames(metric.subset[[metric.names[i]]]) <- gsub(metric.names[i], "FR", colnames(metric.subset[[metric.names[i]]]) )  
}

# make a figure showing distribution of metrics for all cities
pdf(paste0(folder_out, "compare obs vs exp metrics - SMC12.pdf"), height=7, width=7, useDingbats=F)
par(mfrow=c(3,3), mar=c(2.5,2,1,0.5), las=1)
for(k in metric.names) {
  maxVal <- max(max(metric.subset[[k]]), max(data.metrics[[k]][ , "value"]))
  minVal <- min(min(metric.subset[[k]]), min(data.metrics[[k]][ , "value"]))
  hist(as.numeric(as.character(metric.subset[[k]][1,])), col=transp("grey15", 0.1), border=NA, 
       breaks=10, xlim=c(minVal-(maxVal-minVal)/5, maxVal*1.1), ylim=c(0, 10), main="", xlab="")
  
  for(i in 2:50) {
    hist(as.numeric(as.character(metric.subset[[k]][i,])), col=transp("grey15", 0.1), 
         breaks=10, border=NA, add=T)  #seq(0, max(metric.subset[[k]]+1000), 500)
  }
  
  opar <- par(lwd=2)
  par(lwd=3)
  hist(data.metrics[[k]][ , "value"], main="", border=transp("firebrick", 0.99),
       breaks=10, ylab="", bty="n", lwd=4, add=T)
  par(opar)
  legend("topright", legend=k, cex=1.7, col="white", lwd=2, bty="n")
}
dev.off()


}

# # # # # # # # # # # # # # # # # # # # # # # # 
# Make plots of model metric vs observed
boxplot(metric.subset[[1]])
points(data$metrics[grep("mean", data$metrics$name), "value"], pch=19, col="firebrick")

compare.this <- function(metric.name) {
  plot.this <- metric.name
  boxplot(metric.subset[[plot.this]][match(pop.order$Code, colnames(metric.subset[[plot.this]]))], main=metric.name)
  points(data.metrics[[plot.this]]$value[match(pop.order$Code, data.metrics[[plot.this]]$location)], pch=19, col="firebrick")
}

par(mfrow=c(5,1))
for(i in metric.names) {
  compare.this(i)
}

regions <- colnames(metric.subset[["mean"]])
this.region <- regions[1]
plot(x=rep(data.metrics[["mean"]][data.metrics[["mean"]]$location==this.region, "value"], 500),
       y=metric.subset[["mean"]][, which(colnames(metric.subset[["mean"]])==this.region)],
       pch=19, col=transp("dodgerblue3", 0.6), xlim=c(0,8000), ylim=c(0,8000)) 
for(i in 2:19) {
  this.region <- regions[i]
  points(x=rep(data.metrics[["mean"]][data.metrics[["mean"]]$location==this.region, "value"], 500),
       y=metric.subset[["mean"]][, which(colnames(metric.subset[["mean"]])==this.region)],
       pch=19, col=transp("dodgerblue3", 0.6)) 
}

# make a plot of observed vs expected
pdf(paste0(folder_out, "compare obs v exp - SMC12.pdf"), height=9, width=9, useDingbats=F)
par(mfrow=c(3,3), mar=c(3,3.25,0.5,0.5), las=1, cex=0.9, lwd=0.75)
for(j in metric.names) {
  print(j)
  maxVal <- max(max(metric.subset[[j]]), max(data.metrics[[j]][ , "value"]))
  minVal <- min(min(metric.subset[[j]]), min(data.metrics[[j]][ , "value"]))
  expansion <- (maxVal-minVal)/25
  
  boxplot(at=data.metrics[[j]][data.metrics[[j]]$location==this.region, "value"],
          x=metric.subset[[j]][, which(colnames(metric.subset[[j]])==this.region)],
          pch=19, cex=0.5, col=transp("dodgerblue3", 0.6), xlim=c(minVal-(maxVal-minVal)/5, maxVal*1.1), 
          ylim=c(minVal-(maxVal-minVal)/5, maxVal*1.1),
          xlab="", ylab="", boxwex=expansion, axes=F) 
  for(i in 2:19) {
    this.region <- regions[i]
    boxplot(at=data.metrics[[j]][data.metrics[[j]]$location==this.region, "value"],
            x=metric.subset[[j]][, which(colnames(metric.subset[[j]])==this.region)],
            pch=19, cex=0.5, col=transp("dodgerblue3", 0.6), add=T, axes=F,boxwex=expansion) 
  }
  axis(1)
  axis(2)
  abline(0,1, col=transp("grey25", 0.8), lty=3)
  #legend("topright", legend=j, cex=1.7, col="white", lwd=2, bty="n")
}
dev.off()
plot(x=rep(data.metrics[["mean"]][data.metrics[["mean"]]$location=="FR_A", "value"], 500),
     y=metric.subset[["mean"]][, which(colnames(metric.subset[["mean"]])=="FR_A")])
compare.this("q50")

# extract parameter results 
par.names <- c("R0", "CCI", "e_zero", "CJ", "h")
pars <- basic[, par.names]
par.results <- c()
for(i in par.names) {
  par.results = rbind(par.results, 
                      c(mean=mean(pars[, i]), 
                        median=quantile(pars[,i], 0.5), 
                        lower=quantile(pars[,i], 0.025),
                        upper=quantile(pars[,i], 0.975)) )
}
rownames(par.results) <- par.names

par.results["e_zero", ] <- 10^par.results["e_zero", ]
par.results["CJ", ] <- 1/par.results["CJ", ]
# print the means, medians etc.
print(par.results)

### plot the posteriors
priors <- list()
priors[["R0"]] <- c(1,8)
priors[["CCI"]] <- c(0, 1)
priors[["e_zero"]] <- c(0,5)
priors[["CJ"]] <- c(0,1)
priors[["h"]] <- c(0,1)

# plot posterior parameter distributions
pdf(paste0(folder_out, "posterior par dists.pdf"), height=8, width=6.5, useDingbats=F)
par(mfrow=c(5,1), mar=c(2.5,4,0.5,0.5), las=1, cex=0.9)
for(i in par.names) {
  hist(pars[, i], xlim=priors[[i]], breaks=20, main="", xlab="", col=transp("dodgerblue", 0.7))
}
dev.off()

##### plot degree distribution
deg.dist <- c(0, 3, 45, 160, 393, 715, 883, 989, 897, 752,
                      697, 755, 856, 1085, 1224, 1452, 1578, 1622, 1711, 1584,
                      1514, 1355, 1209, 990, 805, 597, 477, 353, 232, 177,
                      126, 90, 69, 54, 36, 29, 21, 16, 10, 5,
                      8, 13, 7, 9, 3, 1, 3, 3, 2, 5,
                      0, 1, 0, 0, 0, 1)
sum(deg.dist)
pdf(paste0(folder_out, "degree dist.pdf"), height=4, width=4, useDingbats=F)
par(mfrow=c(1,1), mar=c(4,4,2,0.5), las=1, cex=0.9)
barplot(height=deg.dist/25622, main="", names.arg=seq(1:56), 
        ylab="frequency", xlab="degree",
     col=transp("dodgerblue", 0.7), bty="n", border="dodgerblue4")
dev.off()

## plot example immunity
pdf(paste0(folder_out, "example inf history.pdf"), height=10, width=6, useDingbats=F)
par(mfrow=c(4,1), mar=c(2.5,5,1,0.5), las=1, cex=0.9)
# no infections occur
strain1 <- c(1,1, 0.6, 0.6, 0.6, 0.36, 0.36, 0.216, 0.216, 0.216, 0.216, 
             0.216, 0.1296, 0.1296)
strain2 <- c(1,1,1, 1, 0.6, 0.36, 0.36,0.36,0.36, 0.216, 0.216,
             0.1296, 0.1296, 0.1296)
plot(strain1, type="l", ylim=c(0,1), col=transp("dodgerblue3", 0.6), 
     lwd=3, ylab="protection \nfrom infection", xlab="", xlim=c(0, 14))
points(strain1, col=transp("dodgerblue3", 0.6), cex=1.4, pch=15)

plot(strain1, type="l", ylim=c(0,1), col=transp("dodgerblue3", 0.6), 
     lwd=3, ylab="protection \nfrom infection", xlab="seasons", xlim=c(0, 14))
lines(strain2, type="l", ylim=c(0,1), col=transp("grey45", 0.6), 
      lwd=3, ylab="protection \nfrom infection")
points(strain2, col=transp("grey45", 0.6), cex=1.4, pch=15)

# infections occur
strain1 <- c(1,1, 0.6, 0.6, 0.6, 0.36, 1, 0.6, 0.6, 0.6, 0.36, 0.36, 1, 1)
strain2 <- c(1,1,1,1, 0.6, 0.36, 0.36, 1, 1, 0.6, 0.6, 
             0.36, 0.36, 0.36)
plot(strain1, type="l", ylim=c(0,1), col=transp("dodgerblue3", 0.6), 
     lwd=3, ylab="protection \nfrom infection", xlab="seasons", xlim=c(0, 14))
points(strain1, col=transp("dodgerblue3", 0.6), cex=1.4, pch=15)

plot(strain1, type="l", ylim=c(0,1), col=transp("dodgerblue3", 0.6), 
     lwd=3, ylab="protection \nfrom infection", xlab="seasons", xlim=c(0, 14))
lines(strain2, type="l", ylim=c(0,1), col=transp("grey45", 0.6), 
      lwd=3, ylab="protection \nfrom infection")
points(strain2, col=transp("grey45", 0.6), cex=1.4, pch=15)

dev.off()





# 
# plot(data$metrics[data$metrics$name=="mean_A", "value"])
# 
# par(mfrow=c(2,1))
# 
# pars <- dbGetQuery(db, "SELECT * from parameters")
# head(pars)
# plot(pars$R0)
# plot(pars$CCI)
# plot(pars$e_zero)
# plot(pars$CJ)
# plot(pars$h)
# 
# 
# dbGetQuery(db , "SELECT * from metrics limit 10")
# 
# dbGetQuery(db , "SELECT * from jobs limit 10")
# dbGetQuery(db , "SELECT * from parameters limit 10")
# 
# dbGetQuery(db, "select *")
#   
# db = dbConnect(drv, "./yucatan_interventions-sm3.sqlite")
# 
# basic <- dbGetQuery(db, 'select mild_EF, severe_EF, vac, catchup, target, catchup_to, M.* 
#                       from parameters P, metrics M, jobs J
#                       where P.serial = M.serial 
#                             and P.serial = J.serial 
#                             and status = \'D\'
#                             and vec_control = 0
#                             and vac_waning = 0
#                             and vec_scenario = 0')