getwd()
setwd("~/Sync/Austin_Work_2/Multi-season flu/")

# read in raw data
flu <- read.csv("20151207 - Sentinelles weekly per 100k.csv", row.names=1,
                na.strings="-")
head(flu)

# print missing weeks for all regions
pdf("Missing weeks in all regions.pdf", height=16, width=10)
par(las=1, mar=c(3,4,0.5,0.5), mfrow=c(11,2))
for(i in 1:22) {
  plot(flu[,i], type="l", xlab="", ylab="inc per 100000")
  mtext(text="week index from 1984-44",side = 1, line = 2, cex=0.6)
  abline(v=which(is.na(flu[,i])==T), col="red", lwd=0.3)
  legend("topright", legend=c(colnames(flu)[i], 
                              length(which(is.na(flu[,i])==T)),
                              paste0(round(length(which(is.na(flu[,i])==T))/nrow(flu)*100,2), "%") ), 
                              bty="n", col="white", lwd=2)
  lines(flu[,i])
}
dev.off()


## read in tom's version which splits year and week
flu <- read.csv("Sentinelles weekly per 100k_tjh.csv", stringsAsFactors=F)
#remove anything between week 20 and week 40
flu.season <- flu[(flu$week >= 40 | flu$week <= 20),]
# cut off 1984
flu.season <- flu.season[!(flu.season$year=="1984"),]
# cut off 1985 & 1986 (more incomplete)
flu.season <- flu.season[!(flu.season$year=="1985"),]
flu.season <- flu.season[!(flu.season$year=="1986" & flu.season$week<40),]

#tidy tail, by removing 2015
flu.season <- flu.season[!(flu.season$year=="2015" & flu.season$week>20),]

# remove provinces with missing data, FR-L (Limousin), FR-T and FR-H (Corse)
flu.season <- flu.season[,!(colnames(flu.season) %in% c("FR.L", "FR.H", "FR.T"))]
  
#label seasons. year is the year with week 40 in it, so 1986, goes week 40 1986 to week 201987
flu.season$season <- 0
for(k in 1986:2015) {
  for(j in 1:nrow(flu.season)) {
    if(flu.season$year[j]==k && flu.season$week[j]>=40){
      flu.season$season[j] <- k 
    }
    if(flu.season$year[j]==k && flu.season$week[j]<=20){
      flu.season$season[j] <- k-1 
    }
  }
}

#print image of missing weeks
pdf("Missing weeks in all regions - data for study.pdf", height=16, width=10)
par(las=1, mar=c(3,4,0.5,0.5), mfrow=c(11,2))
flu <- flu.season
for(i in 3:21) {
  plot(x=seq(1:nrow(flu)), y=flu[,i], type="l", xlab="", ylab="inc per 100000")
  mtext(text="week index from 1986, weeks 40-20 only",side = 1, line = 2, cex=0.6)
  abline(v=which(is.na(flu[,i])==T), col="red", lwd=0.3)
  legend("topright", legend=c(colnames(flu)[i], 
                              paste0("num = ", length(which(is.na(flu[,i])==T))),
                              paste0(round(length(which(is.na(flu[,i])==T))/nrow(flu)*100,2), "%") ), 
         bty="n", col="white", lwd=2)
  lines(flu[,i])
}
dev.off()

# set NAs to 0
flu.season[is.na(flu.season)] <- 0

#write out
write.csv(flu.season, "Sentinelles weekly per 100k clean.csv", row.names=F)

#aggregate by season
library(data.table)
flu.dt <- data.table(flu.season)
flu.dt <- melt(flu.dt, id.vars="season", measure.vars=grep("FR",colnames(flu.dt)))
flu <- flu.dt[, list(total=sum(value)), by=c("season", "variable")]

colnames(flu) <- c("season", "region", "total")
write.csv(flu, "Sentinelles per 100k yearly long.csv", row.names=F)

flu.wide <- dcast(flu, formula=season ~ region)
write.csv(flu.wide, "Sentinelles per 100k yearly wide.csv", row.names=F)
