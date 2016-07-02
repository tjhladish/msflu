source("pairs.panels.R")
require('rjson')
require("RSQLite")

# parse json
abc=fromJSON(file="~/Dropbox/Multi_season_flu_transfer/abc_msflu.json"); 
nmet=length(abc$metrics); 
npar=length(abc$parameters); 
#factors=c(rep(1,nmet), rep(2,npar)); factors %o% factors
abc_metrics_values=rep(0,length(abc$metrics)); for(i in 1:length(abc$metrics)) { abc_metrics_values[i]=abc$metrics[[i]]$value; } 
abc_metrics_names=rep(0,length(abc$metrics)); for(i in 1:length(abc$metrics)) { abc_metrics_names[i]=abc$metrics[[i]]$name; } 
# we want short names if they're defined
for(i in 1:length(abc$metrics)) { if('short_name' %in% names(abc$metrics[[i]])) abc_metrics_names[i]=abc$metrics[[i]]$short_name; } 

# read in db
drv = dbDriver("SQLite")
db = dbConnect(drv, "~/Dropbox/Multi_season_flu_transfer/msmsc_flu.sqlite")
post = dbGetQuery(db, 'select J.*, P.*, M.* from jobs J, parameters P, metrics M where J.serial = P.serial and J.serial = M.serial and smcSet = (select max(smcSet) from jobs) and posterior > -1')
extra_serials = which(names(post)== 'serial')[-1] 
post = post[,-c(extra_serials)]
post = subset(post, select=-c(startTime, duration, attempts, seed, serial, status))

proto = post[1,]
proto[1,] <- NA
proto[1,abc_metrics_names] <- abc_metrics_values
dm = rbind(post, proto)
dm$sim = factor(ifelse(is.na(dm$smcSet), F, T))
dm[dm$sim==F, 4:(3+npar)] = apply(dm[dm$sim==T,4:(3+npar)], 2, median)

#colors <- c(sim = '#00000012', par = 'purple', met = 'orange')
#chars <- c(sim = 20, par = '|', met = '—')
############## end



#names(dm)[4:13] = c('EF', 'Mos move', 'Daily intros', 'Num mos', 'Beta', 'Mean', 'Median', 'Stdev', 'Max', 'Skew')
#names(dm)[14] = "Autocorr"

#pdf("pairs-a.pdf", width=16, height=16)
#png("pairs-a.png", width=1800, height=1800, res=150)
#pairs.panels(dm[,4:22], dm[,23], npar, nmet, points.col='#00000012', box.col='black', box.lwd=0.5, gap=0.5)
#dev.off()

dm$test = dm$R0*(1-dm$CCI)*dm$h
#dm$e_zero = 10**dm$e_zero
#names(dm)[c(4,5,6,7,8)] = c('log10(Daily intros)', 'Mosq/loc', 'Median', 'Std deviation')
png("par-pairs-set9.png", width=1800, height=1340, res=200)
pairs.panels(cbind(dm[,c(4,5,6,7,8)], dm$test), dm$sim, npar=6, nmet=0, points.col='#00000012', box.lwd=1, gap=0.5)
dev.off()

png("par-hist.png", width=2000, height=2500, res=300)
par(mfrow=c(5,1), mar=c(4,4,0.5,1))
for(i in 4:8) {
    hist(dm[,i], nclass = 20, xlim = c(abc$parameters[[i-3]]$par1,abc$parameters[[i-3]]$par2), main="", xlab=abc$parameters[[i-3]]$name)
}
dev.off()