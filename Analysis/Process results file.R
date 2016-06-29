
# set paths for who is working
if(Sys.info()[["user"]] %in% c("Rosalind", "eideregg")) {
  setwd("~/Dropbox/LSHTM/Multi_season_flu_transfer/")
  # folder with outputs in
  folder_data <- "~/Dropbox/LSHTM/Multi_season_flu_transfer/"
  # folder with repo in 
  folder_repo <- "~/Documents/Influenza/msflu/Analysis/"
  # outputs folder
  folder_out <- "~/Sync/LSHTM/Collaboration/Multiseason flu/"  
} else {
  # Tom's paths
}

# # # # # # # # # # # # # # # # # # # # # # # # 
# read in location lookup
lookup <- read.csv(paste0(folder_out, "20160628 - Name_Pop_Lookup.csv"), header=T, na.strings="-")
lookup$Code <- gsub("-", "_", lookup$Code)
# read in epi data
source(paste0(folder_repo, "Read in epidata.R"))


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

par(mfrow=c(3,3), mar=c(3,4,1,0))
for(i in metric.names) {
  compare.xy(i)
}

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
                          and smcSet = 7') 
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
}


# # # # # # # # # # # # # # # # # # # # # # # # 
# Make plots of model metric vs observed
boxplot(metric.subset[[1]])
points(data$metrics[grep("mean", data$metrics$name), "value"], pch=19, col="firebrick")

compare.this <- function(metric.name) {
  plot.this <- metric.name
  boxplot(metric.subset[[plot.this]], main=metric.name)
  points(data$metrics[grep(plot.this, data$metrics$name), "value"], pch=19, col="firebrick")
}

par(mfrow=c(5,1))
for(i in metric.names) {
  compare.this(i)
}
compare.this("q50")



















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