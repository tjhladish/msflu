setwd("~/Dropbox/LSHTM/Multi_season_flu_transfer/")

# folder with outputs in
folder_data <- "~/Dropbox/LSHTM/Multi_season_flu_transfer/"

# folder with repo in 
folder_repo <- "~/Documents/Influenza/msflu/Analysis/"

# outputs folder
folder_out <- "~/Sync/LSHTM/Collaboration/Multiseason flu/"

# # # # # # # # # # # # # # # # # # # # # # # # 
# read in json file
require(jsonlite)
data <- fromJSON(txt=paste0(folder_data, "abc_msflu.json"))
data$metrics
hist(data$metrics[grep("mean", data$metrics$name), "value"])
hist(data$metrics[grep("skew", data$metrics$name), "value"])
hist(data$metrics[grep("q0", data$metrics$name), "value"])
hist(data$metrics[grep("mc", data$metrics$name), "value"])

# # # # # # # # # # # # # # # # # # # # # # # # 
# read in results
require("RSQLite")
drv = dbDriver("SQLite")
db = dbConnect(drv, "./msmsc_flu.sqlite")
dbListTables(db)

# # # # # # # # # # # # # # # # # # # # # # # # 
# read in all metrics
metrics <- dbGetQuery(db, "SELECT * from metrics")

# these are the names of the metrics and the rows we want
metric.names <- c("mean", "q0", "q25", "q50", "q75", "q100", "sd", "skew", "mc")
these.rows <- 60001:70000

# subset the parts of the metrics that we want (last SMC iteration)
metric.subset <- list()
for(i in seq_along(metric.names)) {
  print(i)
  metric.subset[[metric.names[i]]] <- metrics[these.rows, grep(metric.names[i], colnames(metrics))]
}

# # # # # # # # # # # # # # # # # # # # # # # # 
# read in location lookup
lookup <- read.csv(paste0(folder_out, "20160628 - Name_Pop_Lookup.csv"), header=T, na.strings="-")
# read in epi data
source(paste0(folder_repo, "Read in epidata.R"))

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