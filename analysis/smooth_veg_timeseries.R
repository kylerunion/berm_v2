library(tidyverse)
library(dplyr)
library(zoo)
library(readr)
library(ggplot2)
library(ggpmisc)

# this is for setting up the new folder in case the pwd is the project folder
setwd("for_new_github/")

##need "flood" model
#super_model <- readRDS("output/random_forest_flood_model_l8.rds")

#Load the field data
veg<-read.csv("output/core_vegplot_combined.csv", header=T)
veg<-arrange(veg, date)
veg$date<-as.Date(veg$date)
veg$year<-as.numeric(format(veg$date, "%Y"))
veg$mo<-as.numeric(format(veg$date, "%m"))

##smooth veg data series to estimate info on all valid landsat dates
#i<-1; j<-1; k<-1
var<-c("bgm2.areascaled","agm2.core.areascaled", "bgm2.core.stemscaled",
       "agm2.core.stemscaled", "agm2.allom", "bgm2.allometric.rootshoot",
       "bgm2.allometric.rootgreenshoot", "agm2.allom.new", "bgm2.allo.new.rootshoot",
       "bgm2.allo.new.rootgreenshoot",  "rhizm2.allometric.rootshoot",
       "rhizm2.allometric.rootgreenshoot", "rhizm2.allo.new.rootshoot",
       "rhizm2.allo.new.rootgreenshoot","rhizm2.areascaled","rhizm2.core.stemscaled",
       "chl", "LAI", "percentN")
sites <- c("dc", "du", "fb", "fr", "hc", "ns", "sk", "ug")
fluxa <- c("s", "m", "t")
loc <- sort(c(as.vector(outer(sites, seq(1:9), paste0)), as.vector(outer(fluxa, seq(1:8), paste0))))
dates<-seq.Date(as.Date("2013-05-01"), as.Date("2023-03-27"), by= "1 week")
#dates<-dates[as.numeric(format(dates, "%m"))>4&as.numeric(format(dates, "%m"))<11]
out<-data.frame(matrix(ncol=length(var)+2, nrow=length(dates)*length(loc)))
names(out)<-c("date", "plot", var)
loop<-1

for (j in seq_along(loc)){
  for (i in seq_along(var)){
    x<-veg[veg$plot==loc[j], c("date", var[i])]
    x$date<-as.numeric(x$date)
    x<-x[complete.cases(x),]
    
    ifelse(nrow(x)>1,{
      smoothed <- approx(x$date,x[,var[i]], xout=as.numeric(dates),
                         method="linear"
      )
      smoothy<-smoothed
    num<-length(smoothy$x)
    out[loop:(loop+num-1),"plot"]<-rep(loc[j], num)
    out[loop:(loop+num-1),"date"]<-as.Date(smoothy$x, origin = "1970-01-01")
    out[loop:(loop+num-1), var[i]]<-(smoothy$y)
},{
smoothy<-smoothed
num<-length(smoothy$x)
out[loop:(loop+num-1),"plot"]<-rep(loc[j], num)
out[loop:(loop+num-1),"date"]<-as.Date(smoothy$x, origin = "1970-01-01")
out[loop:(loop+num-1), var[i]]<-(NA)
  })}
  loop<-loop+num
}
out$date<-as.Date(out$date, origin="1970-01-01")
out$year<-as.numeric(format(out$date, "%Y"))
out$mo<-as.numeric(format(out$date, "%m"))
out$site<-substr(out$plot, 1,2)
out$site[!(out$site %in% sites)]<- "fa"

##remove unsampled years 
fluxa_all <- c(as.vector(outer(fluxa, seq(1:6), paste0)))
fluxa_tall_extra <- c("t7", "t8")
fluxa_all <- c(fluxa_all, fluxa_tall_extra)
old_berm_sites <- c("sk", "ug", "fb")
old_berm_plots <- c(as.vector(outer(old_berm_sites, seq(1:9), paste0)))
old_berm_all <- c(old_berm_plots, fluxa_all)
rm(old_berm_plots)
new_berm_sites <- c("dc", "du", "fr", "hc", "ns", "sk", "ug", "fb")
new_berm_all <- c(as.vector(outer(new_berm_sites, seq(1:9), paste0)))
rm(new_berm_sites)
fluxa_full_record <- out %>%
  filter(plot %in% fluxa_all)
old_berm_record <- out %>%
  filter(plot %in% old_berm_all & date > as.Date("2016-05-01") & date < as.Date("2016-11-01"))
new_berm_record <- out %>%
  filter(plot %in% new_berm_all & date > as.Date("2021-05-01"))
out <- rbind(fluxa_full_record, old_berm_record, new_berm_record)
out <- out[!(duplicated(out)), ]

#applies mean of core data from each zone to all plots
## this was the method used in O'Connell et al. 2021 but not in Runion et al. 2024
#out$zone <- substr(out$plot, 1, 1)
#cols <- c("bgm2.areascaled", "agm2.core.areascaled")
#out <- out %>%
#  group_by(date, zone) %>%
#  mutate_at(cols, na.aggregate)

names(out)[names(out) == 'chl'] <- "chlorophyll"
out$chlorophyll[out$year<2016]<-NA
out$chlorophyll[out$year==2016&out$mo<5]<-NA
out$LAI[out$year==2016&out$mo<5]<-NA
out$LAI[out$year<2016]<-NA
out<-out[as.Date(out$date)>"2013-06-11",]
out<-out[!is.na(out$date),]

#drops t7 and t8 before they were established
del <- c("t7", "t8")
out <- subset(out, !(plot %in% del & date < as.Date("2017-07-01")))

# removing occurrences of interpolated variables when they shouldn't be
var<-c("rhizm2.allometric.rootshoot",
       "rhizm2.allometric.rootgreenshoot", "rhizm2.allo.new.rootshoot",
       "rhizm2.allo.new.rootgreenshoot","rhizm2.areascaled","rhizm2.core.stemscaled"
)

out[out$date<"2016-05-03",var]<-NA
out[out$date>"2016-10-29"&out$date<"2017-06-30",var]<-NA
out[out$date>"2017-07-18"&out$date<"2017-10-03",var]<-NA

#var<-"bgm2.core.stemscaled"
#out[out$date<"2017-05-22",var]<-NA

var<-"percentN"
out[out$date<"2014-04-24",var]<-NA
out[out$date>"2014-10-03"&out$date<"2015-04-14",var]<-NA
out[out$date>"2015-09-23"&out$date<"2016-05-03",var]<-NA
out[out$date>"2015-09-23"&out$date<"2016-05-03",var]<-NA

#t7 and t8 only have one percentN measurement and this screws up the interpolation. a better solution would be to prevent interpolation among gaps but for now:
out[out$plot %in% c("t7", "t8"),var] <- NA

var <- "LAI"
out[out$date<"2014-04-24",var]<-NA
out[out$date > "2016-10-29" & out$date < "2018-08-22" & out$site == "fa", var] <- NA
out[out$date > "2018-11-06" & out$date < "2019-05-29" & out$site == "fa", var] <- NA
out[out$date > "2019-12-10" & out$date < "2020-09-30" & out$site == "fa", var] <- NA
out[out$date > "2020-09-30" & out$date < "2021-06-23" & out$site == "fa", var] <- NA


var <- "chlorophyll"
out[out$date<"2014-04-24",var]<-NA
out[out$date > "2016-10-25" & out$site == "fa", var] <- NA




write.csv(out, "output/smoothed_veg_timeseries_weekly.csv", row.names = F)
