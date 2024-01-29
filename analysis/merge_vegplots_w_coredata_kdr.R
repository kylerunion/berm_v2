#this script prepares the field data

library(readxl)
library(xlsx)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(lubridate)
library(dplyr)
library(tidyr)
library(readr)
library(tidyverse)
library(zoo)

#data are available at https://gce-lter.marsci.uga.edu/public/app/dataset_details.asp?accession=PLT-GCET-2308

# this is for setting up the new folder in case the pwd is the project folder
setwd("for_new_github/")


#set dates for sample id
#monthly from 2013-06 to present for flux A data
start <- as.Date("2013-06-01")
end <- as.Date("2023-04-01")
start.date <- seq(start, end, by = "month")
rm(start)
rm(end)
#during these dates there were multiple surveys
extra.dates <- as.Date(c("2017-11-15", "2019-01-18", "2019-07-20", "2019-07-29", "2019-08-28", "2019-11-20", "2022-01-30", "2020-03-30"))
start.date <- c(start.date, extra.dates)
start.date <- sort(start.date)

#core_data
#my data
core <- read_csv("data/PLT-GCET-2308_rootcore.csv", skip = 2)
#removes non-data rows
core <- core[-c(1,2),]
#renames columns for ease
names(core) <- c("date", "site", "plot", "species", "live_stem_count", "live_stem_g", "dead_stem_g", "TotalBG_0_10_g", "Rhizomes_0_10_g", "TotalBG_10_30_g", "Rhizomes_10_30_g")
core$date <- as.Date(core$date, format = "%Y-%m-%d")
core <- core %>%
  dplyr::select(date, plot, species, live_stem_g, dead_stem_g, TotalBG_0_10_g, Rhizomes_0_10_g, TotalBG_10_30_g, Rhizomes_10_30_g, site)
core[,4:9] <- sapply(core[4:9], as.numeric)
core$ag.total <- core$live_stem_g + core$dead_stem_g
core <- core[which(core$plot != "dc10"),]
core <- core %>%
  dplyr::arrange(ymd(core$date))
core <- core[!(duplicated(core)), ]
core$date<-as.Date(core$date, format="%m-%d-%Y")
core$survey<-findInterval(core$date, as.Date(start.date))
core$sampleid<-paste(core$survey,core$plot,sep="")

##veg data
#my data
veg <- read_csv("data/PLT-GCET-2308_vegplot.csv", skip = 2)
#removes non-data rows
veg <- veg[-c(1,2),]
#renames columns for ease
names(veg) <- c("date", "site", "plot", "species", "vegplot_area_cm2", "LAI", "chl")
veg$date <- as.Date(veg$date, format = "%Y-%m-%d")
veg[,5:7] <- sapply(veg[5:7], as.numeric)
veg$survey<-findInterval(veg$date, as.Date(start.date))
veg$sampleid<-paste(veg$survey,veg$plot,sep="")
veg$survey <- NULL

#stem height data
#my data
hts <- read_csv("data/PLT-GCET-2308_stemheights.csv", skip = 2)
#removes non-data rows
hts <- hts[-c(1,2),]
#renames columns for ease
names(hts) <- c("date", "site", "plot", "species", "zone", "plot_type", "quadrat_area", "plant_height", "allometric_1", "allometric_2")
hts[,7:10] <- sapply(hts[7:10], as.numeric)
hts$date <- as.Date(hts$date, format = "%Y-%m-%d")

## spartina zone was assigned by height of 75th percentile
# this is already done in the GCE online data, but can be changed here if desired
hts <- hts %>%
  group_by(plot) %>%
  mutate(zone = quantile(plant_height, 0.75))
# short/medium/tall height distinction (cm) based on Rolando et al 2022 https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01187-7
hts$zone[hts$zone > 80] <- "T"
hts$zone[hts$zone <= 80 & hts$zone >= 50] <- "M"
hts$zone[hts$zone < 50] <- "S"
## playing around with what percentile of plant height should be used to assign zone - look at the distribution
table(hts$zone)
#biomass calculations
# this is already done in the GCE online data, but can be changed here if desired
hts$allometric_1 <- ifelse(hts$zone == "T", 
                                  2.718281828459^(-6.011+1.719*log(hts$plant_height)), 
                                  2.718281828459^(-7.311+1.988*log(hts$plant_height)))
hts$allometric_2 <- ifelse(hts$zone == "T", 
                                      2.718281828459^(-4.5987+1.4317*log(hts$plant_height)), 
                                      2.718281828459^(-4.30856+1.24807*log(hts$plant_height)))

hts$survey<-findInterval(hts$date, as.Date(start.date))
hts$sampleid<-paste(hts$survey,hts$plot,sep="")

#leafN data
#my data
leafN <- read_csv("data/PLT-GCET-2308_leafn.csv", skip = 2)
#removes non-data rows
leafN <- leafN[-c(1,2),]
#renames columns for ease
names(leafN) <- c("date", "site", "plot", "percent_n", "percent_c", "cn_ratio")
leafN[,4:6] <- sapply(leafN[4:6], as.numeric)
leafN <- leafN[!duplicated(leafN),]
leafN$date<-as.Date(leafN$date, format="%Y-%m-%d")
leafN$survey<-findInterval(leafN$date, as.Date(start.date))
leafN$sampleid<-paste(leafN$survey,leafN$plot,sep="")
leafN <- leafN %>%
  group_by(sampleid) %>%
  mutate(plot = plot, date = date, percentC = mean(percent_c), percentN = mean(percent_n), site = site, survey = survey) %>%
  distinct(sampleid, .keep_all = T)

#misc data
#this is the same as 'PLT-GCET-2308_plots.csv' except it has plots grouped by Landsat pixel
rtk <- read.csv("data/biomass_plots.csv")
rtk <- rtk %>%
  select(plot, elevation, lat, long)

#count number of stems from quadrat and paste to veg
quadrat_hts <- subset(hts, plot_type == "quadrat")
quadrat_hts <- quadrat_hts %>%
  dplyr::select(date, site, plot, quadrat_area, plant_height)
colnames(quadrat_hts) <- c("date", "site", "plot", "vegplot_area_cm2", "plant_height")
quadrat_counts <- quadrat_hts %>%
  group_by(date, site, plot, vegplot_area_cm2) %>%
  dplyr::summarize(spal_live_stem_count = n())
quadrat_counts$vegplot_area_cm2 <- as.numeric(quadrat_counts$vegplot_area_cm2)

#sampleid
quadrat_counts$date<-as.Date(quadrat_counts$date, format="%Y-%m-%d")
quadrat_counts$survey<-findInterval(quadrat_counts$date, as.Date(start.date))
quadrat_counts$sampleid<-paste(quadrat_counts$survey,quadrat_counts$plot,sep="")
quadrat_counts$survey <- NULL

veg <- merge(veg, quadrat_counts, by = c("date", "sampleid", "site", "plot", "vegplot_area_cm2"))

core_hts <- subset(hts, plot_type == "core")
colnames(core_hts)[4] <- "Core_Area"
core_hts$Core_Area <- pi*(7.62/2)^2 # 3inch core
core_hts <- core_hts %>%
  dplyr::select(date, site, plot, Core_Area, plant_height)
colnames(core_hts) <- c("date", "site", "plot", "core_area_cm2", "plant_height")
core_counts <- core_hts %>%
  group_by(date, site, plot, core_area_cm2) %>%
  dplyr::summarize(core_live_stem_count = n())
core_counts$date<-as.Date(core_counts$date, format="%Y-%m-%d")
core_counts$survey<-findInterval(core_counts$date, as.Date(start.date))
core_counts$sampleid<-paste(core_counts$survey,core_counts$plot,sep="")
core_counts$survey <- NULL
core_counts$date <- NULL

veg <- merge(veg, core_counts, by = c("sampleid", "site", "plot"), all = T)
# this below identifies if you have two separate dates for one sampleid, which fudges everything up. if rows appear here, take the later date and add it to the 'extra.dates' vector above
veg[duplicated(veg$sampleid),]

#sum aboveground_biomass from quadrat (g) and paste to veg
ag <- hts %>% 
  group_by(sampleid, plot_type) %>%
  dplyr::summarize(plot_ag_1 = sum(allometric_1), plot_ag_2 = sum(allometric_2))
quadrat_ag <- subset(ag, plot_type == "quadrat")
quadrat_ag <- quadrat_ag %>%
  dplyr::select(sampleid, plot_ag_1, plot_ag_2)
veg <- merge(veg, quadrat_ag, by = c("sampleid"))
veg[duplicated(veg$sampleid),]

# these are the SPAD to chlorophyll calibration equations, based on the calibrations done in O'Connell et al. 2021 and Runion et al. 2024
# already done in these data
#veg <- veg %>% 
#  mutate(chlorophyll = ifelse(site == "fluxa" | date < as.Date("2021-01-01"),
#    -0.21 + 0.0232 * spad,
#    -0.32 + 0.0378 * spad))

#veg calcs
veg$stemsm2<-(veg$spal_live_stem_count/veg$vegplot_area_cm2)*(100*100)
#core stuff
area <- pi*(7.62/2)^2 # 3inch core
core$area <- rep(area,nrow(core))  
core$bg.total <- core$TotalBG_0_10_g + core$TotalBG_10_30_g
core$ag.total <- ifelse(is.na(core$ag.total) == T, core$live_stem_g + core$dead_stem_g, core$ag.total)
core$rhiz.total <- core$Rhizomes_0_10_g + core$Rhizomes_10_30_g
core$rhiz.shoot <- core$rhiz.total / core$ag.total
core$rhiz.green.shoot <- core$rhiz.total / core$live_stem_g
core$rhizm2.core <- ((core$rhiz.total) / (core$area)) * (100*100)
core$rootgreenshoot <- core$bg.total / core$live_stem_g
core$rootshoot <- core$bg.total / core$ag.total
core$bgm2.core <- (core$bg.total / (core$area)) * (100*100)
core$agm2.core <- (core$ag.total / core$area) * (100*100)

leafN <- leafN %>%
  ungroup() %>%
  dplyr::select(-plot, -date, -site, -survey)

veg<-merge(veg, core, by= c("sampleid", "plot", "site"), all.x=T, all.y=T)
veg$date <- as.Date(ifelse(is.na(veg$date.x), veg$date.y, veg$date.x))
veg <- veg %>%
  ungroup() %>%
  dplyr::select(-date.x, -date.y)
veg<-merge(veg, leafN, by=c("sampleid"), all.x=T, all.y=T)
veg<-merge(veg, rtk, by="plot", all.x=T, all.y=T)
veg <- veg[!duplicated(veg),]

# this below identifies if you have two separate dates for one sampleid, which fudges everything up. if rows appear here, take the later date and add it to the 'extra.dates' vector above
veg[duplicated(veg$sampleid),c(1,2,33)]

veg$bgm2.areascaled <- (veg$bg.total / veg$area) * (100*100)
veg$agm2.core.areascaled <- (veg$ag.total / veg$area) * (100*100)
veg$rhizm2.areascaled <- veg$rhiz.total / veg$area * (100*100)
veg$bgm2.core.stemscaled <- (veg$bg.total / veg$core_live_stem_count) * (veg$spal_live_stem_count / veg$vegplot_area_cm2) * (100*100)
veg$agm2.core.stemscaled <- (veg$ag.total / veg$core_live_stem_count) * (veg$spal_live_stem_count / veg$vegplot_area_cm2) * (100*100)

#applies mean of core data from each zone to all plots
veg$zone <- ifelse(veg$site == "fluxa", substr(veg$plot, 1, 1), NA)
cols <- c("bgm2.core.stemscaled", "agm2.core.areascaled")
test <- veg %>%
  group_by(date, zone) %>%
  mutate_at(cols, na.aggregate)

veg %>% group_by(site) %>% 
  summarise(CHL=mean(chl, na.rm=T), LAI=mean(LAI, na.rm=T),
            bg=mean(bgm2.core.stemscaled, na.rm=T), ag=mean(agm2.core.stemscaled, na.rm=T))

veg$agm2.allom <- veg$plot_ag_1 / veg$vegplot_area_cm2 * 100*100
veg$agm2.allom.new <- veg$plot_ag_2 / veg$vegplot_area_cm2 * 100*100

veg$bgm2.allometric.rootshoot<-veg$agm2.allom*veg$rootshoot
veg$bgm2.allometric.rootgreenshoot<-veg$agm2.allom*veg$rootgreenshoot
veg$bgm2.allo.new.rootshoot<-veg$agm2.allom.new*veg$rootshoot
veg$bgm2.allo.new.rootgreenshoot<-veg$agm2.allom.new*veg$rootgreenshoot
veg$rhizm2.allometric.rootshoot<-veg$agm2.allom*veg$rhiz.shoot
veg$rhizm2.allometric.rootgreenshoot<-veg$agm2.allom*veg$rhiz.green.shoot
veg$rhizm2.allo.new.rootshoot<-veg$agm2.allom.new*veg$rhiz.shoot
veg$rhizm2.allo.new.rootgreenshoot<-veg$agm2.allom.new*veg$rhiz.green.shoot
veg$rhizm2.core.stemscaled<-(veg$rhiz.total/veg$spal_live_stem_count)*(veg$spal_live_stem_count/veg$vegplot_area_cm2)*(100*100)

veg$AGCg<-veg$percentC/100*veg$agm2.allom
veg$wholeplantCg<-veg$percentC/100*(veg$agm2.allom+veg$bgm2.core.stemscaled)
veg$wholeplantCg<-veg$percentC/100*(veg$agm2.core.areascaled+veg$bgm2.areascaled)
plot(veg$date[veg$AGCg>0], veg$wholeplantCg[veg$AGCg>0], xlab="Month", ylab="C in all tissue (g)")

write.table(veg, "output/core_vegplot_combined.csv", row.names=FALSE, sep=",")
