###xgboost belowground biomass
library(tidyverse)
library(dplyr)
library(RcppRoll)
library(sp)
library(lubridate)
library(rgdal)

##process daymet data
##make date equivalent to estimated veg data dates

#2022 data, daymet switched from V4 to V4R1 in 2022, so 2022 data are not in GEE

#berm_daymet is a file with the berm calibration plots joined to the nearest daymet measurement location (1 km x 1 km grid)
#so most berm plots within the same site have the same daymet_joined coordinates
#this joining was done in qgis
berm_daymet <- read_csv("data/berm_daymet.csv")
berm_daymet_coords <- unique(berm_daymet[c("lat", "long")])

#you can get daymet data from google earth engine and summarize monthly there, but daymet switched from V4 to V4R1 in 2022, so 2022 data are not in GEE
# in general, daymet data become available for a calendar year the following spring (~Apr)
library(daymetr)
berm_dm <- data.frame(matrix(ncol = 11, nrow = 0))

for(i in seq(1:nrow(berm_daymet_coords))){
  tmp <- download_daymet(site = "daymet_dl",
                              lat = berm_daymet_coords[i,1],
                              lon = berm_daymet_coords[i,2],
                              start = 2013,
                              end = 2022,
                              internal = TRUE)
  daymet_dl <- tmp$data
  daymet_dl$lat <- as.numeric(berm_daymet_coords[i,1])
  daymet_dl$long <- as.numeric(berm_daymet_coords[i,2])
  berm_dm <- rbind(daymet_dl, berm_dm)
}
berm_dm$date <- as.Date(berm_dm$yday -1, origin = "2022-01-01")
berm_dm$mo <- month(berm_dm$date)
berm_dm <- berm_dm %>% group_by(lat, long, mo) %>%
  summarise(dayl = mean(dayl..s.), prcp = mean(prcp..mm.day.), srad = mean(srad..W.m.2.), swe = mean(swe..kg.m.2.), tmax = mean(tmax..deg.c.), tmin = mean(tmin..deg.c.), vp = mean(vp..Pa.))
berm_dm$yr <- 2022
berm_dm$date <- as.Date(with(berm_dm, paste(yr,mo,"01", sep = "-")), "%Y-%m-%d")
daymet <- inner_join(berm_dm, berm_daymet, by = c("lat", "long"))

daymet$site<-substr(daymet$plot, 1,2)
daymet$site[!(daymet$site %in% c("dc", "du", "fb", "fr","hc", "ns", "sk", "ug"))]<- "fa"

#create UTM coordinates from lat long
cord.dec = SpatialPoints(cbind(daymet$long, daymet$lat), 
                         proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32617"))
# add coordinates to UTM list
daymet$utm_east=cord.UTM@coords[,1]
daymet$utm_north=cord.UTM@coords[,2] 
rm(cord.dec)
rm(cord.UTM)

#GEE data is monthly average, this calculates sum by month
daymet <- rename(daymet, dayl_mean = dayl, prcp_mean = prcp, srad_mean = srad, tmax_mean = tmax, tmin_mean = tmin, vp_mean = vp)
daymet$dayl_sum <- daymet$dayl_mean * days_in_month(daymet$date)
daymet$prcp_sum <- daymet$prcp_mean * days_in_month(daymet$date)
daymet$srad_sum <- daymet$srad_mean * days_in_month(daymet$date)
daymet$tmax_sum <- daymet$tmax_mean * days_in_month(daymet$date)
daymet$tmin_sum <- daymet$tmin_mean * days_in_month(daymet$date)
daymet$vp_sum <- daymet$vp_mean * days_in_month(daymet$date)

##remove variables we don't need and create estimates summarized by month and plot
daymet<-dplyr::select(daymet, -lat, -long, -rtk_m, -date)

##assign the monthly estimate to the middle of the month
daymet$date<-as.Date(paste(daymet$yr, daymet$mo, "15", sep="-"))

#insert pix
plots<-read.csv("data/biomass_plots.csv", header=T)
plots <- plots %>% 
  dplyr::select(plot, pix)
daymet <- merge(daymet, plots, by = "plot")
daymet$key<-paste(daymet$date, daymet$pix)

#average by pix

#test <- daymet %>% 
 
# group_by(key)
 # summarise(across(where(is.numeric), ~ list(summary(.))))


#daymet <- daymet %>%
 # group_by(key) %>%
  #select_if(is.numeric) %>%
  #summarise_all(mean)

daymet<-dplyr::select(daymet, -yr, -mo, -plot,-site)

##calculate PAR from solar radiation via standard formula
daymet$par_sum<-2.114*daymet$srad_sum
daymet$par_tot<-daymet$dayl_sum*daymet$par_sum/1000000

daymet <- daymet %>% distinct()

##create derived variables that we might use in the Machine Learning Model
daymet<-arrange(daymet, pix, date)
daymet<-daymet %>% dplyr::group_by( pix) %>%
  dplyr::mutate(lagg_tmax_1=lag(tmax_mean),
                diff_tmax_1=c(NA, diff(tmax_mean, lag=1)),
                diff_tmax_2=c(rep(NA, 2), diff(tmax_mean, lag=2)),
                diff_tmax_3=c(rep(NA, 3),diff(tmax_mean, lag=3)),
                diff_tmax_4=c(rep(NA, 4), diff(tmax_mean, lag=4)),
                diff_tmax_5=c(rep(NA, 5),diff(tmax_mean, lag=5)),
                diff_tmax_6=c(rep(NA, 6), diff(tmax_mean, lag=6)),
                roll_tmax_2=RcppRoll::roll_meanr(tmax_mean, n=2, fill=NA, na.rm=T),
                roll_tmax_3=RcppRoll::roll_meanr(tmax_mean, n=3, fill=NA, na.rm=T),
                roll_tmax_4=RcppRoll::roll_meanr(tmax_mean, n=4, fill=NA, na.rm=T),
                roll_tmax_5=RcppRoll::roll_meanr(tmax_mean, n=5, fill=NA, na.rm=T),
                roll_tmax_6=RcppRoll::roll_meanr(tmax_mean, n=6, fill=NA, na.rm=T),
                lagg_tmin_1=lag(tmin_mean),
                diff_tmin_1=c(NA, diff(tmin_mean, lag=1)),
                diff_tmin_2=c(rep(NA, 2), diff(tmin_mean, lag=2)),
                diff_tmin_3=c(rep(NA, 3),diff(tmin_mean, lag=3)),
                diff_tmin_4=c(rep(NA, 4), diff(tmin_mean, lag=4)),
                diff_tmin_5=c(rep(NA, 5),diff(tmin_mean, lag=5)),
                diff_tmin_6=c(rep(NA, 6), diff(tmin_mean, lag=6)),
                roll_tmin_2=RcppRoll::roll_meanr(tmin_mean, n=2, fill=NA, na.rm=T),
                roll_tmin_3=RcppRoll::roll_meanr(tmin_mean, n=3, fill=NA, na.rm=T),
                roll_tmin_4=RcppRoll::roll_meanr(tmin_mean, n=4, fill=NA, na.rm=T),
                roll_tmin_5=RcppRoll::roll_meanr(tmin_mean, n=5, fill=NA, na.rm=T),
                roll_tmin_6=RcppRoll::roll_meanr(tmin_mean, n=6, fill=NA, na.rm=T),
                lagg_prcp_1=lag(prcp_mean),
                diff_prcp_1=c(NA, diff(prcp_mean, lag=1)),
                diff_prcp_2=c(rep(NA, 2), diff(prcp_mean, lag=2)),
                diff_prcp_3=c(rep(NA, 3),diff(prcp_mean, lag=3)),
                diff_prcp_4=c(rep(NA, 4), diff(prcp_mean, lag=4)),
                diff_prcp_5=c(rep(NA, 5),diff(prcp_mean, lag=5)),
                diff_prcp_6=c(rep(NA, 6), diff(prcp_mean, lag=6)),
                roll_prcp_2=RcppRoll::roll_meanr(prcp_mean, n=2, fill=NA, na.rm=T),
                roll_prcp_3=RcppRoll::roll_meanr(prcp_mean, n=3, fill=NA, na.rm=T),
                roll_prcp_4=RcppRoll::roll_meanr(prcp_mean, n=4, fill=NA, na.rm=T),
                roll_prcp_5=RcppRoll::roll_meanr(prcp_mean, n=5, fill=NA, na.rm=T),
                roll_prcp_6=RcppRoll::roll_meanr(prcp_mean, n=6, fill=NA, na.rm=T),
                lagg_vp_1=lag(vp_mean),
                diff_vp_1=c(NA, diff(vp_mean, lag=1)),
                diff_vp_2=c(rep(NA, 2), diff(vp_mean, lag=2)),
                diff_vp_3=c(rep(NA, 3),diff(vp_mean, lag=3)),
                diff_vp_4=c(rep(NA, 4), diff(vp_mean, lag=4)),
                diff_vp_5=c(rep(NA, 5),diff(vp_mean, lag=5)),
                diff_vp_6=c(rep(NA, 6), diff(vp_mean, lag=6)),
                roll_vp_2=RcppRoll::roll_meanr(vp_mean, n=2, fill=NA, na.rm=T),
                roll_vp_3=RcppRoll::roll_meanr(vp_mean, n=3, fill=NA, na.rm=T),
                roll_vp_4=RcppRoll::roll_meanr(vp_mean, n=4, fill=NA, na.rm=T),
                roll_vp_5=RcppRoll::roll_meanr(vp_mean, n=5, fill=NA, na.rm=T),
                roll_vp_6=RcppRoll::roll_meanr(vp_mean, n=6, fill=NA, na.rm=T),
                lagg_srad_1=lag(srad_mean),
                diff_srad_1=c(NA, diff(srad_mean, lag=1)),
                diff_srad_2=c(rep(NA, 2), diff(srad_mean, lag=2)),
                diff_srad_3=c(rep(NA, 3),diff(srad_mean, lag=3)),
                diff_srad_4=c(rep(NA, 4), diff(srad_mean, lag=4)),
                diff_srad_5=c(rep(NA, 5),diff(srad_mean, lag=5)),
                diff_srad_6=c(rep(NA, 6), diff(srad_mean, lag=6)),
                roll_srad_2=RcppRoll::roll_meanr(srad_mean, n=2, fill=NA, na.rm=T),
                roll_srad_3=RcppRoll::roll_meanr(srad_mean, n=3, fill=NA, na.rm=T),
                roll_srad_4=RcppRoll::roll_meanr(srad_mean, n=4, fill=NA, na.rm=T),
                roll_srad_5=RcppRoll::roll_meanr(srad_mean, n=5, fill=NA, na.rm=T),
                roll_srad_6=RcppRoll::roll_meanr(srad_mean, n=6, fill=NA, na.rm=T),
                lagg_par_1=lag(par_tot),
                diff_par_1=c(NA, diff(par_tot, lag=1)),
                diff_par_2=c(rep(NA, 2), diff(par_tot, lag=2)),
                diff_par_3=c(rep(NA, 3),diff(par_tot, lag=3)),
                diff_par_4=c(rep(NA, 4), diff(par_tot, lag=4)),
                diff_par_5=c(rep(NA, 5),diff(par_tot, lag=5)),
                diff_par_6=c(rep(NA, 6), diff(par_tot, lag=6)),
                roll_par_2=RcppRoll::roll_sumr(par_tot, n=2, fill=NA, na.rm=T),
                roll_par_3=RcppRoll::roll_sumr(par_tot, n=3, fill=NA, na.rm=T),
                roll_par_4=RcppRoll::roll_sumr(par_tot, n=4, fill=NA, na.rm=T),
                roll_par_5=RcppRoll::roll_sumr(par_tot, n=5, fill=NA, na.rm=T),
                roll_par_6=RcppRoll::roll_sumr(par_tot, n=6, fill=NA, na.rm=T)
  )

daymet$key<-paste(daymet$date, daymet$pix)
daymet<-ungroup(daymet)
daymet<-dplyr::select(daymet, -date,  -pix)
#this is likely too big to push to github so save it locally
write.csv(daymet, "/home/kyle/Documents/Documents/UT/Sapelo/GEE/daymet/berm/output/daymet_bermplots", row.names=F)