##set-up new pixels without ground-truth data (usually site level data) to predict the biophysical model outputs
##This script processes the landsat data (cloud and tide filter, get date/location etc.), estimates pixel flood proportion, 
##calculates atmospherically corrected Land Surface Temperature (LST), 
##calculates soil temperature from LST and Air temperature, and estimate the date of green up
###Version for running on local laptop for an expanded prediction area, based on Hladik et al. 2013 habitat map of spartina alterniflora 

##libraries
library(R.utils); library(caret); library(tidyverse);library(RColorBrewer); library(pls);
library(randomForest); library(data.table); library(raster); library(RcppRoll);

##output and save new models and data sets?? Set to True/False below. 
##Note that True will overwrite the existing models and output datasets, so you'll want a backup if you don't intend this
save_model<-FALSE
save_model<-TRUE


##are we testing a small subset of the data?
small_test<-F

##custom functions and models we need for this, uncomment if not already in memory
source("~/git/berm2/functions/calc_indices_l8.R")
source("~/git/berm2/functions/soil_temp_from_lst_and_air_lags.R")
super_model <- readRDS("output/random_forest_flood_model_l8.rds")
print(super_model)



##read in elevation data as a spatial raster object
#dem<-raster("/home/kyle/Documents/Documents/UT/Sapelo/ModDEM_ascii/mod_dem.txt")



###Read in Landsat observations for site level data to predict; This is created from the Google earth Engine Script
l8_2014 <- read_csv("/home/kyle/Documents/Documents/UT/Sapelo/GEE/landsat/hladik_spal_450/L8_app_plots_hladik_2014_nocloudsmask.csv")
l8_2015 <- read_csv("/home/kyle/Documents/Documents/UT/Sapelo/GEE/landsat/hladik_spal_450/L8_app_plots_hladik_2015_nocloudsmask.csv")
l8_2016 <- read_csv("/home/kyle/Documents/Documents/UT/Sapelo/GEE/landsat/hladik_spal_450/L8_app_plots_hladik_2016_nocloudsmask.csv")
l8_2017 <- read_csv("/home/kyle/Documents/Documents/UT/Sapelo/GEE/landsat/hladik_spal_450/L8_app_plots_hladik_2017_nocloudsmask.csv")
l8_2018 <- read_csv("/home/kyle/Documents/Documents/UT/Sapelo/GEE/landsat/hladik_spal_450/L8_app_plots_hladik_2018_nocloudsmask.csv")
l8_2019 <- read_csv("/home/kyle/Documents/Documents/UT/Sapelo/GEE/landsat/hladik_spal_450/L8_app_plots_hladik_2019_nocloudsmask.csv")
l8_2020 <- read_csv("/home/kyle/Documents/Documents/UT/Sapelo/GEE/landsat/hladik_spal_450/L8_app_plots_hladik_2020_nocloudsmask.csv")
l8_2021 <- read_csv("/home/kyle/Documents/Documents/UT/Sapelo/GEE/landsat/hladik_spal_450/L8_app_plots_hladik_2021_nocloudsmask.csv")
l8_2022 <- read_csv("/home/kyle/Documents/Documents/UT/Sapelo/GEE/landsat/hladik_spal_450/L8_app_plots_hladik_2022_nocloudsmask.csv")
l8_2023 <- read_csv("/home/kyle/Documents/Documents/UT/Sapelo/GEE/landsat/hladik_spal_450/L8_app_plots_hladik_2023_nocloudsmask.csv")
l9_2014_2023 <- read_csv("/home/kyle/Documents/Documents/UT/Sapelo/GEE/landsat/hladik_spal_450/L9_app_plots_hladik_2014_2023_nocloudsmask.csv")

#rbind
landsat<-rbind(l8_2014, l8_2015, l8_2016, l8_2017, l8_2018, l8_2019, l8_2020, l8_2021, l8_2022, l8_2023, l9_2014_2023)
rm(l8_2014, l8_2015, l8_2016, l8_2017, l8_2018, l8_2019, l8_2020, l8_2021, l8_2022, l8_2023, l9_2014_2023)

# some things calculate across years, so it'll be best to combine all data, then split by pixel to do calculations
cut <- T
ceiling(length(unique(landsat$id))/1000)
cut_no <- 16


if(cut==T) {
  newpixels<-landsat[landsat$id %in% unique(landsat$id)[(1*(cut_no*1000)-999):(cut_no*1000)],]
  }
rm(landsat)


if(small_test==T) {
  newpixels<-newpixels[newpixels$fid %in% unique(newpixels$fid)[1:150],]
}


attrib <- read_csv("/home/kyle/Documents/Documents/UT/Sapelo/Mapping/Landsat pixel outline clipped/centroids_Landsat_pixel_hladik_spal_450_mod_dem_geo.csv")
attrib <- attrib %>% 
  dplyr::select(-c(VALUE, VALUE_2, fid_2, '_mean', '_median', area, fid)) %>%
  rename(utm_east = xcoord_2, utm_north = ycoord_2, long = xcoord, lat = ycoord)

newpixels <- inner_join(newpixels, attrib, by = "id")

##process the landsat data to create needed time, date, and location columns
newpixels$time<-substr(newpixels$date,12,19)
newpixels$date<-as.Date(substr(newpixels$date,1,10))
newpixels$year<-as.numeric(format(newpixels$date, "%Y"))
newpixels$doy<-as.numeric(format(newpixels$date, "%j"))
newpixels$mo<-as.numeric(format(newpixels$date, "%m"))
newpixels$pix<-paste(round(newpixels$utm_east,0), round(newpixels$utm_north,0), sep="_")


##filter the landsat 8 data if needed through the pixel_qa landsat mask; 
##note that the current Landsat 8 Google Earth Engine script only returns cloud filtered data
##see https://landsat.usgs.gov/landsat-surface-reflectance-quality-assessment
newpixels$QA_PIXEL<-intToBin(newpixels$QA_PIXEL)
##"000010" from the right means no fill, yes clear, no water, no cloud shadow
#kyle: edited this to reflect QA_PIXEL as opposed to pixel_qa in diff L8 products. ##"01000000" from the right means no fill, no dilated clouds, no cirrus, no cloud, no cloud shadow, no snow, yes clear, no water
## no water above is via Landsat 8's pixel qa mask, and misses a lot of marsh flooding in mixed pixels, we'll handle this below
newpixels$qa_good<-ifelse(str_sub(newpixels$QA_PIXEL,-8,-1)=="01000000",T,F) 
#let's keep in the water yes for the flood_time variable we'll make later
newpixels$qa_good<-ifelse(str_sub(newpixels$QA_PIXEL,-8,-1)=="11000000",T,newpixels$qa_good) 

## subset to only good pixels as indicated by the pixel_qa mask
newpixels<-newpixels[newpixels$qa_good==T,]

##pull lat longs from .geo column
x<-strsplit(newpixels$".geo", "\\[")
x<-sapply(x, "[", 2)
x<-strsplit(x, "\\]")
x<-sapply(x, "[", 1)
x<-strsplit(x, ",")
newpixels$long<-sapply(x, "[", 1)
newpixels$lat<-sapply(x, "[", 2)
rm(x)

##remove columns we don't need now from Landsat 8 Google Earth Engine data
newpixels<-dplyr::select(newpixels,-c('system:index', .geo, QA_RADSAT, SR_QA_AEROSOL, ST_ATRAN, ST_CDIST, ST_DRAD, ST_EMIS, ST_EMSD, ST_QA, ST_TRAD, ST_URAD, QA_PIXEL, qa_good))
head(newpixels)

## Scaling done in GEE (L8L2C2T1)
newpixels$b1<-newpixels$SR_B1
newpixels$b2<-newpixels$SR_B2
newpixels$b3<-newpixels$SR_B3
newpixels$b4<-newpixels$SR_B4
newpixels$b5<-newpixels$SR_B5
newpixels$b6<-newpixels$SR_B6
newpixels$b7<-newpixels$SR_B7
##double-check this, I think GEE takes care of thermal scaling too
##this is the needed scalar for these thermal bands: B10 and B11
newpixels$b10<-newpixels$ST_B10

##clean up columns to remove the ones we now don't need
newpixels<-dplyr::select(newpixels, -c(SR_B1,SR_B2,SR_B3,SR_B4,SR_B5,SR_B6,SR_B7,ST_B10))

##negative reflectance values can occur at scene edges and should be removed
newpixels$b1<-ifelse(newpixels$b1<0,NA, newpixels$b1)
newpixels$b2<-ifelse(newpixels$b2<0,NA, newpixels$b2)
newpixels$b3<-ifelse(newpixels$b3<0,NA, newpixels$b3)
newpixels$b4<-ifelse(newpixels$b4<0,NA, newpixels$b4)
newpixels$b5<-ifelse(newpixels$b5<0,NA, newpixels$b5)
newpixels$b6<-ifelse(newpixels$b6<0,NA, newpixels$b6)
newpixels$b7<-ifelse(newpixels$b7<0,NA, newpixels$b7)
##filter to non-NA values
newpixels<-newpixels[is.na(newpixels$b1)==F,]; newpixels<-newpixels[is.na(newpixels$b7)==F,]
newpixels<-newpixels[is.na(newpixels$b5)==F,];newpixels<-newpixels[is.na(newpixels$b6)==F,]

###calculate a set of standard vegetation and spectral reflectance indices for landsat 8 data using our custom function
indices<-calc_index_l8(newpixels)
newpixels<-cbind(newpixels, indices)
newpixels$pix<-paste(round(newpixels$utm_east,0), round(newpixels$utm_north,0))
newpixels$ndmi<-(newpixels$b1-newpixels$b6)/(newpixels$b1+newpixels$b6) 
newpixels$pheno<-(newpixels$b4-newpixels$b6)/(newpixels$b4+newpixels$b6) ##from TMII paper
newpixels$doy<-as.numeric(format(newpixels$date, "%j"))

##calculate lst with custom function, first we need to merge in the atmospheric correction parameters
#atm$date <- as.Date(atm$date)
#newpixels<-left_join(newpixels, atm, by="date")
#newpixels$lst<-calc_lst(newpixels)
## actually don't because the updated GEE L8 product has surface temp instead of brightness temp
newpixels$lst<-newpixels$b10-273.15


##create rolling pheno mean
##trick to get all dates, just create the full dates dataframe and merge it to the missing dates dataframe
y<-data.frame(date=seq.Date(from=min(newpixels$date, na.rm=T), to=max(newpixels$date, na.rm=T), by=1))
y<-expand.grid(date=y$date, pix=unique(newpixels$pix))

####crashing

newpixels<-full_join(newpixels,y)
newpixels<-arrange(newpixels, pix,date)

##this is an even window roll, because there's an observation for every day
newpixels<-newpixels %>% group_by(pix) %>% mutate(pheno2= roll_mean(pheno, n=30, na.rm=T,align="center", fill=NA))
newpixels<-newpixels[!is.na(newpixels$b2),]
newpixels<-newpixels[!is.na(newpixels$b5),]

newpixels<-arrange(newpixels, date)
newpixels$year<-as.numeric(format(newpixels$date, "%Y"))

##Predict pixelwise flooding with the flood model, set up variables in the right format first
newpixels$mo<-format(newpixels$date, "%m")
newpixels$mo<-ifelse(is.na(newpixels$mo),format(newpixels$landsat_date, "%m"), newpixels$mo)
newpixels$mo<-factor(newpixels$mo, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))

newpixels$flooded <- predict(super_model, newdata=newpixels)
table(newpixels$date, newpixels$flooded)

newpixels$flooded<-as.numeric(newpixels$flooded)
newpixels$flooded<-ifelse(newpixels$flooded==1,0,newpixels$flooded)
newpixels$flooded<-ifelse(newpixels$flooded==2,1,newpixels$flooded)

##summarise flooded observation data by pixel, to get proportion of flooded observations by pixel
flood<- newpixels %>% dplyr::group_by(pix) %>% dplyr::summarise(flood_time= sum(as.numeric(flooded), na.rm=T)
                                                                /length(flooded[!is.na(flooded)]))

newpixels<-merge(newpixels, flood, by="pix", all.x=T)

##now that we know the flood_time for each pixel, we can subset to just dry observations for downstream analysis
newpixels<-newpixels[newpixels$flooded==0,]

##this last gets rid of dates with few satelitte obs.
if(small_test==F) {
  x<-table(newpixels$date);x<-x[x>100]
  newpixels<-newpixels[newpixels$date %in% as.Date(names(x)),]
}

##get elevation data
#plot.locations<-data.frame(newpixels[,c("utm_east", "utm_north")])
#plot.dem<-raster::extract(dem,plot.locations)
#newpixels$elevation<-plot.dem
newpixels$elevation <- newpixels$`_mean`


newpixels<-ungroup(newpixels)
newpixels<-newpixels[!is.na(newpixels$date),]
newpixels$mo<-as.numeric(format(newpixels$date, "%m"))
Sys.time()
newpixels.heat.l8<-setup_heat(newpixels, gce=T)
Sys.time()

if (save_model==T){
  write_csv(newpixels, paste0("/home/kyle/Documents/Documents/UT/Sapelo/berm_output/gce_pixels_landsat_processed_split", cut_no, ".csv"))
  write_csv(newpixels.heat.l8, paste0("output/gce_pixels_landsat8_greenup_via_lst_and_air_lags_split", cut_no, ".csv"))
}
