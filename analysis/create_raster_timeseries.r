###Script to make pretty biomass maps
library(plotKML);library(raster);library(maptools);library(rasterVis)
library(latticeExtra);library(tidyverse); library(rgdal)

####Set the working directory; change path to your working directory####


habitat<-readOGR(dsn = "data/habitatmap/", 
                 layer = "DupClass_MLC")
habitat<-habitat[habitat$CLASS_NAME=="SM"|habitat$CLASS_NAME=="ST"|habitat$CLASS_NAME=="SS",]

habitat2<-readOGR(dsn = "data/habitatmap/", 
                 layer = "centroids_Landsat_pixel_hladik_spal_450_mod_dem_wgs_geo")


####Set up newpixelsetation predictions#####
#newpixels<-read_csv("/home/kyle/Documents/Documents/UT/Sapelo/berm_output/xgb_predicted_biophysical_landsat_allpixels_2014.csv")
split1 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split1.csv")
split2 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split2.csv")
split3 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split3.csv")
split4 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split4.csv")
split5 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split5.csv")
split6 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split6.csv")
split7 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split7.csv")
split8 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split8.csv")
split9 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split9.csv")
split10 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split10.csv")
split11 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split11.csv")
split12 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split12.csv")
split13 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split13.csv")
split14 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split14.csv")
split15 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split15.csv")
split16 <- read_csv("output/xgb2.0_predicted_bgbiomass_landsat8_hladik_split16.csv")



newpixels<-rbind(split1, split2, split3, split4, split5, split6, split7, split8, split9, split10, split11, split12, split13, split14, split15, split16)
rm(split1, split2, split3, split4, split5, split6, split7, split8, split9, split10, split11, split12, split13, split14, split15, split16)

newpixels<-data.frame(newpixels)
newpixels$date<-as.Date(newpixels$date)
newpixels<-data.frame(newpixels)
newpixels<-newpixels[!is.na(newpixels$date),]


##make mask file as a raster file, gives faster masking, use the first masked raster for this
dates<-unique(newpixels$date[newpixels$date>as.Date("2013-09-01")]) 
newpixels1<-newpixels[newpixels$date==dates[1],]
newpixels1<-newpixels1[!is.na(newpixels1$utm_east), ]
newpixels1<-newpixels1[!is.na(newpixels1$predbg.xgb), ]
#newpixels1<-newpixels1[!is.na(newpixels1$predag.allom.l8), ]
proj <- "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
coordinates(newpixels1) <- ~utm_east+utm_north ##converts newpixels to a shape file; cool 
proj4string(newpixels1) <- proj  # define projection system of our data
agarea <- vect2rast(newpixels1, cell.size=30, fname="predbg.xgb", 
                    file.name="output/test6.tif")
agarea<-raster("output/test6.tif")
agarea<-mask(agarea, habitat2)
writeRaster(agarea, file = paste0("results/biomassmaps/mask.tif"), 
            bylayer=TRUE, format="GTiff", overwrite=T)

mask<-raster("results/biomassmaps/mask.tif")
mask<-reclassify(mask, c(0,10000, 1))

###BG BIOMASS####
##create separate BG geotiffs by date
dates<-unique(newpixels$date[newpixels$date>as.Date("2013-09-01")]) 
for (i in 1:length(dates)){
  newpixels1<-newpixels[newpixels$date==dates[i],]
  newpixels1<-newpixels1[!is.na(newpixels1$utm_east), ]
  newpixels1<-newpixels1[!is.na(newpixels1$predbg.xgb), ]
  if(nrow(newpixels1)>1){
    proj <- "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
    coordinates(newpixels1) <- ~utm_east+utm_north ##converts newpixels to a shape file; cool 
    proj4string(newpixels1) <- proj  # define projection system of our data
    newpixels1@bbox <- as.matrix(extent(mask))
    agarea <- vect2rast(newpixels1, cell.size=30, fname="predbg.xgb", bbox = as.matrix(extent(mask)),
                        file.name=paste0("results/biomassmaps/bg_rasters/predbg_", format(dates[i], "%Y%m%d"), ".tif"))
  }
}

##read in all the files as rasters, stack them together and mask all in one step
setwd("results/biomassmaps/bg_rasters")
files<-list.files(pattern = ".tif$")
ag<-stack( files)
agarea<-crop(ag, mask)

#agarea<-mask(ag, mask)

plot(agarea[[1]])

##write out the masked files
names<-names(agarea)
#writeRaster(x = agarea, filename = paste0("0_", names), bylayer=T, format = "GTiff", overwrite = T)


##plot the masked files
#setwd("results/biomassmaps/final_bg_rasters")
files<-list.files(pattern = ".tif$")
ag<-stack( files)
rasterNames<-as.Date(str_split(names(ag), "_", simplify = T)[,2], format="%Y%m%d")
x<-which(rasterNames>as.Date("2013-12-31") &rasterNames<as.Date("2023-12-15")&(format(rasterNames, "%m") %in% c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")))
x<-which((format(rasterNames, "%Y") %in% c("2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023"))&(format(rasterNames, "%m") %in% c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")))
rasterNames<-rasterNames[x]
ag<-ag[[x]]
dat.label<-expression(paste("belowground biomass (g ", m^-2,")"))
my.at=seq(0,max(newpixels$predbg.xgb, na.rm=T)+100, by=50)
label.at=seq(0,max(newpixels$predbg.xgb, na.rm=T)+100, by=200)
mycolorkey <- list(at=my.at, ## where the colors change
                   labels=list(
                     at=label.at ## where to print labels
                   ))
label.at=seq(0,max(newpixels$predbg.xgb, na.rm=T)+100, by=200)
mycolorkey <- list(at=my.at, ## where the colors change
                   labels=list(
                     at=label.at, ## where to print labels
                     cex=1.25
                   ))


rasterVis::levelplot(ag, col.regions = rev(terrain.colors(255)), cuts=254, 
                     margin=FALSE, at=my.at, colorkey=mycolorkey ,main=list(dat.label,cex=2),
                     names.attr=as.character(rasterNames),layout=c(4,5),
                     scales=list(draw=FALSE ), # remove axes labels & ticks
                     xlim=c(mask@extent[1], mask@extent[2]),ylim=c(mask@extent[3],mask@extent[4]),
                     panel=function(...){
                       panel.levelplot(...)
                     })

dev.print(file="../../../results/biomassmaps/bg_xgb_stack.png",
          device=png, bg="white", width=720, height=730) 

#create gif
trellis.device('png', file = 'for_gif/myPic%03d.png')
rasterVis::levelplot(ag, col.regions = rev(terrain.colors(255)), cuts=254, 
                     margin=FALSE, at=my.at, colorkey=mycolorkey ,main=list(dat.label,cex=1.5),
                     names.attr=as.character(rasterNames),layout=c(1,1),
                     scales=list(draw=FALSE ), # remove axes labels & ticks
                     xlim=c(mask@extent[1], mask@extent[2]),ylim=c(mask@extent[3],mask@extent[4]),
                     panel=function(...){
                       panel.levelplot(...)
                     })

## list file names and read in
library(magick)
setwd("results/biomassmaps/bg_rasters")
imgs <- list.files("for_gif/", full.names = TRUE, pattern = ".png")
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 5)

## view animated image
img_animated

## save to disk
image_write(image = img_animated,
            path = "for_gif/gif/pred_bg.gif")

###AG BIOMASS####
##create separate BG geotiffs by date
setwd(paste(root_path,"/UGA/Projects/seagrant", sep=""))
dates<-unique(newpixels$date[newpixels$date>as.Date("2013-05-01")]) 
for (i in 1:length(dates)){
  newpixels1<-newpixels[newpixels$date==dates[i],]
  newpixels1<-newpixels1[!is.na(newpixels1$utm_east), ]
  newpixels1<-newpixels1[!is.na(newpixels1$predag.allom.l8), ]
  if(nrow(newpixels1)>1){
    proj <- "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
    coordinates(newpixels1) <- ~utm_east+utm_north ##converts newpixels to a shape file; cool 
    proj4string(newpixels1) <- proj  # define projection system of our data
    bgarea <- vect2rast(newpixels1, cell.size=30, fname="predag.allom.l8", 
                        file.name=paste0("results/biomassmaps/final_ag_rasters/predag_", format(dates[i], "%Y%m%d"), ".tif"))
  }
}
##read in all the files as rasters, stack them together and mask all in one step
setwd("/home/jessica/Documents/UGA/Projects/seagrant/results/biomassmaps/final_ag_rasters")
files<-list.files(pattern = ".tif$")
bg<-stack( files)
bgarea<-mask(bg, mask)

##write out the masked files
names<-names(bgarea)
writeRaster(bgarea, names, bylayer=T, format="GTiff", overwrite=T)


##plot the masked files
setwd("/home/jessica/Documents/UGA/Projects/seagrant/results/biomassmaps/final_ag_rasters")
files<-list.files(pattern = ".tif$")
bg<-stack( files)
rasterNames<-as.Date(str_split(names(bg), "_", simplify = T)[,2], format="%Y%m%d")
x<-which(rasterNames>as.Date("2013-12-31") &rasterNames<as.Date("2019-01-01")&(format(rasterNames, "%m") %in% c("03", "09", "09", "12")))
x<-which((format(rasterNames, "%Y") %in% c("2014", "2016", "2018"))&(format(rasterNames, "%m") %in% c("03", "06", "09", "12")))
rasterNames<-rasterNames[x]
bg<-bg[[x]]
dat.label<-expression(paste("aboveground biomass (g ", m^-2,")"))
my.at=seq(0,max(newpixels$predag.allom.l8, na.rm=T)+25, by=25)
label.at=seq(0,max(newpixels$predag.allom.l8, na.rm=T)+25, by=100)
mycolorkey <- list(at=my.at, ## where the colors change
                   labels=list(
                     at=label.at, ## where to print labels
                     cex=1.25
                   ))

levelplot(bg, col.regions = rev(terrain.colors(255)), cuts=254, 
          margin=FALSE, at=my.at, colorkey=mycolorkey ,main=list(dat.label,cex=2),
          names.attr=as.character(rasterNames),layout=c(4,3),
          scales=list(draw=FALSE ), # remove axes labels & ticks
          xlim=c(mask@extent[1], mask@extent[2]),ylim=c(mask@extent[3],mask@extent[4]),
          panel=function(...){
            panel.levelplot(...)
            #panel.text(471500, 3480200, cex=1.5,
            #            month)
          })
dev.print(file="../../../results/biomassmaps/ag_xgb_stack.png",
          device=png, bg="white", width=720, height=730) 

####Foliar N####
##create separate geotiffs by date
setwd(paste(root_path,"/UGA/Projects/seagrant", sep=""))
dates<-unique(newpixels$date[newpixels$date>as.Date("2013-09-01")]) 
for (i in 1:length(dates)){
  newpixels1<-newpixels[newpixels$date==dates[i],]
  newpixels1<-newpixels1[!is.na(newpixels1$utm_east), ]
  newpixels1<-newpixels1[!is.na(newpixels1$predn.l8), ]
  if(nrow(newpixels1)>1){
    proj <- "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
    coordinates(newpixels1) <- ~utm_east+utm_north ##converts newpixels to a shape file; cool 
    proj4string(newpixels1) <- proj  # define projection system of our data
    bgarea <- vect2rast(newpixels1, cell.size=30, fname="predn.l8", 
                        file.name=paste0("results/biomassmaps/final_N_rasters/predn_", format(dates[i], "%Y%m%d"), ".tif"))
  }
}
##read in all the files as rasters, stack them together and mask all in one step
setwd("/home/jessica/Documents/UGA/Projects/seagrant/results/biomassmaps/final_N_rasters")
files<-list.files(pattern = ".tif$")
bg<-stack( files)
bgarea<-mask(bg, mask)

##write out the masked files
names<-names(bgarea)
writeRaster(bgarea, names, bylayer=T, format="GTiff", overwrite=T)


##plot the masked files
setwd("/home/jessica/Documents/UGA/Projects/seagrant/results/biomassmaps/final_N_rasters")
files<-list.files(pattern = ".tif$")
bg<-stack( files)
rasterNames<-as.Date(str_split(names(bg), "_", simplify = T)[,2], format="%Y%m%d")
x<-which(rasterNames>as.Date("2013-12-31") &rasterNames<as.Date("2019-01-01")&(format(rasterNames, "%m") %in% c("03", "06", "09", "12")))
x<-which((format(rasterNames, "%Y") %in% c("2014", "2016", "2018"))&(format(rasterNames, "%m") %in% c("03", "06", "09", "12")))
rasterNames<-rasterNames[x]
bg<-bg[[x]]
dat.label<-"% foliar N"
my.at=c(seq(0.4, max(newpixels$predn.l8, na.rm=T)+0.1,by=0.05))
label.at=c(seq(0.4, max(newpixels$predn.l8, na.rm=T)+0.1,by=0.1))
mycolorkey <- list(at=my.at, ## where the colors change
                   labels=list(
                     at=label.at, ## where to print labels
                     cex=1.25
                   ))

levelplot(bg, col.regions = rev(terrain.colors(255)), cuts=254, 
          margin=FALSE, at=my.at, colorkey=mycolorkey ,main=list(dat.label,cex=2),
          names.attr=as.character(rasterNames),layout=c(4,3),
          scales=list(draw=FALSE ), # remove axes labels & ticks
          xlim=c(471000, 474100),ylim=c(3477000,3481000),
          panel=function(...){
            panel.levelplot(...)
            #panel.text(471500, 3480200, cex=1.5,
            #            month)
          })
dev.print(file="../../../results/biomassmaps/N_xgb_stack.png",
          device=png, bg="white", width=720, height=730) 

###CHL####
##create separate geotiffs by date
setwd(paste(root_path,"/UGA/Projects/seagrant", sep=""))
dates<-unique(newpixels$date[newpixels$date>as.Date("2013-09-01")]) 
for (i in 1:length(dates)){
  newpixels1<-newpixels[newpixels$date==dates[i],]
  newpixels1<-newpixels1[!is.na(newpixels1$utm_east), ]
  newpixels1<-newpixels1[!is.na(newpixels1$predchl.l8), ]
  if(nrow(newpixels1)>1){
    proj <- "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
    coordinates(newpixels1) <- ~utm_east+utm_north ##converts newpixels to a shape file; cool 
    proj4string(newpixels1) <- proj  # define projection system of our data
    bgarea <- vect2rast(newpixels1, cell.size=30, fname="predchl.l8", 
                        file.name=paste0("results/biomassmaps/final_CHL_rasters/predchl_", format(dates[i], "%Y%m%d"), ".tif"))
  }
}
##read in all the files as rasters, stack them together and mask all in one step
setwd("/home/jessica/Documents/UGA/Projects/seagrant/results/biomassmaps/final_CHL_rasters")
files<-list.files(pattern = ".tif$")
bg<-stack( files)
bgarea<-mask(bg, mask)

##write out the masked files
names<-names(bgarea)
writeRaster(bgarea, names, bylayer=T, format="GTiff", overwrite=T)


##plot the masked files
setwd("/home/jessica/Documents/UGA/Projects/seagrant/results/biomassmaps/final_CHL_rasters")
files<-list.files(pattern = ".tif$")
bg<-stack( files)
rasterNames<-as.Date(str_split(names(bg), "_", simplify = T)[,2], format="%Y%m%d")
x<-which(rasterNames>as.Date("2013-12-31") &rasterNames<as.Date("2019-01-01")&(format(rasterNames, "%m") %in% c("03", "06", "09", "12")))
x<-which((format(rasterNames, "%Y") %in% c("2014", "2016", "2018"))&(format(rasterNames, "%m") %in% c("03", "06", "09", "12")))
rasterNames<-rasterNames[x]
bg<-bg[[x]]
dat.label<-"% Chlorophyll"
my.at=seq(0.4, max(newpixels$predchl.l8, na.rm=T)+0.05,by=0.05)
mycolorkey <- list(at=my.at, ## where the colors change
                   labels=list(
                     at=my.at, ## where to print labels
                     cex=1.25
                   ))

levelplot(bg, col.regions = rev(terrain.colors(255)), cuts=254, 
          margin=FALSE, at=my.at, colorkey=mycolorkey ,main=list(dat.label,cex=2),
          names.attr=as.character(rasterNames),layout=c(4,3),
          scales=list(draw=FALSE ), # remove axes labels & ticks
          xlim=c(471000, 474100),ylim=c(3477000,3481000),
          panel=function(...){
            panel.levelplot(...)
            #panel.text(471500, 3480200, cex=1.5,
            #            month)
          })
dev.print(file="../../../results/biomassmaps/CHL_xgb_stack.png",
          device=png, bg="white", width=720, height=730) 

####LAI####
##create separate geotiffs by date
setwd(paste(root_path,"/UGA/Projects/seagrant", sep=""))
dates<-unique(newpixels$date[newpixels$date>as.Date("2013-09-01")]) 
for (i in 1:length(dates)){
  newpixels1<-newpixels[newpixels$date==dates[i],]
  newpixels1<-newpixels1[!is.na(newpixels1$utm_east), ]
  newpixels1<-newpixels1[!is.na(newpixels1$predlai.l8), ]
  if(nrow(newpixels1)>1){
    proj <- "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
    coordinates(newpixels1) <- ~utm_east+utm_north ##converts newpixels to a shape file; cool 
    proj4string(newpixels1) <- proj  # define projection system of our data
    bgarea <- vect2rast(newpixels1, cell.size=30, fname="predlai.l8", 
                        file.name=paste0("results/biomassmaps/final_lai_rasters/predlai_", format(dates[i], "%Y%m%d"), ".tif"))
  }
}
##read in all the files as rasters, stack them together and mask all in one step
setwd("/home/jessica/Documents/UGA/Projects/seagrant/results/biomassmaps/final_lai_rasters")
files<-list.files(pattern = ".tif$")
bg<-stack( files)
bgarea<-mask(bg, mask)

##write out the masked files
names<-names(bgarea)
writeRaster(bgarea, names, bylayer=T, format="GTiff", overwrite=T)


##plot the masked files
setwd("/home/jessica/Documents/UGA/Projects/seagrant/results/biomassmaps/final_lai_rasters")
files<-list.files(pattern = ".tif$")
bg<-stack( files)
rasterNames<-as.Date(str_split(names(bg), "_", simplify = T)[,2], format="%Y%m%d")
x<-which(rasterNames>as.Date("2013-12-31") &rasterNames<as.Date("2019-01-01")&(format(rasterNames, "%m") %in% c("03", "06", "09", "12")))
x<-which((format(rasterNames, "%Y") %in% c("2014", "2016", "2018"))&(format(rasterNames, "%m") %in% c("03", "06", "09", "12")))
rasterNames<-rasterNames[x]
bg<-bg[[x]]
dat.label<-"Leaf Area Index"
my.at=c(seq(0.4, 1.5,by=0.05))
label.at=c(seq(0.4, 1.5,by=0.1))
mycolorkey <- list(at=my.at, ## where the colors change
                   labels=list(
                     at=label.at, ## where to print labels
                     cex=1.25
                   ))

levelplot(bg, col.regions = rev(terrain.colors(255)), cuts=254, 
          margin=FALSE, at=my.at, colorkey=mycolorkey ,main=list(dat.label,cex=2),
          names.attr=as.character(rasterNames),layout=c(4,3),
          scales=list(draw=FALSE ), # remove axes labels & ticks
          xlim=c(471000, 474100),ylim=c(3477000,3481000),
          panel=function(...){
            panel.levelplot(...)
            #panel.text(471500, 3480200, cex=1.5,
            #            month)
          })
dev.print(file="../../../results/biomassmaps/lai_xgb_stack.png",
          device=png, bg="white", width=720, height=730) 


####Inundation####
##create separate geotiffs by date
setwd(paste(root_path,"/UGA/Projects/seagrant", sep=""))
dates<-unique(newpixels$date[newpixels$date>as.Date("2013-09-01")]) 
for (i in 1:length(dates)){
  newpixels1<-newpixels[newpixels$date==dates[i],]
  newpixels1<-newpixels1[!is.na(newpixels1$utm_east), ]
  newpixels1<-newpixels1[!is.na(newpixels1$local_hiwater), ]
  if(nrow(newpixels1)>1){
    proj <- "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
    coordinates(newpixels1) <- ~utm_east+utm_north ##converts newpixels to a shape file; cool 
    proj4string(newpixels1) <- proj  # define projection system of our data
    bgarea <- vect2rast(newpixels1, cell.size=30, fname="local_hiwater", 
                        file.name=paste0("results/biomassmaps/final_inundation_rasters/inundation_", format(dates[i], "%Y%m%d"), ".tif"))
  }
}
##read in all the files as rasters, stack them together and mask all in one step
setwd("/home/jessica/Documents/UGA/Projects/seagrant/results/biomassmaps/final_inundation_rasters")
files<-list.files(pattern = ".tif$")
bg<-stack( files)
bgarea<-mask(bg, mask)

##write out the masked files
names<-names(bgarea)
writeRaster(bgarea, names, bylayer=T, format="GTiff", overwrite=T)


##plot the masked files
setwd("/home/jessica/Documents/UGA/Projects/seagrant/results/biomassmaps/final_inundation_rasters")
files<-list.files(pattern = ".tif$")
bg<-stack( files)
rasterNames<-as.Date(str_split(names(bg), "_", simplify = T)[,2], format="%Y%m%d")
x<-which(rasterNames>as.Date("2013-12-31") &rasterNames<as.Date("2019-01-01")&(format(rasterNames, "%m") %in% c("02", "05", "09", "12")))
x<-which((format(rasterNames, "%Y") %in% c("2014", "2016", "2018"))&(format(rasterNames, "%m") %in% c("02", "05", "09", "12")))
rasterNames<-rasterNames[x]
bg<-bg[[x]]
dat.label<-"Inundation intensity"
my.at=c(seq(0, 0.7,by=0.05), 0.8,0.9,1.0,2.0)
label.at=c(seq(0, 1,by=0.1),2)
mycolorkey <- list(at=my.at, ## where the colors change
                   labels=list(
                     at=label.at, ## where to print labels
                     cex=1.25
                   ))

levelplot(bg, col.regions = hcl.colors(255, palette = "viridis", alpha = NULL, rev = T, fixup = TRUE), 
          cuts=254, 
          margin=FALSE, at=my.at, colorkey=mycolorkey ,main=list(dat.label,cex=2),
          names.attr=as.character(rasterNames),layout=c(4,3),
          scales=list(draw=FALSE ), # remove axes labels & ticks
          xlim=c(471000, 474100),ylim=c(3477000,3481000),
          panel=function(...){
            panel.levelplot(...)
            #panel.text(471500, 3480200, cex=1.5,
            #            month)
          })
dev.print(file="../../../results/biomassmaps/inundation_xgb_stack.png",
          device=png, bg="white", width=720, height=730) 


####Green Up####
##create separate geotiffs by date
setwd(paste(root_path,"/UGA/Projects/seagrant", sep=""))
dates<-unique(newpixels$date[newpixels$date>as.Date("2013-09-01")]) 
for (i in 1:length(dates)){
  newpixels1<-newpixels[newpixels$date==dates[i],]
  newpixels1<-newpixels1[!is.na(newpixels1$utm_east), ]
  newpixels1<-newpixels1[!is.na(newpixels1$greendoy), ]
  if(nrow(newpixels1)>1){
    proj <- "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
    coordinates(newpixels1) <- ~utm_east+utm_north ##converts newpixels to a shape file; cool 
    proj4string(newpixels1) <- proj  # define projection system of our data
    bgarea <- vect2rast(newpixels1, cell.size=30, fname="greendoy", 
                        file.name=paste0("results/biomassmaps/final_greenup_rasters/greenup_", format(dates[i], "%Y%m%d"), ".tif"))
  }
}
##read in all the files as rasters, stack them together and mask all in one step
setwd("/home/jessica/Documents/UGA/Projects/seagrant/results/biomassmaps/final_greenup_rasters")
files<-list.files(pattern = ".tif$")
bg<-stack( files)
bgarea<-mask(bg, mask)

##write out the masked files
names<-names(bgarea)
writeRaster(bgarea, names, bylayer=T, format="GTiff", overwrite=T)


##plot the masked files
setwd("/home/jessica/Documents/UGA/Projects/seagrant/results/biomassmaps/final_greenup_rasters")
files<-list.files(pattern = ".tif$")
bg<-stack( files)
rasterNames<-as.Date(str_split(names(bg), "_", simplify = T)[,2], format="%Y%m%d")
x<-which((format(rasterNames, "%Y") %in% c("2014", "2015", "2016", "2017", "2018","2019")&(format(rasterNames, "%m")=="09")))
rasterNames<-rasterNames[x]
bg<-bg[[x]]
dat.label<-"Green-up DOY"
panel.label<-"Green-up DOY"
my.at=c(seq(20, 60,by=1))
label.at=c(seq(10, 75,by=5))
mycolorkey <- list(at=my.at, ## where the colors change
                   labels=list(
                     at=label.at, ## where to print labels
                     cex=1.25
                   ))
#y<-summary(bg)
#range<-(y[5,]-y[1,])

levelplot(bg, col.regions = hcl.colors(255, palette = "viridis", alpha = NULL, rev = T, fixup = TRUE), 
          cuts=254, 
          margin=FALSE, at=my.at, colorkey=mycolorkey ,main=list(dat.label,cex=2),
          names.attr=as.character(format(rasterNames, "%Y")),layout=c(3,2),
          scales=list(draw=FALSE ), # remove axes labels & ticks
          xlim=c(471000, 474100),ylim=c(3477000,3481000),
          panel=function(...){
            panel.levelplot(...)
            #panel.text(471500, 3480200, cex=1.5,
            #            range[x])
          })
dev.print(file="../../../results/biomassmaps/greenup_xgb_stack.png",
          device=png, bg="white", width=720, height=730) 



####Land Surface Temp####
##create separate geotiffs by date
setwd(paste(root_path,"/UGA/Projects/seagrant", sep=""))
dates<-unique(newpixels$date[newpixels$date>as.Date("2013-09-01")]) 
for (i in 1:length(dates)){
  newpixels1<-newpixels[newpixels$date==dates[i],]
  newpixels1<-newpixels1[!is.na(newpixels1$utm_east), ]
  newpixels1<-newpixels1[!is.na(newpixels1$lst), ]
  if(nrow(newpixels1)>1){
    proj <- "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
    coordinates(newpixels1) <- ~utm_east+utm_north ##converts newpixels to a shape file; cool 
    proj4string(newpixels1) <- proj  # define projection system of our data
    bgarea <- vect2rast(newpixels1, cell.size=30, fname="lst", 
                        file.name=paste0("results/biomassmaps/final_lst_rasters/lst_", format(dates[i], "%Y%m%d"), ".tif"))
  }
}
##read in all the files as rasters, stack them together and mask all in one step
setwd("/home/jessica/Documents/UGA/Projects/seagrant/results/biomassmaps/final_lst_rasters")
files<-list.files(pattern = ".tif$")
bg<-stack( files)
bgarea<-raster::mask(bg, mask)

##write out the masked files
names<-names(bgarea)
writeRaster(bgarea, names, bylayer=T, format="GTiff", overwrite=T)


##plot the masked files
setwd("/home/jessica/Documents/UGA/Projects/seagrant/results/biomassmaps/final_lst_rasters")
files<-list.files(pattern = ".tif$")
bg<-stack( files)
rasterNames<-as.Date(str_split(names(bg), "_", simplify = T)[,2], format="%Y%m%d")
x<-which((format(rasterNames, "%Y") %in% c("2014", "2015", "2016", "2017", "2018","2019")&(format(rasterNames, "%m")=="02")))
rasterNames<-rasterNames[x]
bg<-bg[[x]]
dat.label<-"Land Surface Temperature (°C)"
my.at=c(seq(5, 25,by=0.5))
label.at=c(seq(5, 25,by=2))
mycolorkey <- list(at=my.at, ## where the colors change
                   labels=list(
                     at=label.at, ## where to print labels
                     cex=1.25
                   ))

levelplot(bg, #col.regions = rev(soc_pal(255)),
          #col.regions = rev(hcl.colors(255, palette = "viridis", alpha = NULL, rev = T, fixup = TRUE)),
          #col.regions = colorRampPalette(rev(SAGA_pal[["SG_COLORS_RED_BLACK"]]))(255),
          #col.regions = colorRampPalette((SAGA_pal[["SG_COLORS_DEFAULT_BRIGHT"]]))(255),
          col.regions = rev(terrain.colors(255)),
          cuts=254, 
          margin=FALSE, at=my.at, colorkey=mycolorkey ,main=list(dat.label,cex=2),
          names.attr=as.character(rasterNames),layout=c(3,2),
          scales=list(draw=FALSE ), # remove axes labels & ticks
          xlim=c(471000, 474100),ylim=c(3477000,3481000),
          panel=function(...){
            panel.levelplot(...)
            #panel.text(471500, 3480200, cex=1.5,
            #            month)
          })
dev.print(file="../../../results/biomassmaps/lst_stack.png",
          device=png, bg="white", width=720, height=730) 