###This script applies already created aboveground biophysical models to new 
###prediction data sets, usually broad scale site level data
###Version for running on local laptop for an expanded prediction area, based on Hladik et al. 2013 habitat map of spartina alterniflora 

library(R.utils); library(tidyverse);library(RColorBrewer)
library(data.table); library(xgboost); library(randomForest)
library("paradox");  library("mlr3");library("mlr3learners"); 
library("mlr3tuning"); library("randomizr")

ls <- list.files(path = "/home/kyle/Documents/Documents/UT/Sapelo/berm_output", pattern = "processed_split")
for(i in seq(1,length(ls))){
  
#reading in splits of dataset
cut_no <- i

#####Site level prediction set data####
###Read in processed Landsat 8 data for the site level prediction set ####
newpixels<-read_csv(paste0("/home/kyle/Documents/Documents/UT/Sapelo/berm_output/gce_pixels_landsat_processed_split", cut_no, ".csv"))
newpixels$date<-as.Date(newpixels$date)


##Read in the greenup and surface temperature estimates for the site level prediction set
newpixels.heat.l8<-read_csv(paste0("output/gce_pixels_landsat8_greenup_via_lst_and_air_lags_split", cut_no, ".csv"))
newpixels.heat.l8$greenup<-as.Date(newpixels.heat.l8$greenup.soil)

##variable names for use in model prediction for aboveground biophysical models
feature.names<-c("b1","b2","b3","b4","b5","b6","b7","ndvi", "pheno", "vari")

##note for each of the models below, the following column names are expected as predictor features:
#feature.names<-c("b1","b2","b3","b4","b5","b6","b7","ndvi", "pheno", "vari")
##the above was handled by the newpixels processing script

####Predict aboveground biomass####
##first load the AGB model
agmod1<-xgb.load("output/xgb2.0_agb_model")

##now predict the data
newpixels<-data.frame(newpixels)
newpixels$predag.allom.l8 <-(predict(agmod1, data.matrix(newpixels[,feature.names])))

####Predict foliar N####
##first load the foliar N model
Nmod1<-xgb.load("output/xgb2.0_foliarN_model")
newpixels<-data.frame(newpixels)

##now predict the data
newpixels$predn.l8<- (predict(Nmod1, data.matrix(newpixels[,feature.names])))


####Predict LAI####
##first load the model
laimod1<-xgb.load("output/xgb2.0_lai_model")

##now predict the data
newpixels<-data.frame(newpixels)
newpixels$predlai.l8<- (predict(laimod1, data.matrix(newpixels[,feature.names])))


####Predict CHL####
##first load the model
feature.names<-c("b1","b2", "b3","b4","b5","b6","b7", "gari", "ndmi", "pheno", "vari", "ndvi")
chlmod1<-xgb.load("output/xgb2.0_chl_model")

##now predict the data
newpixels<-data.frame(newpixels)
newpixels$predchl.l8<- (predict(chlmod1, data.matrix(newpixels[,feature.names])))


##save a copy of the prediction results for later use
write_csv(newpixels, path = paste0("/home/kyle/Documents/Documents/UT/Sapelo/berm_output/xgb_predicted_biophysical_landsat_hladik_split", cut_no, ".csv"))
}
