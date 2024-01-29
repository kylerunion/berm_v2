###This script creates the aboveground biophysical models from the ground-truth field data
###It also estimates pixel flood status and percent flooded, date of pixel greenup and LST calculations

library(R.utils); library(caret); library(tidyverse);library(RColorBrewer)
library(data.table); library(xgboost); library(randomForest)
library("paradox");  library("mlr3");library("mlr3learners"); 
library("mlr3tuning"); library("randomizr")
library(devtools)
install_github('theislab/kBET')
library(kBET)


##output and save new models and data sets?? Set to True/False below. 
##Note that True will overwrite the existing models and output datasets, so you'll want a backup if you don't intend this
#save_model<-FALSE
save_model<- T

## WHICH VERSION ARE YOU CREATING
#v <- 1
v <- 2

color<-colorRampPalette(brewer.pal(name="Dark2", n = 8))(9)
trans.color<-addalpha(color,0.9)
all_sites <- c("deanc", "dupli", "fluxa", "fluxb", "folly", "huntc", "north", "skida", "ugami")
old_sites <- c("fa", "fb", "sk", "ug")
old_sites_hc <- c("fluxa", "fluxb", "huntc", "skida", "ugami")
cols_all <- trans.color[match(all_sites, all_sites)]
cols_old <- trans.color[match(old_sites, all_sites)]
cols_old_hc <- trans.color[match(old_sites_hc, all_sites)]


####Ground-truth data for field observations####
##Read in location data for vegetation ground truth plots
plots<-read.csv("data/biomass_plots.csv", header=T)
##not all dates have all pixels, so we need to make the averages here in order to avoid duplicate elevations later
pixels<- plots %>% group_by(pix) %>% summarise(elevation=mean(elevation), 
                                               utm_east=mean(utm_east), 
                                               utm_north=mean(utm_north), 
                                               lat=mean(lat), long=mean(long))
plots<-dplyr::select(plots, plot, pix, utm_east, utm_north, elevation)

###Read in Landsat observations for ground-truth plots, collected from Google Earth Engine Script####
### read in the landsat8 data for the main vegetation data plots
landsat8<-read.csv("data/L8_bermplots_2013_2023_noclouds.csv",
                  header = T)
landsat9 <- read.csv("data/L9_bermplots_2021_2023.csv", 
                     header = T)
landsat8$platform <- 8
landsat9$platform <- 9
landsat <- rbind(landsat8, landsat9)

##process the landsat data to create needed time, date, and location columns
landsat$time<-substr(landsat$date,12,19)
landsat$date<-as.Date(substr(landsat$date,1,10))
landsat$year<-as.numeric(format(landsat$date, "%Y"))
landsat$doy<-as.numeric(format(landsat$date, "%j"))
landsat$mo<-as.numeric(format(landsat$date, "%m"))
landsat<-merge(landsat, plots, by="plot")
##filter the landsat 8 data if needed through the pixel_qa landsat mask; 
##note that the current Landsat 8 Google Earth Engine script only returns cloud filtered data
##see https://landsat.usgs.gov/landsat-surface-reflectance-quality-assessment
landsat$QA_PIXEL<-intToBin(landsat$QA_PIXEL)
##"01000000" from the right means no fill, no dilated clouds, no cirrus, no cloud, no cloud shadow, no snow, yes clear, no water
## no water above is via Landsat 8's pixel qa mask, and misses a lot of marsh flooding in mixed pixels, we'll handle this below
landsat$qa_good<-ifelse(str_sub(landsat$QA_PIXEL,-8,-1)=="01000000",T,F) 
#let's keep in the water yes for the flood_time variable we'll make later
landsat$qa_good<-ifelse(str_sub(landsat$QA_PIXEL,-8,-1)=="11000000",T,landsat$qa_good) 

## subset to only good pixels as indicated by the pixel_qa mask
landsat<-landsat[landsat$qa_good==T,]

##remove columns we don't need now from Landsat 8 Google Earth Engine data
landsat<-dplyr::select(landsat,-c(system.index, .geo, QA_RADSAT, SR_QA_AEROSOL, ST_ATRAN, ST_CDIST, ST_DRAD, ST_EMIS, ST_EMSD, ST_QA, ST_TRAD, ST_URAD, QA_PIXEL, qa_good))
landsat$year<-as.numeric(format(landsat$date, "%Y"))

## Scaling done in GEE (L8L2C2T1)
landsat$b1<-landsat$SR_B1
landsat$b2<-landsat$SR_B2
landsat$b3<-landsat$SR_B3
landsat$b4<-landsat$SR_B4
landsat$b5<-landsat$SR_B5
landsat$b6<-landsat$SR_B6
landsat$b7<-landsat$SR_B7

##GEE takes care of thermal scaling too
landsat$b10<-landsat$ST_B10

##clean up columns to remove the ones we now don't need
landsat<-dplyr::select(landsat, -c(SR_B1,SR_B2,SR_B3,SR_B4,SR_B5,SR_B6,SR_B7,ST_B10))

##negative reflectance values can occur at scene edges and should be removed
landsat$b1<-ifelse(landsat$b1<0,NA, landsat$b1)
landsat$b2<-ifelse(landsat$b2<0,NA, landsat$b2)
landsat$b3<-ifelse(landsat$b3<0,NA, landsat$b3)
landsat$b4<-ifelse(landsat$b4<0,NA, landsat$b4)
landsat$b5<-ifelse(landsat$b5<0,NA, landsat$b5)
landsat$b6<-ifelse(landsat$b6<0,NA, landsat$b6)
landsat$b7<-ifelse(landsat$b7<0,NA, landsat$b7)

##filter to non-NA values
landsat<-landsat[is.na(landsat$b1)==F,]; landsat<-landsat[is.na(landsat$b7)==F,]; 
landsat<-landsat[is.na(landsat$b5)==F,];landsat<-landsat[is.na(landsat$b6)==F,]

##filter out date that looks weird
landsat <- landsat %>% filter(date != as.Date("2022-11-23"))
#and events that had really high b3 values
landsat <- landsat %>% filter(b3 < 0.15)


##add in a site variable for grouping the data
landsat$site<-"fluxa"
landsat$site<-ifelse(str_sub(landsat$plot,1,2)=="ug", "ugami", landsat$site)
landsat$site<-ifelse(str_sub(landsat$plot,1,2)=="sk", "skida", landsat$site)
landsat$site<-ifelse(str_sub(landsat$plot,1,2)=="fb", "fluxb", landsat$site)
landsat$site<-ifelse(str_sub(landsat$plot,1,2)=="dc", "deanc", landsat$site)
landsat$site<-ifelse(str_sub(landsat$plot,1,2)=="du", "dupli", landsat$site)
landsat$site<-ifelse(str_sub(landsat$plot,1,2)=="fr", "folly", landsat$site)
landsat$site<-ifelse(str_sub(landsat$plot,1,2)=="hc", "huntc", landsat$site)
landsat$site<-ifelse(str_sub(landsat$plot,1,2)=="ns", "north", landsat$site)
table(landsat$site)

plot(landsat$date, landsat$b3, col = landsat$platform)

##load functions
source("functions/calc_indices_l8.R")
source("functions/soil_temp_from_lst_and_air_lags.R")

###calculate a set of standard vegetation and spectral reflectance indices for landsat 8 data using our custom function
indices<-calc_index_l8(landsat)
landsat<-cbind(landsat, indices)

##calculate lst with custom function
## this represents an update from O'Connell et al 2021; the updated GEE L8 product has surface temp instead of brightness temp
landsat$lst<-landsat$b10-273.15
landsat<-landsat[!is.na(landsat$lst),]

##create rolling pheno mean
##trick to get all dates, just create the full dates dataframe and merge it to the missing dates dataframe
y<-data.frame(date=seq.Date(from=min(landsat$date, na.rm=T), to=max(landsat$date, na.rm=T), by=1))
y<-expand.grid(date=y$date, plot=unique(landsat$plot))
landsat<-full_join(landsat,y)
landsat<-arrange(landsat, plot,date)
##this is an even window roll, because there's an observation for every day
landsat<-landsat %>% dplyr::group_by(plot) %>%
  dplyr::mutate(pheno2= RcppRoll::roll_mean(pheno, n=30, na.rm=T,align="center", fill=NA))
landsat<-landsat[!is.na(landsat$b2),]

##Use the tidal flooding prediction random forest model (super_model); 
##We'll predict tidal flooding and filter the data to dry observations
super_model <- readRDS("output/random_forest_flood_model_l8.rds")
print(super_model)

## create some more of the needed columns the flood model expects
landsat$mo<-factor(landsat$mo, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))

##predict flooding and then process the prediction column to interpretable labels
landsat$flooded <- predict(super_model, newdata=landsat)
prop.table(table(landsat$flooded))
landsat$est.flood<-ifelse(landsat$flooded==0, "dry", "wet")
landsat$est.flood<-factor(landsat$est.flood)
length(table(landsat$date, landsat$flooded))

##read in the smoothed vegetation data, to merge with the landsat data
data<-read.csv("output/smoothed_veg_timeseries_weekly.csv", header=T, stringsAsFactors = F)
data<-arrange(data, date)
data$date<-as.Date(data$date)
data$field_lai<-data$LAI
data<-dplyr::select(data, -LAI)

#BERMv1.0
if(v == 1){
data <- data[data$site %in% old_sites,]
data <- data[data$date < as.Date("2019-09-01"),]
}

##quality control on pixels
x<-table(landsat$date, landsat$est.flood);
##this gets rid of dates with many wet obs
#x<-x[x[,2]<v*10,] 
##this last gets rid of dates with few dry satelitte obs.
#x<-x[x[,1]>v*5,] 
#comment this out for the percent flooded. added it in later
#landsat<-landsat[landsat$date %in% as.Date(row.names(x)),]
#landsat<-landsat[landsat$est.flood=="dry",]

landsat<-landsat[is.na(landsat$b5)==F,]
landsat<-arrange(landsat, date)


heat_data <- inner_join(data, plots, by = "plot")
  
  
#create a better heat interpolation
heat.l8<-setup_heat(landsat, gce=F)
if(save_model == T){
  write_csv(heat.l8, paste0("output/berm", v, ".0_pixels_landsat8_heat_processed.csv"))
}

###Merge lsat and field landsat < 10 days apart
##find the week in landsat landsat that is closest to a week in field survey  to create a merge key
data$start<-data$date-10; data$end<-data$date+11
data<-data.table(data)
setkey(data, plot, start, end)

##add on empty rows with the date and plot of early landsat observations so that these make it into the merged landsat
earlydates<-unique(landsat[,c("date", "plot")])
earlydates<-earlydates[earlydates$date<min(data$start),]
var<-names(data)
tmp<-data[1,]
tmp2<-tmp
for (i in 2:nrow(earlydates)){
  tmp<-rbind(tmp, tmp2)
}
tmp[,1:ncol(tmp)]<-NA
#tmp$date<-earlydates$date
#tmp$plot<-earlydates$plot
#data<-rbind(data, tmp)

landsat<-landsat[!(landsat$plot %in% c("tb1", "tb2")),]
landsat<-dplyr::select(landsat, -mo)
data<-dplyr::select(data, -mo, -site)
data<-data[!is.na(data$date),]

data$start<-as.Date(data$date)-10; data$end<-as.Date(data$date)+11
data<-data.table(data)
setkey(data, plot, start, end)

##now create an index in the lsat landsat, here we'll just repeat the date column
##twice so that if the landsat date is within the field interval, we'll get a match
landsat$lsat_date<-as.Date(landsat$date)
landsat$start<-landsat$date; landsat$end<-landsat$date
landsat<-dplyr::select(landsat, -date)

###this removes dates before the field landsat; but we want to predict those dates to get rolling means
#landsat<-landsat[landsat$lsat_date>=min(data$start),]
landsat<-data.table(landsat)
setkey(landsat, plot, start, end)


#join if the date of landsat is within the end and start interval of field,
##return nothing for rows without matches, return two rows for multiple matches (we'll pick the best one next)
data<-foverlaps(landsat, data, type="any", by.x=c("plot", "start", "end"),
                by.y=c("plot", "start", "end"), nomatch=NA, mult="all")

#return veg data and landsat data back to tibbles;
#data.table objects have some different behavior than data.frames and tibbles and require care when using
landsat<-as_tibble(landsat); data<-as_tibble(data)

##subset field obs to those that are close to a landsat date,
##eliminates duplicate field estimates for the same landsat date
data$date<-as.Date(ifelse(is.na(data$date),data$lsat_date, data$date), origin="1970-01-01")
data$diff<-abs(data$date-data$lsat_date)

data<-data %>% dplyr::group_by(lsat_date, plot) %>%
  filter(diff==min(diff, na.rm=T))
##we don't need this, already decided above the allowed date range
#data<-data %>%  filter(abs(diff)< 10)
data<-data %>% dplyr::select(-c(start,end, i.start, i.end, diff))

data$landsat_doy<-data$doy
data$mo<-as.numeric(format(data$date, "%m"))

##create pixel averages and pixel average flood status
data<-dplyr::mutate(data, pixel=paste(year,landsat_doy,pix, site,sep="_"))
data<-ungroup(data)

##summarize flood status, to write out freq of flood observations by pixel for later use
#took care of flooding data while commenting out things that removed flooded pixels. so don't save this again unless you address that
flood<-data %>%
  dplyr::group_by(date, pix,site, flooded) %>% dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)
jberm <- subset(flood, date < as.Date("2016-11-01"))
prop.table(table(flood$flooded))
flood$flooded<-as.numeric(flood$flooded)
flood$flooded<-ifelse(flood$flooded==1,0,flood$flooded)
flood$flooded<-ifelse(flood$flooded==2,1,flood$flooded)
flood_by_year <- flood %>% group_by(pix, year) %>% summarise(flood_time = sum(as.numeric(flooded), na.rm = T)/length(flooded[!is.na(flooded)]), count = length(flooded[!is.na(flooded)]))
flood_by_year <- flood_by_year %>%
  na.omit() %>%
  group_by(year) %>%
  summarise(mean_flood_time = mean(flood_time), sum_count = sum(count))
write.csv(flood_by_year, "output/flood_by_year.csv")
flood<- flood %>% dplyr::group_by(pix) %>% dplyr::summarise(flood_time= sum(as.numeric(flooded), na.rm=T)
                                                            /length(flooded[!is.na(flooded)]), count = length(flooded[!is.na(flooded)]))

jberm <- subset(jberm, site %in% c("fluxa", "skida", "fluxb", "ugami"))
prop.table(table(jberm$flooded))
jberm$flooded<-as.numeric(jberm$flooded)
jberm$flooded<-ifelse(jberm$flooded==1,0,jberm$flooded)
jberm$flooded<-ifelse(jberm$flooded==2,1,jberm$flooded)
jberm<- jberm %>% dplyr::group_by(pix) %>% dplyr::summarise(flood_time= sum(as.numeric(flooded), na.rm=T)
                                                            /length(flooded[!is.na(flooded)]))
#don't need to be resaving this all the time, particularly when you're filtering out water observations
#if(save_model==T){
#  write_csv(flood, "output/flooded_percent_from_landsat_berm_pixels.csv")
#}

##get rid of variables we don't need now to prepare for pixel averaging of only dry pixels
data<-data %>% dplyr::select(-c( plot,est.flood,  time, flooded,pixel,mo, 
                                 year,lsat_date,landsat_doy, elevation, utm_east, utm_north, lat, long
)) 
data<-data %>% dplyr::group_by(date, pix,site) %>% dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)
data<-ungroup(data)
data$mo<-as.numeric(format(data$date, "%m"))
data$year<-as.numeric(format(data$date, "%Y"))
data$doy<-as.numeric(format(data$date, "%j"))
data$group<-paste(data$site, data$mo)
data<-data[!is.na(data$site),]
data<-data[!is.na(data$b1),]
data$id<-as.numeric(row.names(data))
data$spad<-data$chlorophyll


##XGBoost
##Set up the resampling scheme for nested resampling
outer<-1
inner<-3

clusters <-data$group
blocks <-  data$site
rep<-1:outer
vars<-paste0("cv",rep)
cvs<-matrix(data=NA,nrow=nrow(data),ncol=length(vars))

###create the cross valdidated folds####
###29382
set.seed(999) 
for( i in 1: outer){
  cvs[,i]<- block_and_cluster_ra(blocks = blocks,
                                 clusters = clusters,
                                 prob_each = c(.35, .65))
}       

##require tidbit lai to go into the training data, as we only have one month for that location
cvs<-data.frame(cvs); names(cvs)<-vars
head(table(clusters, cvs$cv1));head(table(blocks, cvs$cv1))

#clusters; blocks
table(cvs$cv1)
head(table(blocks, cvs$cv1))

##set up training/testing index for landsat data for later access
for (i in 1:outer){
  data[[paste0("train.", i)]] <- cvs[,i]}

table(data$train.1)
table(cvs$cv1)

#analysis_table <- data.frame(matrix(ncol = 8, nrow = 0))
#colnames(analysis_table) <- c("n_test", "n_train", "mae", "rmse", "nrmse", "cov_rmse", "cor", "treatment")

##create a list of parameters to focus on tuning, must specify type (Int=Integer, Dbl=continuous number, Fct=factor/character)
tune_type <- "new"
       tune_ps = ParamSet$new(list(
  ParamDbl$new("eta", lower = 0.01, upper = 0.1),
  ParamInt$new("max_depth", lower = 2, upper = 8),
  ParamInt$new("min_child_weight", lower = 1, upper = 5),
  ParamDbl$new("subsample", lower = 0.15, upper = 1),
  ParamDbl$new("colsample_bytree", lower = 0.25, upper = 1)
))

tune_ps

run <- 1

###AGBIOMASS#####
tr<-data

tr<-dplyr::filter(tr, agm2.allom>10)
feature.names<-c("b2","b3","b4","b5","b6","b7","ndvi", "pheno", "ndmi", "vari", "gci", "mtvi1")
#feature.names<-c("b3","ndvi", "pheno", "vari")
resp<-"agm2.allom"
tr<-tr[tr$train.1==1,c(resp, feature.names)]
tr<-tr[complete.cases(tr),]

traintask = TaskRegr$new(id = "agbiomass", backend = tr[,c(resp, feature.names)], target = resp)
traintask$select(feature.names)

##set up learner
##select a resampling stategy for the datasets and preformance measure to use to tune the paramters
##set up how the tuner should iterate through the hyperparamters to try; 
##set up how the tuning should stop: here it stops after 20 iterations
##see other options at: https://mlr3book.mlr-org.com/tuning.html
learner = lrn("regr.xgboost", predict_type = "response",  booster="gbtree",  objective   = "reg:squarederror", nrounds=500)
#resamp=rsmp("holdout", ratio=0.85)
resamp=rsmp("cv", folds=inner)
measure = msr("regr.mse")
tuner = mlr3tuning::tnr("grid_search", resolution = 25)
terminator = trm("stagnation", iters=4, threshold=0.1)

instance = mlr3tuning::TuningInstanceSingleCrit$new(
  task = traintask,
  learner = learner,
  resampling = resamp,
  measure = measure,
  search_space = tune_ps,
  terminator = terminator
)
instance

##run the optimization
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
#instance$result_learner_param_vals

##store the best parameters back in the learner
learner$param_set$values = instance$result_learner_param_vals
learner$train(traintask)

##see the result, eg, the best configuration of parameters found during the iterations
param <- instance$result_x_domain

resp<-tr$agm2.allom
clf <- xgboost(data        = data.matrix(tr[,c(feature.names)]),
               label       = resp,
               booster="gbtree",
               nrounds     = 700,
               params = param,
               verbose=0)

mat <- xgb.importance (feature_names = feature.names,model = clf)
xgb.plot.importance (importance_matrix = mat[1:11]) 
p <-predict(clf, data.matrix(tr[, feature.names]))
plot(resp, p, ylim=c(0,max(c(p, resp), na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T))); abline(0,1)
sqrt((sum((resp-p)^2, na.rm=T))/length(p))

tr_results<- data.frame(na.omit(data$agm2.allom[data$train.1==1]), p, na.omit(data$site[data$train.1==1 & data$agm2.allom > 0]), na.omit(data$year[data$train.1==1 & data$agm2.allom > 0]))
colnames(tr_results) <- c("resp", "p", "site", "year")
tr_results$year <- as.character(tr_results$year)
tr_metrics <- tr_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")

test <- data

test<-dplyr::filter(test, agm2.allom>10)
resp<-test$agm2.allom[test$train.1==0]
resp <- na.omit(resp)
p <-predict(clf, data.matrix(test[test$train.1==0, feature.names]))
rmse<-sqrt((sum((resp-p)^2, na.rm=T))/length(p))
print(rmse)
rmse<-sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)]))
print(rmse)
plot(resp, p, ylim=c(0,max(c(p, resp),na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T)), pch=19,
     col=trans.color[c(1:9)][factor(test$site[test$train.1==0])],
     cex.lab=1.5, cex.axis=1.2)
abline(0,1)
mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(sprintf("%.1f",round(rmse,1))) ~ g ~ m^-2))
legend("bottomright", legend=unique(test$site), 
       title=expression(underline(site)),title.adj = 0.25,
       col=trans.color[c(1:9)][factor(unique(test$site))], pch=15, bty="n")

##aboveground biophysical results
#rmse[n]
rmse <- sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)]))
#nRMSE
nrmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))/(max(resp, na.rm=T)-min(resp, na.rm=T))
#COV RMSE (alternate nRMSE)
cov_rmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))/(mean(resp, na.rm=T))
#mae[n]
mae <- (sum((p-resp), na.rm=T))/length(resp[!is.na(resp)])
#cor
cor <- cor(resp, p, use="pairwise.complete.obs")
#N
n_test <- nrow(data[!is.na(data$agm2.allom)&data$train.1==0,])
n_train <- nrow(data[!is.na(data$agm2.allom)&data$train.1==1,])
#table
agb_tab <- as.data.frame(rbind(c(mae, rmse, nrmse, cov_rmse, cor, n_test, n_train)))

agmod1<-clf 
#clf<-agmod1
data<-data.frame(data)
#data <- old
data$predag.allom.l8<-(predict(agmod1, data.matrix(data[,feature.names])))
maxx<-max(c(data$agm2.allom, data$predag.allom.l8), na.rm=T)
minn<-min(c(data$agm2.allom, data$predag.allom.l8), na.rm=T)
sqrt((sum((data$agm2.allom-data$predag.allom.l8)^2, na.rm=T))/length(data$predag.allom.l8[!is.na(data$predag.allom.l8)]))

data <- data[order(data$site),]


if(save_model==TRUE){
  ##save a copy of the model for future use
  saveRDS(agmod1, paste0("output/xgb", v, ".0_agb_model.rds"))
  xgb.save(agmod1,  paste0("output/xgb", v, ".0_agb_model"))
  agmod1<-xgb.load(paste0("output/xgb", v, ".0_agb_model"))
  
  ##create a graph of measured vs predicted for training and testing data
  jpeg(file=paste0("results/xgb", v, ".0_l8_ag_allometeric.jpg"),
       height=5, width=5, res=200, units="in")
  par(mar=c(4,5,0.5,0.75), oma=c(0,0.2,0.5,0.5), mfrow=c(1,1))
  plot(data$agm2.allom, data$predag.allom.l8,
       xlab=expression(paste("measured AG biomass g ", m^-2, sep="")),
       ylab =expression(paste("predicted AG biomass g ", m^-2, sep="")) ,
       ylim=c(minn,maxx),xlim=c(minn,maxx),  pch=c(1,19)[factor(data$train)],
       cex=c(1,0.8)[factor(data$train)], col=cols_all[c(1:9)][factor(data$site)],
       cex.lab=1.5, cex.axis=1.2)
  abline(0,1)
  mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(round(rmse,1)) ~ g~m^-2), cex=1.2)
  legend("bottomright",inset=c(0,0.07),legend=c("train", "test"), pch=c(19,1), bty="n", title=expression(underline(data)), title.adj = 0.3)
  legend("bottomright", inset=c(0.2,0),legend=unique(data$site[!is.na(data$site)&data$site!="xtidb"]),
         title=expression(underline(site)),title.adj = 0.25,
         col=cols_all[c(1:9)][factor(unique(data$site[!is.na(data$site)&data$site!="xtidb"]))], pch=15, bty="n")
  dev.off()
}

##results by site/year
test_results<- data.frame(na.omit(test$agm2.allom[test$train.1==0]), p, na.omit(test$site[test$train.1==0]), na.omit(test$year[test$train.1==0]))
colnames(test_results) <- c("resp", "p", "site", "year")
test_results$year <- as.character(test_results$year)
test_metrics <- test_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")
uniques <- unique(test_metrics$value)
ag_tab <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(ag_tab) <- c("value", "rmse", "nrmse", "cor", "model", "version", "run")
model <- "Aboveground Biomass"
version <- paste0(v, ".0")
for (i in uniques) {
  resp <- test_metrics$resp[test_metrics$value == i]
  p <- test_metrics$p[test_metrics$value == i]
  rmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))
  nrmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))/(max(resp, na.rm=T)-min(resp, na.rm=T))
  cor <- cor(resp, p)
  n_test <- length(resp)
  n_train <- length(tr_metrics$resp[tr_metrics$value == i])
  met_tab <- as.data.frame(rbind(c(i, rmse, nrmse, cor, n_test, n_train, model, version, run)))
  colnames(met_tab) <- c("value", "rmse", "nrmse", "cor", "n_test","n_train", "model", "version", "run")
  ag_tab <- rbind(ag_tab,met_tab)
}
#all tr & test data
all_tr_test_results_ag <- rbind(tr_results, test_results)
all_tr_test_results_ag$model <- "Aboveground Biomass"
all_tr_test_results_ag$split <- c(rep("train", length(tr_results$resp)), rep("test", length(test_results$resp)))

## commonality analysis
library(yhat)
if(v == 2){
CC_agb=yhat::commonalityCoefficients(test_results, "p", list("site", "year"))
CC_agb_tab_per <- as.data.frame(rbind(c("agb", CC_agb$CC[5], CC_agb$CC[6], CC_agb$CC[7])))
CC_agb_tab_r2 <- as.data.frame(rbind(c("agb", CC_agb$CC[1], CC_agb$CC[2], CC_agb$CC[3])))
}

###PERCENT N#####

tune_ps = ParamSet$new(list(
  ParamDbl$new("eta", lower = 0.01, upper = 0.3),
  ParamInt$new("max_depth", lower = 2, upper = 6),
  ParamInt$new("min_child_weight", lower = 1, upper = 13),
  ParamDbl$new("subsample", lower = 0.15, upper = 0.8),
  ParamDbl$new("colsample_bytree", lower = 0.25, upper = 1)
))
#####
###for some reason this doesn't create a prediction with the new MLR3 package, so fitting model with old MLR package
tr<-data
tr<-dplyr::filter(tr, agm2.allom>10)
feature.names<-c("b1","b2","b3","b4","b5","b6","b7","ndvi", "pheno", "vari")
resp<-"percentN"
tr<-tr[tr$train.1==1,c(resp, feature.names)]
tr<-tr[complete.cases(tr),]

traintask = TaskRegr$new(id = "id", backend = tr[,c(resp, feature.names)], target = resp)
traintask$select(feature.names)

##set up learner
learner = lrn("regr.xgboost", predict_type = "response",  booster="gbtree",  objective   = "reg:squarederror", nrounds=500)
#resamp=rsmp("holdout", ratio=0.85)
resamp=rsmp("cv", folds=inner)
measure = msr("regr.mse")
tuner = mlr3tuning::tnr("grid_search", resolution = 25)
terminator = trm("stagnation", iters=4, threshold=0.1)

instance = mlr3tuning::TuningInstanceSingleCrit$new(
  task = traintask,
  learner = learner,
  resampling = resamp,
  measure = measure,
  search_space = tune_ps,
  terminator = terminator
)
instance

##run the optimization
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
#instance$result_learner_param_vals

##store the best parameters back in the learner
learner$param_set$values = instance$result_learner_param_vals
learner$train(traintask)

##see the result, eg, the best configuration of parameters found during the iterations
param <- instance$result_x_domain

resp<-tr$percentN
clf <- xgboost(data        = data.matrix(tr[,c(feature.names)]),
               label       = resp,
               booster="gbtree",
               nrounds     = 400,
               params = param,
               verbose=0)

mat <- xgb.importance (feature_names = feature.names,model = clf)
xgb.plot.importance (importance_matrix = mat[1:18]) 
p <-predict(clf, data.matrix(tr[, feature.names]))
plot(resp, p, ylim=c(0,max(c(p, resp), na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T))); abline(0,1)
sqrt((sum((resp-p)^2, na.rm=T))/length(p))

tr_results<- data.frame(na.omit(data$percentN[data$train.1==1]), p, na.omit(data$site[data$train.1==1 & data$percentN > 0]), na.omit(data$year[data$train.1==1 & data$percentN > 0]))
colnames(tr_results) <- c("resp", "p", "site", "year")
tr_results$year <- as.character(tr_results$year)
tr_metrics <- tr_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")

test <- data
test <- test %>% drop_na(percentN)
resp<-test$percentN[test$train.1==0]
resp <- na.omit(resp)
p <-predict(clf, data.matrix(test[test$train.1==0, feature.names]))
rmse<-sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)]))
print(rmse)
plot(resp, p, ylim=c(0,max(c(p, resp),na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T)), pch=19,
     col=trans.color[c(1:9)][factor(test$site[test$train.1==0])], 
     cex.lab=1.5, cex.axis=1.2)
abline(0,1)
mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(sprintf("%.2f",round(rmse,2))) ~ g ~ m^-2), cex=1.2)
legend("bottomright", legend=unique(test$site), 
       title=expression(underline(site)),title.adj = 0.25,
       col=trans.color[c(1:9)][factor(unique(test$site))], pch=15, bty="n")

##aboveground biophysical results
#rmse[n]
rmse <- sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)]))
#nRMSE
nrmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))/(max(resp, na.rm=T)-min(resp, na.rm=T))
#COV RMSE (alternate nRMSE)
cov_rmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))/(mean(resp, na.rm=T))
#mae[n]
mae <- (sum((p-resp), na.rm=T))/length(resp[!is.na(resp)])
#cor
cor <- cor(resp, p, use="pairwise.complete.obs")
#N
n_test <- nrow(data[!is.na(data$percentN)&data$train.1==0,])
n_train <- nrow(data[!is.na(data$percentN)&data$train.1==1,])
#table
n_tab <- as.data.frame(rbind(c(mae, rmse, nrmse, cov_rmse, cor, n_test, n_train)))
tab <- rbind(agb_tab, n_tab)

Nmod1<-clf 
#clf<-Nmod1
data<-data.frame(data)
data$predn.l8<-(predict(Nmod1, data.matrix(data[,feature.names])))
maxx<-max(c(data$percentN, data$predn.l8), na.rm=T)
minn<-min(c(data$percentN, data$predn.l8), na.rm=T)
sqrt((sum((data$percentN-data$predn.l8)^2, na.rm=T))/length(data$predn.l8[!is.na(data$predn.l8)]))

if(save_model==TRUE){
  ##save a copy of the model for future use
  saveRDS(Nmod1, paste0("output/xgb", v, ".0_foliarN_model.rds"))
  xgb.save(Nmod1, paste0("output/xgb", v, ".0_foliarN_model"))
  Nmod1<-xgb.load(paste0("output/xgb", v, ".0_foliarN_model"))
  
  ##create a graph of measured vs predicted for training and testing data
  jpeg(file=paste0("results/xgb", v, ".0_l8_N.jpg"), res=200, units="in",
       height=5, width=5)
  par(mar=c(4,5,0.5,0.75), oma=c(0,0.2,0.5,0.5), mfrow=c(1,1))
  plot(data$percentN, data$predn.l8, 
       xlab="measured % foliar N", 
       ylab ="predicted % foliar N",
       ylim=c(0,maxx),xlim=c(0,maxx),  pch=c(1,19)[factor(data$train)], 
       cex=c(1,0.8)[factor(data$train)], col=cols_all[c(1:9)][factor(data$site)],
       cex.lab=1.5, cex.axis=1.2)
  abline(0,1)
  mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(sprintf("%.2f",round(rmse,2))) ~ "%"), cex=1.2)
  legend("bottomright",inset=c(0,0.07),legend=c("train", "test"), pch=c(19,1), bty="n", title=expression(underline(data)), title.adj = 0.3)
  legend("bottomright", inset=c(0.2,0),legend=unique(data$site[!is.na(data$site)&data$site!="xtidb"]), 
         title=expression(underline(site)),title.adj = 0.25,
         col=trans.color[c(1:9)][factor(unique(data$site[!is.na(data$site)&data$site!="xtidb"]))], pch=15, bty="n")
  dev.off()
}

##results by site/year
test_results<- data.frame(na.omit(test$percentN[test$train.1==0]), p, na.omit(test$site[test$train.1==0]), na.omit(test$year[test$train.1==0]))
colnames(test_results) <- c("resp", "p", "site", "year")
test_results$year <- as.character(test_results$year)
test_metrics <- test_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")
uniques <- unique(test_metrics$value)
model <- "Foliar N"
version <- paste0(v, ".0")
for (i in uniques) {
  resp <- test_metrics$resp[test_metrics$value == i]
  p <- test_metrics$p[test_metrics$value == i]
  rmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))
  nrmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))/(max(resp, na.rm=T)-min(resp, na.rm=T))
  cor <- cor(resp, p)
  n_test <- length(resp)
  n_train <- length(tr_metrics$resp[tr_metrics$value == i])
  met_tab <- as.data.frame(rbind(c(i, rmse, nrmse, cor, n_test, n_train, model, version, run)))
  colnames(met_tab) <- c("value", "rmse", "nrmse", "cor", "n_test","n_train", "model", "version", "run")
  ag_tab <- rbind(ag_tab,met_tab)
}
#all tr & test data
all_tr_test_results_fn <- rbind(tr_results, test_results)
all_tr_test_results_fn$model <- "Foliar Nitrogen"
all_tr_test_results_fn$split <- c(rep("train", length(tr_results$resp)), rep("test", length(test_results$resp)))

## commonality analysis
if(v == 2){
CC_n=yhat::commonalityCoefficients(test_results, "p", list("site", "year"))
CC_n_tab_per <- as.data.frame(rbind(c("n", CC_n$CC[5], CC_n$CC[6], CC_n$CC[7])))
CC_n_tab_r2 <- as.data.frame(rbind(c("n", CC_n$CC[1], CC_n$CC[2], CC_n$CC[3])))
}

####LAI#####
tune_type <- "new"
tune_ps = ParamSet$new(list(
  ParamDbl$new("eta", lower = 0.01, upper = 0.15),
  ParamInt$new("max_depth", lower = 4, upper = 8),
  ParamInt$new("min_child_weight", lower = 1, upper = 5),
  ParamDbl$new("subsample", lower = 0.5, upper = 0.8),
  ParamDbl$new("colsample_bytree", lower = 0.5, upper = 0.9)
))
#####
data$train.1[data$train.1==0&data$site=="fluxa"][which.max(data$field_lai[data$train.1==0&data$site=="fluxa"])]<-1

tr<-data
tr<-dplyr::filter(tr, agm2.allom>10)
feature.names<-c("b1","b2","b3","b4","b5","b6","b7", "ndvi", "pheno", "vari", "ndmi", "gari", "rdvi", "mtvi2")
resp<-"field_lai"
tr<-tr[tr$train.1==1,c(resp, feature.names)]
tr<-tr[complete.cases(tr),]

traintask = TaskRegr$new(id = "id", backend = tr[,c(resp, feature.names)], target = resp)
traintask$select(feature.names)

##set up learner
learner = lrn("regr.xgboost", predict_type = "response",  booster="gbtree",  objective   = "reg:squarederror", nrounds=500)
#resamp=rsmp("holdout", ratio=0.85)
resamp=rsmp("cv", folds=inner)
measure = msr("regr.mse")
tuner = mlr3tuning::tnr("grid_search", resolution = 25)
terminator = trm("stagnation", iters=4, threshold=0.1)

instance = mlr3tuning::TuningInstanceSingleCrit$new(
  task = traintask,
  learner = learner,
  resampling = resamp,
  measure = measure,
  search_space = tune_ps,
  terminator = terminator
)
instance

##run the optimization
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
#instance$result_learner_param_vals

##store the best parameters back in the learner
learner$param_set$values = instance$result_learner_param_vals
learner$train(traintask)

##see the result, eg, the best configuration of parameters found during the iterations
param <- instance$result_x_domain

resp<-tr$field_lai
clf <- xgboost(data        = data.matrix(tr[,c(feature.names)]),
               label       = resp,
               booster="gbtree",
               nrounds     = 100,
               params = param,
               verbose=0)

mat <- xgb.importance (feature_names = feature.names,model = clf)
xgb.plot.importance (importance_matrix = mat[1:15]) 
p <-predict(clf, data.matrix(tr[, feature.names]))
plot(resp, p, ylim=c(0,max(c(p, resp), na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T))); abline(0,1)
sqrt((sum((resp-p)^2, na.rm=T))/length(p))

tr_results<- data.frame(na.omit(data$field_lai[data$train.1==1]), p, na.omit(data$site[data$train.1==1 & data$field_lai > 0]), na.omit(data$year[data$train.1==1 & data$field_lai > 0]))
colnames(tr_results) <- c("resp", "p", "site", "year")
tr_results$year <- as.character(tr_results$year)
tr_metrics <- tr_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")

test <- data
test <- test %>% drop_na(field_lai)
resp<-test$field_lai[test$train.1==0]
resp <- na.omit(resp)
p <-predict(clf, data.matrix(test[test$train.1==0, feature.names]))
rmse<-sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)]))
print(rmse)
plot(resp, p, ylim=c(0,max(c(p, resp),na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T)), pch=19,
     col=trans.color[c(1:9)][factor(test$site[test$train.1==0])],
     cex.lab=1.5, cex.axis=1.2)
abline(0,1)
mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(sprintf("%.1f",round(rmse,1))) ~ g ~ m^-2), cex=1.2)
legend("bottomright", legend=unique(test$site), 
       title=expression(underline(site)),title.adj = 0.25,
       col=trans.color[c(1:9)][factor(unique(test$site))], pch=15, bty="n")

##aboveground biophysical results
#rmse[n]
rmse <- sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)]))
#nRMSE
nrmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))/(max(resp, na.rm=T)-min(resp, na.rm=T))
#COV RMSE (alternate nRMSE)
cov_rmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))/(mean(resp, na.rm=T))
#mae[n]
mae <- (sum((p-resp), na.rm=T))/length(resp[!is.na(resp)])
#cor
cor <- cor(resp, p, use="pairwise.complete.obs")
#N
n_test <- nrow(data[!is.na(data$field_lai)&data$train.1==0,])
n_train <- nrow(data[!is.na(data$field_lai)&data$train.1==1,])
#table
lai_tab <- as.data.frame(rbind(c(mae, rmse, nrmse, cov_rmse, cor, n_test, n_train)))
tab <- rbind(tab, lai_tab)

laimod1<-clf 
#clf<-laimod1

data<-data.frame(data)
data$predlai.l8<-(predict(laimod1, data.matrix(data[,feature.names])))
maxx<-max(c(data$field_lai, data$predlai.l8), na.rm=T)
minn<-min(c(data$field_lai, data$predlai.l8), na.rm=T)
sqrt((sum((data$field_lai-data$predlai.l8)^2, na.rm=T))/length(data$predlai.l8[!is.na(data$predlai.l8)]))

if(save_model==TRUE){
  ##save a copy of the model for future use
  saveRDS(laimod1, paste0("output/xgb", v, ".0_lai_model.rds"))
  xgb.save(laimod1, paste0("output/xgb", v, ".0_lai_model"))
  laimod1<-xgb.load(paste0("output/xgb", v, ".0_lai_model"))
  
  ##create a graph of measured vs predicted for training and testing data
  jpeg(file=paste0("results/xgb", v, ".0_l8_lai.jpg"), res=200, units="in",
       height=5, width=5)
  par(mar=c(4,5,0.5,0.75), oma=c(0,0.2,0.5,0.5), mfrow=c(1,1))
  plot(data$field_lai, data$predlai.l8, 
       xlab="measured LAI", 
       ylab ="predicted LAI",
       ylim=c(0,maxx),xlim=c(0,maxx),  pch=c(1,19)[factor(data$train)], 
       cex=c(1,0.8)[factor(data$train)], col=trans.color[c(1:9)][factor(data$site)],
       cex.lab=1.5, cex.axis=1.2)
  abline(0,1)
  mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(sprintf("%.2f",round(rmse,2)))), cex=1.2)
  legend("bottomright",inset=c(0,0.07),legend=c("train", "test"), pch=c(19,1), bty="n", title=expression(underline(data)), title.adj = 0.3)
  legend("bottomright", inset=c(0.2,0),legend=unique(data$site[!is.na(data$site)&data$site!="xtidb"]), 
         title=expression(underline(site)),title.adj = 0.25,
         col=trans.color[c(1:9)][factor(unique(data$site[!is.na(data$site)&data$site!="xtidb"]))], pch=15, bty="n")
  dev.off()
}
print(rmse)

##results by site/year
test_results<- data.frame(na.omit(test$field_lai[test$train.1==0]), p, na.omit(test$site[test$train.1==0]), na.omit(test$year[test$train.1==0]))
colnames(test_results) <- c("resp", "p", "site", "year")
test_results$year <- as.character(test_results$year)
test_metrics <- test_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")
uniques <- unique(test_metrics$value)
model <- "LAI"
version <- paste0(v, ".0")
for (i in uniques) {
  resp <- test_metrics$resp[test_metrics$value == i]
  p <- test_metrics$p[test_metrics$value == i]
  rmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))
  nrmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))/(max(resp, na.rm=T)-min(resp, na.rm=T))
  cor <- cor(resp, p)
  n_test <- length(resp)
  n_train <- length(tr_metrics$resp[tr_metrics$value == i])
  met_tab <- as.data.frame(rbind(c(i, rmse, nrmse, cor, n_test, n_train, model, version, run)))
  colnames(met_tab) <- c("value", "rmse", "nrmse", "cor", "n_test","n_train", "model", "version", "run")
  ag_tab <- rbind(ag_tab,met_tab)
}
#all tr & test data
all_tr_test_results_lai <- rbind(tr_results, test_results)
all_tr_test_results_lai$model <- "Leaf Area Index"
all_tr_test_results_lai$split <- c(rep("train", length(tr_results$resp)), rep("test", length(test_results$resp)))

## commonality analysis
if(v == 2){
CC_lai=yhat::commonalityCoefficients(test_results, "p", list("site", "year"))
CC_lai_tab_per <- as.data.frame(rbind(c("lai", CC_lai$CC[5], CC_lai$CC[6], CC_lai$CC[7])))
CC_lai_tab_r2 <- as.data.frame(rbind(c("lai", CC_lai$CC[1], CC_lai$CC[2], CC_lai$CC[3])))

}else{}


####CHL#####
tune_type <- "old"
       tune_ps = ParamSet$new(list(
  ParamDbl$new("eta", lower = 0.05, upper = 0.3),
  ParamInt$new("max_depth", lower = 2, upper = 8),
  ParamInt$new("min_child_weight", lower = 1, upper = 10),
  ParamDbl$new("subsample", lower = 0.25, upper = 0.8),
  ParamDbl$new("colsample_bytree", lower = 0.15, upper = 0.8)
))
#####
tr<-data
tr<-dplyr::filter(tr, agm2.allom>10)
feature.names<-c("b1","b2", "b3","b4","b5","b6","b7", "gari", "ndmi", "pheno", "vari", "gndvi", "ndvi", "mtvi1")
resp<-"chlorophyll"
tr<-tr[tr$train.1==1,c(resp, feature.names)]
tr<-tr[complete.cases(tr),]

traintask = TaskRegr$new(id = "id", backend = tr[,c(resp, feature.names)], target = resp)
traintask$select(feature.names)

##set up learner
learner = lrn("regr.xgboost", predict_type = "response",  booster="gbtree",  objective   = "reg:squarederror", nrounds=500)
#resamp=rsmp("holdout", ratio=0.85)
resamp=rsmp("cv", folds=inner)
measure = msr("regr.mse")
tuner = mlr3tuning::tnr("grid_search", resolution = 25)
terminator = trm("stagnation", iters=4, threshold=0.1)

instance = mlr3tuning::TuningInstanceSingleCrit$new(
  task = traintask,
  learner = learner,
  resampling = resamp,
  measure = measure,
  search_space = tune_ps,
  terminator = terminator
)
instance

##run the optimization
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
#instance$result_learner_param_vals

##store the best parameters back in the learner
learner$param_set$values = instance$result_learner_param_vals
learner$train(traintask)

##see the result, eg, the best configuration of parameters found during the iterations
param <- instance$result_x_domain

resp<-tr$chlorophyll
clf <- xgboost(data        = data.matrix(tr[,c(feature.names)]),
               label       = resp,
               booster="gbtree",
               nrounds     = 100,
               params = param,
               verbose=0)

mat <- xgb.importance (feature_names = feature.names,model = clf)
xgb.plot.importance (importance_matrix = mat[1:20]) 
p <-predict(clf, data.matrix(tr[, feature.names]))
plot(resp, p, ylim=c(0,max(c(p, resp), na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T))); abline(0,1)
sqrt((sum((resp-p)^2, na.rm=T))/length(p))

tr_results<- data.frame(na.omit(data$chlorophyll[data$train.1==1]), p, na.omit(data$site[data$train.1==1 & data$chlorophyll > 0]), na.omit(data$year[data$train.1==1 & data$chlorophyll > 0]))
colnames(tr_results) <- c("resp", "p", "site", "year")
tr_results$year <- as.character(tr_results$year)
tr_metrics <- tr_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")

test <- data
test <- test %>% drop_na(chlorophyll)
resp<-test$chlorophyll[test$train.1==0]
resp <- na.omit(resp)
p <-predict(clf, data.matrix(test[test$train.1==0, feature.names]))
rmse<-sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)]))
print(rmse)
plot(resp, p, ylim=c(0,max(c(p, resp),na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T)), pch=19,
     col=trans.color[c(1:9)][factor(test$site[test$train.1==0])],
     cex.lab=1.5, cex.axis=1.2)
abline(0,1)
mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(sprintf("%.2f",round(rmse,2))) ), cex=1.2)
legend("bottomright", legend=unique(test$site), 
       title=expression(underline(site)),title.adj = 0.25,
       col=trans.color[c(1:9)][factor(unique(test$site))], pch=15, bty="n")


##aboveground biophysical results
#rmse[n]
rmse <- sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)]))
#nRMSE
nrmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))/(max(resp, na.rm=T)-min(resp, na.rm=T))
#COV RMSE (alternate nRMSE)
cov_rmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))/(mean(resp, na.rm=T))
#mae[n]
mae <- (sum((p-resp), na.rm=T))/length(resp[!is.na(resp)])
#cor
cor <- cor(resp, p, use="pairwise.complete.obs")
#N
n_test <- nrow(data[!is.na(data$chlorophyll)&data$train.1==0,])
n_train <- nrow(data[!is.na(data$chlorophyll)&data$train.1==1,])
#table
chl_tab <- as.data.frame(rbind(c(mae, rmse, nrmse, cov_rmse, cor, n_test, n_train)))

chlmod1<-clf 
#clf<-chlmod1

data<-data.frame(data)
data$predchl.l8<-(predict(chlmod1, data.matrix(data[,feature.names])))
maxx<-max(c(data$chlorophyll, data$predchl.l8), na.rm=T)
minn<-min(c(data$chlorophyll, data$predchl.l8), na.rm=T)
sqrt((sum((data$chlorophyll-data$predchl.l8)^2, na.rm=T))/length(data$predchl.l8[!is.na(data$predchl.l8)]))

if(save_model==TRUE){
  ##save a copy of the model for future use
  
  saveRDS(chlmod1, paste0("output/xgb", v, ".0_chl_model.rds"))
  xgb.save(chlmod1, paste0("output/xgb", v, ".0_chl_model"))
  chlmod1<-xgb.load(paste0("output/xgb", v, ".0_chl_model"))
  
  ##create a graph of measured vs predicted for training and testing data
  jpeg(file=paste0("results/xgb", v, ".0_l8_chl.jpg"), res=200, units="in",
       height=5, width=5)
  par(mar=c(4,5,0.5,0.75), oma=c(0,0.2,0.5,0.5), mfrow=c(1,1))
  plot(data$chlorophyll, data$predchl.l8, 
       xlab =expression(paste("measured CHL mg ", g^-1, sep="")) , 
       ylab =expression(paste("predicted CHL mg ", g^-1, sep="")) , 
       ylim=c(0,maxx),xlim=c(0,maxx),  pch=c(1,19)[factor(data$train)], 
       cex=c(1,0.8)[factor(data$train)], col=trans.color[c(1:9)][factor(data$site)],
       cex.lab=1.5, cex.axis=1.2)
  abline(0,1)
  mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(sprintf("%.2f",round(rmse,2))) ~ mg ~ g^-1), cex=1.2)
  legend("bottomright",inset=c(0,0.07),legend=c("train", "test"), pch=c(19,1), bty="n", title=expression(underline(data)), title.adj = 0.3)
  legend("bottomright", inset=c(0.2,0),legend=unique(data$site[!is.na(data$site)&data$site!="xtidb"]), 
         title=expression(underline(site)),title.adj = 0.25,
         col=trans.color[c(1:9)][factor(unique(data$site[!is.na(data$site)&data$site!="xtidb"]))], pch=15, bty="n")
  dev.off()
}

##results by site/year
test_results<- data.frame(na.omit(test$chlorophyll[test$train.1==0]), p, na.omit(test$site[test$train.1==0]), na.omit(test$year[test$train.1==0]))
colnames(test_results) <- c("resp", "p", "site", "year")
test_results$year <- as.character(test_results$year)
test_metrics <- test_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")
uniques <- unique(test_metrics$value)
model <- "Chlorophyll"
version <- paste0(v, ".0")
for (i in uniques) {
  resp <- test_metrics$resp[test_metrics$value == i]
  p <- test_metrics$p[test_metrics$value == i]
  rmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))
  nrmse <- (sqrt((sum((resp-p)^2, na.rm=T))/length(resp[!is.na(resp)])))/(max(resp, na.rm=T)-min(resp, na.rm=T))
  cor <- cor(resp, p)
  n_test <- length(resp)
  n_train <- length(tr_metrics$resp[tr_metrics$value == i])
  met_tab <- as.data.frame(rbind(c(i, rmse, nrmse, cor, n_test, n_train, model, version, run)))
  colnames(met_tab) <- c("value", "rmse", "nrmse", "cor", "n_test","n_train", "model", "version", "run")
  ag_tab <- rbind(ag_tab,met_tab)
}
#all tr & test data
all_tr_test_results_chl <- rbind(tr_results, test_results)
all_tr_test_results_chl$model <- "Chlorophyll"
all_tr_test_results_chl$split <- c(rep("train", length(tr_results$resp)), rep("test", length(test_results$resp)))

## commonality analysis
if(v == 2){
CC_chl=yhat::commonalityCoefficients(test_results, "p", list("site", "year"))
CC_chl_tab_per <- as.data.frame(rbind(c("chl", CC_chl$CC[5], CC_chl$CC[6], CC_chl$CC[7])))
CC_chl_tab_r2 <- as.data.frame(rbind(c("chl", CC_chl$CC[1], CC_chl$CC[2], CC_chl$CC[3])))
}else{}

#####

CC_data_per <- data.frame(matrix(ncol = 4, nrow = 0))
CC_data_per <- rbind(CC_data_per, CC_agb_tab_per, CC_n_tab_per, CC_lai_tab_per, CC_chl_tab_per)
colnames(CC_data_per) <- c("variable", "site", "year", "site:year")
CC_data_r2 <- data.frame(matrix(ncol = 4, nrow = 0))
CC_data_r2 <- rbind(CC_data_r2, CC_agb_tab_r2, CC_n_tab_r2, CC_lai_tab_r2, CC_chl_tab_r2)
colnames(CC_data_r2) <- c("variable", "site", "year", "site:year")

#create results table
tabl <- rbind(tab, chl_tab)
colnames(tabl) <- c("MAE", "RMSE", "nRMSE", "COV RMSE", "Correlation", "Ntest", "Ntrain")
tabl$Variable <- c("Aboveground Biomass g m-2", "% Foliar N", "Leaf Area Index", "% Foliar Chlorophyll")
tabl$nRMSE <- tabl$nRMSE*100
tabl$`COV RMSE` <- tabl$`COV RMSE`*100

#rounding
tabl <- tabl %>% mutate_at(vars(MAE, RMSE, Correlation), list(~ round(., 2)))
tabl <- tabl %>% mutate_at(vars(nRMSE, 'COV RMSE'), list(~ round(., 1)))
tabl <- tabl %>% mutate_at(vars(Ntest, Ntrain), list(~ round(., 0)))

#saving tables
if(save_model==TRUE){
  write.table(tabl, paste0("results/ag_table_berm", v, ".0.csv"), row.names = F, sep = ",")}
if(save_model==TRUE){
  write.table(CC_data_per, paste0("results/cc_per_ag_table_by_site_year_v", v, ".0.csv"), row.names = F, sep = ",")
  }
if(save_model==TRUE){
  write.table(CC_data_r2, paste0("results/cc_r2_ag_table_by_site_year_v", v, ".0.csv"), row.names = F, sep = ",")
  }
if(save_model==TRUE){
  write.table(ag_tab, paste0("results/ag_table_by_site_year_v", v, ".0.csv"), row.names = F, sep = ",")}
all_tr_test_results <- rbind(all_tr_test_results_ag, all_tr_test_results_fn, all_tr_test_results_lai, all_tr_test_results_chl)
if(save_model==TRUE){
  write.table(all_tr_test_results, paste0("results/ag_all_tr_test_results_v", v, ".0.csv"), row.names = F, sep = ",")}


#####


##save a copy of the prediction results for later use
if(save_model==T){
  write_csv (data, file =  paste0("output/xgb", v, "_predicted_biophysical_seagrant_landsat8.csv"))
}



