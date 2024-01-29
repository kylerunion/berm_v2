##after running the process groundtruth dataset for xgboost belowround biomass.r script, 
##run this script to actually build the model
##this model uses a nested resampling workflow with spatial crossvalidation to train, tune, build, and then test the model
##note that predicting novel pixels code (eg site level data) is at the end of this script

library(xgboost); library("paradox");  library("mlr3");library("mlr3learners"); 
library("data.table");  library("caret"); library("tidyverse"); library("mlr3tuning"); library("randomizr")
library(RColorBrewer)
library(kBET)

## WHICH VERSION ARE YOU CREATING
#v <- 1
v <- 2

##output and save new models and data sets?? Set to True/False below. 
##Note that True will overwrite the existing models and output datasets, so you'll want a backup if you don't intend this
save_model<-T
#save_model<-TRUE

##set up some colors for plots
color<-colorRampPalette(brewer.pal(name="Dark2", n = 8))(9)
trans.color<-addalpha(color,0.9)
all_sites <- c("deanc", "dupli", "fluxa", "fluxb", "folly", "huntc", "north", "skida", "ugami")
old_sites <- c("fluxa", "fluxb", "skida", "ugami")
old_sites_hc <- c("fluxa", "fluxb", "huntc", "skida", "ugami")
cols_all <- trans.color[match(all_sites, all_sites)]
cols_old <- trans.color[match(old_sites, all_sites)]
cols_old_hc <- trans.color[match(old_sites_hc, all_sites)]

data<-read_csv(paste0("/home/kyle/Documents/Documents/UT/Sapelo/berm_output/processed_bg", v, ".0.csv"))

data$id<-1:nrow(data)
##subset the data to dates where we have at least a year of antecedant conditions, 
##and remove dates without groundtruth data for this model
data<-data[data$date>"2014-06-01"&data$date<"2023-03-15",]
#data<-data[!(data$site == "fluxa"),]
with(data, pix[which.max(bgm2.core.stemscaled)])

#data<-data[!(data$pix == "faf"),]
#data<-data[!(data$pix == "faa"),]
#data<-data[!(data$pix == "fac"),]
#data<-data[!(data$pix == "fad"),]
#data<-data[!(data$pix == "fae"),]
#data<-data[!(data$pix == "fag"),]


data$group<-paste(data$site, data$date)

##subset to summer months/same as jberm
#data$mo <- as.numeric(data$mo)
#data <- subset(data, mo >= 6 & mo <= 10)

##we'll be creating training and testing datasets as part of that process, 
##so let's ensure we save a copy of the data in memory by creating a new data set called train
train <- data
#train <- subset(data, date < as.Date("2016-11-01"))
train$mo<-as.numeric(format(train$date, "%m"))  

##select to just the vars we want to consider as predictors for training the model
train<-dplyr::select(train, id,date,site, year, group,bgm2.core.stemscaled,
                     elevation, greendoy, growingday,
                     predag.allom.l8, predchl.l8, predlai.l8, photo, predn.l8,
                     lagg_ag_1, lagg_chl_1, lagg_lai_1, lagg_photo_1, lagg_N_1, lagg_ndvi_1,
                     roll_ag_2, roll_chl_2, roll_lai_2, roll_photo_2, roll_N_2, roll_ndvi_2,
                     roll_ag_3, roll_chl_3, roll_lai_3, roll_photo_3, roll_N_3, roll_ndvi_3,
                     roll_ag_4, roll_chl_4, roll_lai_4, roll_photo_4, roll_N_4, roll_ndvi_4,
                     roll_ag_5, roll_chl_5, roll_lai_5, roll_photo_5, roll_N_5, roll_ndvi_5,
                     diff_ag_1, diff_chl_1, diff_lai_1, diff_photo_1, diff_N_1, diff_ndvi_1,
                     diff_ag_2, diff_chl_2, diff_lai_2, diff_photo_2, diff_N_2, diff_ndvi_2,
                     diff_ag_3, diff_chl_3, diff_lai_3, diff_photo_3, diff_N_3, diff_ndvi_3,
                     diff_ag_4, diff_chl_4, diff_lai_4, diff_photo_4, diff_N_4, diff_ndvi_4,
                     diff_ag_5, diff_chl_5, diff_lai_5, diff_photo_5, diff_N_5, diff_ndvi_5,
                     diff_N_growing, diff_chl_growing, diff_ag_growing, diff_lai_growing,
                     growag, growchl, growlai, growphoto, growN,
                     local_hiwater, local_lowater, flood_time, 
                     lagg_lochi_1, roll_lochi_2, roll_lochi_3, roll_lochi_4,  roll_lochi_5,
                     lagg_loclo_1, roll_loclo_2, roll_loclo_3, roll_loclo_4,  roll_loclo_5,
                     lst, lagg_lst_1, roll_lst_2, roll_lst_3, roll_lst_4, roll_lst_5,
                     diff_lst_2,diff_lst_3,diff_lst_4,diff_lst_5,
                     dayl_mean, par_tot, prcp_mean, tmax_mean, tmin_mean,  vp_mean,
                     lagg_par_1, roll_par_2, roll_par_3, roll_par_4, roll_par_5,
                     lagg_prcp_1, roll_prcp_2, roll_prcp_3, roll_prcp_4, roll_prcp_5,
                     lagg_tmax_1, roll_tmax_2, roll_tmax_3,roll_tmax_4,roll_tmax_5,
                     lagg_tmin_1, roll_tmin_2, roll_tmin_3,roll_tmin_4, roll_tmin_5,
                     lagg_vp_1, roll_vp_2, roll_vp_3,roll_vp_4,roll_vp_5
)

##remove rows with missing data
train<-train[complete.cases(train),]
train$rowid<-as.numeric(rownames(train))

##feature names are the variables that will be predictor candidates, field data, dates, and row descriptors should be excluded
feature.names<-names(train)
feature.names<-feature.names[!(feature.names %in% c("date", "site","year","id", "rowid", "set","val","bgm2.core.stemscaled", 
                                                    "bgm2.allometric.rootshoot", "bgm2.allometric.rootgreenshoot",
                                                    "bgm2.allo.new.rootshoot", "bgm2.allo.new.rootgreenshoot", "mo",
                                                    "group"))]
feature.names
if(write_tx == T) {
  feature.names <- read_csv("/home/kyle/git/berm_tx/output/potential_features_bgb_tx_202310.csv")
feature.names <- feature.names[["feature.names"]]
feature.names <- feature.names[1:106]
  }
feature.names
  
##Set up the resampling scheme for nested resampling
##outer is the number of outer crossvalidations
outer<-5
##note that for this workflow, inner models are just 1 model, 
##eg we divide the outer model training data just once into a single inner training and testing model 
##the inner training model is for hyperparameter tunning and feature selection, 
##which can be tested against an inner model testing set; 
##This preserves the outer testing set for final model fitting

##set up clusters and blocks, used to keep certain observations together during nested resampling
clusters <-train$group ##note, group is paste(site, date)
blocks <-  train$site
rep<-1:outer
vars<-paste0("cv",rep)
cvs<-matrix(data=NA,nrow=nrow(train),ncol=length(vars))
###create the cross valdidated folds, which will keep obsers from same site and date together####
set.seed(29382) 
#set.seed(54564)
for( i in 1: outer){
  cvs[,i]<- block_and_cluster_ra(blocks = blocks,
                                 clusters = clusters,
                                 prob_each = c(.3, .7))
}                          
cvs[1:4,]
cvs<-data.frame(cvs); names(cvs)<-vars
head(table(clusters, cvs$cv1));head(table(blocks, cvs$cv1))
head(table(clusters, cvs$cv2));head(table(blocks, cvs$cv2))
#take a look at the clusters; blocks
table(cvs$cv1)
table(blocks, cvs$cv1)
table(blocks, cvs$cv2)
table(blocks, cvs$cv3)

##create columns in the larger dataset that will be named "train.1" or "train.2" etc, 
##and that will have 0 or 1 if the row is part of that training set, where 0 means it's in the testing set; 
##this will provide a record for later access
for (i in 1:outer){
  train[[paste0("train.", i)]] <- cvs[,i]
}

for (i in 1:outer) {
  data[[paste0("train.", i)]]<-NA
  data[[paste0("train.", i)]] [data$id %in% train$id[train[[paste0("train.", i)]] ==1 ] ]<-1
  data[[paste0("train.", i)]] [data$id %in% train$id[train[[paste0("train.", i)]] ==0 ] ]<-0
}
table(data$train.1);table(data$train.2)
table(cvs$cv1);table(cvs$cv2)
train[train$site=="ugami"&train$train.2==0,"date"]
train[train$site=="fluxb"&train$train.2==0,"date"]
train[train$site=="fluxa"&train$train.2==1,"date"]

####Select the features to use####
###XGBoost####
resp<-"bgm2.core.stemscaled" ## the y variable

##create a cut off for feature importance across the cross validations; if the feature is less important it will be discarded
cut.off<-0.005
rep<-1:outer
vars<-paste0("train.",rep)
feats.all<-matrix(data=NA, nrow=length(feature.names), ncol=length(vars))
feats.all<-data.frame(feats.all); names(feats.all)<-vars
feats.all$features<-feature.names

##create a list of parameters to focus on tuning, must specify type (Int=Integer, Dbl=continuous number, Fct=factor/character)
tune_ps = ParamSet$new(list(
  #ParamFct$new("booster",  levels = c("gbtree","gblinear")),
  ParamInt$new("gamma", lower = 3, upper = 10),
  ParamDbl$new("eta", lower = 0.01, upper = 0.3),
  ParamInt$new("max_depth", lower = 4, upper = 10),
  ParamInt$new("min_child_weight", lower = 1, upper = 10),
  ParamDbl$new("subsample", lower = 0.25, upper = 1),
  ParamDbl$new("colsample_bytree", lower = 0.25, upper = 0.8)
))
tune_ps

####go through all CVs to get average features #####

for (i in 1:outer) {
  tr<-train [train[[paste0("train.", i)]] ==1, ]
  traintask = TaskRegr$new(id = "bgbiomass", backend = tr[,c(resp, "group", feature.names)], target = resp)
  traintask$select(feature.names)
  ##set up learner
  ##select a resampling stategy for the datasets and preformance measure to use to tune the paramters
  ##set up how the tuner should iterate through the hyperparamters to try; 
  ##set up how the tuning should stop: here it stops after 20 iterations
  ##see other options at: https://mlr3book.mlr-org.com/tuning.html
  learner = lrn("regr.xgboost", predict_type = "response",  booster="gbtree",  
                objective   = "reg:squarederror", nrounds=300, verbose=0)
  #resamp = rsmp("holdout")
  resamp=rsmp("holdout", ratio=0.85)
  measure = msr("regr.mse")
  tuner = mlr3tuning::tnr("grid_search", resolution = 20)
  #tuner = tnr("random_search")
  #terminator = trm("evals", n_evals = 10)
  terminator = trm("stagnation", iters=3, threshold=0.1)
  
  instance = mlr3tuning::TuningInstanceSingleCrit$new(
    task = traintask,
    learner = learner,
    resampling = resamp,
    measure = measure,
    search_space = tune_ps,
    terminator = terminator
  )
  instance
  tuner = mlr3tuning::tnr("grid_search", resolution = 20)
  #tuner = tnr("random_search")
  
  ##run the optimization
  tuner$optimize(instance)
  
  ##see the result, eg, the best configuration of parameters found during the iterations
  instance$result_learner_param_vals
  #see the benchmark results of all evaluations
  instance$archive$data
  ##store the best parameters back in the learner
  learner$param_set$values = instance$result_learner_param_vals
  learner$train(traintask)
  
  ##gather the features to use
  y=learner$importance(); x=names(y)
  mar=c(3,15,2,2)
  features.all<-names(y[y>cut.off])
  features.all
  
  if(length(features.all)<length(feature.names)){
    feats.all[,i]<- c(features.all, rep(NA, length(feature.names)-length(features.all)))
  } else{
    feats.all[,i]<- features.all
  }
}

##summarize the features from the gathered features
features.all<-data.frame(var=feature.names, train.1=NA, train.2=NA, train.3=NA)
for( i in 1: outer){
  features.all [[paste0("train.", i)]]<-ifelse(features.all$var %in% feats.all[[paste0("train.", i)]], 1,0)
}

##create objects that record which features to retain and summarizes their importance across the crossvalidated trials
features.all$vote<-rowSums(features.all[,2:(outer+1)])
features.all.use<-features.all$var[features.all$vote>0.5*outer]

features<-features.all.use
features


###looking at all the data
library(ggplot2)
#ggplot(data = data, aes_string(x = var[19], y = as.character("bgm2.core.stemscaled"))) + 
 # geom_point(aes(color = site))


#####Train the models with the selected features#####
###train model 1#####
run <- 1
resp<-"bgm2.core.stemscaled" ## the y variable
tr<-train[train$train.1==1,]
traintask = TaskRegr$new(id = "bgbiomass", backend = tr[,c(resp,  features)], target = resp)
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
instance$result_learner_param_vals
param <- instance$result_x_domain

resp<-tr$bgm2.core.stemscaled
clf <- xgboost(data        = data.matrix(tr[,features]),
               label       = resp,
               booster="gbtree",
               nrounds     = 125,
               params = param,
               verbose=0)

mat <- xgb.importance (feature_names = features,model = clf)
xgb.plot.importance (importance_matrix = mat[1:30]) 

##plot the training data predicted vs measured
p <-(as.integer(round(predict(clf, data.matrix(train[train$train.1==1,features])))))
plot(train$bgm2.core.stemscaled[train$train.1==1], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.1==1]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.1==1]))), col=factor(train$site[train$train.1==1])); abline(0,1)
sqrt((sum((train$bgm2.core.stemscaled[train$train.1==1]-p)^2, na.rm=T))/length(p))

tr_results<- data.frame(na.omit(train$bgm2.core.stemscaled[train$train.1==1]), p, na.omit(train$site[train$train.1==1 & train$bgm2.core.stemscaled > 0]), na.omit(train$year[train$train.1==1 & train$bgm2.core.stemscaled > 0]))
colnames(tr_results) <- c("resp", "p", "site", "year")
tr_results$year <- as.character(tr_results$year)
tr_metrics <- tr_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")

##plot the testing data predicted vs measured
p <-(as.integer(round(predict(clf, data.matrix(train[train$train.1==0,features])))))
plot(train$bgm2.core.stemscaled[train$train.1==0], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.1==0]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.1==0]))), col=factor(train$site[train$train.1==0])); abline(0,1)
rmset<-sqrt((sum((train$bgm2.core.stemscaled[train$train.1==0]-p)^2, na.rm=T))/length(p))
print(rmset)

## results
resp <- train$bgm2.core.stemscaled[train$train.1==0]
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
n_test <- nrow(train[!is.na(train$bgm2.core.stemscaled)&train$train.1==0,])
n_train <- nrow(train[!is.na(train$bgm2.core.stemscaled)&train$train.1==1,])

bg_1_tab <- as.data.frame(rbind(c(mae, rmse, nrmse, cov_rmse, cor, n_test, n_train, run)))


##results by site/year
test_results<- data.frame(na.omit(train$bgm2.core.stemscaled[train$train.1==0]), p, na.omit(train$site[train$train.1==0]), na.omit(train$year[train$train.1==0]))
colnames(test_results) <- c("resp", "p", "site", "year")
test_results$year <- as.character(test_results$year)
test_metrics <- test_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")
uniques <- unique(test_metrics$value)
bg_tab <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(bg_tab) <- c("value", "rmse", "nrmse", "corr", "model", "version", "run")
model <- "Belowground Biomass"
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
  bg_tab <- rbind(bg_tab,met_tab)
}
#all tr & test data
all_tr_test_results_bg_r1 <- rbind(tr_results, test_results)
all_tr_test_results_bg_r1$model <- "Belowground Biomass"
all_tr_test_results_bg_r1$split <- c(rep("train", length(tr_results$resp)), rep("test", length(test_results$resp)))
all_tr_test_results_bg_r1$run <- run

## commonality analysis
CC_data <- as.data.frame(matrix(ncol = 5, nrow = 0))
CC_bg1=yhat::commonalityCoefficients(test_results, "p", list("site", "year"))
CC_bg1_per_tab <- as.data.frame(rbind(c("bg", CC_bg1$CC[5], CC_bg1$CC[6], CC_bg1$CC[7], paste0(run))))
CC_bg1_r2_tab <- as.data.frame(rbind(c("bg", CC_bg1$CC[1], CC_bg1$CC[2], CC_bg1$CC[3], paste0(run))))
CC_data_per <- rbind(CC_data, CC_bg1_per_tab)
CC_data_r2 <- rbind(CC_data, CC_bg1_r2_tab)


##save the model as it's own object
xgmod1<-clf

#clf<-xgmod1

#### Train set 2 #####
run <- 2
tr<-train [train$train.2==1, ]
resp<-"bgm2.core.stemscaled" ## the y variable
traintask = TaskRegr$new(id = "bgbiomass", backend = tr[,c(resp, features)], target = resp)
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
instance$result_learner_param_vals
param <- instance$result_x_domain


##Code with early stopping rounds based on hold out set performance
resp<-tr$bgm2.core.stemscaled
clf <- xgboost(data        = data.matrix(tr[,features]),
               label       = resp,
               booster="gbtree",
               nrounds     = 125,
               params = param, 
               verbose=0)

mat <- xgb.importance (feature_names = features,model = clf)
xgb.plot.importance (importance_matrix = mat[1:30]) 
p <-(as.integer(round(predict(clf, data.matrix(train[train$train.2==1,features])))))
plot(train$bgm2.core.stemscaled[train$train.2==1], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.2==1]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.2==1]))), col=factor(train$site[train$train.2==1])); abline(0,1)
sqrt((sum((train$bgm2.core.stemscaled[train$train.2==1]-p)^2, na.rm=T))/length(p))

tr_results<- data.frame(na.omit(train$bgm2.core.stemscaled[train$train.2==1]), p, na.omit(train$site[train$train.2==1 & train$bgm2.core.stemscaled > 0]), na.omit(train$year[train$train.2==1 & train$bgm2.core.stemscaled > 0]))
colnames(tr_results) <- c("resp", "p", "site", "year")
tr_results$year <- as.character(tr_results$year)
tr_metrics <- tr_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")

p <-(as.integer(round(predict(clf, data.matrix(train[train$train.2==0,features])))))
plot(train$bgm2.core.stemscaled[train$train.2==0], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.2==0]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.2==0]))), col=factor(train$site[train$train.2==0])); abline(0,1)
rmset<-sqrt((sum((train$bgm2.core.stemscaled[train$train.2==0]-p)^2, na.rm=T))/length(p))
print(rmset)

## results
resp <- train$bgm2.core.stemscaled[train$train.2==0]
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
n_test <- nrow(train[!is.na(train$bgm2.core.stemscaled)&train$train.2==0,])
n_train <- nrow(train[!is.na(train$bgm2.core.stemscaled)&train$train.2==1,])

bg_2_tab <- as.data.frame(rbind(c(mae, rmse, nrmse, cov_rmse, cor, n_test, n_train, run)))


##results by site/year
test_results<- data.frame(na.omit(train$bgm2.core.stemscaled[train$train.2==0]), p, na.omit(train$site[train$train.2==0]), na.omit(train$year[train$train.2==0]))
colnames(test_results) <- c("resp", "p", "site", "year")
test_results$year <- as.character(test_results$year)
test_metrics <- test_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")
uniques <- unique(test_metrics$value)
model <- "Belowground Biomass"
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
  bg_tab <- rbind(bg_tab,met_tab)
}
#all tr & test data
all_tr_test_results_bg_r2 <- rbind(tr_results, test_results)
all_tr_test_results_bg_r2$model <- "Belowground Biomass"
all_tr_test_results_bg_r2$split <- c(rep("train", length(tr_results$resp)), rep("test", length(test_results$resp)))
all_tr_test_results_bg_r2$run <- run

## commonality analysis
CC_bg2=yhat::commonalityCoefficients(test_results, "p", list("site", "year"))
CC_bg2_per_tab <- as.data.frame(rbind(c("bg", CC_bg2$CC[5], CC_bg2$CC[6], CC_bg2$CC[7], paste0(run))))
CC_bg2_r2_tab <- as.data.frame(rbind(c("bg", CC_bg2$CC[1], CC_bg2$CC[2], CC_bg2$CC[3], paste0(run))))
CC_data_per <- rbind(CC_data_per, CC_bg2_per_tab)
CC_data_r2 <- rbind(CC_data_r2, CC_bg2_r2_tab)

xgmod2<-clf
#clf<-xgmod2

####Train set 3 ####
run = 3
## this moves the largest ugami bgb from the testing set into the training set
train$train.3[train$train.3==0&train$site=="ugami"][which.max(train$bgm2.core.stemscaled[train$train.3==0&train$site=="ugami"])]<-1
## so let's move the min over the the testing set because otherwise there's only 1 testing data
train$train.3[train$train.3==1&train$site=="ugami"][which.min(train$bgm2.core.stemscaled[train$train.3==1&train$site=="ugami"])]<-0
tr<-train [train$train.3==1, ]
resp<-"bgm2.core.stemscaled" ## the y variable
traintask = TaskRegr$new(id = "bgbiomass", backend = tr[,c(resp, features)], target = resp)
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
instance$result_learner_param_vals
param <- instance$result_x_domain


##Code with early stopping rounds based on hold out set performance
resp<-tr$bgm2.core.stemscaled
clf <- xgboost(data        = data.matrix(tr[,features]),
               label       = resp,
               booster="gbtree",
               nrounds     = 125,
               params = param, 
               verbose=0)

mat <- xgb.importance (feature_names = features,model = clf)
xgb.plot.importance (importance_matrix = mat[1:40]) 
p <-(as.integer(round(predict(clf, data.matrix(train[train$train.3==1,features])))))
plot(train$bgm2.core.stemscaled[train$train.3==1], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.3==1]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.3==1]))), col=factor(train$site[train$train.3==1])); abline(0,1)
sqrt((sum((train$bgm2.core.stemscaled[train$train.3==1]-p)^2, na.rm=T))/length(p))

tr_results<- data.frame(na.omit(train$bgm2.core.stemscaled[train$train.3==1]), p, na.omit(train$site[train$train.3==1 & train$bgm2.core.stemscaled > 0]), na.omit(train$year[train$train.3==1 & train$bgm2.core.stemscaled > 0]))
colnames(tr_results) <- c("resp", "p", "site", "year")
tr_results$year <- as.character(tr_results$year)
tr_metrics <- tr_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")

p <-(as.integer(round(predict(clf, data.matrix(train[train$train.3==0,features])))))
plot(train$bgm2.core.stemscaled[train$train.3==0], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.3==0]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.3==0]))), col=factor(train$site[train$train.3==0])); abline(0,1)
rmset<-sqrt((sum((train$bgm2.core.stemscaled[train$train.3==0]-p)^2, na.rm=T))/length(p))
print(rmset)

## results
resp <- train$bgm2.core.stemscaled[train$train.3==0]
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
n_test <- nrow(train[!is.na(train$bgm2.core.stemscaled)&train$train.3==0,])
n_train <- nrow(train[!is.na(train$bgm2.core.stemscaled)&train$train.3==1,])

bg_3_tab <- as.data.frame(rbind(c(mae, rmse, nrmse, cov_rmse, cor, n_test, n_train, run)))


##results by site/year
test_results<- data.frame(na.omit(train$bgm2.core.stemscaled[train$train.3==0]), p, na.omit(train$site[train$train.3==0]), na.omit(train$year[train$train.3==0]))
colnames(test_results) <- c("resp", "p", "site", "year")
test_results$year <- as.character(test_results$year)
test_metrics <- test_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")
uniques <- unique(test_metrics$value)
model <- "Belowground Biomass"
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
  bg_tab <- rbind(bg_tab,met_tab)
}
#all tr & test data
all_tr_test_results_bg_r3 <- rbind(tr_results, test_results)
all_tr_test_results_bg_r3$model <- "Belowground Biomass"
all_tr_test_results_bg_r3$split <- c(rep("train", length(tr_results$resp)), rep("test", length(test_results$resp)))
all_tr_test_results_bg_r3$run <- run

## commonality analysis
CC_bg3=yhat::commonalityCoefficients(test_results, "p", list("site", "year"))
CC_bg3_per_tab <- as.data.frame(rbind(c("bg", CC_bg3$CC[5], CC_bg3$CC[6], CC_bg3$CC[7], paste0(run))))
CC_bg3_r2_tab <- as.data.frame(rbind(c("bg", CC_bg3$CC[1], CC_bg3$CC[2], CC_bg3$CC[3], paste0(run))))
CC_data_per <- rbind(CC_data_per, CC_bg3_per_tab)
CC_data_r2 <- rbind(CC_data_r2, CC_bg3_r2_tab)

xgmod3<-clf
#clf<-xgmod3

####Train set 4 ####
run <- 4
tr<-train [train$train.4==1, ]
resp<-"bgm2.core.stemscaled" ## the y variable
traintask = TaskRegr$new(id = "bgbiomass", backend = tr[,c(resp, features)], target = resp)
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
instance$result_learner_param_vals
param <- instance$result_x_domain


##Code with early stopping rounds based on hold out set performance
resp<-tr$bgm2.core.stemscaled
clf <- xgboost(data        = data.matrix(tr[,features]),
               label       = resp,
               booster="gbtree",
               nrounds     = 125,
               params = param, 
               verbose=0)

mat <- xgb.importance (feature_names = features,model = clf)
xgb.plot.importance (importance_matrix = mat[1:30]) 
p <-(as.integer(round(predict(clf, data.matrix(train[train$train.4==1,features])))))
plot(train$bgm2.core.stemscaled[train$train.4==1], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.4==1]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.4==1]))), col=factor(train$site[train$train.4==1])); abline(0,1)
sqrt((sum((train$bgm2.core.stemscaled[train$train.4==1]-p)^2, na.rm=T))/length(p))

tr_results<- data.frame(na.omit(train$bgm2.core.stemscaled[train$train.4==1]), p, na.omit(train$site[train$train.4==1 & train$bgm2.core.stemscaled > 0]), na.omit(train$year[train$train.4==1 & train$bgm2.core.stemscaled > 0]))
colnames(tr_results) <- c("resp", "p", "site", "year")
tr_results$year <- as.character(tr_results$year)
tr_metrics <- tr_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")

p <-(as.integer(round(predict(clf, data.matrix(train[train$train.4==0,features])))))
plot(train$bgm2.core.stemscaled[train$train.4==0], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.4==0]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.4==0]))), col=factor(train$site[train$train.4==0])); abline(0,1)
rmset<-sqrt((sum((train$bgm2.core.stemscaled[train$train.4==0]-p)^2, na.rm=T))/length(p))
print(rmset)

## results
resp <- train$bgm2.core.stemscaled[train$train.4==0]
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
n_test <- nrow(train[!is.na(train$bgm2.core.stemscaled)&train$train.4==0,])
n_train <- nrow(train[!is.na(train$bgm2.core.stemscaled)&train$train.4==1,])

bg_4_tab <- as.data.frame(rbind(c(mae, rmse, nrmse, cov_rmse, cor, n_test, n_train, run)))


##results by site/year
test_results<- data.frame(na.omit(train$bgm2.core.stemscaled[train$train.4==0]), p, na.omit(train$site[train$train.4==0]), na.omit(train$year[train$train.4==0]))
colnames(test_results) <- c("resp", "p", "site", "year")
test_results$year <- as.character(test_results$year)
test_metrics <- test_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")
uniques <- unique(test_metrics$value)
model <- "Belowground Biomass"
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
  bg_tab <- rbind(bg_tab,met_tab)
}
#all tr & test data
all_tr_test_results_bg_r4 <- rbind(tr_results, test_results)
all_tr_test_results_bg_r4$model <- "Belowground Biomass"
all_tr_test_results_bg_r4$split <- c(rep("train", length(tr_results$resp)), rep("test", length(test_results$resp)))
all_tr_test_results_bg_r4$run <- run

## commonality analysis
CC_bg4=yhat::commonalityCoefficients(test_results, "p", list("site", "year"))
CC_bg4_per_tab <- as.data.frame(rbind(c("bg", CC_bg4$CC[5], CC_bg4$CC[6], CC_bg4$CC[7], paste0(run))))
CC_bg4_r2_tab <- as.data.frame(rbind(c("bg", CC_bg4$CC[1], CC_bg4$CC[2], CC_bg4$CC[3], paste0(run))))
CC_data_per <- rbind(CC_data_per, CC_bg4_per_tab)
CC_data_r2 <- rbind(CC_data_r2, CC_bg4_r2_tab)

xgmod4<-clf
#clf<-xgmod4

####Train set 5####
run <- "5"
tr<-train [train$train.5==1, ]
resp<-"bgm2.core.stemscaled" ## the y variable
traintask = TaskRegr$new(id = "bgbiomass", backend = tr[,c(resp, features)], target = resp)
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
instance$result_learner_param_vals
param <- instance$result_x_domain


##Code with early stopping rounds based on hold out set performance
resp<-tr$bgm2.core.stemscaled
clf <- xgboost(data        = data.matrix(tr[,features]),
               label       = resp,
               booster="gbtree",
               nrounds     = 125,
               params = param, 
               verbose=0)

mat <- xgb.importance (feature_names = features,model = clf)
xgb.plot.importance (importance_matrix = mat[1:30]) 
p <-(as.integer(round(predict(clf, data.matrix(train[train$train.5==1,features])))))
plot(train$bgm2.core.stemscaled[train$train.5==1], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.5==1]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.5==1]))), col=factor(train$site[train$train.5==1])); abline(0,1)
sqrt((sum((train$bgm2.core.stemscaled[train$train.5==1]-p)^2, na.rm=T))/length(p))

tr_results<- data.frame(na.omit(train$bgm2.core.stemscaled[train$train.5==1]), p, na.omit(train$site[train$train.5==1 & train$bgm2.core.stemscaled > 0]), na.omit(train$year[train$train.5==1 & train$bgm2.core.stemscaled > 0]))
colnames(tr_results) <- c("resp", "p", "site", "year")
tr_results$year <- as.character(tr_results$year)
tr_metrics <- tr_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")

p <-(as.integer(round(predict(clf, data.matrix(train[train$train.5==0,features])))))
plot(train$bgm2.core.stemscaled[train$train.5==0], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.5==0]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.5==0]))), col=factor(train$site[train$train.5==0])); abline(0,1)
rmset<-sqrt((sum((train$bgm2.core.stemscaled[train$train.5==0]-p)^2, na.rm=T))/length(p))
print(rmset)

## results
resp <- train$bgm2.core.stemscaled[train$train.5==0]
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
n_test <- nrow(train[!is.na(train$bgm2.core.stemscaled)&train$train.5==0,])
n_train <- nrow(train[!is.na(train$bgm2.core.stemscaled)&train$train.5==1,])

bg_5_tab <- as.data.frame(rbind(c(mae, rmse, nrmse, cov_rmse, cor, n_test, n_train, run)))

##results by site/year
test_results<- data.frame(na.omit(train$bgm2.core.stemscaled[train$train.5==0]), p, na.omit(train$site[train$train.5==0]), na.omit(train$year[train$train.5==0]))
colnames(test_results) <- c("resp", "p", "site", "year")
test_results$year <- as.character(test_results$year)
test_metrics <- test_results %>%
  pivot_longer(cols = c("site", "year"), names_to = "variable", values_to = "value")
uniques <- unique(test_metrics$value)
model <- "Belowground Biomass"
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
  bg_tab <- rbind(bg_tab,met_tab)
}
#all tr & test data
all_tr_test_results_bg_r5 <- rbind(tr_results, test_results)
all_tr_test_results_bg_r5$model <- "Belowground Biomass"
all_tr_test_results_bg_r5$split <- c(rep("train", length(tr_results$resp)), rep("test", length(test_results$resp)))
all_tr_test_results_bg_r5$run <- run

## commonality analysis
CC_bg5=yhat::commonalityCoefficients(test_results, "p", list("site", "year"))
CC_bg5_per_tab <- as.data.frame(rbind(c("bg", CC_bg5$CC[5], CC_bg5$CC[6], CC_bg5$CC[7], paste0(run))))
CC_bg5_r2_tab <- as.data.frame(rbind(c("bg", CC_bg5$CC[1], CC_bg5$CC[2], CC_bg5$CC[3], paste0(run))))
CC_data_per <- rbind(CC_data_per, CC_bg5_per_tab)
CC_data_r2 <- rbind(CC_data_r2, CC_bg5_r2_tab)


xgmod5<-clf
#clf<-xgmod5

#create results table
tabl <- rbind(bg_1_tab, bg_2_tab, bg_3_tab, bg_4_tab, bg_5_tab)
colnames(tabl) <- c("MAE", "RMSE", "nRMSE", "COV RMSE", "Correlation", "Ntest", "Ntrain", "Run")
if(save_model==TRUE){
  write.table(tabl, paste0("results/bg_table_berm", v, ".0.csv"), row.names = F, sep = ",")
}

#create table of metrics
if(save_model==TRUE){
  write.table(bg_tab, paste0("results/bg_table_by_site_year_v", v, ".0.csv"), row.names = F, sep = ",")
}

#commonality analysis table
colnames(CC_data_per) <- c("variable", "site", "year", "site:year", "run")
colnames(CC_data_r2) <- c("variable", "site", "year", "site:year", "run")
if(save_model==TRUE){
  write.table(CC_data_per, paste0("results/cc_per_bg_table_by_site_year_v", v, ".0.csv"), row.names = F, sep = ",")
  write.table(CC_data_r2, paste0("results/cc_r2_bg_table_by_site_year_v", v, ".0.csv"), row.names = F, sep = ",")
}

# all train test results table
all_tr_test_results <- rbind(all_tr_test_results_bg_r1, all_tr_test_results_bg_r2, all_tr_test_results_bg_r3, all_tr_test_results_bg_r4, all_tr_test_results_bg_r5)
if(save_model==TRUE){
  write.table(all_tr_test_results, paste0("results/bg_all_tr_test_results_v", v, ".0.csv"), row.names = F, sep = ",")
}



###gather features importance from the final fitted models from the outer model fitting####
rep<-1:outer
vars<-paste0("train.", rep)
feats.all.import<-matrix(data=NA, nrow=length(features), ncol=length(vars))
feats.all.import<-data.frame(feats.all.import); names(feats.all.import)<-vars
feats.all.import$features<-features

mat <- xgb.importance (feature_names = features,model = xgmod1)
x=mat$Feature; y=mat$Gain
import<-data.frame(features=x, import=as.numeric(y))
feats.all.import[match(import$features, feats.all.import$features),1]<- import$import

mat <- xgb.importance (feature_names = features,model = xgmod2)
x=mat$Feature; y=mat$Gain
import<-data.frame(features=x, import=as.numeric(y))
feats.all.import[match(import$features, feats.all.import$features),2]<- import$import

mat <- xgb.importance (feature_names = features,model = xgmod3)
x=mat$Feature; y=mat$Gain
import<-data.frame(features=x, import=as.numeric(y))
feats.all.import[match(import$features, feats.all.import$features),3]<- import$import

mat <- xgb.importance (feature_names = features,model = xgmod4)
x=mat$Feature; y=mat$Gain
import<-data.frame(features=x, import=as.numeric(y))
feats.all.import[match(import$features, feats.all.import$features),4]<- import$import

mat <- xgb.importance (feature_names = features,model = xgmod5)
x=mat$Feature; y=mat$Gain
import<-data.frame(features=x, import=as.numeric(y))
feats.all.import[match(import$features, feats.all.import$features),5]<- import$import

data$pred.1[!is.na(data$train.1)]<-(as.integer(round(predict(xgmod1, data.matrix(data[!is.na(data$train.1),features])))))
data$pred.2[!is.na(data$train.2)]<-(as.integer(round(predict(xgmod2, data.matrix(data[!is.na(data$train.2),features])))))
data$pred.3[!is.na(data$train.3)]<-(as.integer(round(predict(xgmod3, data.matrix(data[!is.na(data$train.3),features])))))
data$pred.4[!is.na(data$train.4)]<-(as.integer(round(predict(xgmod4, data.matrix(data[!is.na(data$train.4),features])))))
data$pred.5[!is.na(data$train.5)]<-(as.integer(round(predict(xgmod5, data.matrix(data[!is.na(data$train.5),features])))))

y<-dplyr::select(data, train.1, train.2, train.3,  train.4, train.5,
                 pred.1, pred.2, pred.3, pred.4, pred.5, bgm2.core.stemscaled)

rmse<-data.frame(set=1:(outer), train=rep(NA, (outer)), test=rep(NA, (outer)), 
                 trainmean=rep(NA, (outer)), trainmed=rep(NA, (outer)))
for (i in 1:outer){
  tr<-y[,c(i,i+outer,ncol(y))]
  colnames(tr)<-c("train", "pred", "resp")
  trainset<-tr[which(tr$train==1), ]
  testset<-tr[which(tr$train==0),]
  
  rmse$train[i]<-sqrt((sum((trainset$resp-trainset$pred)^2, na.rm=T))/length(trainset$resp[!is.na(trainset$resp)]))
  rmse$test[i]<-sqrt((sum((testset$resp-testset$pred)^2, na.rm=T))/length(testset$resp[!is.na(testset$resp)]))
  rmse$trainmean[i]<-mean(trainset$pred, na.rm=T)
  rmse$trainmed[i]<-median(trainset$pred, na.rm=T)
}
summary(rmse)
rmse<-arrange(rmse,set)
rmse

data$predbg.xgb<-rowMeans(data[,c("pred.1","pred.2", "pred.3", "pred.4", "pred.5")])
maxx<-max(c(data$bgm2.core.stemscaled, data$predbg.xgb), na.rm=T)
minn<-min(c(data$bgm2.core.stemscaled, data$predbg.xgb), na.rm=T)
sqrt((sum((data$bgm2.core.stemscaled-data$predbg.xgb)^2, na.rm=T))/length(data$predbg.xgb[!is.na(data$predbg.xgb)]))



write_csv(rmse, file = paste0("results/bg_rmse_table", v, ".0.csv"))



#jpeg(file="results/xgb_l8_bg_all_months.jpg", units="in",
 #    height=5.25, width=7.25, res=200)
#par(mar=c(4,5,1,0.75), oma=c(0,0.2,0.5,0.5), mfrow=c(1,2))
#layout(matrix(c(1,2),1,2,byrow = TRUE), c(1,2), TRUE)
#boxplot(rmse[,   c("train","test")], 
#        ylab =expression(paste("RMSE g ", m^-2, sep="")), 
#        cex.lab=1.5, cex.axis=1.2)
#mtext(side=3, adj=0.95, line=0.25, "a", cex=1.2)

jpeg(file=paste0("results/xgb", v,".0_l8_bg_all_months.jpg"), units="in",
     height=6, width=6, res=300)
par(mar=c(4,5,1,0.75), oma=c(0,0.2,0.5,0.5), mfrow=c(1,1))
plot( data$predbg.xgb, data$bgm2.core.stemscaled,
      xlab=expression(paste("predicted BG biomass g ", m^-2, sep="")), 
      ylab =expression(paste("measured BG biomass g ", m^-2, sep="")) , 
      ylim=c(minn,maxx),xlim=c(minn,maxx),  pch=c(19), 
      cex=c(1), col=trans.color[c(1:9)][factor(data$site)],
      cex.lab=1.5, cex.axis=1.2)
abline(0,1)
#mtext(side=3, adj=0.05, line=-2, bquote(RMSE[train]: ~ .(round(mean(rmse$train, na.rm=T),1)) ~ g~m^-2), cex=1.2)
#mtext(side=3, adj=0.05, line=-4, bquote(RMSE[test]: ~ .(round(mean(rmse$test, na.rm=T),1)) ~ g~m^-2), cex=1.2)
mtext(side=3, adj=0.05, line=-4, bquote(RMSE[test]: ~ .(round(mean(rmse$test, na.rm=T),1)) ~ g~m^-2), cex=2)
legend("bottomright", inset=c(0.2,0),legend=unique(data$site[!is.na(data$site)]), title=expression(underline(site)),title.adj = 0.25,
       col=trans.color[c(1:9)][factor(unique(data$site[!is.na(data$site)]))], pch=15, bty="n")
#mtext(side=3, adj=0.97, line=0.25, "b", cex=1.2)
dev.off()

write_csv(data, "output/bg_model_output.csv")


feat<-data.frame(all=c(features, rep(NA,length(feature.names)-length(features))), 
                 potential=feature.names)

##save the results
if(save_model==T){
  write_csv (data, paste0("output/xgb", v, ".0_predicted_bgbiomass_seagrant_plots_landsat8.csv"))
  write_csv (feat, paste0("output/xgb", v, ".0_features.csv"))
  write_csv (feats.all.import, paste0("output/xgb", v, ".0_all_import_features.csv"))
  write_csv (rmse, paste0("output/xgb", v, ".0_model_rmse_seagrant_plots.csv"))
  
  saveRDS(xgmod1, file = paste0("output/xgboost", v, ".0_bgbiomass_1.rda"))
  saveRDS(xgmod2, file = paste0("output/xgboost", v, ".0_bgbiomass_2.rda"))
  saveRDS(xgmod3, file = paste0("output/xgboost", v, ".0_bgbiomass_3.rda"))
  saveRDS(xgmod4, file = paste0("output/xgboost", v, ".0_bgbiomass_4.rda"))
  saveRDS(xgmod5, file = paste0("output/xgboost", v, ".0_bgbiomass_5.rda"))
  xgb.save(xgmod1,  paste0("output/xgboost", v, ".0_bgbiomass_1"))
  xgb.save(xgmod2,  paste0("output/xgboost", v, ".0_bgbiomass_2"))
  xgb.save(xgmod3,  paste0("output/xgboost", v, ".0_bgbiomass_3"))
  xgb.save(xgmod4,  paste0("output/xgboost", v, ".0_bgbiomass_4"))
  xgb.save(xgmod5,  paste0("output/xgboost", v, ".0_bgbiomass_5"))
  
}

data <- read_csv(paste0("output/xgb", v, ".0_predicted_bgbiomass_seagrant_plots_landsat8.csv"))
bgb_goodness_of_fit <- data %>%
  dplyr::select(c(14, 255:266))
bgb_goodness_of_fit <- bgb_goodness_of_fit[complete.cases(bgb_goodness_of_fit),]

sspe.1 <- sum(bgb_goodness_of_fit$bgm2.core.stemscaled - bgb_goodness_of_fit$pred.1)^2
sspe.2 <- sum(bgb_goodness_of_fit$bgm2.core.stemscaled - bgb_goodness_of_fit$pred.2)^2
sspe.3 <- sum(bgb_goodness_of_fit$bgm2.core.stemscaled - bgb_goodness_of_fit$pred.3)^2
sspe.4 <- sum(bgb_goodness_of_fit$bgm2.core.stemscaled - bgb_goodness_of_fit$pred.4)^2
sspe.5 <- sum(bgb_goodness_of_fit$bgm2.core.stemscaled - bgb_goodness_of_fit$pred.5)^2

u_bias.1 <- (length(bgb_goodness_of_fit) * (mean(bgb_goodness_of_fit$bgm2.core.stemscaled) - mean(bgb_goodness_of_fit$pred.1))^2) / sspe.1
u_bias.2 <- (length(bgb_goodness_of_fit) * (mean(bgb_goodness_of_fit$bgm2.core.stemscaled) - mean(bgb_goodness_of_fit$pred.2))^2) / sspe.2
u_bias.3 <- (length(bgb_goodness_of_fit) * (mean(bgb_goodness_of_fit$bgm2.core.stemscaled) - mean(bgb_goodness_of_fit$pred.3))^2) / sspe.3
u_bias.4 <- (length(bgb_goodness_of_fit) * (mean(bgb_goodness_of_fit$bgm2.core.stemscaled) - mean(bgb_goodness_of_fit$pred.4))^2) / sspe.4
u_bias.5 <- (length(bgb_goodness_of_fit) * (mean(bgb_goodness_of_fit$bgm2.core.stemscaled) - mean(bgb_goodness_of_fit$pred.5))^2) / sspe.5

bg_1_lm <- lm(pred.1 ~ bgm2.core.stemscaled, data = bgb_goodness_of_fit)
bg_2_lm <- lm(pred.2 ~ bgm2.core.stemscaled, data = bgb_goodness_of_fit)
bg_3_lm <- lm(pred.3 ~ bgm2.core.stemscaled, data = bgb_goodness_of_fit)
bg_4_lm <- lm(pred.4 ~ bgm2.core.stemscaled, data = bgb_goodness_of_fit)
bg_5_lm <- lm(pred.5 ~ bgm2.core.stemscaled, data = bgb_goodness_of_fit)

p1 <- predict(bg_1_lm, bgb_goodness_of_fit)
p2 <- predict(bg_2_lm, bgb_goodness_of_fit)
p3 <- predict(bg_3_lm, bgb_goodness_of_fit)
p4 <- predict(bg_4_lm, bgb_goodness_of_fit)
p5 <- predict(bg_5_lm, bgb_goodness_of_fit)

u_b1_1.1 <- ((summary(bg_1_lm)$coefficients[2] - 1)^2 * sum(bgb_goodness_of_fit$pred.1 - mean(bgb_goodness_of_fit$pred.1))^2) / sspe.1
u_b1_1.2 <- ((summary(bg_2_lm)$coefficients[2] - 1)^2 * sum(bgb_goodness_of_fit$pred.2 - mean(bgb_goodness_of_fit$pred.1))^2) / sspe.2
u_b1_1.3 <- ((summary(bg_3_lm)$coefficients[2] - 1)^2 * sum(bgb_goodness_of_fit$pred.3 - mean(bgb_goodness_of_fit$pred.1))^2) / sspe.3
u_b1_1.4 <- ((summary(bg_4_lm)$coefficients[2] - 1)^2 * sum(bgb_goodness_of_fit$pred.4 - mean(bgb_goodness_of_fit$pred.1))^2) / sspe.4
u_b1_1.5 <- ((summary(bg_5_lm)$coefficients[2] - 1)^2 * sum(bgb_goodness_of_fit$pred.5 - mean(bgb_goodness_of_fit$pred.1))^2) / sspe.5

u_e.1 <- sum(p1 - bgb_goodness_of_fit$bgm2.core.stemscaled)^2 / sspe.1
u_e.2 <- sum(p2 - bgb_goodness_of_fit$bgm2.core.stemscaled)^2 / sspe.2
u_e.3 <- sum(p3 - bgb_goodness_of_fit$bgm2.core.stemscaled)^2 / sspe.3
u_e.4 <- sum(p4 - bgb_goodness_of_fit$bgm2.core.stemscaled)^2 / sspe.4
u_e.5 <- sum(p5 - bgb_goodness_of_fit$bgm2.core.stemscaled)^2 / sspe.5





ls <- list.files(path = "/home/kyle/Documents/Documents/UT/Sapelo/berm_output", pattern = "landsat8_hladik_split")
Sys.time()
for(i in seq(1,length(ls))){
  
  #reading in splits of dataset
  cut_no <- i


v <- 2
outer<-5

##read in results previously saved, if needed
features<-read_csv(paste0("output/xgb", v, ".0_features.csv"))
features<-features$all[!is.na(features$all)]
xgmod1<-xgb.load(paste0("output/xgboost", v, ".0_bgbiomass_1"))
xgmod2<-xgb.load(paste0("output/xgboost", v, ".0_bgbiomass_2"))
xgmod3<-xgb.load(paste0("output/xgboost", v, ".0_bgbiomass_3"))
xgmod4<-xgb.load(paste0("output/xgboost", v, ".0_bgbiomass_4"))
xgmod5<-xgb.load(paste0("output/xgboost", v, ".0_bgbiomass_5"))

##Now we can use the models to predict novel data, typically at the site level
##make sure the novel data have run through the process_new_pixels_for_xgb.r script
newpixels<-read.csv(paste0("/home/kyle/Documents/Documents/UT/Sapelo/berm_output/xgb_processed_landsat8_hladik_split", cut_no, ".csv"), stringsAsFactors = F)
newpixels$date<-as.Date(newpixels$date)

newout<-matrix(data=NA, nrow=nrow(newpixels), ncol=outer)
newout[,1]<-(as.integer(round(predict(xgmod1, data.matrix(newpixels[,features])))))
newout[,2]<-(as.integer(round(predict(xgmod2, data.matrix(newpixels[,features])))))
newout[,3]<-(as.integer(round(predict(xgmod3, data.matrix(newpixels[,features])))))
newout[,4]<-(as.integer(round(predict(xgmod4, data.matrix(newpixels[,features])))))
newout[,5]<-(as.integer(round(predict(xgmod5, data.matrix(newpixels[,features])))))

##save the predicted belowground biomass as the mean of the predictions from all 5 models
newpixels$predbg.xgb<-rowMeans(newout)
##set a floor for the predictions, though probably not even needed with XGB
newpixels$predbg.xgb<-ifelse(newpixels$predbg.xgb<0,50,newpixels$predbg.xgb)
newpixels$predag.allom.l8 <-ifelse(newpixels$predag.allom.l8<0,50,newpixels$predag.allom.l8)

##newpixels is too huge to keep all the variables, lets parse it down just to interesting ones for plotting results
newpixels<-dplyr::select(newpixels, pix,  year,  utm_east, utm_north, greenup, date, obs,  predn.l8,predchl.l8,
                         predlai.l8, predag.allom.l8, lst, ndvi, mo, greendoy,doy, growingday, photo, growag,  growN, 
                         diff_N_growing, diff_lai_growing, diff_chl_growing, diff_ag_growing,
                         growchl, growlai, growphoto, utm_east,utm_north,  elevation, water,hiwater, maxwater,  
                         lowater, local_hiwater, local_lowater,flood_time, dayl_mean,  prcp_mean,  srad_mean, par_sum,par_tot,
                         tmax_mean,  tmin_mean,  vp_mean, dayl_sum, prcp_sum,srad_sum,  tmax_sum,tmin_sum,vp_sum, 
                         predbg.xgb , flood_time)

##save the results
if(save_model==T){
  write_csv(newpixels, paste0("output/xgb", v, ".0_predicted_bgbiomass_landsat8_hladik_split", cut_no, ".csv"))
}
}
