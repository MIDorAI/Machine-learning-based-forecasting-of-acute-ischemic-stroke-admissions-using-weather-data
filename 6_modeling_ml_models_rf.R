###############clear Rstudio environment variables ####################
rm(list=ls())
gc()
#.rs.restartR()

start.time <- Sys.time()
##################libraries#######################################
library(dplyr)
library(lubridate)
library(caret)
library(xgboost)
library(zoo)
library(ranger)  

#read  config
if(file.exists("config.yml")) {
  conf <- config::get(file = "config.yml")
} else {
  conf <- config::get(file = "config_default.yml")
}


site.name <- conf$site
######################## Import prepared data  #######
daily <- read.csv(file = file.path(getwd(), "data/daily_level.csv"))


###################### Source customRF function ##################################

source(file.path(getwd(), "customRF.R"))

######################################## Modelling #################################
wb <-openxlsx:::createWorkbook()
openxlsx:::addWorksheet(wb, "Daily")


	#1-b.################################################# Models for Ischemic count################################### 
	print("split to train and test set")
	##split train and test
	test_index <- which(year(daily$admission_date) ==  max(year(daily$admission_date)))
	train_features_daily <- features_daily_ischemic[-test_index,]
	test_features_daily <- features_daily_ischemic[test_index,]
	

	 train_ischemic_count_daily <- Ischemic_count_daily[-test_index]
	 test_ischemic_count_daily <- Ischemic_count_daily[test_index]
	 train_daily <- as.data.frame(cbind(train_features_daily,train_ischemic_count_daily))

	
	### RF
	  # hyper parameter
	try({
	 set.seed(1492)
	 tunegrid <- expand.grid(.mtry = seq(17, 20, 1),
	                          .ntree = seq(1000, 2000, 500), 
	                          .nodesize = 5)
	 
	  #fit model
	  forest <- train(train_ischemic_count_daily~.,
	                  data = train_daily, 
	                  method = customRF, 
	                  trControl = timecontrol_cv, 
					          preProcess = c("center","scale"),
	                  metric = 'RMSE', 
	                  tuneGrid = tunegrid)
	  
	  # final grid 
	    final_grid <- expand.grid(mtry = forest$bestTune$mtry)
	 
	 # final rf model with chosen hyper parameter
	  X_train = train_features_daily 
	  y_train = train_ischemic_count_daily
	  X_test = test_features_daily 
	  y_test = test_ischemic_count_daily
	  print("fitting forest based on chosen hyperparameter")
	  forest_daily_ischmeic_count = train(X_train, y_train,
	                                   method = 'rf',
	                                   metric = 'RMSE',
	                                   preProcess = c("center","scale"),
	                                   tuneGrid = final_grid,
	                                   ntree = forest$bestTune$ntree,
	                                   verbosity = TRUE)
	  #predict
	  preds <- predict(forest_daily_ischmeic_count, newdata = as.data.frame(test_features_daily), type = "raw")
	  #metrics
	  rmse <- paste("RMSE of random forest daily model for ischmeic count", RMSE(pred = preds, obs = test_ischemic_count_daily))
	  mae <- paste("MAE of Poisson daily model for ischmeic count", MAE(pred = preds, obs = test_ischemic_count_daily))
	  openxlsx:::writeData(wb = wb, x = rmse, sheet = "Daily", withFilter = FALSE, startRow = 13)
	  openxlsx:::writeData(wb = wb, x = mae, sheet = "Daily", withFilter = FALSE, startRow = 14)
	  #save model
	  saveRDS(object = forest_daily_ischmeic_count, file = paste("./results/forest_daily_ischemic_count_", site.name, ".rda", sep = ""))
	  rm(forest_daily_ischmeic_count, forest, preds, rmse, mae)
	}, silent=TRUE)
	
	
	

end.time <- Sys.time()

print(end.time - start.time)
print("Machine learning model - RF done")
