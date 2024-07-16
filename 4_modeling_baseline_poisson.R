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

######################################## Modelling #################################
wb <-openxlsx:::createWorkbook()
openxlsx:::addWorksheet(wb, "Daily")

#1.################################### Daily ########################################
	print("Models for daily count")
	daily <- na.locf(daily,na.rm = FALSE)
	daily <- na.locf(daily,fromLast = TRUE)
	####features
	# weather features standardize
	print("Standardize the input features")
	features_daily_weather <- subset(daily, select = -c(admission_date, total_count, count_I63,count_I61,
	                                                    count_I60, Ischemic_count, Bleeding_count, 
	                                                    mean_prior_week_ischemic, median_prior_week_ischemic,
	                                                    mean_prior_week_bleeding, median_prior_week_bleeding,
	                                                    mean_prior_week_total, median_prior_week_total,
	                                                    day_of_month,day_of_year, month, wday, year, week_num))
	
	#features_daily_weather <- scale(features_daily_weather) # default: center = TRUE, scale = TRUE
	
	# calendar and mean prior week feature 
	features_daily_calendar_total_count <- subset(daily, select = c(day_of_month,day_of_year, month, wday, year, week_num, mean_prior_week_total, median_prior_week_total))
	features_daily_calendar_ischemic <- subset(daily, select = c(day_of_month, day_of_year, month, wday, year, week_num, mean_prior_week_ischemic, median_prior_week_ischemic))
	features_daily_calendar_bleeding <- subset(daily, select = c(day_of_month, day_of_year, month, wday, year, week_num, mean_prior_week_bleeding, median_prior_week_bleeding))
	
	# bind weather and calendar features for all three outputs
	features_daily_total <- cbind(features_daily_weather, features_daily_calendar_total_count)
	features_daily_ischemic <- cbind(features_daily_weather, features_daily_calendar_ischemic)
	features_daily_bleeding <- cbind(features_daily_weather, features_daily_calendar_bleeding)
	
	rm(features_daily_weather,features_daily_calendar_total_count,features_daily_calendar_ischemic,features_daily_calendar_bleeding)
	

	
	#1-a.################################################# Models for Ischemic count################################### 
	print("split to train and test set")
	##split train and test
	test_index <- which(year(daily$admission_date) ==  max(year(daily$admission_date)))
	train_features_daily <- features_daily_ischemic[-test_index,]
	test_features_daily <- features_daily_ischemic[test_index,]
	

	 train_ischemic_count_daily <- Ischemic_count_daily[-test_index]
	 test_ischemic_count_daily <- Ischemic_count_daily[test_index]
	 train_daily <- as.data.frame(cbind(train_features_daily,train_ischemic_count_daily))

	print("Poisson")
	### Poisson
	  #fit the model
	try({
	  set.seed(1492)
	  poisson_daily_ischmeic_count <- glm(train_ischemic_count_daily ~ ., data = train_daily, family = poisson(link = "log"))
	  #print summary
	  summary(poisson_daily_ischmeic_count)
	  #predict
	  preds <- predict(poisson_daily_ischmeic_count, newdata = as.data.frame(test_features_daily), type = "response")
	  #metrics
	  rmse <- paste("RMSE of Poisson daily model for ischmeic count",RMSE(pred = preds,obs = test_ischemic_count_daily))
	  mae <- paste("MAE of Poisson daily model for ischmeic count",MAE(pred = preds,obs = test_ischemic_count_daily))
	  openxlsx:::writeData(wb = wb,x = rmse,sheet = "Daily", withFilter = FALSE,startRow = 11)
	  openxlsx:::writeData(wb = wb,x = mae,sheet = "Daily", withFilter = FALSE,startRow = 12)
	  ##save model
	  saveRDS(object = poisson_daily_ischmeic_count,file = paste("./results/poisson_daily_ischemic_count_",site.name,".rda",sep = ""))
	  rm(poisson_daily_ischmeic_count, preds, rmse, mae)
	}, silent=TRUE)
	print("Random Forest")
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
print("Baseline poisson modelling done")
