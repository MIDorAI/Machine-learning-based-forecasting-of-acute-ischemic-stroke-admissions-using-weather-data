
###############clear Rstudio environment variables ####################
rm(list=ls())
gc()

##################libraries#######################################
library(geosphere)
library(dplyr)
library(lubridate)
library(stringr)

#read  config
if(file.exists("config.yml")){
  conf <- config::get(file = "config.yml")
}else{
  conf <- config::get(file = "config_default.yml")
}


print("loading stroke cohort")
###################load stroke cohort #############################
df.stroke.cohort <- read.csv(file  = file.path(getwd(),
                                               "data_delivered/p21_stroke_cohort_20231003.csv")
                             , sep = ";")

colnames(df.stroke.cohort) <- c("encounter_id", "patient_id", "diagnoseart",
                                "icd_code", "admission_date", "discharge_date",
                                "discharge_reason", "patient_zip", "wohnort", 
                                "patient_sex", "verweildauer_intensiv",
                                "beatmungsstunden" , "alter_in_jahren_am_aufnahmetag")


df.stroke.cohort$admission_datetime <- df.stroke.cohort$admission_date

df.stroke.cohort$icd_family <- substr(df.stroke.cohort$icd_code, 1,3)

################## checking for missing dates in admission, discharge and recorded date ##############
print(" Handling missing date")
df.stroke.cohort$admission_date <- coalesce(df.stroke.cohort$admission_date, 
                                            df.stroke.cohort$discharge_date)
print(paste("Missing dates count after coalesce:", length(which(is.na(df.stroke.cohort$admission_date)))))

if(length(which(is.na(df.stroke.cohort$admission_date))) > 0){
  df.stroke.cohort <- df.stroke.cohort[-c(which(is.na(df.stroke.cohort$admission_date))),]
}

print(paste("Stroke cohort rows after removing missing dates:", nrow(df.stroke.cohort)))
df.stroke.cohort$admission_date <- as.character(as.Date(df.stroke.cohort$admission_date))


###################extract latitude and longitude ###################
print("extracting lat and long based on PLZ")
df.plz <- read.csv(file  = file.path(getwd(), "data/plz.csv"))
df.plz$plz <- str_pad(df.plz$plz, 5, pad = "0")
df.stroke.cohort$patient_zip <- as.character(df.stroke.cohort$patient_zip)

#replace rows with NA as postal code with hospital location
df.stroke.cohort$patient_zip[which(is.na(df.stroke.cohort$patient_zip))]<- conf$PLZ
#replace rows that have postal code less than 5 digits with hospital location 
df.stroke.cohort$patient_zip[which(nchar(df.stroke.cohort$patient_zip)!=5)] <- conf$PLZ

df.stroke.cohort <- left_join(df.stroke.cohort,df.plz,c("patient_zip" = "plz"))
if(length(which(is.na(df.stroke.cohort$latitude) | is.na(df.stroke.cohort$longitude) ))>0){
  df.stroke.cohort <- df.stroke.cohort[-c(which(is.na(df.stroke.cohort$latitude) | is.na(df.stroke.cohort$longitude) )),]
  
}
############################stroke daily counts#################################
print("stroke daily counts - three different outputs")
count.info <- df.stroke.cohort%>%
  group_by(Date = as.Date(admission_date))%>%
  summarise(total_count = length(unique(encounter_id))
            ,count_I63 = length(unique(encounter_id[ c(which(startsWith(icd_family,"I63")))]))
            ,count_I61 = length(unique(encounter_id[ c(which(startsWith(icd_family,"I61")))]))
            ,count_I60 = length(unique(encounter_id[ c(which(startsWith(icd_family,"I60")))]))
            ,Ischemic_count = count_I63
            ,Bleeding_count = count_I61 + count_I60
            
  )
count.info.all <- data.table::data.table()
count.info.all$Date <- seq(ymd(min(count.info$Date)),ymd(max(count.info$Date)),by='days')
count.info.all <- left_join(count.info.all,count.info,"Date")
count.info.all[is.na(count.info.all)] <- 0

###################load appropriate conf files#############################
variables <- config::get(file = file.path(getwd(),"conf/WeatherVariables.yml"))
sites <- config::get(file = file.path(getwd(),"conf/sites.yml"))

###############load the weather data downloaded from DWD  and meta index dataframe#######
print("Load weather data and find the appropriate stations for patients")
load(file = file.path(getwd(),"data/weather_data_all_updated.RData"))
weather.data$MESS_DATUM =as.character(weather.data$MESS_DATUM)
load(file = file.path(getwd(),"data/metaIndex.RData"))

#################calculate the distance bin and angle based on clinic location and weather station ######
metaIndex$distance_from_clinic <- apply(X = data.frame(1:nrow(metaIndex)),MARGIN = 1,FUN = function(X){(distHaversine(p1 = c(sites$UMM$clinic.lon,sites$UMM$clinic.lat), p2  = c(metaIndex$lon[X],metaIndex$lat[X]))/1000)})
metaIndex$distance_bin<- cut(metaIndex$distance_from_clinic, breaks=c(0,10,  50, 100, 150, 200, 300,400,500,max(metaIndex$distance_from_clinic)))

metaIndex$angle <- apply(X = data.frame(1:nrow(metaIndex)),MARGIN = 1,FUN = function(X){bearing(p1 = c(sites$UMM$clinic.lat,sites$UMM$clinic.lon),p2 = c(metaIndex$lat[X],metaIndex$lon[X]))})
metaIndex$angle<-(metaIndex$angle + 360) %% 360

arr_4=c("N","NE","E","SE","S","SW","W","NW")
metaIndex$direction <-arr_4[ceiling(metaIndex$angle/45)]

#########################################station which contains all the variables##################################
weather.data <- subset(weather.data,select = -c(snow_pack_height , fresh_snow_depth,precipitation_form,
  precipitation_form,snow_depth,min_precipitation_height ))

#weather.data <- na.locf(weather.data)
station.agg.all.vars <- metaIndex%>%
  group_by(Stations_id)%>%
  mutate(unique_variable_count = length(unique(res_var)))%>%
  filter(unique_variable_count == length(variables))%>%
  dplyr::select(Stations_id,lat,lon)%>%
  distinct()

########find closest station from the list of station that has all the variables
df.stroke.cohort$closest.from.plz <- station.agg.all.vars$Stations_id[as.data.frame(RANN::nn2(station.agg.all.vars[,c("lat","lon")] ,df.stroke.cohort[,c("latitude","longitude")], k=1))$nn.idx]
df.stroke.cohort$distance <- apply(X = data.frame(1:nrow(df.stroke.cohort)),MARGIN = 1,FUN = function(X){(distHaversine(p1 = c(metaIndex[which(metaIndex$Stations_id ==df.stroke.cohort$closest.from.plz[X])[1],"lon"]
                                                                                                                                       ,metaIndex[which(metaIndex$Stations_id ==df.stroke.cohort$closest.from.plz[X])[1],"lat"])                                                                                                                             , p2  = c(df.stroke.cohort$longitude[X],df.stroke.cohort$latitude[X]))/1000)})
##############################fill for zero dates based on clinic location################################

print("fill for zero dates based on clinic location")
type_1 <- df.stroke.cohort[,c("admission_date","closest.from.plz")]
closest.to.site <- station.agg.all.vars$Stations_id[as.data.frame(RANN::nn2(station.agg.all.vars[,c("lat","lon")] ,data.frame(sites$UMM$clinic.lat,sites$UMM$clinic.lon), k=1))$nn.idx]
zero_dates <- count.info.all$Date[c(which(count.info.all$total_count == 0))]
type_1_zero <- data.table::data.table()
type_1_zero$admission_date <- zero_dates
type_1_zero$admission_date <- as.character(type_1_zero$admission_date)
type_1_zero$closest.from.plz <- closest.to.site

type_1 <- rbind(type_1,as.data.frame(type_1_zero))
type_1 <- type_1[c(order(type_1$admission_date)),]

####################################adding weather#############################################
  print("adding weather and lagged values")
  #####################add current days weather based on the closest station############
  type_1 <- left_join(type_1,weather.data,by = c("admission_date" = "MESS_DATUM" , "closest.from.plz" ="STATIONS_ID"))
  #####################add PT temperature based on site location - station id from sites.yml############
  # read PT temp weather data
  PT.data<-  read.csv(file  = file.path(getwd(),"data/PT_Dwd_Agg.csv"))
  type_1$PT.station <- sites[[conf$site]]$PT.Station.id
  type_1 <- left_join(type_1,PT.data,by = c("admission_date" = "MESS_DATUM" , "PT.station" ="STATIONS_ID"))
  type_1[sapply(type_1, is.infinite)] <- NA
 ################## add lagged weather data based on clinic site ################
  lag_names  <- paste("lag_",seq(1,7),sep = "")
  old_name <- names(weather.data)
  for(i in 1:length(lag_names)){
    type_1$closest.to.site <- closest.to.site
    type_1$admission_date_lag <- as.character(as.Date(type_1$admission_date) - i)
    names(weather.data) <- c(names(weather.data)[1:2],paste("lag_",i,old_name[3:ncol(weather.data)],sep = ""))
    type_1 <- left_join(type_1,weather.data,by = c("admission_date_lag" = "MESS_DATUM" , "closest.to.site" ="STATIONS_ID"))
    #print(i)
  }
  names(weather.data) <- old_name
  
  
#####################Transform patient level to different resolution#################################
 type_1 <- subset(type_1,select = -c(closest.from.plz,PT.station,closest.to.site,admission_date_lag))


##1.###################################Daily################################################ 
  print("Creating daily level data")
  
  daily_level <- type_1%>%
    group_by(admission_date)%>%
    summarise_all(funs(mean), na.rm = TRUE)
  
  daily_level$day_of_month <- day(daily_level$admission_date)
  daily_level$day_of_year <- yday(daily_level$admission_date)
  daily_level$month <- month(daily_level$admission_date)
  daily_level$wday <- as.factor(wday(daily_level$admission_date,week_start =1)) 
  
  
  daily_level$year <- year(daily_level$admission_date)
  daily_level$week_num <-as.integer(format(as.Date(daily_level$admission_date),"%W"))
  
  
  
  week.agg <- count.info.all%>%
                    group_by(year = year(count.info.all$Date)
                            # ,week_num = week(count.info.all$Date)
                             ,week_num = as.integer(format(as.Date(count.info.all$Date),"%W"))
                             )%>%
                    summarize(mean_prior_week_total = mean(total_count,na.rm=TRUE)
                           ,median_prior_week_total = median(total_count,na.rm=TRUE)
                           ,mean_prior_week_ischemic = mean(Ischemic_count,na.rm=TRUE)
                           ,median_prior_week_ischemic = median(Ischemic_count,na.rm=TRUE)
                           ,mean_prior_week_bleeding = mean(Bleeding_count,na.rm=TRUE)
                           ,median_prior_week_bleeding = median(Bleeding_count,na.rm=TRUE)
                           )
  
  
  #join to get mean prior week and median prior week variables
  week.agg$week_num<- week.agg$week_num + 1
  daily_level <- left_join(daily_level,week.agg,c("year" = "year", "week_num"="week_num"))
  week.agg$week_num<- week.agg$week_num - 1
  
  #join to get output variable
  count.info.all$Date <- as.character(count.info.all$Date)
  daily_level <- right_join(daily_level,count.info.all,c("admission_date" = "Date"))
  
  #drop week number which is equal to minimum and min year
  # daily_level<- daily_level[-c(which(daily_level$week_num == min(daily_level$week_num) & 
  #                                      daily_level$year == min(daily_level$year))),] 
  
  
  
  write.csv(daily_level,file = file.path(getwd(),"data_delivered/daily_level_p21.csv"),row.names = FALSE)

