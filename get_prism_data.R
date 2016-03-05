# Download and install PRISM data into PostgreSQL and values for BBS routes
#Taken from Ethan White. https://github.com/weecology/bbs-forecasting

library(prism)
library(raster)
library(dplyr)
library(tidyr)
dataFolder='~/data/'
options(prism.path = paste(dataFolder,'prism', sep=''))
routes=read.csv(paste(dataFolder,'bbs/BBS_routes.csv', sep=''))

downloadPrismData=function(){
  months <- c(1:12)
  clim_vars <- c("ppt", "tmin", "tmean", "tmax")
  for (month in months){
    for (clim_var in clim_vars){
      get_prism_monthlys(type=clim_var, year = 1966:1981, mon = month, keepZip=F)
    }
  }
}

#downloadPrismData()

locations = routes %>%
  mutate(siteID=paste(countrynum, statenum, route,sep='-')) %>%
  dplyr::select(siteID, loni, lati) %>%
  rename(long=loni, lat=lati)

coordinates(locations) <- c("long", "lat")
#prism_files <- grep("_[0-9]{4}[0]",ls_prism_data()[,1],value=T)
prism_stacked <- prism_stack(ls_prism_data())
extracted <- raster::extract(prism_stacked, locations)
prism_bbs_data <- data.frame(siteID = locations$siteID, coordinates(locations), extracted)
prism_bbs_data <- prism_bbs_data %>%
                  gather(date, value, 4:ncol(prism_bbs_data)) %>%
                  tidyr::extract(date, c("clim_var", "year", "month"),
                                 "PRISM_([:alpha:]*)_stable_[:alnum:]*_([:digit:]{4})([:digit:]{2})_") %>%
  spread(clim_var, value)
prism_bbs_data$year <- as.numeric(prism_bbs_data$year)
prism_bbs_data$month <- as.numeric(prism_bbs_data$month)
write.csv(prism_bbs_data, './prism_bbs_data.csv', row.names = FALSE)
rm(extracted, locations, prism_stacked)
#database <- src_sqlite("./data/bbsforecasting.sqlite", create = TRUE)
#mydata <- copy_to(database, prism_bbs_data, temporary = FALSE,
#                  indexes = list(c("site_id", "year", "month")))
