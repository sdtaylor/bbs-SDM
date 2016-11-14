library(prism)
library(raster)
library(magrittr)
library(dplyr)
library(sp)
source('config.R')
options(prism.path = "~/data/prismdata")

#TODO: put in download functionality to make this stand alone

########################################################
#Download prism data if not already there
###
download_prism = function(){
  
}

###################################################################
#Helper function to calculate some of the bioclim variables,
#like "precip in coldest month"
###################################################################
max_min_combo=function(vec1,vec2,max=TRUE){
  #Return the value in vec1 in the position where
  #vec2 is either highest or lowest. But 1st check for na 
  #values. 
  if(any(is.na(vec1)) | any(is.na(vec2))){
    return(NA)
  } else  if(max){
    return(vec1[which.max(vec2)])
  } else {
    return(vec1[which.min(vec2)])
  }
}


########################################################
#Create raster stack of 19 bioclim vars from monthly prism data
#Function used in raster::overlay() which accepts
make_bioclim_from_prism = function(prism_year_stacked, year){
  
  df=as.data.frame(prism_year_stacked) 
  df$cell_id = 1:nrow(df)
  df = df %>%
    tidyr::gather(layer_name, value, -cell_id) %>%
    tidyr::extract(layer_name, c("clim_var", "year", "month"),
                   "PRISM_([:alpha:]*)_stable_[:alnum:]*_([:digit:]{4})([:digit:]{2})_",
                   convert=TRUE) %>%
    tidyr::spread(clim_var, value)

  #Process the quarter ones first.
  quarter_info=data.frame(month=1:12, quarter=c(1,1,1,2,2,2,3,3,3,4,4,4))
  bioclim_quarter_data= df %>%
    dplyr::left_join(quarter_info, by='month') %>%
    dplyr::group_by(cell_id, quarter) %>%
    dplyr::summarize(precip=sum(ppt), temp=mean(tmean)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarize(bio8=max_min_combo(temp, precip, max=TRUE),
              bio9=max_min_combo(temp, precip, max=FALSE),
              bio10=max(temp),
              bio11=min(temp),
              bio16=max(precip),
              bio17=min(precip),
              bio18=max_min_combo(precip, temp, max=TRUE),
              bio19=max_min_combo(precip, temp, max=FALSE)) %>%
    dplyr::ungroup()  
 
  #Next the yearly ones, joining the quartely ones  back in at the end. 
  bioclim_data=df %>%
    dplyr::group_by(cell_id) %>%
    dplyr::mutate(monthly_temp_diff=tmax-tmin) %>%
    dplyr::summarize(bio1=mean(tmean),
              bio2=mean(monthly_temp_diff),
              bio4=sd(tmean)*100,
              bio5=max_min_combo(tmax,tmean,max=TRUE),
              bio6=max_min_combo(tmin,tmean,max=FALSE),
              bio12=sum(ppt),
              bio13=max(ppt),
              bio14=min(ppt),
              bio15=cv(ppt)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(bio7=bio5-bio6,
                  bio3=(bio2/bio7)*100) %>%
    dplyr::full_join(bioclim_quarter_data, by=c('cell_id'))
  
     
  
  template_raster=raster::raster(nrows = prism_year_stacked@nrows, ncols = prism_year_stacked@ncols, 
                            crs = prism_year_stacked@crs, ext = prism_year_stacked@extent)
  
  bioclim_names = paste('bio',1:19, sep = '')
  for(this_bioclim_var in bioclim_names){
    raster_obj = raster::setValues(template_raster, magrittr::extract2(bioclim_data, this_bioclim_var))
    raster_filename=paste0(bioclim_data_folder,this_bioclim_var,'_',year,'.tif')
    writeRaster(raster_obj, raster_filename, formate='GTiff')
  }

}



########################################################
#Convert all availableyears of prism data to bioclim vars
###
convert_prism_to_raster=function(){
  for(this_year in timeRange){
    this_year_files = ls_prism_data() %>% 
      dplyr::filter(grepl(as.character(this_year),files))

    this_year_stack = prism_stack(this_year_files)
  
    make_bioclim_from_prism(this_year_stack, this_year)
      
    
  }
}

#######################################################
#Apply model to bioclim verification data and return result as a raster stack
#######################################################

#Make a list of RasterStacks of all the testing data. One rasterstack is all the bioclim vars
#for that year. Labels reflect the name so it can be used inside raster::predict()
bioclim_test_data=list()
for(test_year in testing_years){
  bioclim_test_files = list.files(bioclim_data_folder, pattern=as.character(test_year), full.names = TRUE)
  bioclim_test_files =  bioclim_test_files[!grepl('xml', bioclim_test_files)]
    
  test_year_stack=raster::stack(bioclim_test_files)
  
  for(layer_num in 1:length(test_year_stack@layers)){
    layer_name = test_year_stack[[layer_num]]@data@names
    layer_name = stringr::word(layer_name, 1, 1, sep='_')
    test_year_stack[[layer_num]]@data@names = layer_name
  }
  
  bioclim_test_data = append(bioclim_test_data, test_year_stack)
}

#Take a distribution model output and make probability maps
#of all testing years, returned in a raster stack. 
apply_model_to_bioclim = function(m){
  predictions_stack = raster::stack()
  for(i in seq_along(testing_years)){
    predictions_stack=stack(predictions_stack, predict(bioclim_test_data[[i]], m, type='response'))
    predictions_stack[[i]]@data@names=as.character(testing_years[i])
  }
  
  return(predictions_stack)
}

###################################################################
#Extract bioclim data for all routes in the training years.
#Used for initial model fit
####################################################################
get_bioclim_data=function(){
  routes=read.csv(paste(dataFolder, 'BBS_routes.csv', sep='')) %>%
    dplyr::mutate(siteID=paste(countrynum, statenum, route,sep='-')) %>%
    dplyr::select(siteID, long=loni, lat=lati)
  
  routes = SpatialPointsDataFrame(cbind(routes$long, routes$lat), data=routes, 
                                  proj4string = CRS('+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0'))
  
  #Create a regex expression to use only bioclim data from the years being trained and tested on. 
  #Saves on processing time
  relevant_years = paste0(c(training_years, testing_years), '|', collapse='')
  relevant_years = substr(relevant_years, 1, nchar(relevant_years)-1)
  bioclim_files = list.files(bioclim_data_folder, pattern=relevant_years, full.names = TRUE)
  
  #Filter out any xml files that sneak in there. 
  bioclim_files =  bioclim_files[!grepl('xml', bioclim_files)]
  
  bioclim_stack=raster::stack(bioclim_files)
  
  route_data = as.data.frame(raster::extract(bioclim_stack, routes))
  
  route_data$siteID = routes$siteID
  
  #Take raster names (bio12_1995, bio13_1995, etc) and convert year to a column and
  #bioX values to individual columns
  route_data = route_data %>%
    tidyr::gather(var_year, value, -siteID) %>%
    tidyr::separate(var_year, c('var','year'), '_', convert = TRUE) %>%
    tidyr::spread(var, value)
            
  return(route_data)       
}
