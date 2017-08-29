library(sp)
library(raster)
library(dplyr)

#######################################################
#Create a continuous grid of a certain cell size across the extent of the template raster
#######################################################
templateRaster=raster::raster('~/data/yearly_bioclim/bio1_1971.tif')
create_grid=function(cellsize, templateRaster){
  minX=bbox(templateRaster)[1,1]
  maxX=bbox(templateRaster)[1,2]
  minY=bbox(templateRaster)[2,1]
  maxY=bbox(templateRaster)[2,2]
  
  lowerX=maxX-minX
  lowerY=maxY-minY
  
  numCellsX=round((lowerX+1)/cellsize)
  numCellsY=round((lowerY+1)/cellsize)
  g=GridTopology(c(minX, minY), c(cellsize,cellsize), c(numCellsX, numCellsY))
  g=SpatialGrid(g, proj4string = CRS(crs(templateRaster)@projargs))
}


#######################################################
#Takes a grid cell and site spatial df.
#Returns a df of colnames('cellID','siteID')
#######################################################
assign_sites_to_grid=function(g, sites){
  g=as(g, 'SpatialPolygons')
  if(!identicalCRS(g, sites)){
    sites = spTransform(sites, CRSobj = crs(g))
  }
  x=data.frame(spatial_cell_id=over(sites, g), siteID=sites@data$siteID) %>%
    filter(!is.na(spatial_cell_id))
  return(x)
}

###################################################################
#Extract bioclim data for all routes.
#Used for initial model fit
####################################################################
get_bioclim_data=function(){
  #Create a regex expression to use only bioclim data from the years being trained and tested on. 
  #Saves on processing time
  relevant_years = paste0(c(training_years, testing_years), '|', collapse='')
  relevant_years = substr(relevant_years, 1, nchar(relevant_years)-1)
  bioclim_files = list.files(bioclim_data_folder, pattern=relevant_years, full.names = TRUE)
  
  #Filter out any xml files that sneak in there. 
  bioclim_files =  bioclim_files[!grepl('xml', bioclim_files)]
  
  bioclim_stack=raster::stack(bioclim_files)
  
  route_data = as.data.frame(raster::extract(bioclim_stack, routes_spatial))
  
  route_data$siteID = routes_spatial$siteID
  
  #Take raster names (bio12_1995, bio13_1995, etc) and convert year to a column and
  #bioX values to individual columns
  route_data = route_data %>%
    tidyr::gather(var_year, value, -siteID) %>%
    tidyr::separate(var_year, c('var','year'), '_', convert = TRUE) %>%
    tidyr::spread(var, value)
  
  return(route_data)       
}

###########################################################
#Fractions skill score from Roberts and Lean 2008
########################################################
fractions_skill_score = function(y_obs, y_pred){
  numerator   = mean((y_pred - y_obs)^2)
  denominator = sum(y_pred^2 + y_obs^2) / length(y_obs)
  return(1 - (numerator/denominator))
}

###########################################################
#Choose a probability threshold by maximizing specificity + sensitivity
###########################################################

sensitivity = function(observed, predicted){
  correctly_predicted = observed==predicted
  
  tp=sum(correctly_predicted & predicted==1)
  fn=sum((!correctly_predicted) & predicted==0)
  return( sum(tp) / (sum(tp)+sum(fn)) )
}

specificity = function(observed, predicted){
  correctly_predicted = observed==predicted
  
  tn=sum(correctly_predicted & predicted==0)
  fp=sum((!correctly_predicted) & predicted==1)
  return( sum(tn) / (sum(tn)+sum(fp)) )
}

get_threshold = function(observed, predicted){
  max_sss=0
  threshold=0
  for(t in seq(0.01, 0.99, 0.01)){
    predicted_binary = ifelse(predicted>t, 1,0)
    
    sens=sensitivity(observed, predicted_binary)
    spec=specificity(observed, predicted_binary)
    if(sens + spec > max_sss){
      max_sss=(sens+spec)
      threshold=t
    }
  }
  return(threshold)
}
