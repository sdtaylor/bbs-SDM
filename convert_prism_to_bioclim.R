library(prism)
library(raster)
library(magrittr)

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
create_bioclim_stack = function(prism_year_stacked){
  
  #Get location of the different variables for each month
  stack_info=data.frame()
  for(layer in 1:length(prism_year_stacked@layers)){
    layer_name = prism_year_stacked[[layer]]@data@names
    layer_var  = stringr::word(layer_name, 2,2, '_')
    layer_month= stringr::word(layer_name, 5,5, '_')
    layer_month= as.integer(stringr::str_sub(layer_month, 5, 6))
    
    stack_info = stack_info %>%
      bind_rows(data.frame(layer=layer, var=layer_var, month=layer_month))
  }
   
  
  y
  
  #Assign each month to a quarter 
  stack_info = stack_info %>%
    left_join(data.frame(month=1:12, quarter=c(1,1,1,2,2,2,3,3,3,4,4,4)), by='month')
  
  bio1_layers = stack_info %>%
    filter(var == 'tmean') %>%
    extract2('layer')
  
  bio1_stack = prism_year_stacked[[bio1_layers]]
  
  bio1_raster = overlay(bio1_stack, fun=mean)
  
  
  bio2_layers = stack_info %>%
    
  
    
}



########################################################
#
###
