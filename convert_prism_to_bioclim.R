library(prism)
library(raster)
library(magrittr)

options(prism.path = "~/data/prismdata")

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
  
     
}



########################################################
#
###
