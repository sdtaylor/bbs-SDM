source('config.R')


####################################################################
#Load the bioclim variables from sqlite.
get_bioclim_data=function(){
  #Query sqlite database for the bioclim_bbs_data table. If it exists, return it.
  #Otherwise create it from the raw prism data.
  if('bioclim_bbs_data' %in% src_tbls(database)){
    return(collect(tbl(database, sql('SELECT * from bioclim_bbs_data')), n=Inf))
  } else { 
    print("bioclim data table not found. need to run load_scaled_bioclim.R")
  } 
}

####################################################################
#Load the data frame that tells which site is in which cell at all the 
#different spatial scales
get_spatial_grid_info=function(){
  #Query sqlite database for thetable. If it exists, return it.
  #Otherwise throw an error
  if('spatial_grid_info' %in% src_tbls(database)){
    return(collect(tbl(database, sql('SELECT * from spatial_grid_info')), n=Inf))
  } else { 
    stop("spatial_grid_info not found, need to run load_scaled_bioclim.R")
  } 
}

###########################################################
#Fractions skill score from Roberts and Lean 2008
########################################################
fractions_skill_score = function(y_obs, y_pred){
  numerator   = mean((y_pred - y_obs)^2)
  denominator = sum(y_pred^2 + y_obs^2) / length(y_obs)
  return(1 - (numerator/denominator))
}
