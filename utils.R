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



specificity=function(observed, predicted){
  correctly_predicted = observed==predicted
  
  tn=sum(correctly_predicted & predicted==0)
  fp=sum((!correctly_predicted) & predicted==1)
  return( sum(tn) / (sum(tn)+sum(fp)) )
}
sensitivity=function(observed, predicted){
  correctly_predicted = observed==predicted
  
  tp=sum(correctly_predicted & predicted==1)
  fn=sum((!correctly_predicted) & predicted==0)
  return( sum(tp) / (sum(tp)+sum(fn)) )
}


#The probability threshold which maximizes sensitivity+specificity
get_threshold = function(observations, probabilites){
  max_sss=0
  threshold=0
  for(proposed_threshold in seq(0.001, 1, 0.001)){
    predicted_binary = (probabilites > proposed_threshold)*1
    
    sens=sensitivity(observations, predicted_binary)
    spec=specificity(observations, predicted_binary)
    if(sens + spec > max_sss){
      max_sss=(sens+spec)
      threshold=proposed_threshold
    }
  }
  return(threshold)
}
