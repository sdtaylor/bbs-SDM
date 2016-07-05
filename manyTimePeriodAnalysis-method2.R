#This method models at the smallest temporal scale (1 year) and upscales by temporal averaging instead of training on the larger and larger scales.

library(dplyr)
library(tidyr)
library(doParallel)
library(magrittr)
library(stringr)
library(DBI)
library(ggplot2)
library(RColorBrewer)
library(Metrics)

sppNames=read.csv('~/data/bbs/BBS_species.csv') %>%
  rename(name=english_common_name, Aou=AOU) %>%
  select(name, Aou)


#Focal species from Waitling et al. 2013
#Wild turkey, greater sage-grouse, Greater prairie chicken, Gambel's quail, Red-headed woodpecker, red-bellied woodpecker, fish crow
#carolina chickadee, carolina wren, bewick's wren, backmans sparrow, tricolored blackbird
focal_spp=c(3100, 3090, 3050, 2950, 4060, 4090, 4900, 7360, 7180, 7190, 5750, 5000)


source('config.R')

results=tbl(database, 'modelResults')


############################################
#Build temporal scaling info

model_set_matrix=data.frame()
setID=1

testing_years=sort(timeRange[!timeRange %in% training_years])

for(this_temporal_scale in temporal_scales){
  #if(this_temporal_scale==1){next} #the base predictions are all temporal scale of 1
  
  numSets=floor((length(testing_years)/this_temporal_scale))
  
  temporal_replicate=c()
  for(i in 1:this_temporal_scale){
    temporal_replicate=c(temporal_replicate, 1:numSets)
  }
  temporal_replicate=sort(temporal_replicate)
  
  while(length(temporal_replicate)<(length(testing_years))){
    temporal_replicate=c(-1, temporal_replicate)
  }
  
  for(this_spatial_scale in spatial_cell_sizes){
    this_set_df=data.frame(temporal_scale=this_temporal_scale, spatial_scale=this_spatial_scale, Year=as.factor(testing_years), temporal_replicate=temporal_replicate, setID=setID)
    setID=setID+1
    model_set_matrix=bind_rows(model_set_matrix, this_set_df)
  }
  
}



#############################################
#Tur
apply_temporal_scale=function(df, this_temporal_scale){
  df=df %>%
    select(-setID)
  y=model_set_matrix %>%
    filter(temporal_scale==this_temporal_scale) %>%
    left_join(df, by=c('Year','spatial_scale'='cellSize')) %>%
    group_by(temporal_replicate, spatial_scale, cellID, temporal_scale) %>%
    summarize(presence=max(presence), prediction=(1-prod(1-prediction)))
  return(y)
}

#############################################


x=data.frame()
for(this_sp in focal_spp){
  results_this_sp=results %>%
    filter(Aou==this_sp) %>%
    collect() 
  if(nrow(results_this_sp)==0){next}
  
  for(this_temporal_scale in temporal_scales){
    temp=results_this_sp %>%
      apply_temporal_scale(this_temporal_scale) %>%
      group_by(spatial_scale,temporal_replicate) %>%
      summarize(auc=auc(presence, prediction), logloss=logLoss(as.integer(presence), prediction)) %>%
      ungroup()
    
    temp$Aou=this_sp
    temp$temporal_scale=this_temporal_scale
    x=x %>%
      bind_rows(temp)
  }
}


x=x %>%
  filter(temporal_replicate>0) %>%
  mutate(year_lag=(temporal_replicate * temporal_scale)-(temporal_scale/2) )

for(this_sp in focal_spp){
  this_sp_data=x %>% filter(Aou==this_sp)
  this_sp_name=sppNames %>%
    filter(Aou==this_sp) %>%
    extract2('name')
  if(nrow(this_sp_data)==0){print(paste('skipping:',this_sp)); next}
  this_plot=ggplot(this_sp_data, aes(y=auc, x=year_lag)) + 
            geom_point() + geom_line() +
            facet_grid(temporal_scale ~ spatial_scale) +
            ggtitle(paste(this_sp,this_sp_name))
  print(this_plot)

}


