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
  dplyr::select(name, Aou)


#Focal species from Waitling et al. 2013
#Wild turkey, greater sage-grouse, Greater prairie chicken, Gambel's quail, Red-headed woodpecker, red-bellied woodpecker, fish crow
#carolina chickadee, carolina wren, bewick's wren, backmans sparrow, tricolored blackbird
#focal_spp=c(3100, 3090, 3050, 2950, 4060, 4090, 4900, 7360, 7180, 7190, 5750, 5000)
focal_spp=c(7360, #Carolina chickadee
            6010, #painted bunting
            3100, #wild turky
            4100 #Golden-fronted Woodpecker
)

source('config.R')

#results=tbl(database, 'modelResults')
results=read.csv('results/bbsSDMResults.csv')

#############################################

brier=function(obs, pred, reliability=FALSE, bins=c(0.2,0.5,0.7)){
  return(mean((pred-obs)^2))
}

brier_reliability=function(obs, pred){
  x=verification::verify(obs, pred, bins=T, thresholds = seq(0,1,0.1), show=FALSE)$bs.reliability
  return(x)
}

############################################
#Forecast skill over time. all species in 1 graph
temp=results %>%
  dplyr::filter(windowID>1) %>% #Don't use the training data
  dplyr::group_by(temporal_scale, cellSize, windowID, modelName, offset, Aou) %>%
  dplyr::summarize(brier=1-brier(presence, prediction)) %>%
  dplyr::ungroup() %>%
  tidyr::spread(modelName, brier) %>%
  dplyr::mutate(skill=1-(gbm/naive)) %>%
  dplyr::filter(!is.infinite(skill)) %>%
  #dplyr::mutate(skill=ifelse(is.nan(skill), 0, ifelse(skill<(-1), -1, skill))) %>%
  dplyr::group_by(temporal_scale, cellSize, windowID, Aou) %>%
  dplyr::summarize(skill=mean(skill), model_brier=mean(gbm), naive_brier=mean(naive)) %>%
  dplyr::mutate(time_lag=(windowID*temporal_scale) - (temporal_scale/2)) %>%
  dplyr::left_join(sppNames)

ggplot(temp, aes(x=time_lag, y=skill, colour=name, group=name)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0) +
  theme_bw() +
  facet_grid(cellSize~temporal_scale, labeller=label_both) 

#################################################
#Forecast skill over time. grid of temporal_scale~spatial_scale
#1 graph for each species
for(this_sp in focal_spp){
  results_this_sp=results %>%
    filter(Aou==this_sp) 
  if(nrow(results_this_sp)==0){next}
  
  temp=results_this_sp %>%
    dplyr::filter(windowID>1) %>% #Don't use the training data
    dplyr::group_by(temporal_scale, cellSize, windowID, modelName, offset) %>%
    dplyr::summarize(brier=1-brier(presence, prediction)) %>%
    dplyr::ungroup() %>%
    tidyr::spread(modelName, brier) %>%
    #mutate(skill=(gbm-naive)/(1-naive)) %>%
    #mutate(skill=1-(gbm/naive)) %>%
    #mutate(skill=ifelse(is.nan(skill), 0, ifelse(skill<(-1), -1, skill))) %>%
    dplyr::group_by(temporal_scale, cellSize, windowID) %>%
    dplyr::summarize(model_brier=mean(gbm), naive_brier=mean(naive)) %>%
    dplyr::mutate(time_lag=(windowID*temporal_scale) - (temporal_scale/2))
  
  sp_name=sppNames %>%
    dplyr::filter(Aou==this_sp) %>%
    extract2('name')
  
  
  this_sp_plot=ggplot(temp, aes(x=time_lag)) + 
    geom_point(aes(y=model_brier), color='red') + 
    geom_line(aes(y=model_brier), color='red')+
    geom_point(aes(y=naive_brier), color='black') + 
    geom_line(aes(y=naive_brier), color='black')+
    facet_grid(cellSize~temporal_scale, labeller=label_both) +
    ggtitle(sp_name)
    
  print(this_sp_plot)
  
  
}


