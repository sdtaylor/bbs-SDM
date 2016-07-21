#This method models at the smallest temporal scale (1 year) and upscales by temporal averaging instead of training on the larger and larger scales.

library(dplyr)
library(tidyr)
library(magrittr)
library(stringr)
library(ggplot2)
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

#######################################################3
#fractional skill score graph

results=read.csv('results/bbsSDMResults.csv') %>%
  filter(!is.na(Aou))

for(this_sp in unique(results$Aou)){
  for(this_model in c('gbm','naive')){
  #results_this_sp=results %>%
  #  filter(Aou==this_sp) %>%
  #  spread(modelName, mse) %>%
  #  mutate(skill=1-(gbm/naive)) %>%
  #  group_by(temporal_scale, cellSize) %>%
  # summarize(skill=mean(skill)) %>%
  #  ungroup()
  results_this_sp =results %>%
    filter(Aou==this_sp, modelName==this_model) %>%
    group_by(temporal_scale, cellSize) %>%
    summarize(skill=mean(fss), skill_sd=sd(fss)) %>%
    ungroup()
  
  
  sp_name=sppNames %>%
    dplyr::filter(Aou==this_sp) %>%
    extract2('name')
  
  this_sp_plot=ggplot(results_this_sp, aes(as.factor(cellSize), as.factor(temporal_scale), fill=skill, label=paste(round(skill, 2), round(skill_sd, 2), sep='-'))) +
    geom_raster() +
    scale_fill_gradient(low='grey100', high='grey40') + 
    geom_text() +
    ggtitle(paste(sp_name, this_model, sep=' - '))
    
  print(this_sp_plot)
  }
}


#######################################################3
#fractional skill score graph with time lags

focal_spp=c(7350, #Black capped chickadee. no shift
            3100, #wild turkey. northern shift (only in La Sorte & Thompson 2007 )
            6100, #Summer Tanager. northern shift
            7190, #bewicks wren. southern shift 
            5170 #Purple finch. no shift
)


sppNames=read.csv('~/data/bbs/BBS_species.csv') %>%
  rename(name=english_common_name, Aou=AOU) %>% 
  #filter(Aou %in% focal_spp) %>%
  dplyr::select(name, Aou)

results=read.csv('results/bbsSDMResults_with_timelags_all.csv') %>%
  filter(!is.na(Aou))


#A few spp in one graph

results_summarized =results %>%
  filter(modelName=='gbm', Aou %in% focal_spp) %>%
  group_by(temporal_scale, cellSize, windowID, modelName, Aou) %>%
  summarize(skill=mean(fss), skill_sd=sd(fss)) %>%
  ungroup() %>%
  mutate(time_lag=(temporal_scale * windowID) - (temporal_scale/2)) %>%
  left_join(sppNames, by='Aou')

ggplot(results_summarized, aes(x=time_lag, y=skill, color=name, group=name)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8) +
  theme_bw() +
  facet_grid(cellSize~temporal_scale, labeller=label_both)


##Iterate over each spp

for(this_sp in unique(results$Aou)[50:80]){
  #for(this_model in c('gbm','naive')){
    #results_this_sp=results %>%
    #  filter(Aou==this_sp) %>%
    #  spread(modelName, mse) %>%
    #  mutate(skill=1-(gbm/naive)) %>%
    #  group_by(temporal_scale, cellSize) %>%
    # summarize(skill=mean(skill)) %>%
    #  ungroup()
    results_this_sp =results %>%
      filter(Aou==this_sp,temporal_scale %in% c(5,10), cellSize %in% c(0.1, 1)) %>%
      #filter(Aou==this_sp) %>%
      group_by(temporal_scale, cellSize, windowID, modelName) %>%
      summarize(skill=mean(fss), skill_sd=sd(fss)) %>%
      ungroup() %>%
      mutate(time_lag=(temporal_scale * windowID) - (temporal_scale/2))

    
    results_this_sp$temporal_scale = factor(results_this_sp$temporal_scale, levels=c(5,10), labels = c('Temporal Scale: 5yrs','Temporal Scale: 10yrs'), ordered = TRUE)
    results_this_sp$cellSize = factor(results_this_sp$cellSize, levels=c(0.1,1.0), labels = c('Spatial Scale: 0.1°','Spatial Scale: 1.0°'), ordered = TRUE)
    
    sp_name=sppNames %>%
      dplyr::filter(Aou==this_sp) %>%
      extract2('name')
    
    this_sp_plot=ggplot(results_this_sp, aes(x=time_lag, y=skill, colour=modelName, group=modelName)) +
      geom_point(size=3) +
      geom_line(linetype='dotted', size=2) +
      geom_smooth(method='lm', se=FALSE, linetype='solid', size=1.3, aes(group=modelName)) +
      scale_colour_manual(values=c('#0072B2','#E69F00')) +
      #geom_hline(yintercept = 0.75) +
      theme_bw() +
      facet_grid(cellSize~temporal_scale, labeller = label_value) +
      ggtitle(sp_name) + 
      theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 20), 
        plot.title = element_text(size = 21), 
        legend.text = element_text(size = 23), 
        legend.key = element_rect(fill = "white"), 
        legend.background = element_rect(fill = "white"), 
        strip.text.x=element_text(size=22),
        strip.text.y=element_text(size=22),
        legend.position = "bottom", legend.direction = "horizontal") +labs(x = "Years Into Future", 
        y = "Fractions Skill Score", colour = NULL) + 
      theme(axis.text = element_text(size = 12))
    print(this_sp_plot)
  #}
}



this_sp_plot=ggplot(filter(results_this_sp,temporal_scale %in% c(5,10), cellSize %in% c(0.1, 1)), aes(x=time_lag, y=skill, colour=modelName, group=modelName)) +
  filter(temporal_scale %in% c(5,10), cellSize %in% c(0.1, 1))+
  geom_point() +
  geom_line(linetype='dotted') +
  geom_smooth(method='lm', se=FALSE, linetype='solid', size=0.9, aes(group=modelName)) +
  geom_hline(yintercept = 0.75) +
  theme_bw() +
  facet_grid(cellSize~temporal_scale, labeller=label_both) +
  ggtitle(sp_name)
print(this_sp_plot)


