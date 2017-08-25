library(dplyr)
library(tidyr)
library(doParallel)
library(magrittr)
library(stringr)
library(gbm)

####################################################################################
#Adjustments based on where this script is being run.
#################################################################################
#If running this on hipergator, use "Rscript <thisScript> hipergator" 
#if running locally use 'Rscript <thisScript> local'
#if running in rstudio nothing special is needed. 

args=commandArgs(trailingOnly = TRUE)

#If the 1st argument is na (ie, no argument), then this script is being run inside rstudio
if(is.na(args[1])){
  print('Running locally (probably rstudio)')
  numProcs=2
} else if(args[1]=='local') {
  print('Running locally (probably cli)')
  dataFolder='~/data/bbs/'
  numProcs=2
} else if(args[1]=='hipergator') {
  print('Running on hipergator')
  numProcs=16
}

####################################################################
#Configuration
####################################################################
source('config.R')
source('utils.R')

#################################################################
#Data loading and pre-processing
################################################################

counts=read.csv(paste(dataFolder, 'BBS_counts.csv', sep=''))
routes=read.csv(paste(dataFolder, 'BBS_routes.csv', sep='')) %>%
  dplyr::mutate(siteID=paste(countrynum, statenum, route,sep='-')) %>%
  dplyr::select(siteID, long=loni, lat=lati)
species=read.csv(paste(dataFolder, 'BBS_species.csv', sep=''))
weather=read.csv(paste(dataFolder, 'BBS_weather.csv', sep='')) %>%
  mutate(siteID=paste(countrynum, statenum, route,sep='-')) %>%
  dplyr::select(siteID,Year=year, RPID=rpid,runtype)
  
#Spatial route object to use where needed
routes_spatial = SpatialPointsDataFrame(cbind(routes$long, routes$lat), data=routes, 
                                proj4string = CRS('+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0'))

#Some records are of genus only and "unidentified". Get rid of those.
unidSpp=species %>%
  dplyr::select(Aou=AOU, name=english_common_name) %>%
  mutate(unID=ifelse(str_sub(name, 1,5)=='unid.', 1,0)  ) %>%
  filter(unID==1) %>%
  extract2('Aou')

#Filter weather to years of study so it can be used to calculate occData next.
weather=weather %>%
  filter(Year %in% c(testing_years, training_years))

occData=counts %>%
  filter(Year %in% c(testing_years, training_years)) %>%
  mutate(siteID=paste(countrynum, statenum, Route,sep='-')) %>%
  dplyr::select(Aou, siteID,year=Year, RPID)

#Remove anything not identified to species
occData=occData %>%
  filter(!Aou %in% unidSpp)

#Remove waterbirds and nocturnal birds.
occData=occData %>%
  filter(!Aou <=2880, !(Aou >=3650 & Aou<=3810), !(Aou>=3900 & Aou <=3910), !(Aou>=4160 & Aou <=4210), Aou!=7010)

#Remove any data where weather was not suitable or did not follow standard
#bbs protocal.
occData=occData %>%
  left_join(weather, by=c('siteID','year'='Year','RPID')) %>%
  filter(runtype==1, RPID==101) %>%
  dplyr::select(-runtype, -RPID)

#A list of sites and the years they were sampled. To be used in making absences for 
#species occurances. Keep only sites with the minium specified number of years in the
#10 year period.
all_routes_surveyed = occData %>%
  dplyr::select(year, siteID) %>%
  distinct() %>%
  mutate(timeframe = ifelse(year %in% training_years, 'training','testing')) %>%
  group_by(timeframe, siteID) %>%
  summarize(num_years_sampled = n()) %>%
  ungroup() %>%
  filter(num_years_sampled>=minimum_years) %>%
  select(-num_years_sampled)

#Aggregate the  occurances down to presence in training/testing periods only
occData = occData %>%
  mutate(timeframe = ifelse(year %in% training_years, 'training','testing')) %>%
  select(Aou, siteID, timeframe) %>%
  distinct()
  
#bioclim values for all bbs routes from PRISM data.
bioclim_data=get_bioclim_data()

#Some sites have na bioclim values. Mostly canada sites and the occasional one on water. 
bioclim_data = bioclim_data %>%
  filter(!is.na(bio1))

#Only keep data for sites that were actually surveyed
bioclim_data = bioclim_data %>%
  filter(siteID %in% all_routes_surveyed$siteID)


  # dplyr::select(year, siteID) %>%
  # dplyr::mutate(keep='yes') %>%
  # dplyr::right_join(bioclim_data, by=c('siteID','year')) %>%
  # dplyr::filter(keep=='yes') %>%
  # dplyr::select(-keep)

#Split into train and test and aggregate
bioclim_data_training = bioclim_data %>%
  dplyr::filter(year %in% training_years) %>%
  select(-year) %>%
  group_by(siteID) %>%
  summarise_all(mean)

bioclim_data_testing = bioclim_data %>%
  dplyr::filter(year %in% testing_years) %>%
  select(-year) %>%
  group_by(siteID) %>%
  summarise_all(mean)

######################################################################
#Prepare spatial and temporal aggregation info for all sites
#######################################################################

#Assign each site to a cell for each spatial scale. 
spatial_grid_info = data.frame()

for(this_spatial_scale in spatial_scales){
  grid = create_grid(this_spatial_scale, templateRaster)
  grid_info = assign_sites_to_grid(grid, routes_spatial)
  grid_info$spatial_scale = this_spatial_scale 
  
  spatial_grid_info = spatial_grid_info %>%
    dplyr::bind_rows(grid_info)
}

#Assign each year to a temporal cell.
# temporal_grid_info = data.frame()
# for(this_temporal_scale in temporal_scales){
#   num_temporal_cells = floor(length(testing_years)/this_temporal_scale)
#   
#   cell_identifier=c()
#   for(i in 1:this_temporal_scale){
#     cell_identifier=c(cell_identifier, 1:num_temporal_cells)
#   }
#   cell_identifier=sort(cell_identifier)
#   
#   while(length(cell_identifier) < length(testing_years)){
#     cell_identifier = c(cell_identifier,NA)
#   }
#   
#   temporal_grid_info = temporal_grid_info %>%
#     dplyr::bind_rows(data.frame(temporal_cell_id = cell_identifier, year=testing_years, temporal_scale = this_temporal_scale))
# }

all_routes_surveyed = all_routes_surveyed %>%
  left_join(spatial_grid_info, by='siteID')
  #left_join(temporal_grid_info, by='year')

#Remove some NA values due to some temporal scales not dividing equally into number of testing years
all_routes_surveyed = all_routes_surveyed[complete.cases(all_routes_surveyed),]

#Restrict upscaled cells in the verification to ones that have a minimum number of sites within them.
#Minimum sites are the median number of sites for that spatial scale
median_cell_counts = all_routes_surveyed %>%
  dplyr::filter(timeframe == 'testing') %>%
  dplyr::group_by(spatial_scale, spatial_cell_id) %>%
  dplyr::summarize(num_sites = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(spatial_scale) %>%
  dplyr::summarise(median_site_count = median(num_sites)) %>%
  dplyr::ungroup()

#this data.frame stores how the sites are aggregated within each spatial grain size.
#Also filters out any cells that don't meet the minimum site requirement.  
testing_cell_info = all_routes_surveyed %>%
  dplyr::group_by(spatial_scale, spatial_cell_id) %>%
  dplyr::tally() %>%
  dplyr::ungroup() %>%
  dplyr::left_join(median_cell_counts, by='spatial_scale') %>%
  dplyr::filter(n>=median_site_count) %>%
  dplyr::select(-median_site_count, -n) %>%
  dplyr::left_join(all_routes_surveyed, by=c('spatial_scale', 'spatial_cell_id'))

#Final tallys of effective sample size. for information only. 
#final_counts = testing_cell_info %>%
#  dplyr::group_by(spatial_scale, temporal_scale) %>%
#  dplyr::summarize(n=n())

rm(grid, grid_info, num_temporal_cells, cell_identifier, spatial_grid_info, median_cell_counts)

###################################################################
#Setup parallel processing
####################################################################
cl=makeCluster(numProcs)
registerDoParallel(cl)

#####################################################################
#Iterate thru spp, building SDM's for each windowsize, offset, and model.
#Parallel processing happens over the ~250 species
####################################################################
focal_spp=c(7360, #Carolina chickadee
           6010, #painted bunting
            #3100, #wild turky
            4100 #Golden-fronted Woodpecker
)

return_site_level_predictions=TRUE

#finalDF=foreach(thisSpp=unique(occData$Aou)[1:3], .combine=rbind, .packages=c('dplyr','tidyr','magrittr','DBI')) %do% {
finalDF=foreach(thisSpp=focal_spp, .combine=rbind, .packages=c('dplyr','tidyr','gbm')) %dopar% {
#finalDF=foreach(thisSpp=focal_spp, .combine=rbind, .packages=c('dplyr','tidyr','magrittr','DBI','RPostgreSQL')) %dopar% {
  this_spp_results=data.frame()
  
  thisSpp_occurances = occData %>%
    dplyr::filter(Aou==thisSpp) %>%
    dplyr::mutate(presence=1)
  
  thisSpp_data_training = thisSpp_occurances %>%
    filter(timeframe == 'training') %>%
    dplyr::right_join(bioclim_data_training, by = 'siteID') %>%
    dplyr::mutate(presence = ifelse(is.na(presence), 0, presence))
  
  #Skip rare species that end up with a low sample size after all the  filtering
  if(sum(thisSpp_data_training$presence)<100){
    return(data.frame())
  }
  
  model=gbm(modelFormula, n.trees=5000, distribution = 'bernoulli', interaction.depth = 4, shrinkage=0.001, 
            data= thisSpp_data_training)
  perf=gbm.perf(model, plot.it=FALSE)
  #model = glm(modelFormula, family='binomial', data=thisSpp_data_training)
  
  thisSpp_data_testing = thisSpp_occurances %>%
    filter(timeframe == 'testing') %>%
    dplyr::right_join(bioclim_data_testing, by = 'siteID') %>%
    dplyr::mutate(presence = ifelse(is.na(presence), 0, presence))

  predictions = thisSpp_data_testing %>%
    dplyr::select(siteID, presence) 
  
  predictions$prediction = predict(model, n.trees=perf, newdata=thisSpp_data_testing, type='response')
  
    for(this_spatial_scale in spatial_scales){
      
    scaled_prediction = testing_cell_info %>%
      dplyr::filter(spatial_scale == this_spatial_scale) %>%
      dplyr::left_join(predictions, by='siteID') %>%
      dplyr::filter(!is.na(presence))
      #group_by(spatial_cell_id, temporal_cell_id) %>%
      #summarize(presence = max(presence), prediction = max((prediction>0.5)*1) ) %>%
      #ungroup()
    
    scaled_prediction$spatial_scale = this_spatial_scale

    score=fractions_skill_score(scaled_prediction$presence, scaled_prediction$prediction)

    if(return_site_level_predictions){
      this_spp_results = this_spp_results %>%
        dplyr::bind_rows(scaled_prediction)
    } else{
      this_spp_results = this_spp_results %>%
        dplyr::bind_rows(data.frame(Aou=thisSpp, spatial_scale=this_spatial_scale, score=score))
    }
  }

  this_spp_results$Aou=thisSpp
  return(this_spp_results)
 
}

finalDF$dataset='bbs_method1'
write_csv(finalDF, resultsFile)


