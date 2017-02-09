library(dplyr)
library(tidyr)
library(doParallel)
library(magrittr)
library(stringr)
library(DBI)

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
  dataFolder='~/data/bbs/'
  numProcs=2
  resultsFile=paste('./results/bbsSDMResults.csv',sep='')
  #rawResultsFile='./results/bbsSDMResultsRaw.csv'
  rawResultsFile='./results/bbsSDMResultsRaw.sqlite'

} else if(args[1]=='local') {
  print('Running locally (probably cli)')
  dataFolder='~/data/bbs/'
  numProcs=2
  resultsFile=paste('./results/bbsSDMResults.csv',sep='')
  rawResultsFile='./results/bbsSDMResultsRaw.csv'
  
} else if(args[1]=='hipergator') {
  print('Running on hipergator')
  #dataFolder='/scratch/lfs/shawntaylor/data/bbs/'
  dataFolder='/ufrc/ewhite/shawntaylor/data/bbs/'
  numProcs=16
  resultsFile='./results/bbsSDMResults.csv'
  rawResultsFile='/scratch/lfs/shawntaylor/data/bbs/bbsSDMResults_ManyTimePeriods_Raw.csv'
  
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
routes=read.csv(paste(dataFolder, 'BBS_routes.csv', sep=''))
species=read.csv(paste(dataFolder, 'BBS_species.csv', sep=''))
weather=read.csv(paste(dataFolder, 'BBS_weather.csv', sep='')) %>%
  mutate(siteID=paste(countrynum, statenum, route,sep='-')) %>%
  dplyr::select(siteID,Year=year, RPID=rpid,runtype)
  
#Some records are of genus only and "unidentified". Get rid of those.
#TODO: Maybe find someway to incorperate them, because now some sites of false-negs. 
unidSpp=species %>%
  dplyr::select(Aou=AOU, name=english_common_name) %>%
  mutate(unID=ifelse(str_sub(name, 1,5)=='unid.', 1,0)  ) %>%
  filter(unID==1) %>%
  extract2('Aou')

#Filter weather to years of study so it can be used to calculate occData next.
weather=weather %>%
  filter(Year %in% c(testing_years, training_years))

occData=counts %>%
  filter(Year %in% c(training_years, testing_years)) %>%
  mutate(siteID=paste(countrynum, statenum, Route,sep='-')) %>%
  dplyr::select(Aou, siteID,Year, RPID)

#Remove anything not identified to species
occData=occData %>%
  filter(!Aou %in% unidSpp)

#Remove waterbirds and nocturnal birds.
occData=occData %>%
  filter(!Aou <=2880, !(Aou >=3650 & Aou<=3810), !(Aou>=3900 & Aou <=3910), !(Aou>=4160 & Aou <=4210), Aou!=7010)

#Remove any data where weather was not suitable or did not follow standard
#bbs protocal.
occData=occData %>%
  left_join(weather, by=c('siteID','Year','RPID')) %>%
  filter(runtype==1, RPID==101) %>%
  dplyr::select(-runtype, -RPID)

#A list of sites and the years they were sampled
siteList= occData %>%
  dplyr::select(Year, siteID) %>%
  distinct() %>%
  mutate(Year=as.factor(Year))

#bioclim values for all bbs routes from PRISM data.
bioclimData=get_bioclim_data()

#Some sites have na bioclim values. Mostly canada sites and the occasional one on water. 
bioclimData = bioclimData %>%
  filter(!is.na(value))


#####################################################################
# Define how years will be aggregated in all the temporal scales
# Also define "sets", sets being each unique spatiotemporal grain size. 
define_year_aggregations=function(year_list){
  model_sets=data.frame()
  set_id=1
  num_years = length(year_list)
  for(this_spatial_scale in spatial_scales){
    for(this_temporal_scale in temporal_scales){
      #The number of yearly sets given this window size and the number of years in study set.
      #ie for 10 years of data with 3 year window size = 3 sets (and 1 year leftover)
      numSets=floor(num_years/this_temporal_scale)
      
      cell_id=c()
      for(i in 1:this_temporal_scale){
        cell_id=c(cell_id, 1:numSets)
      }
      cell_id=sort(cell_id)
      
      #For window sizes not a multiple of the number of years in study, add -1 to remaining years.
      while(length(cell_id)<num_years){
        cell_id=c(cell_id,-1)
      }
  
      thisSetDF=data.frame(temporal_scale=this_temporal_scale, spatial_scale=this_spatial_scale, Year=year_list, temporal_cell_id=cell_id, set_id=set_id)
      set_id=set_id+1
      model_sets=bind_rows(model_sets, thisSetDF)
    }
  }
  return(model_sets)
}

model_sets = define_year_aggregations(testing_years) %>%
  bind_rows(define_year_aggregations(training_years))

#####################################################################
#Convert sites to cells
#Enforce a minimum number of sites within each cell.
#######################################################################
#What sites are within what cells at all scales.
#combine with what years each of those sites is sampled. 
spatial_grid_info=get_spatial_grid_info() %>%
  left_join(siteList, by='siteID') %>%
  group_by(spatial_cell_id, spatial_scale, Year) %>%
  summarize(num_sites=n()) %>%
  ungroup() %>%
  filter(!is.na(Year))

cell_size_median_sites=spatial_grid_info %>%
  group_by(spatial_scale, Year) %>%
  summarize(median_num_sites=median(num_sites))%>%
  ungroup()

#A list of cell's that will be used within each spatial scale
spatial_grid_info = spatial_grid_info %>%
  left_join(cell_size_median_sites, by=c('spatial_scale','Year')) %>%
  filter(num_sites>=median_num_sites) %>%
  dplyr::select(-num_sites, -median_num_sites)

#Where all the sites are actually located. 
site_id_cell_id=get_spatial_grid_info()

#####################################################################
#Build site information for each set_id
#Also ensure a minimum number of years covered in the temporal scales
#######################################################################
#Filter sites based on coverage within a particular temporal_cell_id
#calculate weather data for all those sites. 
aggregated_bioclim_data=data.frame()

for(this_set_id in unique(model_sets$set_id)){
  #Get the yearly aggregations to use (ie. set1: 80-84, set2: 85-89, etc) and other info about this set
  this_set_years=model_sets %>% filter(set_id==this_set_id) 
  this_spatial_scale=unique(this_set_years$spatial_scale)
  this_temporal_scale=unique(this_set_years$temporal_scale)
  this_set_years=this_set_years %>% dplyr::select(Year, temporal_cell_id)
  
  #Temporal aggregation
  this_set_bioclim=bioclimData %>%
    filter(year %in% c(training_years, testing_years), spatial_scale==this_spatial_scale) %>%
    rename(Year=year) %>%
    left_join(this_set_years, by='Year') %>%
    group_by(temporal_cell_id, bioclim_var, spatial_cell_id) %>%
    summarize(value=mean(value)) %>%
    ungroup() 

  this_set_bioclim$set_id=this_set_id
  
  aggregated_bioclim_data = bind_rows(aggregated_bioclim_data, this_set_bioclim)
  }
rm(this_set_years, this_temporal_scale, this_set_id, bioclimData)

###################################################################
#BBS occurance data is for individual sites. Need to convert all those to presences
#in cells at all the different spatial scales.
#Doing the same thing for the different temporal scales will eat up 10's of GB of RAM,
#so it's done on a species by species basis in process_sp_observations()
occData=occData %>%
  left_join(site_id_cell_id, by='siteID') %>%
  dplyr::select(-siteID) %>%
  distinct()

####################################################################
#Process data for a single species to a specific window size in prep for modeling.
#-Inserts absences; -splits years into specified window sizes (temporal_cell_id's);
#
###################################################################
process_sp_observations=function(sp_observations, this_set_id){
  this_temporal_scale=model_sets %>% filter(set_id==this_set_id) %>% extract2('temporal_scale')
  this_spatial_scale=model_sets %>% filter(set_id==this_set_id) %>% extract2('spatial_scale')
  
  sp_observations$presence=1
  
  #Combine each observations with it's assocated bioclim data
  sp_data = sp_observations %>%
    right_join(dplyr::filter(aggregated_bioclim_data, set_id==this_set_id), by='spatial_cell_id') %>%
    mutate(presence=ifelse(is.na(presence), 0, 1)) 
  
  
  
  
  
  #Get the window ID for all the years, which assigns years for temporal averaging according to the
  #temporal scale of this set ID
  thisSetYears=model_sets %>% 
    filter(set_id==this_set_id) %>% 
    gather(Year, temporal_cell_id, -temporal_scale, -set_id, -spatial_scale, -offset) %>%
    dplyr::select(Year, temporal_cell_id)
  
  sp_observations$Year=as.factor(sp_observations$Year) #change to factor to work in join
  
  #Summarize presence across the window ID's of this temporal scale
  sp_observations=sp_observations %>%
    left_join(thisSetYears, by='Year') %>%
    dplyr::select(Aou, cellID, temporal_cell_id) %>%
    distinct() %>%
    mutate(presence=1)
  
  #Merge with the site data matrix to get absences & bioclim data at the same time
  #this is a left_join here because siteDataMatrix only includes sites that have been
  #filtered for adequate coverage inside the temporal and spatial scales. 
  x=  siteDataMatrix %>%
    filter(set_id==this_set_id) %>%
    dplyr::select(-cellSize) %>%
    left_join(sp_observations, by=c('cellID','temporal_cell_id')) %>%  
    mutate(presence=ifelse(is.na(presence), 0, 1)) 
  
  return(x)
}

##################################################################
#Clear out the modelResults table from the database
#####################################################################
#Not working at the moment. Just need to do this manually whenever I run this.
#dbSendQuery(database$con, "TRUNCATE modelResults")

###################################################################
#Setup parallel processing
####################################################################
cl=makeCluster(numProcs)
registerDoParallel(cl)

#Sending results to the external database needs to be in a seperate function due to
#sql connection issues with multiple processes. 
updateResults=function(results){
  source('databaseConfig.R')
  database=src_postgres(dbname = dbName, host = dbHost, user = dbUser, password = dbPw)
  rm(dbName, dbHost, dbUser, dbPw)
  db_insert_into(con=database$con, table='modelResults', values=modelResults)
  #dbWriteTable(conn=database$con, name='modelResults', value=as.data.frame(results), append=TRUE, row.names=FALSE)
  dbDisconnect(database$con)
}

#####################################################################
#Iterate thru spp, building SDM's for each windowsize, offset, and model.
#Parallel processing happens over the ~250 species
####################################################################
#Whether to write results to the postgres DB or keep in a local DF.
#Keeping it in a local DF is used for testing on a small number of species. 
writeToDB=TRUE

#focal_spp=c(7360, #Carolina chickadee
#           6010, #painted bunting
#            3100, #wild turky
#            4100 #Golden-fronted Woodpecker
#)

#finalDF=foreach(thisSpp=unique(occData$Aou)[1:3], .combine=rbind, .packages=c('dplyr','tidyr','magrittr','DBI')) %do% {
finalDF=foreach(thisSpp=unique(occData$Aou), .combine=rbind, .packages=c('dplyr','tidyr','magrittr','DBI','RPostgreSQL')) %dopar% {
#finalDF=foreach(thisSpp=focal_spp, .combine=rbind, .packages=c('dplyr','tidyr','magrittr','DBI','RPostgreSQL')) %dopar% {
  thisSppResults=data.frame()
  for(this_set_id in model_sets$set_id){
    this_spatial_scale=model_sets %>% filter(set_id==this_set_id) %>% extract2('spatial_scale') %>% unique()
    this_temporal_scale=model_sets %>% filter(set_id==this_set_id) %>% extract2('temporal_scale') %>% unique()

    #Process the data. excluding sites with low coverage, add bioclim variables, aggregating years into single widow size 
    #occurance, labeling those occurances, etc. 
    this_sp_observations=dplyr::filter(occData, Aou==thisSpp, spatial_scale==this_spatial_scale)
    if(nrow(this_sp_observations)==0){next}
    thisSppData=process_sp_observations(sp_observations=this_sp_observations,this_set_id=this_set_id)
    thisSppData$Aou=thisSpp
    
    #Only use species that have >20 sites with at least 1 occurance
    #in the training period
    if( thisSppData %>%
            filter(temporal_cell_id==1, presence==1) %>%
            dplyr::select(cellID) %>%
            distinct() %>%
            nrow() < 20 ) { 
      print('Skipping, too few occurances')
      next }
    
    #Most models need presence/absence treated as factors.
    thisSppData$presence = as.factor(thisSppData$presence)

    # -1 temporal_cell_id was just years that didn't eventy divide into this window size. Don't need them any more
    #thisSppData = thisSppData %>%
    #  filter(temporal_cell_id != (-1))

    #A template to create results for each model in modelsToUse
    modelResultsTemplate=thisSppData %>%
      dplyr::select(presence, cellID, set_id, temporal_cell_id)
    
    #All the different model results for this set will be pasted together in here
    modelResults=data.frame()
    
    #Source the model script here inside the parallel loop so the packages get loaded in all the parallel threads.
    source('bbsSDMModels-ManyTimePeriods.R')
    
    #Iterate thru all the models. They are trained on T1 and tested/validated on T2-Tn. Results are returned for *all* timeperiods.
    #Don't record stuff that cause model errrors
    for(thisModel in modelsToUse){
      
      predictions=try(sdmModels(data=thisSppData, modelName=thisModel, modelFormula=modelFormula))
      if(class(predictions)!='try-error'){

        modelResultsTemplate$modelName=thisModel
        modelResultsTemplate$prediction=predictions
        
        modelResults=modelResults %>%
          bind_rows(modelResultsTemplate)
      }
    }
    
    #Setup a results dataframe for TV validation plot accuracy. 
    #for each site in each window ID this creates a T1_actual and T2_actual. 
    #Note there are many T2's to compare (from ever increasing gap between time), but only a single T1. 
    #This pulls out the T1 (always temporal_cell_id 1) and copies it over all the T2 (every other temporal_cell_id)
    #with left_join, and renames things accordingly
    
    #modelResults = modelResults %>%
    #  filter(temporal_cell_id!=1) %>%
    #  rename(T2_actual=presence, T2_prob=prediction) %>%
    #  left_join(  filter(modelResults, temporal_cell_id==1) %>% 
    #                dplyr::select(presence, cellID, modelName, prediction) %>% 
    #                rename(T1_actual=presence, T1_prob=prediction), by=c('cellID','modelName'))
    

    
    #Get temporal validation plot accuracy from Rapacciuolo et al. 2014
    #Not using this at the moment. So model results will be the predictions and 
    #observations for every species/site/model/windowsize. 
    #source('bbsSDM-temporal_validaton_accuracy.R')
    #modelResults=getAccuracyTV(modelResults)

    #Species and window size for this set of models. 
    modelResults = modelResults %>%
      mutate(Aou=thisSpp,set_id=this_set_id, cellSize=this_spatial_scale, temporal_scale=this_temporal_scale, offset=this_offset) 
    

    
    #Append results to the database results table or the results DF
    if(writeToDB){
      updateResults(modelResults)
    } else {
      #add results to final dataframe to be written as a csv - CSV Output
      thisSppResults=bind_rows(thisSppResults, modelResults)
    }
    
    
  } 
  #This gets returned to be added to the finalDF dataframe. 
  #Not needed when writing results to DB - CSV OutPut
  if(!writeToDB){
    return(thisSppResults)
  }
}

#Not needed when writing results to DB
#write.csv(finalDF,resultsFile,row.names = FALSE)
