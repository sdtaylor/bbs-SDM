library(dplyr)
library(tidyr)
library(doParallel)
library(magrittr)
library(stringr)
library(DBI)
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
#Define all the different model sets. A single set is a specific spatial/temporal grain combination.
#ie. Set 1 is a spatial grain size of 40km and temporal grain size of 1 year
#set 45 is (potentially) a spatial grain size of 100km and temporal grain size of 3 years.
#The model_sets data.frame also holds info on how the adjacent individual years will be aggregated into larger grains 
define_year_aggregations=function(year_list, data_type){
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
  model_sets$data_type=data_type
  return(model_sets)
}

model_sets = define_year_aggregations(training_years, data_type = 'training') %>%
   bind_rows(define_year_aggregations(testing_years, data_type = 'testing'))

#####################################################################
#Information on the spatial aggregation of indivudal sites.
#Also enforce a minimum number of sites within each spatiotemporal cell.
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
#Build site bioclim data for each grain size (defined by set id's)
#######################################################################
aggregated_bioclim_data=data.frame()

for(this_set_id in unique(model_sets$set_id)){
  this_set_years=model_sets %>% filter(set_id==this_set_id) 
  this_spatial_scale=unique(this_set_years$spatial_scale)
  this_temporal_scale=unique(this_set_years$temporal_scale)
  this_set_years=this_set_years %>% dplyr::select(Year, temporal_cell_id, data_type)
  
  #Temporal aggregation
  this_set_bioclim=bioclimData %>%
    filter(year %in% c(training_years, testing_years), spatial_scale==this_spatial_scale) %>%
    rename(Year=year) %>%
    left_join(this_set_years, by='Year') %>%
    group_by(temporal_cell_id, bioclim_var, spatial_cell_id, data_type) %>%
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
process_sp_observations=function(sp_observations, this_set_id, this_data_type){
  this_temporal_scale=model_sets %>% filter(set_id==this_set_id) %>% extract2('temporal_scale') %>% unique()
  this_spatial_scale=model_sets %>% filter(set_id==this_set_id) %>% extract2('spatial_scale') %>% unique()
  
  sp_observations$presence=1
  
  #Convert the years in the observations to temporal_cell_ids
  this_set_years=model_sets %>% 
    filter(set_id==this_set_id, data_type==this_data_type) %>% 
    dplyr::select(Year, temporal_cell_id)
  
  sp_observations = sp_observations %>%
    left_join(this_set_years, by='Year') %>%
    dplyr::select(-Year) %>%
    distinct()
  
  sp_observations$presence=1
  
  #Combine each observations with it's assocated bioclim data
  bioclim_data = aggregated_bioclim_data %>%
    dplyr::filter(set_id == this_set_id, data_type==this_data_type)
  
  sp_data = bioclim_data %>%
    left_join(sp_observations, by=c('spatial_cell_id','temporal_cell_id')) %>%
    mutate(presence=ifelse(is.na(presence), 0, 1)) %>%
    dplyr::filter(temporal_cell_id>0)
  
  sp_data = sp_data %>%
    spread(bioclim_var, value)
  
  return(sp_data)
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
#Iterate thru spp, building SDM's for each windowsize, and model.
#Parallel processing happens over the ~250 species
####################################################################
#Whether to write results to the postgres DB or keep in a local DF.
#Keeping it in a local DF is used for testing on a small number of species. 
writeToDB=FALSE

#focal_spp=c(7360, #Carolina chickadee
#           6010, #painted bunting
#            3100, #wild turky
#            4100 #Golden-fronted Woodpecker
#)

#finalDF=foreach(thisSpp=unique(occData$Aou)[1:3], .combine=rbind, .packages=c('dplyr','tidyr','magrittr','DBI')) %do% {
finalDF=foreach(thisSpp=c(3320, 4080, 5110), .combine=rbind, .packages=c('dplyr','tidyr','magrittr','DBI','RPostgreSQL','gbm')) %dopar% {
#finalDF=foreach(thisSpp=focal_spp, .combine=rbind, .packages=c('dplyr','tidyr','magrittr','DBI','RPostgreSQL')) %dopar% {
  thisSppResults=data.frame()
  for(this_set_id in unique(model_sets$set_id)){
    this_spatial_scale=model_sets %>% filter(set_id==this_set_id) %>% extract2('spatial_scale') %>% unique()
    this_temporal_scale=model_sets %>% filter(set_id==this_set_id) %>% extract2('temporal_scale') %>% unique()

    #Process the species observations to to presence/absence in aggregated cells & add in bioclim data
    this_sp_observations=dplyr::filter(occData, Aou==thisSpp, spatial_scale==this_spatial_scale)
    
    #Skip rare species that end up with a low sample size after all the  filtering
    training_sample_size = this_sp_observations %>%
      filter(Year %in% training_years) %>%
      nrow()
    if(training_sample_size < 100){next}
    
    thisSpp_training_data= this_sp_observations %>%
      filter(Year %in% training_years) %>%
      process_sp_observations(this_set_id=this_set_id, this_data_type = 'training')
    
    thisSpp_testing_data= this_sp_observations %>%
      filter(Year %in% testing_years) %>%
      process_sp_observations(this_set_id=this_set_id, this_data_type = 'testing')
    
    #Only use species that have >20 sites with at least 1 occurance
    #in the training period
  #  if( thisSppData %>%
  #          filter(temporal_cell_id==1, presence==1) %>%
  #          dplyr::select(cellID) %>%
  #          distinct() %>%
  #          nrow() < 20 ) { 
  #    print('Skipping, too few occurances')
   #   next }
    

    model=gbm(modelFormula, n.trees=5000, distribution = 'bernoulli', interaction.depth = 4, shrinkage=0.001, 
              data= thisSpp_training_data)
    perf=gbm.perf(model, plot.it=FALSE)
    #model=glm(modelFormula, family='binomial', data= thisSpp_training_data)
    thisSpp_testing_data$prediction = predict(model, n.trees=perf, newdata=thisSpp_testing_data, type='response')
    thisSpp_testing_data$Aou = thisSpp
    thisSpp_testing_data$temporal_scale = this_temporal_scale
    thisSpp_testing_data$spatial_scale  = this_spatial_scale
    
    
    score=fractions_skill_score(thisSpp_testing_data$presence, thisSpp_testing_data$prediction)

    if(save_verbose_results){
      thisSppResults = thisSppResults %>%
        bind_rows(select(thisSpp_testing_data, temporal_cell_id, spatial_cell_id, set_id, Aou, spatial_scale, temporal_scale, presence, prediction))
    } else {
      #Species and window size for this set of models. 
      thisSppResults = thisSppResults %>%
        bind_rows(data.frame(Aou=thisSpp,set_id=this_set_id, spatial_scale=this_spatial_scale, temporal_scale=this_temporal_scale, score=score) )
    }

  } 
  #Append results to the database results table or return to foreach()
  if(writeToDB){
    updateResults(thisSppResults)
  } else {
    return(thisSppResults)
  }
}

if(!writeToDB){
  write.csv(finalDF,resultsFile,row.names = FALSE)
}