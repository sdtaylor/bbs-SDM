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
  dataFolder='/scratch/lfs/shawntaylor/data/bbs/'
  #dataFolder='/ufrc/ewhite/shawntaylor/data/bbs/'
  numProcs=32
  resultsFile='./results/bbsSDMResults_with_timelags_prob_response.csv'
  rawResultsFile='/scratch/lfs/shawntaylor/data/bbs/bbsSDMResults_ManyTimePeriods_Raw.csv'
  
}

####################################################################
#Configuration
####################################################################
source('config.R')

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
  filter(Year %in% timeRange)

occData=counts %>%
  filter(Year %in% timeRange) %>%
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
source('get_prism_data.R')
bioclimData=get_bioclim_data()

#Some sites have na bioclim values. Mostly canada sites and the occasional one on water. 
bioclimData = bioclimData %>%
  filter(!is.na(bio1))


#####################################################################
#Build list of training/validation year sets
###################################################################
# Using the defined range of years and the window sizes above, define
# training and testing years for all combinations of window sizes.
#Then repeat for every available yearly offset. 
modelSetMatrix=data.frame()
setID=1

for(this_spatial_scale in spatial_cell_sizes){
  for(thisOffset in yearlyOffsets){
    for(this_temporal_scale in temporal_scales){
      #The number of yearly sets given this window size and the number of years in study set.
      #ie for 10 years of data with 3 year window size = 3 sets (and 1 year leftover)
      numSets=floor((length(timeRange)-thisOffset)/this_temporal_scale)
      
      #Produce an array partitioning each of the study years into a set for this window size.
      #the 1st set (identified by 1's) is added on at the end, so that excess years due to window
      #size not being a multiple of total years can be used as padding to reduce temporal autocorrelation.
      #Hopefully that makes sense. ask shawn to explain if not. 
      windowIdentifier=c()
      ones=c()
      for(i in 1:this_temporal_scale){
        windowIdentifier=c(windowIdentifier, 2:numSets)
        ones=c(ones, 1)
      }
      windowIdentifier=sort(windowIdentifier)
      
      #For window sizes not a multiple of the number of years in study, add -1 to remaining years.
      while(length(c(ones,windowIdentifier))<length(timeRange)-thisOffset){
        windowIdentifier=c(-1,windowIdentifier)
      }
      windowIdentifier=c(ones, windowIdentifier)
      
      #Pad the beginning offset years with -1
      while(length(windowIdentifier)<length(timeRange)){
        windowIdentifier=c(-1,windowIdentifier)
      }
      
      thisSetDF=data.frame(temporal_scale=this_temporal_scale, spatial_scale=this_spatial_scale, Year=timeRange, offset=thisOffset, windowID=windowIdentifier, setID=setID)
      setID=setID+1
      modelSetMatrix=bind_rows(modelSetMatrix, thisSetDF)
    }
  }
}

modelSetMatrix = modelSetMatrix %>%
  spread(Year, windowID)
rm(windowIdentifier, setID, thisSetDF)

#####################################################################
#Convert sites to cells
#Enforce a minimum number of sites within each cell.
#######################################################################
#What sites are within what cells at all scales.
#combine with what years each of those sites is sampled. 
spatial_grid_info=get_spatial_grid_info() %>%
  left_join(siteList, by='siteID') %>%
  group_by(cellID, cellSize, Year) %>%
  summarize(num_sites_in_cell=n()) %>%
  ungroup() %>%
  filter(!is.na(Year))

cell_size_median_sites=spatial_grid_info %>%
  group_by(cellSize, Year) %>%
  summarize(median_num_sites=median(num_sites_in_cell))%>%
  ungroup()

#A list of cell's that will be used within each spatial scale
spatial_grid_info = spatial_grid_info %>%
  left_join(cell_size_median_sites, by=c('cellSize','Year')) %>%
  filter(num_sites_in_cell>=median_num_sites) %>%
  dplyr::select(-median_num_sites)

#Where all the sites are actually located. 
site_id_cell_id=get_spatial_grid_info()

#####################################################################
#Build site information for each setID
#Also ensure a minimum number of years covered in the temporal scales
#######################################################################
#Filter sites based on coverage within a particular windowID

#Holds the bioclim variables for each cell with each temporal/spatial scale
siteDataMatrix=data.frame()

#This df holds how many sites are within each cell/windowID
total_sites_per_cell_window_id=data.frame()
  

for(thisSetID in unique(modelSetMatrix$setID)){
  #Get the yearly sets to use (ie. set1: 80-84, set2: 85-89, etc) and other info about this set
  thisSetYears=modelSetMatrix %>% filter(setID==thisSetID) %>% gather(Year, windowID, -temporal_scale, -setID, -spatial_scale, -offset)
  #thisWindowSize=modelSetMatrix %>% filter(setID==thisSetID) %>% extract2('windowSize')
  this_spatial_scale=modelSetMatrix %>% filter(setID==thisSetID) %>% extract2('spatial_scale')
  this_temporal_scale=modelSetMatrix %>% filter(setID==thisSetID) %>% extract2('temporal_scale')
  thisSetYears=thisSetYears %>% dplyr::select(Year, windowID)
  
  #List of sites and the number of years coverage they have within a particular windowID within this setID
  thisSetSiteInfo= spatial_grid_info %>%
    filter(cellSize==this_spatial_scale) %>%
    left_join(thisSetYears, by='Year') %>%
    group_by(windowID, cellID, cellSize) %>%
    summarize(n_years=n(), total_sites=sum(num_sites_in_cell)) %>% #total number of sites in this spatial cell AND within this temporal window of >= 1 years
    ungroup() %>%
    #This full_join fills in sites that are not covered *period* in the 1st training window. They need to have 0 value
    #to be excluded from the analysis in the next step.
    full_join( data.frame(windowID=1, cellID=unique(spatial_grid_info$cellID)), by=c('cellID','windowID')) 
  
  thisSetSiteInfo$n_years[is.na(thisSetSiteInfo$n_years)]=0
  
  #List of sites to drop completely because they don't have adequate coverage
  #in the training set
  dropSites= thisSetSiteInfo %>%
    filter(windowID==1 & n_years < this_temporal_scale*0.8) %>%
    extract2('cellID') %>%
    unique()

  thisSetSiteInfo=thisSetSiteInfo %>%
    filter(!cellID %in% dropSites) %>%
    #Of those remaining, drop any sites on a per set basis if they don't have
    #adequate coverage in a particule set
    filter((windowID > 1 & n_years>= this_temporal_scale*0.66) | (windowID==1)) %>%
    dplyr::select(-n_years)

  #Save this info to use the total site samples per cell and temporal window for later. 
  total_sites_per_cell_window_id = thisSetSiteInfo %>%
    mutate(setID=thisSetID) %>%
    bind_rows(total_sites_per_cell_window_id)
  
  #Average all bioclim over this spatial scale. subset to sites
  #with adequate temporal coverage. 
  thisSetWeather=bioclimData %>%
    filter(year %in% timeRange, cellSize==this_spatial_scale) %>%
    rename(Year=year) %>%
    mutate(Year=as.factor(Year)) %>%
    left_join(thisSetYears, by='Year') %>%
    group_by(windowID, cellID) %>%
    summarize_each(funs(mean), -Year,-cellID, -cellSize) %>%
    ungroup() %>%
    filter(cellID %in% thisSetSiteInfo$cellID)
  
  thisSetSiteInfo = thisSetSiteInfo %>%
    left_join(thisSetWeather, by=c('windowID','cellID')) %>%
    filter(!is.na(bio1))
  
  thisSetSiteInfo$setID=thisSetID
  
  siteDataMatrix = bind_rows(siteDataMatrix, thisSetSiteInfo)
  }
rm(thisSetSiteInfo, thisSetWeather, thisSetYears, this_temporal_scale, dropSites, thisSetID, bioclimData)

###################################################################
#BBS occurance data is for individual sites. Need to convert all those to presences
#in cells at all the different spatial scales while also accounting for the different 
#Temporal scales.
occData=occData %>%
  left_join(site_id_cell_id, by='siteID') %>%
  dplyr::select(-siteID) 
  #distinct()

####################################################################
#Process data for a single species to a specific window size in prep for modeling.
#-Inserts absences; -splits years into specified window sizes (windowID's);
#
###################################################################
processSpDataToWindowSize=function(spData, thisSetID){
  this_temporal_scale=modelSetMatrix %>% filter(setID==thisSetID) %>% extract2('temporal_scale')
  this_spatial_scale=modelSetMatrix %>% filter(setID==thisSetID) %>% extract2('spatial_scale')
  
  #Get the window ID for all the years, which assigns years for temporal averaging according to the
  #temporal scale of this set ID
  thisSetYears=modelSetMatrix %>% 
    filter(setID==thisSetID) %>% 
    gather(Year, windowID, -temporal_scale, -setID, -spatial_scale, -offset) %>%
    dplyr::select(Year, windowID)
  
  spData$Year=as.factor(spData$Year) #change to factor to work in join
  
  #Summarize presence across the window ID's of this temporal scale
  spData=spData %>%
    left_join(thisSetYears, by='Year') %>%
    dplyr::select(Aou, cellID, windowID) %>%
    #distinct() %>%
    mutate(presence=1)
  
  #Fill in absences by subtracting the total number of sites sampled in `total_sites_per_cell_window_id`
  #with the presences in spData for each cell/windowID. 
  num_sampled_sites=total_sites_per_cell_window_id %>%
    filter(setID==thisSetID)
  
  for(i in 1:nrow(num_sampled_sites)){
    this_window_id=num_sampled_sites$windowID[i]
    this_cell_id  =num_sampled_sites$cellID[i]
    total_sites   =num_sampled_sites$total_sites[i]
    
    current_samples=spData %>%
      filter(windowID==this_window_id, cellID==this_cell_id) %>%
      nrow()
    
    difference=total_sites-current_samples
    if(difference < 0){ 
      warning(paste('total sites < current samples for some reason', difference))
    } else {
    absences= data.frame(windowID=this_window_id, cellID=this_cell_id, presence=0)
    absences= do.call('rbind', replicate(difference, absences, simplify = FALSE ))
  
    spData = spData %>%
      bind_rows(absences)
    } 
  }
  
  spData = spData %>%
    group_by(cellID, windowID) %>%
    summarise(presence=mean(presence)) %>%
    ungroup()
  
  #Merge with the site data matrix to get absences & bioclim data at the same time
  #this is a left_join here because siteDataMatrix only includes sites that have been
  #filtered for adequate coverage inside the temporal and spatial scales. 
  x=  siteDataMatrix %>%
    filter(setID==thisSetID) %>%
    dplyr::select(-cellSize) %>%
    left_join(spData, by=c('cellID','windowID')) 
    #mutate(presence=ifelse(is.na(presence), 0, 1)) 
  
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
  dbWriteTable(conn=database$con, name='modelResults', value=as.data.frame(results), append=TRUE, row.names=FALSE)
  dbDisconnect(database$con)
}

#####################################################################
#Iterate thru spp, building SDM's for each windowsize, offset, and model.
#Parallel processing happens over the ~250 species
####################################################################
#Whether to write results to the postgres DB or keep in a local DF.
#Keeping it in a local DF is used for testing on a small number of species. 
writeToDB=FALSE

focal_spp=c(5840, #swamp sparrow. no shift
            7350, #Black capped chickadee. no shift
            3100, #wild turkey. northern shift (only in La Sorte & Thompson 2007 )
            6100, #Summer Tanager. northern shift
            7190, #bewicks wren. southern shift 
            5170 #Purple finch. no shift
)

parallel_process_iteration=expand.grid(Aou=focal_spp, setID=unique(modelSetMatrix$setID))

#finalDF=foreach(thisSpp=unique(occData$Aou)[1:3], .combine=rbind, .packages=c('dplyr','tidyr','magrittr','DBI')) %do% {
#finalDF=foreach(thisSpp=unique(occData$Aou)[1:2], .combine=rbind, .packages=c('dplyr','tidyr','magrittr','DBI','RPostgreSQL')) %dopar% {
finalDF=foreach(i=1:nrow(parallel_process_iteration), .combine=rbind, .packages=c('dplyr','tidyr','magrittr','DBI','RPostgreSQL')) %dopar% {

    thisSpp=parallel_process_iteration$Aou[i]
    thisSetID=parallel_process_iteration$setID[i]
    
    this_spatial_scale=modelSetMatrix %>% filter(setID==thisSetID) %>% extract2('spatial_scale')
    this_temporal_scale=modelSetMatrix %>% filter(setID==thisSetID) %>% extract2('temporal_scale')
    this_offset=modelSetMatrix %>% filter(setID==thisSetID) %>% extract2('offset')
    
    
    #Process the data. excluding sites with low coverage, add bioclim variables, aggregating years into single widow size 
    #occurance, labeling those occurances, etc. 
    thisSppData=dplyr::filter(occData, Aou==thisSpp, cellSize==this_spatial_scale)
    if(nrow(thisSppData)==0){return(NA)}
    thisSppData=processSpDataToWindowSize(spData=thisSppData,thisSetID=thisSetID)
    thisSppData$Aou=thisSpp
    
    #Only use species that have >20 sites with at least 1 occurance
    #in the training period
    if( thisSppData %>%
            filter(windowID==1, presence==1) %>%
            dplyr::select(cellID) %>%
            distinct() %>%
            nrow() < 20 ) { 
      print('Skipping, too few occurances')
      return(NA) }
    
    #Most models need presence/absence treated as factors.
    #thisSppData$presence = as.numeric(thisSppData$presence)

    # -1 windowID was just years that didn't eventy divide into this window size. Don't need them any more
    #thisSppData = thisSppData %>%
    #  filter(windowID != (-1))

    #A template to create results for each model in modelsToUse
    modelResultsTemplate=thisSppData %>%
      dplyr::select(presence, cellID, setID, windowID)
    
    #All the different model results for this set will be pasted together in here
    modelResults=data.frame()
    
    #Source the model script here inside the parallel loop so the packages get loaded in all the parallel threads.
    source('bbsSDMModels-ManyTimePeriods.R')
    
    #Iterate thru all the models. They are trained on T1 and tested/validated on T2-Tn. Results are returned for *all* timeperiods.
    for(thisModel in modelsToUse){
      predictions=sdmModels(data=thisSppData, modelName=thisModel, modelFormula=modelFormula)
      #colnames(predictions)=c(paste(thisModel,'prob',sep='-'))
      
      modelResultsTemplate$modelName=thisModel
      modelResultsTemplate$prediction=predictions
      
      modelResults=modelResults %>%
        bind_rows(modelResultsTemplate)
    }

    #site_totals=thisSppData %>%
    #  dplyr::select(windowID, cellID, total_sites) %>%
    #  dplyr::distinct
    
    
    
    modelResults = modelResults %>%
      group_by(cellID, windowID, modelName) %>%
      summarise(presence=mean(as.numeric(as.character(presence))), prediction_mean=mean(prediction), predicion_sd=sd(prediction), n=n()) %>%
      ungroup() %>%
      mutate(error=(prediction_mean-presence)^2) %>%
      filter(windowID>1) %>%
      group_by(modelName, windowID) %>%
      summarize(mse=mean(error), num_cubes=n(), sum_of_sq=mean(prediction_mean^2 + presence^2)) %>%
      ungroup() %>%
      mutate(fss = 1-(mse/sum_of_sq))

    #Setup a results dataframe for TV validation plot accuracy. 
    #for each site in each window ID this creates a T1_actual and T2_actual. 
    #Note there are many T2's to compare (from ever increasing gap between time), but only a single T1. 
    #This pulls out the T1 (always windowID 1) and copies it over all the T2 (every other windowID)
    #with left_join, and renames things accordingly
    
    #modelResults = modelResults %>%
    #  filter(windowID!=1) %>%
    #  rename(T2_actual=presence, T2_prob=prediction) %>%
    #  left_join(  filter(modelResults, windowID==1) %>% 
    #                dplyr::select(presence, cellID, modelName, prediction) %>% 
    #                rename(T1_actual=presence, T1_prob=prediction), by=c('cellID','modelName'))
    

    
    #Get temporal validation plot accuracy from Rapacciuolo et al. 2014
    #Not using this at the moment. So model results will be the predictions and 
    #observations for every species/site/model/windowsize. 
    #source('bbsSDM-temporal_validaton_accuracy.R')
    #modelResults=getAccuracyTV(modelResults)

    #Species and window size for this set of models. 
    modelResults = modelResults %>%
      mutate(Aou=thisSpp,setID=thisSetID, cellSize=this_spatial_scale, temporal_scale=this_temporal_scale, offset=this_offset) 
    

    
    #Append results to the database results table or return to finalDF if testing
    if(writeToDB){
      updateResults(modelResults)
    } else {
      return(modelResults)
    }
    
}

#Not needed when writing results to DB
write.csv(finalDF,resultsFile,row.names = FALSE)
