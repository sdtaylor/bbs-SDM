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
  dplyr::select(Aou, siteID,year=Year, -RPID)

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
source('bioclim_utils.R')
bioclimData=get_train_bioclim_data()

#Some sites have na bioclim values. Mostly canada sites and the occasional one on water. 
bioclimData = bioclimData %>%
  filter(!is.na(bio1))


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
  
  thisSppData=dplyr::filter(occData, Aou==thisSpp)
  
  thisSppData = bioclimData %>%
    left_join(dplyr::filter(occData, Aou==thisSpp), by=c('siteID','year')) %>%
    filter(!is.na(bio1))
  
  
  for(thisSetID in modelSetMatrix$setID){
    this_spatial_scale=modelSetMatrix %>% filter(setID==thisSetID) %>% extract2('spatial_scale')
    this_temporal_scale=modelSetMatrix %>% filter(setID==thisSetID) %>% extract2('temporal_scale')
    this_offset=modelSetMatrix %>% filter(setID==thisSetID) %>% extract2('offset')
    
    
    #Process the data. excluding sites with low coverage, add bioclim variables, aggregating years into single widow size 
    #occurance, labeling those occurances, etc. 
    if(nrow(thisSppData)==0){next}
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
      next }
    
    #Most models need presence/absence treated as factors.
    thisSppData$presence = as.factor(thisSppData$presence)

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
    

    #Species and window size for this set of models. 
    modelResults = modelResults %>%
      mutate(Aou=thisSpp,setID=thisSetID, cellSize=this_spatial_scale, temporal_scale=this_temporal_scale, offset=this_offset) 
    

    
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
