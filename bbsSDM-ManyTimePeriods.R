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
  numProcs=32
  resultsFile=paste('/scratch/lfs/shawntaylor/data/bbs/bbsSDMResults_ManyTimePeriods_Raw.csv',sep='')
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
if(length(timeRange)/max(windowSizes) < 2 ){ stop(paste('Window size',max(windowSizes),'too long for',length(timeRange),'year time range',sep=' ')) }

for(thisOffset in yearlyOffsets){
  for(thisWindowSize in windowSizes){
    #The number of yearly sets given this window size and the number of years in study set.
    #ie for 10 years of data with 3 year window size = 3 sets (and 1 year leftover)
    numSets=floor((length(timeRange)-thisOffset)/thisWindowSize)
    
    #Produce an array partitioning each of the study years into a set for this window size.
    #the 1st set (identified by 1's) is added on at the end, so that excess years due to window
    #size not being a multiple of total years can be used as padding to reduce temporal autocorrelation.
    #Hopefully that makes sense. ask shawn to explain if not. 
    windowIdentifier=c()
    ones=c()
    for(i in 1:thisWindowSize){
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
    
    thisSetDF=data.frame(windowSize=thisWindowSize, Year=timeRange, windowID=windowIdentifier, setID=setID)
    setID=setID+1
    modelSetMatrix=bind_rows(modelSetMatrix, thisSetDF)
  }
}

modelSetMatrix = modelSetMatrix %>%
  spread(Year, windowID)
rm(numSets, windowIdentifier, thisSetDF, ones)

#####################################################################
#Build site information for each setID
#######################################################################
#Filter sites based on coverage within a particular windowID
#calculate weather data for all those sites. 
siteDataMatrix=data.frame()
for(thisSetID in unique(modelSetMatrix$setID)){
  #Get the yearly sets to use (ie. set1: 80-84, set2: 85-89, etc) and the sets window size
  thisSetYears=modelSetMatrix %>% filter(setID==thisSetID) %>% gather(Year, windowID, -windowSize, -setID)
  thisWindowSize=thisSetYears %>% dplyr::select(windowSize) %>% distinct() %>% extract2('windowSize')
  thisSetYears=thisSetYears %>% dplyr::select(Year, windowID)
  
  #List of sites and the number of years coverage they have within a particular windowID within this setID
  thisSetSiteInfo= siteList %>%
    left_join(thisSetYears, by='Year') %>%
    group_by(windowID, siteID) %>%
    summarize(nYears=n()) %>%
    #This fills in sites that are not covered *period* in the 1st training window. They need to have 0 value
    #to be excluded from the analysis in the next step.
    full_join( data.frame(windowID=1, siteID=unique(siteList$siteID)), by=c('siteID','windowID')) 
  
  thisSetSiteInfo$nYears[is.na(thisSetSiteInfo$nYears)]=0
  
  #List of sites to drop completely because they don't have adequate coverage
  #in the training set( set 1)
  dropSites= thisSetSiteInfo %>%
    filter(windowID==1 & nYears < thisWindowSize*0.66) %>%
    extract2('siteID') %>%
    unique()

  thisSetSiteInfo=thisSetSiteInfo %>%
    filter(!siteID %in% dropSites) %>%
    #Of those remaining, drop any sites on a per set basis if they don't have
    #adequate coverage in a particule set
    filter((windowID > 1 & nYears >= thisWindowSize*0.66) | (windowID<=1))
  
  thisSetWeather=bioclimData %>%
    filter(year %in% timeRange) %>%
    rename(Year=year) %>%
    mutate(Year=as.factor(Year)) %>%
    left_join(thisSetYears, by='Year') %>%
    group_by(windowID, siteID) %>%
    summarize_each(funs(mean), -Year,-siteID)
  
  thisSetSiteInfo = thisSetSiteInfo %>%
    left_join(thisSetWeather, by=c('windowID','siteID')) %>%
    filter(!is.na(bio1))
  
  thisSetSiteInfo$setID=thisSetID
  
  siteDataMatrix = bind_rows(siteDataMatrix, thisSetSiteInfo)
  }
rm(thisSetSiteInfo, thisSetWeather, thisSetYears, thisWindowSize, dropSites)

####################################################################
#Process data for a single species to a specific window size in prep for modeling.
#-Inserts absences; -splits years into specified window sizes (windowID's);
#
###################################################################
processSpDataToWindowSize=function(spData, thisSetID){
  #Get the window ID for all the years
  thisSetYears=modelSetMatrix %>% filter(setID==thisSetID) %>% gather(Year, windowID, -windowSize, -setID)
  thisWindowSize=thisSetYears %>% dplyr::select(windowSize) %>% distinct() %>% extract2('windowSize')
  thisSetYears=thisSetYears %>% dplyr::select(Year, windowID)
  
  spData$Year=as.factor(spData$Year) #change to factor to work in join
  
  #Get presence data for this species summarized by the windowIDs
  x=spData %>%
    left_join(thisSetYears, by='Year') %>%
    dplyr::select(Aou, siteID, windowID) %>%
    distinct() %>%
    mutate(presence=1)
  
  #Merge with the site data matrix to get absences & bioclim data at the same time
  #this is a left_join here because siteDataMatrix only includes sites that have been
  #filtered for adequate coverage inside this window size. 
  x=  siteDataMatrix %>%
    filter(setID==thisSetID) %>%
    left_join(x, by=c('siteID','windowID')) %>%  
    mutate(presence=ifelse(is.na(presence), 0, 1)) %>%
    dplyr::select(-nYears)
  
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
  dbWriteTable(conn=database$con, name='modelResults', value=as.data.frame(results), append=TRUE)
  dbDisconnect(database$con)
}

#####################################################################
#Iterate thru spp, building SDM's for each windowsize, offset, and model.
#Parallel processing happens over the ~250 species
####################################################################
#finalDF=foreach(thisSpp=unique(occData$Aou)[1:3], .combine=rbind, .packages=c('dplyr','tidyr','magrittr','DBI')) %do% {
finalDF=foreach(thisSpp=unique(occData$Aou)[1:10], .combine=rbind, .packages=c('dplyr','tidyr','magrittr','DBI','RPostgreSQL')) %dopar% {
  thisSppResults=data.frame()
  for(thisSetID in modelSetMatrix$setID){

    thisWindowSize=modelSetMatrix$windowSize[modelSetMatrix$setID==thisSetID]
    #Process the data. excluding sites with low coverage, add bioclim variables, aggregating years into single widow size 
    #occurance, labeling those occurances, etc. 
    thisSppData=filter(occData, Aou==thisSpp)
    thisSppData=processSpDataToWindowSize(spData=thisSppData,thisSetID=thisSetID)
    thisSppData$Aou=thisSpp
    
    #Only use species that have >20 occurances in time period 1
    if( sum(as.integer(as.character(thisSppData %>%
            filter(windowID==1) %>%
            extract2('presence'))))  < 20 ) { 
      print('Skipping, too few occurances')
      next }
    
    #Most models need presence/absence treated as factors.
    thisSppData$presence = as.factor(thisSppData$presence)

    # -1 windowID was just years that didn't eventy divide into this window size. Don't need them any more
    thisSppData = thisSppData %>%
      filter(windowID != (-1))

    #A template to create results for each model in modelsToUse
    modelResultsTemplate=thisSppData %>%
      dplyr::select(windowID, presence, siteID, setID)
    
    #All the different model results will be pasted together in here
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
    
    #Setup a results dataframe for TV validation plot accuracy. 
    #for each site in each window ID this creates a T1_actual and T2_actual. 
    #Note there are many T2's to compare (from ever increasing gap between time), but only a single T1. 
    #This pulls out the T1 (always windowID 1) and copies it over all the T2 (every other windowID)
    #with left_join, and renames things accordingly
    modelResults = modelResults %>%
      filter(windowID!=1) %>%
      rename(T2_actual=presence, T2_prob=prediction) %>%
      left_join(  filter(modelResults, windowID==1) %>% 
                    dplyr::select(presence, siteID, modelName, prediction) %>% 
                    rename(T1_actual=presence, T1_prob=prediction), by=c('siteID','modelName'))
    

    
    #Get temporal validation plot accuracy from Rapacciuolo et al. 2014
    #Not using this at the moment. So model results will be the predictions and 
    #observations for every species/site/model/windowsize. 
    #source('bbsSDM-temporal_validaton_accuracy.R')
    #modelResults=getAccuracyTV(modelResults)

    #Species and window size for this set of models. 
    modelResults = modelResults %>%
      mutate(Aou=thisSpp, windowSize=thisWindowSize, setID=thisSetID) 
    
    #add results to final dataframe to be written as a csv - CSV Output
    #thisSppResults=bind_rows(thisSppResults, modelResults)
    
    #Append results to the database results table. - PostGRES output
    updateResults(modelResults)
  } 
  #This gets returned to be added to the finalDF dataframe. 
  #Not needed when writing results to DB - CSV OutPut
  #return(thisSppResults)
}

#Not needed when writing results to DB
#write.csv(finalDF,resultsFile,row.names = FALSE)
