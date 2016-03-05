library(dplyr)
library(tidyr)
library(doParallel)
library(magrittr)
library(stringr)

#If running this on hipergator, use "Rscript <thisScript> hipergator" 
args=commandArgs(trailingOnly = TRUE)

#If the 1st argument is na (ie, no argument), then this script is being run inside rstudo
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
  
#If not na, then check to see if the 1st arg means running in the ufl hpc.  
} else if(args[1]=='hipergator') {
  print('Running on hipergator')
  dataFolder='/scratch/lfs/shawntaylor/data/bbs/'
  numProcs=32
  resultsFile=paste('/scratch/lfs/shawntaylor/data/bbs/bbsSDMResults_ManyTimePeriods_Raw.csv',sep='')
  rawResultsFile='/scratch/lfs/shawntaylor/data/bbs/bbsSDMResults_ManyTimePeriods_Raw.csv'
  
}


counts=read.csv(paste(dataFolder, 'BBS_counts.csv', sep=''))
routes=read.csv(paste(dataFolder, 'BBS_routes.csv', sep=''))
species=read.csv(paste(dataFolder, 'BBS_species.csv', sep=''))
weather=read.csv(paste(dataFolder, 'BBS_weather.csv', sep='')) %>%
  mutate(siteID=paste(countrynum, statenum, route,sep='-')) %>%
  dplyr::select(siteID,Year=year, RPID=rpid,runtype)
  

timeRange=c(1971:2014)
windowSizes=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)

#These are the offsets which scooch the analysis forward n years at a time to average out climatic variability
#0 means no offset. 
yearlyOffsets=0:10 #There is no check for setting this so high there is no room for the large window sizes.


#Some records are of genus only and "unidentified". Get rid of those.
#Or maybe find someway to incorperate them, because now some sites of false-negs. 
unidSpp=species %>%
  dplyr::select(Aou=AOU, name=english_common_name) %>%
  mutate(unID=ifelse(str_sub(name, 1,5)=='unid.', 1,0)  ) %>%
  filter(unID==1) %>%
  extract2('Aou')

#Filter weather to years of study so it can be used to calculate occData next.
weather=weather %>%
  filter(Year %in% timeRange)

#data frame for data just in the study years and without unidentified spp. 
#and only from runs where the BBS standards =1
occData=counts %>%
  filter(Year %in% timeRange) %>%
  mutate(siteID=paste(countrynum, statenum, Route,sep='-')) %>%
  filter(!Aou %in% unidSpp) %>%
  dplyr::select(Aou, siteID,Year, RPID) %>%
  left_join(weather, by=c('siteID','Year','RPID')) %>%
  filter(runtype==1) %>%
  dplyr::select(-runtype, RPID)


#A list of sites and the years they were sampled
siteList= occData %>%
  dplyr::select(Year, siteID) %>%
  distinct() %>%
  mutate(Year=as.factor(Year))

#bioclim values for all bbs routes from PRISM data.
source('get_prism_data.R')
bioclimData=get_bioclim_data()

#Some sites have na values.
bioclimData = bioclimData %>%
  filter(!is.na(bio1))


###################################################################
#Setup parallel processing
cl=makeCluster(numProcs)
registerDoParallel(cl)

#####################################################################
#Build list of training/validation year sets
###################################################################
# Using the defined range of years and the window sizes above, define
# training and testing years for all combinations of window sizes
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
#modelSetMatrix$setID=1:nrow(modelSetMatrix)
rm(numSets, windowIdentifier, thisSetDF, ones)

###################################################################
#Define site change types. These are used in analysis to look at accuracy
#of changed vs unchaged sites.
###################################################################
siteChanges=data.frame(changeType=c('absent-absent','present-present','absent-present','present-absent'),
                       T1_actual=as.factor(c(0,1,0,1)),
                       T2_actual=as.factor(c(0,1,1,0)))

#####################################################################
#Build site information for each setID
#######################################################################
#If your thought the last section was confusing, don't even bother
#reading this one. 
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

####################################################################
#Process model temporal validation accuracy from Rapacciuolo et al. 2014
processAccTV=function(results){

  #Calculate gain and loss curves
  Gains_curve=glm(gains ~ ns(deltaMweighted, df = 2), weights=rep(1, nrow(filter(results, changeType!='present-present'))), family=binomial,
                  data=filter(results, changeType!='present-present'))
  Losses_curve=glm(losses ~ ns(deltaMweighted, df = 2), weights=rep(1, nrow(filter(results, changeType!='absent-absent'))), family=binomial,
                  data=filter(results, changeType!='absent-absent'))
  
  #This fits the curve to the actual site data. ah la fig. 3a
  Simulated_gain <- predict(Gains_curve, newdata = results, se.fit = FALSE, type = "response")
  Simulated_loss <- predict(Losses_curve, newdata = results, se.fit = FALSE, type = "response")
  
  #the model y
  yModel=Simulated_gain-Simulated_loss
  
  #Pull out deltaMweighted for clarity
  deltaMweighted=results$deltaMweighted
  
  #Now the actual equation #2
  #In eq. 2 there is a yIdeal. That is a 1:1 line with deltaMweighted. So using deltaMweighted in place
  #of that is correct. 
  numerator = sum( abs(yModel - deltaMweighted) * deltaMweighted )
  divisor   = sum( deltaMweighted)
  #return(1 - (numerator/ divisor) )
  return( 1 - weighted.mean(abs(yModel - deltaMweighted), abs(deltaMweighted) ) )
}

###################################################################
#Process correct classification rate on changed or stable sites separtely
#( of sites that are absent-present, or present-absent, true-positives + true-negatives / total)
processCCR=function(results, type){
  if(type=='changed'){
    x=results %>%
      filter(changeType %in% c('absent-present','present-absent'))
  } else if(type=='stable'){
    x=results %>%
      filter(changeType %in% c('present-present','absent-absent'))
  }
  ccr=sum(with(x, T2_actual==T2_prediction)) / nrow(x)
}



#####################################################################
#Model to use. Details for each one are in bbsDMModels.R
modelsToUse=c('gbm')

#####################################################################
# Formula to pass to models
modelFormula=as.formula('presence ~ bio1+bio2+bio4+bio5+bio6+bio7+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio16+bio17+bio18+bio19')


#####################################################################
#Iterate thru spp

finalDF=foreach(thisSpp=unique(occData$Aou)[1:3], .combine=rbind, .packages=c('dplyr','tidyr','magrittr')) %do% {
#finalDF=foreach(thisSpp=unique(occData$Aou), .combine=rbind, .packages=c('dplyr','tidyr','magrittr')) %dopar% {
#for(thisSpp in unique(occData$Aou)[1:2]){
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
    
    
    #The raw results of all models/species/sites
    rawModelResults = modelResults %>%
      mutate(Aou=thisSpp, windowSize=thisWindowSize) 
    
    #####################
    #Now start calculation temporal validation curve accuracy
    #Start with the 1st parameters which can be calculated rowwise
    #Calculate deltaM and deltaMweighted
    modelResults = modelResults %>%
      mutate(deltaM = T2_prob-T1_prob) %>%
      mutate(deltaMweighted = ifelse( deltaM<0, deltaM/T1_prob, 
                                      ifelse(deltaM==0, 0, 
                                             deltaM/(1-T1_prob))))
    
    #Define stable and changed sites
    modelResults = modelResults %>%
      left_join(siteChanges, by=c('T1_actual','T2_actual'))
    
    #Isolate gains (absent-presnt) and losses (present-absent)
    modelResults = modelResults %>%
      mutate(gains= ifelse(changeType=='absent-present', 1,0), 
             losses= ifelse(changeType=='present-absent',1,0))
    
    #T2 binary prediction for ccr calculation
    modelResults = modelResults %>%
      mutate(T2_prediction= ifelse(T2_prob>0.5, 1, 0))
    
    accuracyTV=data.frame()
    #Now get the accuracy for each windowID and model
    #This requires fitting some glm's, so it needs to be inside for loops
    for(thisWindowID in unique(modelResults$windowID)){
      for(thisModelName in unique(modelResults$modelName)){
        
        #Processing of some results fail because of low variation in prediction scores, or low gains/losses. If that happens
        #mark it as -1
        thisAccuracy=try( processAccTV( filter(modelResults, windowID==thisWindowID, modelName==thisModelName)) , silent = TRUE)
        if(class(thisAccuracy)=='try-error'){ thisAccuracy=-1}
        
        thisCCRstable=processCCR(filter(modelResults, windowID==thisWindowID, modelName==thisModelName), type='stable' )
        thisCCRchanged=processCCR(filter(modelResults, windowID==thisWindowID, modelName==thisModelName), type='changed' )
        
        thisAccuracy=data.frame(accTV=thisAccuracy,
                                modelName=thisModelName,
                                thisWindowID=thisWindowID,
                                ccrStable=thisCCRstable,
                                ccrChanged=thisCCRchanged)
        
        accuracyTV = accuracyTV %>%
          bind_rows(thisAccuracy)
      }
    }

    #Include the window size and setID in the final accuracy results and add to this spp results
    accuracyTV$windowSize=thisWindowSize
    accuracyTV$setID=thisSetID
    accuracyTV$Aou=thisSpp

    #Add this window size iteration to results for this species.
    thisSppResults=bind_rows(thisSppResults, rawModelResults)
  } 
  #This gets returned to be added to the finalDF dataframe. 
  return(thisSppResults)
}

write.csv(finalDF,resultsFile,row.names = FALSE)
#write.csv(rawModelResults, rawResultsFile, row.names=FALSE)
