#Essentially a scratch pad for me to run various anallysis and graphs inside rstudio. 

library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(broom)
library(Metrics)
source('multiplot.R')
results=read.csv('results/bbsSDMResults_ManyTimePeriods.csv') %>%
  filter(modelName=='rf')
  
sppNames=read.csv('~/data/bbs/BBS_species.csv') %>%
  rename(name=english_common_name, Aou=AOU) %>%
  select(name, Aou)


with(results %>% filter(windowSize==2), plot(accTV~thisWindowID))

for(thisSp in unique(results$Aou)[1:10]){
  with(results %>% filter(Aou==thisSp), plot(accTV~thisWindowID, main=thisSp, ylim=c(0,1)))
  
}

windowSizeOneResults=results %>% filter(windowSize==1)

slopes = results %>%
  group_by(Aou, windowSize) %>%
  do(model=lm(accTV~thisWindowID, data=.)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(slope=unlist(model)[2]$'coefficients.thisWindowID') %>%
  select(-model)

for(thisSp in unique(slopes$Aou)[1:100]){
  with(slopes %>% filter(Aou==thisSp), plot(slope~windowSize, main=thisSp, ylim=c(-0.6,0.5)))
  abline(0,0)
}


for(thisWindowSize in unique(results$windowSize)){
  for(thisSp in unique(resu))
}

results= results %>%
  mutate(year=(thisWindowID * windowSize) - (windowSize*0.5))

######################3
#Scale each species results around 0
results$accTVscaled=0
for(thisSp in unique(results$Aou)){
  results$accTVscaled[results$Aou==thisSp]=scale(results$accTV[results$Aou==thisSp], scale=FALSE)
}


#A single graph for each species modeling decay rates of different window sizes
for(thisSp in unique(results$Aou)[50:70]){
  speciesName=sppNames %>% filter(Aou==thisSp) %>% extract2('name')
  thisPlot=ggplot(filter(results, Aou==thisSp, windowSize %in% c(3,6,9)), aes(x=year, y=ccrChanged, group=windowSize, colour=as.factor(windowSize))) + 
    geom_point(size=4) + geom_line(size=2) + 
    #stat_smooth(se=FALSE) +
    ylim(0,1) +
    ggtitle(speciesName)+
    theme_grey(base_size = 25)
  print(thisPlot)
}

##########################################################
#Repeat original 2 tiem period analysis


results=read.csv('results/bbsSDMResults_ManyTimePeriods.csv') %>%
  filter(windowSize==5, modelName=='rf')

ggplot(results, aes(as.factor(thisWindowID), accTVscaled)) + geom_boxplot()

thisPlot=ggplot(results, aes(x=year, y=accTVscaled, group=Aou, colour=as.factor(Aou))) + 
  geom_point() + geom_line()
  #stat_smooth(se=FALSE) +
  ylim(-0.5,0.5)
print(thisPlot)


#########################################################3
##########################################################
#Raw results analysis
#Create probability vs probability graphs

siteChanges=data.frame(changeType=c('absent-absent','present-present','absent-present','present-absent'),
                       T1_actual=c(0,1,0,1),
                       T2_actual=c(0,1,1,0))

results=read.csv('results/bbsSDMResults_ManyTimePeriods_Raw.csv') %>%
  left_join(siteChanges, by=c('T1_actual','T2_actual'))


results= results %>%
  mutate(yearLag=(windowID * windowSize) - (windowSize*0.5))

#Create prop-prop graphs one species at a time. 

resultsCompactAou=results %>%
  mutate(T2_prob=round(T2_prob, digits=1)) %>%
  filter(!changeType %in% c('absent-present','present-absent')) %>%
  group_by(windowSize, yearLag, T2_prob, Aou) %>%
  summarize(actualPct=sum(T2_actual)/n(), n=n())

for(thisSp in unique(resultsCompactAou$Aou)[1:90]){
  thisPlot=ggplot(filter(resultsCompactAou, windowSize==3, Aou==thisSp), aes(x=T2_prob, y=actualPct, colour=as.factor(yearLag), group=as.factor(yearLag)))+
    geom_smooth(se=FALSE)+
    geom_point(aes(size=n))+
    geom_abline(intercept=0, slope=1) +
    ylim(0,1)+xlim(0,1)+
    ggtitle(paste(thisSp,sep=' '))
  print(thisPlot)
}


#Create prob-prob graphs with all species and sites together
resultsCompact=results %>%
  mutate(T2_prob=round(T2_prob, digits=2)) %>%
  #filter(changeType %in% c('absent-present','present-absent')) %>%
  filter(changeType %in% c('present-present','absent-absent')) %>%
  group_by(windowSize, yearLag, T2_prob) %>%
  summarize(actualPct=sum(T2_actual)/n(), n=n())

for(thisWindowSize in unique(resultsCompact$windowSize)){
  
  thisPlot=ggplot(filter(resultsCompact, windowSize==thisWindowSize), aes(x=T2_prob, y=actualPct, colour=as.factor(yearLag), group=as.factor(yearLag)))+
    #geom_smooth(se=FALSE)+
    geom_point(size=3)+
    geom_abline(intercept=0, slope=1) +
    ylim(0,1)+xlim(0,1)+
    ggtitle(paste(thisWindowSize,'Year Window Size',sep=' ')) +
    theme_grey(base_size = 25)
  print(thisPlot)

}

oneOff=resultsCompact %>%
  filter(windowSize==3, yearLag %in% c(7.5,16.5,31.5))
thisPlot=ggplot(oneOff, aes(x=T2_prob, y=actualPct, colour=as.factor(yearLag), group=as.factor(yearLag)))+
  geom_smooth(se=FALSE)+
  geom_point()+
  geom_abline(intercept=0, slope=1) +
  ylim(0,1)+xlim(0,1)+
  ggtitle(paste('3 Year Window Size',sep=' '))+
  ylab('Percentage across all sites/species of actual occurances')+
  xlab('Probability of occurance across all sites/species')+
  labs(linetype='Time Lag in Years')
print(thisPlot)

#############################################################
#cohens kappa
kappaBinary=function(df){
  
}


#########################################################3
##########################################################
#Raw results analysis
#per spp time lags with proportion of changed sites in there as well. 

siteChanges=data.frame(changeType=c('absent-absent','present-present','absent-present','present-absent'),
                       T1_actual=c(0,1,0,1),
                       T2_actual=c(0,1,1,0))

results=read.csv('results/bbsSDMResults_ManyTimePeriods_WindowSizes5s-14_Raw.csv') %>%
  mutate(yearLag=(windowID * windowSize) - (windowSize*0.5)) %>%
  left_join(siteChanges, by=c('T1_actual','T2_actual')) %>%
  mutate(changed=ifelse(changeType %in% c('absent-present','present-absent'), 'yes','no'))

migrantStatus=read.csv('~/data/bbs/bbsMigrationStatus.csv')

#Shawns metric. 0 is a perfect prediction for 2nd time period. -1 is a perfect false-positive, +1 is a perfect false negative
results=results %>%
  mutate(shawnScore=T2_actual - T2_prob)

#For each spp/windowSize/timelag: 
#rangeChange: number of occurances in T1 / number of occurance in T2
#siteChangePercent: percent of all surveyed sites that were absent-present or present-absent between T1 and T2
#T2sites: total occurances in T2
#T1sites: total occurances in T1
#sitesSurveyed: total sites surveyed in T2
rangeChanges=results %>%
  group_by(yearLag, windowSize, windowID, Aou) %>%
  summarize(rangeChange=sum(T1_actual)/sum(T2_actual), siteChangePercent=sum(changed=='yes')/n(), 
            T2sites=sum(T2_actual), T1sites=sum(T1_actual), sitesSurveyed=n(),
            numChangedSites=sum(changed=='yes'),numSitesDropped=sum(changeType=='present-absent'), numSitesAdded=sum(changeType=='absent-present')) %>%
  ungroup()

resultsCompact=results %>%
  mutate(T2_prediction= ifelse(T2_prob>0.5, 1, 0), probDiff=T2_prob-T1_prob) %>%
  group_by(yearLag, windowSize, windowID, Aou, changeType) %>%
  summarize(ccr=sum(T2_actual==T2_prediction)/n(), auc=Metrics::auc(T2_actual, T2_prob), probDiffMean=mean(probDiff), probDiffsd=sd(probDiff)) %>%
  ungroup()

resultsCompact = resultsCompact %>%
  left_join(rangeChanges, by=c('Aou', 'windowSize', 'windowID', 'yearLag'))


#get slope of observed change shift. 
slopes = resultsCompact %>%
  group_by(Aou, windowSize) %>%
  do(slope=tidy(lm(siteChangePercent~yearLag, data=.))$estimate[2]) %>%
  filter(slope!='NaN') %>%
  ungroup() %>%
  mutate(slope=as.numeric(slope))

slopes = resultsCompact %>%
  group_by(Aou, windowSize) %>%
  mutate(pValue=glance(lm(siteChangePercent~yearLag, data=.)$p.value[1]))  %>%
  filter(pValue!='NaN') %>%
  ungroup() %>%
  mutate(pValue=as.numeric(pValue), slope=as.numeric(pValue))

#Birds in watling et al 2013
watlingBirds=c(3100, 3090, 3050, 2950, 4060, 4090, 4900, 7360, 7180, 7190, 5750, 5000)

sppToGraph = slopes %>%
  filter(slope>0.001 & windowSize==7) %>%
  select(Aou) %>%
  distinct() %>%
  #left_join(migrantStatus, by='Aou') %>%
  #filter(migrant=='nomadic') %>%
  extract2('Aou')

for(thisSp in sppToGraph){
#for(thisSp in unique(resultsCompact$Aou)[100:150]){
  speciesName=sppNames %>% filter(Aou==thisSp) %>% extract2('name')
  thisPlot=ggplot(filter(resultsCompact, Aou==thisSp, windowSize %in% c(5,9,14)), aes(x=yearLag)) + 
    #CCR lines
    geom_point(size=4, aes(y=auc, group=windowSize, colour=as.factor(windowSize))) +
    geom_line(size=2, aes(y=auc, group=windowSize, colour=as.factor(windowSize))) + 
    #Proportion lines
    geom_point(size=4, aes(y=siteChangePercent, group=windowSize, colour=as.factor(windowSize))) +
    geom_line(size=2,linetype=2, aes(y=siteChangePercent, group=windowSize, colour=as.factor(windowSize))) + 
    #prob diff lines
    #geom_point(size=4, aes(y=probDiffMean, group=windowSize, colour=as.factor(windowSize))) +
    #geom_line(size=2,linetype=1, aes(y=probDiffMean, group=windowSize, colour=as.factor(windowSize))) + 
    #stat_smooth(se=FALSE) +
    ylim(-0.5,1) +
    ggtitle(speciesName)+
    theme_grey(base_size = 25)
  print(thisPlot)
}

####################################################
ggplot(filter(resultsCompact, yearLag %in% c(35.0, 37.5), changeType=='absent-present')) +
  geom_point(aes(y=numSitesAddedProportion, x=probDiffMean, group=windowSize, colour=as.factor(windowSize))) +
  #ylim(0,1) + xlim(-0.75, 0.75) + ggtitle('Number of sites added per species after 38 years') + 
  xlab('Mean difference in probability for those sites (Probability T2 - Probability T1)')

ggplot(filter(resultsCompact, yearLag %in% c(21.0, 22.5), changeType=='absent-present', windowSize!=9)) +
  geom_point(aes(y=numSitesAddedProportion, x=probDiffMean, group=windowSize, colour=as.factor(windowSize)))+
  ylim(0,1) + xlim(-0.75, 0.75) + ggtitle('Number of sites added per species after 21 years') + 
  xlab('Mean difference in probability for those sites (Probability T2 - Probability T1)')

ggplot(filter(resultsCompact, yearLag %in% c(35.0, 37.5), changeType=='present-absent')) +
  geom_point(aes(y=numSitesDropped, x=probDiffMean, group=windowSize, colour=as.factor(windowSize))) +
  ylim(0,500) + xlim(-0.75, 0.75) + ggtitle('Number of sites dropped per species after 38 years') + 
  xlab('Mean difference in probability for those sites (Probability T2 - Probability T1) ')

ggplot(filter(resultsCompact, yearLag %in% c(21.0, 22.5), changeType=='present-absent', windowSize!=9)) +
  geom_point(aes(y=numSitesDropped, x=probDiffMean, group=windowSize, colour=as.factor(windowSize)))+
  ylim(0,500) + xlim(-0.75, 0.75) + ggtitle('Number of sites dropped per species after 21 years') + 
  xlab('Mean difference in probability for those sites (Probability T2 - Probability T1)')


###################################################
#Trajectory of distribution (site) growth vs trajectory of forecast accuracy

forecastTrajectory=resultsCompact %>%
  filter(siteChangePercent>0.1, windowSize %in% c(3,5,7,9), changed=='yes') %>%
  group_by(Aou,windowSize) %>%
  do(aucModel=lm(auc~yearLag, data=.), percentChangeModel=lm(siteChangePercent~yearLag, data=.)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(aucSlope=unlist(aucModel)[2]$'coefficients.yearLag', changeSlope=unlist(percentChangeModel)[2]$'coefficients.yearLag')


####################################################
#Just histograms
siteChanges=data.frame(changeType=c('absent-absent','present-present','absent-present','present-absent'),
                       T1_actual=c(0,1,0,1),
                       T2_actual=c(0,1,1,0))


results=read.csv('results/bbsSDMResults_ManyTimePeriods_WindowSizes5-14_Raw.csv') %>%
  mutate(yearLag=(windowID * windowSize) - (windowSize*0.5)) %>%
  left_join(siteChanges, by=c('T1_actual','T2_actual')) %>%
  mutate(changed=ifelse(changeType %in% c('absent-present','present-absent'), 'yes','no'))

#For each spp/windowSize/timelag: 
#rangeChange: number of occurances in T1 / number of occurance in T2
#siteChangePercent: percent of all surveyed sites that were absent-present or present-absent between T1 and T2
#T2sites: total occurances in T2
#T1sites: total occurances in T1
#sitesSurveyed: total sites surveyed in T2
rangeChanges=results %>%
  filter(windowSize==5) %>%
  group_by(yearLag, windowSize, windowID, Aou) %>%
  summarize(rangeChange=sum(T1_actual)/sum(T2_actual), siteChangePercent=sum(changed=='yes')/n(), 
            T2sites=sum(T2_actual), T1sites=sum(T1_actual), sitesSurveyed=n(),
            numChangedSites=sum(changed=='yes'),numSitesDropped=sum(changeType=='present-absent'), numSitesAdded=sum(changeType=='absent-present')) %>%
  ungroup()

speciesRangeBins= rangeChanges %>%
  filter(windowSize==5, yearLag==12.5) %>%
  mutate(sppRangeBin= ifelse(T1sites<=50, 'low',
                             ifelse(T1sites<=600, 'med', 'high')))%>%
  dplyr::select(Aou, T1sites, sppRangeBin)

resultsCompact=results %>%
  filter(windowSize==5, yearLag == 12.5) %>% #JUST 5!!!!!!!
  mutate(T2_prediction= ifelse(T2_prob>0.5, 1, 0), probDiff=T2_prob-T1_prob) %>%
  group_by(yearLag, windowSize, windowID, Aou, changeType) %>%
  summarize(ccr=sum(T2_actual==T2_prediction)/n(), auc=Metrics::auc(T2_actual, T2_prob), probDiffMean=mean(probDiff), probDiffsd=sd(probDiff)) %>%
  ungroup() %>%
  left_join(speciesRangeBins, by='Aou')
  




thisPlotData=filter(resultsCompact, changeType=='present-absent', yearLagBin == 12.5)
ggplot(thisPlotData, aes(x=probDiffMean)) +


ggplot(thisPlotData) +
  geom_freqpoly(aes(x=probDiffMean, colour=as.factor(windowSize)))



for(thisRangeSize in unique(resultsCompact$sppRangeBin)){
  for(thisChangeType in unique(resultsCompact$changeType)) {
    thisPlotTitle=paste(thisChangeType, thisRangeSize, '1971-1975 occurance', sep=' ')
    
    thisPlot=ggplot(filter(resultsCompact, sppRangeBin==thisRangeSize, changeType == thisChangeType))+
      geom_histogram(aes(x=probDiffMean), bins=50)+
      ggtitle(thisPlotTitle) +
      #scale_y_continuous(limit=c(0,100)) +
      ylim(c(0,90))+
      xlim(c(-1,1)) +
      geom_vline(xintercept = 0, size=1)
      theme_bw()
    
    print(thisPlot)
    
  }
}

######################################################################
#Extreme deltas
siteChanges=data.frame(changeType=c('absent-absent','present-present','absent-present','present-absent'),
                       T1_actual=c(0,1,0,1),
                       T2_actual=c(0,1,1,0))


results=read.csv('results/bbsSDMResults_ManyTimePeriods_WindowSizes5-14_Raw.csv') %>%
  mutate(yearLag=(windowID * windowSize) - (windowSize*0.5), deltaProb=T2_prob-T1_prob) %>%
  left_join(siteChanges, by=c('T1_actual','T2_actual')) %>%
  mutate(changed=ifelse(changeType %in% c('absent-present','present-absent'), 'yes','no'))

#For each spp/windowSize/timelag: 
#rangeChange: number of occurances in T1 / number of occurance in T2
#siteChangePercent: percent of all surveyed sites that were absent-present or present-absent between T1 and T2
#T2sites: total occurances in T2
#T1sites: total occurances in T1
#sitesSurveyed: total sites surveyed in T2
rangeChanges=results %>%
  filter(windowSize==5) %>%
  group_by(yearLag, windowSize, windowID, Aou) %>%
  summarize(rangeChange=sum(T1_actual)/sum(T2_actual), siteChangePercent=sum(changed=='yes')/n(), 
            T2sites=sum(T2_actual), T1sites=sum(T1_actual), sitesSurveyed=n(),
            numChangedSites=sum(changed=='yes'),numSitesDropped=sum(changeType=='present-absent'), numSitesAdded=sum(changeType=='absent-present')) %>%
  ungroup()

resultsCompact=results %>%
  filter(windowSize==5, yearLag == 12.5) %>% #JUST 5!!!!!!!
  mutate(T2_prediction= ifelse(T2_prob>0.5, 1, 0), probDiff=T2_prob-T1_prob) %>%
  group_by(yearLag, windowSize, windowID, Aou, changeType) %>%
  summarize(ccr=sum(T2_actual==T2_prediction)/n(), auc=Metrics::auc(T2_actual, T2_prob), probDiffMean=mean(probDiff), probDiffsd=sd(probDiff)) %>%
  ungroup() %>%
  left_join(speciesRangeBins, by='Aou')



extremeDeltas=results %>%
  filter(windowSize==5) %>%
  filter((T1_prob <=0.2 & T2_prob >= 0.8) | (T1_prob>=0.8 & T2_prob<=0.2)) %>%
  mutate(highProbOf=ifelse(T2_prob>=0.8, 'absent-present','present-absent'))
  
#Sites that have a high probability of gaining or losing species
x=extremeDeltas %>%
  filter(yearLag==22.5) %>%
  group_by(siteID) %>%
  summarize(numSpp=n_distinct(Aou))

  
ggplot(extremeDeltas, aes(highProbOf, fill=changeType)) + geom_bar() + facet_grid(. ~ yearLag)

hist(extremeDeltas$deltaProb)

for(thisLag in unique(extremeDeltas$yearLag)){
  thisPlotTitle=paste(thisLag, 'year lag', sep=' ')
  thisPlot=ggplot(filter(extremeDeltas, yearLag==thisLag), aes(highProbOf, fill=changeType)) + 
    geom_bar() +
    ggtitle(thisPlotTitle)
  print(thisPlot)
  
}

#############################################################
#Ethans graph

resultsCompact=results %>%
  filter(windowSize==5) %>% #JUST 5!!!!!!!
  group_by(yearLag, windowID, Aou, changeType) %>%
  summarize(ccr=sum(T2_actual==T2_prediction)/n(), auc=Metrics::auc(T2_actual, T2_prob), probDiffMean=mean(probDiff), probDiffsd=sd(probDiff)) %>%
  ungroup() %>%
  left_join(speciesRangeBins, by='Aou')





