library(dplyr)
library(tidyr)
library(Metrics)
library(ggplot2)
setwd('~/projects/bbs')
source('multiplot.R')
library(stringr)


results=read.csv('~/data/bbs/bbsSDMResults.csv')
sppNames=read.csv('~/data/bbs/BBS_species.csv') %>%
  select(Aou=AOU, name=english_common_name)

results$T2_actual=as.numeric(as.character(results$T2_actual))

#calculate a binary (not probability) prediction
results = results %>%
  mutate(T2_prediction = ifelse(T2_prob>0.5, 1, 0))

changeCatagories=data.frame(changeType=c('absent-absent','present-present','absent-present','present-absent'),
                            changeCat=c('stable','stable','changed','changed'))

#get accuracy for each species in sites that changed and sites that did not change
speciesCCR=results %>%
  left_join(changeCatagories, by='changeType') %>%
  group_by(Aou, modelName, changeCat,featureSet) %>%
  summarize(ccr= auc(T2_actual, T2_prob), n=n() ) %>%
  ungroup()%>%
  group_by(Aou, changeCat,n) %>%
  summarize(meanAcc=mean(ccr, na.rm=TRUE), sdAcc=sd(ccr, na.rm=TRUE)) %>%
  left_join(sppNames, by='Aou') %>%
  filter(changeCat=='stable') %>%
  arrange(meanAcc)

#get accuracy for each species in sites that changed and sites that did not change

#Add in % site change between the two time periods
speciesCCR = speciesCCR %>%
  left_join(speciesRangeChange, by='Aou')

speciesCCR = speciesCCR %>%
  tidyr::spread(changeCat,ccr)
#######################################################
#Stacked plot like the Rapacciuoio paper
#Using a smoothing spline
rangeHistogram=ggplot(speciesCCR %>% dplyr::select(change) %>% distinct(), aes(change)) + geom_histogram()
stablePlot=ggplot(speciesCCR, aes(x=change, y=stable, colour=modelName, group=modelName))+stat_smooth(se=FALSE, method='gam', formula=y~ns(x,4))+ylim(0,1)+theme(legend.text=element_text(size=25)) 
changedPlot=ggplot(speciesCCR, aes(x=change, y=changed, colour=modelName, group=modelName))+stat_smooth(se=FALSE, method='gam', formula=y~ns(x,4))+ylim(0,1)+theme(legend.text=element_text(size=25)) 
multiplot(rangeHistogram, stablePlot, changedPlot, cols=1)

#without the smoother
rangeHistogram=ggplot(speciesCCR %>% dplyr::select(change) %>% distinct(), aes(change)) + geom_histogram()
stablePlot=ggplot(speciesCCR, aes(x=change, y=ccrStable, colour=modelName, group=modelName))+geom_point()
changedPlot=ggplot(speciesCCR, aes(x=change, y=ccrChanged, colour=modelName, group=modelName))+geom_point()
multiplot(rangeHistogram, stablePlot, changedPlot, cols=1)
##################################################
#Use ggplot to look at the same graphs but also point size reflect sample size
#ggplot(speciesCCR, aes(x=change, y=ccrChanged,size=changedN))+geom_point() +theme_bw()
#ggplot(speciesCCR, aes(x=change, y=ccrStable,size=stableN))+geom_point() +theme_bw()

######################################
#Get diversity
#siteRichness=results %>% 
#  group_by(siteID) %>%
#  summarize(actualRichness=sum(T2_actual), predictedRichness=sum(T2_prediction), probabilityRichness=sum(T2_prob))
