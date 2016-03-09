library(splines)

####################################################################
#The actual temporal validation plot method from Rapacciuolo et al. 2014
###################################################################
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

###################################################################
#Define site change types. These are used in analysis to look at accuracy
#of changed vs unchaged sites.
###################################################################
siteChanges=data.frame(changeType=c('absent-absent','present-present','absent-present','present-absent'),
                       T1_actual=as.factor(c(0,1,0,1)),
                       T2_actual=as.factor(c(0,1,1,0)))

###################################################################
#Process temporal validation plot accuracy. A little preproccing first, then processAccTV
#for each model
###################################################################
getAccuracyTV=function(modelResults){
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

  return(accuracyTV)  
}




