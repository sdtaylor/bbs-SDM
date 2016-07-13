library(gbm)
library(rpart)
library(randomForest)
#library(neuralnet)
#library(party)
library(dplyr)

#Working so far. rf, gbm, glm, cta

#Function wrapper for many different SDM models. Takes the data frame of Presence,TimePeriod, predictor1, predictor2,....
#Trains on TimePeriod==T1 and returns prediction probabilities on TimePeriod==T2
#Note to self. Use as.formula() here instead of putting presence~. everywhere
sdmModels=function(data, modelName, modelFormula){
  
  #In the main script the presence/absence gets set to factors. Some models need it to be numeric though. 
  #if(modelName %in% c('gbm','ann','naive')) {
  #  data$presence=as.integer(as.character(data$presence))
  #}
  
  #Train on a single time period (always set 1), but test on all of them. including the training data. 
  #probability scores from training data are needed for temporal validation plots. 
  trainData=data %>% filter(windowID==1)
  testData=data
  
  if(modelName=='rf'){ #presence/absence needs to be factor
    model=randomForest(modelFormula, data= trainData)
    return(predict(model, type='prob', newdata = testData)[,2])
    
  } else if(modelName=='gbm'){ #presence/absence needs to be numeric
    model=gbm(modelFormula, n.trees=5000, distribution = 'gaussian', interaction.depth = 4, shrinkage=0.001, 
              data= trainData)
    perf=gbm.perf(model, plot.it=FALSE)
    x=predict(model, n.trees=perf, newdata=testData, type='response')
    return(x)
    
  } else if(modelName=='glm'){ #Would like also to put polynomials for all terms in here before the AIC model selection.
    model=glm(modelFormula, family='binomial', data=trainData)
    #model=step(model, direction='both') #In initial testing, AIC parameter selection shows very little improvement overall. 
    x=predict(model, type='response',newdata=testData)
    return(x)
    
  } else if(modelName=='ann'){ #prediction here not working yet
    model=neuralnet(modelFormula, hidden=7, learningrate=0.03, data=trainData)
    x=compute(model, testData)
    x=compute(model, as.matrix(dplyr::select(testData, -presence)))
    
  } else if(modelName=='naive'){
    #Mean probability of occurance over all training years for each cell.
    mean_presence_per_site=trainData %>%
      group_by(cellID) %>% 
      #summarize(predict=mean(as.integer(as.character(presence)))) %>%
      summarize(predict=mean(presence)) %>%
      ungroup()
    
    testData = testData %>%
      left_join(mean_presence_per_site, by='cellID')
    
    return(testData$predict)
    
  } else if(modelName=='cta'){
    model=rpart(modelFormula, method='class', data=trainData)
    return(predict(model, type='prob', newdata=testData)[,2])
    
  } else if(modelName=='gam'){
    
  } else if(modelName=='rfLogit') { #random forest of logistic regression trees
    ctrl=mob_control(verbose=FALSE, minsplit = 25, bonferroni=FALSE, alpha=0.001)
    model=party::mob(presence~bio1 + bio2+bio3 + bio4 + bio5 + bio6 + bio7| bio1+bio2+bio3+bio4+bio5+bio6+bio7, data=trainData, model=glinearModel, family=binomial(), control=ctrl)
    
  } else if(modelName=='mistnet') { #Harris 2015
    
  } else {
    stop(paste('Model',modelName,'not found',sep=' '))
  }
}


##########################
#Testing

#data$presence=as.integer(as.character(data$presence))

#data$presence=as.factor(data$presence)

#ctrl=mob_control(verbose=TRUE, minsplit = 25, bonferroni=FALSE, alpha=0.001)
#model=party::mob(presence~bio1 + bio2+bio3 + bio4 + bio5 + bio6 + bio7| bio1+bio2+bio3+bio4+bio5+bio6+bio7, data=trainData, model=glinearModel, family=binomial(), control=ctrl)
#model=glmtree(presence~bio1 + bio2+bio3 + bio4 + bio5 + bio6 + bio7| bio1+bio2+bio3+bio4+bio5+bio6+bio7, data=trainData, model=glinearModel, family=binomial())

#x=predict(model, newdata=testData)

##################
