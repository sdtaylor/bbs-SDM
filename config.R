
database=src_sqlite('./bbs_method2_data.sqlite', create = T)

resultsFile='./results/bbs_results_method2.csv'

bioclim_data_folder='~/data/yearly_bioclim/'

dataFolder='~/data/bbs/'

#Window sizes to use. This is the temporal scale
#windowSizes=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
temporal_scales=c(1,3,5,10)

training_years=1971:1980
testing_years=2000:2014

#Spatial scales to use. 
#in km on a side. 
#spatial_scales=c(40, 80, 160, 320)
#in Deg lat/long. Calculated as km*0.01 (very rough approximiation but suitable for this analysis)
spatial_scales=c(0.4, 0.8, 1.6, 3.2)

# Formula to pass to models
modelFormula=as.formula('presence ~ bio1+bio2+bio4+bio5+bio6+bio7+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio16+bio17+bio18+bio19')
