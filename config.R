#bioclim data folder. Used to store yearly rasters of bioclim vars made from prism data
bioclim_data_folder='~/data/yearly_bioclim/'

dataFolder='~/data/bbs/'

resultsFile=paste('./results/bbs_results_method1.csv',sep='')

training_years=1971:1980
testing_years=2001:2010

#Mininum years in either the training or testing period to keep a site
minimum_years = 5

#Grain sizes in time and space
temporal_scales=c(1,3,5,10)

#Spatial scales to use. 
#in km on a side. 
#spatial_scales=c(40, 80, 160, 320)
#in Deg lat/long. Calculated as km*0.01 (very rough approximiation but suitable for this analysis)
spatial_scales=c(0.4, 0.8, 1.6, 3.2)
spatial_scales_km=c(40, 80, 160, 320)

# Formula to pass to models
modelFormula=as.formula('presence ~ bio1+bio2+bio4+bio5+bio6+bio7+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio16+bio17+bio18+bio19')


options(prism.path = "~/data/prism")
