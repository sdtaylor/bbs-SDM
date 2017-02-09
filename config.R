#bioclim data folder. Used to store yearly rasters of bioclim vars made from prism data
bioclim_data_folder='~/data/yearly_bioclim/'

dataFolder='~/data/bbs/'

resultsFile=paste('./results/bbs_results_method1.csv',sep='')

training_years=1971:1980
testing_years=2000:2014

#Grain sizes in time and space
temporal_scales=c(1,3,5,10)

#Spatial scale to use. 
#in Deg lat/long
spatial_scales=c(0.1, 1.0, 2.0, 4.0)
#spatial_scales=c(0.1, 0.25, 0.5, 1.0, 2.0)
#in km
#spatial_scales=c(11,  28,   56,  111, 212)
#in number of 4km grid cells
#spatial_scales=c(3,    7,    14,  28,  53)


# Formula to pass to models
modelFormula=as.formula('presence ~ bio1+bio2+bio4+bio5+bio6+bio7+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio16+bio17+bio18+bio19')


options(prism.path = "~/data/prism")
