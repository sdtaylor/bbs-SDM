#database details are in a seperate file that is ignore by git so it's not published on github. 
source('databaseConfig.R')
database=src_postgres(dbname = dbName, host = dbHost, user = dbUser, password = dbPw)
rm(dbName, dbHost, dbUser, dbPw)

#Years in the data to consider.
timeRange=c(1971:2014)

#bioclim data folder. Used to store yearly rasters of bioclim vars made from prism data
bioclim_data_folder='~/data/yearly_bioclim/'

#Window sizes to use. This is the temporal scale
#windowSizes=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
temporal_scales=c(1,3,5,10)

training_years=1971:1980
testing_years=2000:2014

#Spatial scale to use. 
#in Deg lat/long
spatial_scales=c(0.1, 1.0, 2.0)
#spatial_scales=c(0.1, 0.25, 0.5, 1.0, 2.0)
#in km
#spatial_scales=c(11,  28,   56,  111, 212)
#in number of 4km grid cells
#spatial_scales=c(3,    7,    14,  28,  53)


#These are the offsets which scooch the analysis forward n years at a time to average out climatic variability
#0 means no offset. 
#There is no check for setting this so high there is no room for the large window sizes.
yearlyOffsets=0:5

#Model to use. Details for each one are in bbsDMModels.R
modelsToUse=c('gbm','naive')

# Formula to pass to models
modelFormula=as.formula('presence ~ bio1+bio2+bio4+bio5+bio6+bio7+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio16+bio17+bio18+bio19')
