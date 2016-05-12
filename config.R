#database details are in a seperate file that is ignore by git so it's not published on github. 
source('databaseConfig.R')
database=src_postgres(dbname = dbName, host = dbHost, user = dbUser, password = dbPw)
rm(dbName, dbHost, dbUser, dbPw)

#Years in the data to consider.
timeRange=c(1971:2014)

#Window sizes to use. This is the temporal scale
#windowSizes=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
windowSizes=c(5,10)

training_years=1971:1980

#Spatial scale to use. These are width & height of grid cells in degrees. 
spatial_cell_sizes=c(0.1, 0.25, 0.5, 1.0, 2.0)

#These are the offsets which scooch the analysis forward n years at a time to average out climatic variability
#0 means no offset. 
#There is no check for setting this so high there is no room for the large window sizes.
yearlyOffsets=0:0

#Model to use. Details for each one are in bbsDMModels.R
modelsToUse=c('gbm')

# Formula to pass to models
modelFormula=as.formula('presence ~ bio1+bio2+bio4+bio5+bio6+bio7+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio16+bio17+bio18+bio19')
