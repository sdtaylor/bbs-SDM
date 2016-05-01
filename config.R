#database details are in a seperate file that is ignore by git so it's not published on github. 
source('databaseConfig.R')
database=src_postgres(dbname = dbName, host = dbHost, user = dbUser, password = dbPw)
rm(dbName, dbHost, dbUser, dbPw)


spatial_cell_sizes=c(0.1, 0.25, 0.5, 1.0, 2.0)
