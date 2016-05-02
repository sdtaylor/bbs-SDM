# Download and install PRISM data into PostgreSQL and values for BBS routes

library(prism)
library(raster)
library(dplyr)
library(tidyr)
library(stringr)
library(rgeos)
library(magrittr)
options(prism.path = "~/data/prism")

#Years of prism data to 
years_to_use=1966:2014
#Load the DB. create = TRUE does not seem 


source('config.R')

#######################################################
#Create a continuous grid of a certain cell size across the extent of the template raster
#######################################################
create_grid=function(cellsize, templateRaster){
  minX=bbox(templateRaster)[1,1]
  maxX=bbox(templateRaster)[1,2]
  minY=bbox(templateRaster)[2,1]
  maxY=bbox(templateRaster)[2,2]
  
  lowerX=maxX-minX
  lowerY=maxY-minY
  
  numCellsX=round((lowerX+1)/cellsize)
  numCellsY=round((lowerY+1)/cellsize)
  g=GridTopology(c(minX, minY), c(cellsize,cellsize), c(numCellsX, numCellsY))
  g=SpatialGrid(g, proj4string = CRS(crs(templateRaster)@projargs))
}

#######################################################
#Subset a map of grid cells to those cells that have at least
#one point in them. 
#bbs median number of sites in each cell size:
#cell size:        0.05, 0.1, 0.25, 0.5, 1.0, 2.0
#median num sites: 1   ,  1,  1   ,  1 ,  4 ,  16
#######################################################
subset_grid=function(g, sp){
  g=as(g, 'SpatialPolygons')
  cells_to_keep = data.frame(grid_id=over(sp, g)) %>%
    filter(!is.na(grid_id)) %>%
    distinct() %>%
    extract2('grid_id')
  g=as(g, 'SpatialPolygonsDataFrame')
  g=g[cells_to_keep,]
  return(g)
}

#for(size in c(0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0)){
#  prism_grid=create_grid(size, p)
#  prism_grid=subset_grid(prism_grid, locations)
#  plot(p)
#  plot(prism_grid, add=T)

#}

#######################################################
#Takes a grid cell and site spatial df.
#Returns a df of colnames('cellID','cellSize','siteID')
#######################################################
assign_sites_to_grid=function(g, sites, cellSize){
  g=as(g, 'SpatialPolygons')
  x=data.frame(cellID=over(sites, g), siteID=sites@data$siteID) %>%
    filter(!is.na(cellID))
  x$cellSize=cellSize
  return(x)
}

#######################################################
#Much quicker method for extracting data from a rasterstack w/ polygons. 
#######################################################
extract_polygon=function(raster_stack, sp, radius=NULL){
  if(class(raster_stack)!='RasterStack'){stop('This function only provides value when using a raster stack. For a single raster it is no faster than raster::extract()')}
  #Need check to make sure CRS is teh same
  ######
  #If sp == spatialPointsDataFrame, convert it to a circular polygon w/ size=radius
  ######
  sp=as(sp, 'SpatialPolygons')
  cellNumbers=cellFromPolygon(raster_stack, sp)
  
  idVec=c()
  for(i in seq_along(cellNumbers)){
    idVec=c(idVec,rep(i, length(cellNumbers[[i]])))
  }
  cellNumbers=unlist(cellNumbers)
  cellData=extract(raster_stack, cellNumbers)
  df=data.frame(id=idVec, cellData)
}

#######################################################
#Downloads the raw prism rasters into the folder specified above. 
#######################################################
download_prism=function(){
  months <- c(1:12)
  clim_vars <- c("ppt", "tmin", "tmean", "tmax")
  for (month in months){
    for (clim_var in clim_vars){
      get_prism_monthlys(type=clim_var, year = years_to_use, mon = month, keepZip=F)
    }
  }
}


########################################################
#Takes output of ls_prism_data() and makes sure all files
#within a year range are there. ie. for all 12 months and 
#all 4 variables. 
#########################################################
check_if_prism_files_present=function(prism_ls, years){
  #Extract out all the variable names and dates
  prism_ls = prism_ls %>%
    rowwise() %>%
    mutate(var=strsplit(files, '_')[[1]][2], yearMonth=strsplit(files, '_')[[1]][5]) %>%
    mutate(year=as.integer(substr(yearMonth,1,4)), month=as.integer(substr(yearMonth, 5,6))) %>%
    dplyr::select(-files, -yearMonth) %>%
    filter(year %in% years) %>%
    arrange(var, year, month)
  
  #Setup a list of what should be there
  to_check=expand.grid(var=c('ppt','tmax','tmean','tmin'), year=years, month=1:12, stringsAsFactors = FALSE) %>%
    arrange(var,year,month)
  
  #If these two data frames are equal, then all prism raster data in the years specified is present.
  #all.equal returns TRUE if they are equal, and a description of the discrepency if they are not equal,
  #hence the isTRUE wrapper. 
  return(isTRUE(base::all.equal(to_check, prism_ls, check.attributes=FALSE)))
}

########################################################
#Return prism data for each site from the sqlite DB. If it's not loaded in
#the DB, then extract values from the raw prism rasters, store in DB, and return them.
########################################################
get_prism_data=function(){
  #Query sqlite database for the prism_bbs_data table. If it exists, return it.
  #Otherwise create it from the raw prism data.
  if('prism_bbs_data' %in% src_tbls(database)){
    return(collect(tbl(database, sql('SELECT * from prism_bbs_data'))))
  } else { 
    
    print('PRISM data table not found, processing raw data')
    
    #Load the bbs data locations and convert them to a spatial object.
    #Stop here if bbs data isn't available. Could also make this query the DB as well.
    bbs_routes <- try(read.csv("~/data/bbs/BBS_routes.csv"))
    if(class(bbs_routes)=='try-error'){stop("Can't load bbs_data.csv inside get_prism_data()")}
    locations = bbs_routes %>%
      mutate(siteID=paste(countrynum, statenum, route,sep='-')) %>%
      dplyr::select(siteID, long=loni, lat=lati)
    #locations=locations[3818:4010,] #just use a few points when testing
    coordinates(locations) <- c("long", "lat")
    crs(locations)=CRS('+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0')
    
    #Check to see if all the raw data in the years specified are downloaded,
    #download everything again if not. (Getting only what is needed, say if a 
    #previous download failed, might be overly complicated)
    if(!check_if_prism_files_present(ls_prism_data(), years_to_use)){
      download_prism()
    }
    
    #Load the prism data
    prism_stacked <- prism_stack(ls_prism_data())
    
    #Build the different grids at each spatial scale
    #And list of which site is within each cell at each scale
    spatial_grids=list()
    spatial_grid_info=data.frame()
    for(cell_size in spatial_cell_sizes){
      this_grid=create_grid(cell_size, prism_stacked)
      this_grid=subset_grid(this_grid, locations)
      spatial_grids[paste('size',cell_size, sep='-')]=this_grid
      
      spatial_grid_info = spatial_grid_info %>%
        bind_rows(assign_sites_to_grid(this_grid, locations, cell_size))
    }
    

    #Extract the prism values for each cell size
    prism_bbs_data=data.frame()
    for(cell_size in spatial_cell_sizes){
      extracted <- extract_polygon(prism_stacked, spatial_grids[[paste('size',cell_size,sep='-')]]) %>%
        group_by(id) %>%
        summarize_each(funs(mean(.,na.rm=T))) %>%
        rename(cellID=id)

      extracted$cellSize=cell_size
    
    #prism_bbs_data <- data.frame(siteID = locations$siteID, coordinates(locations), extracted)
    prism_bbs_data <- extracted %>%
      gather(date, value, -cellID, -cellSize) %>%
      tidyr::extract(date, c("clim_var", "year", "month"),
                     "PRISM_([:alpha:]*)_stable_[:alnum:]*_([:digit:]{4})([:digit:]{2})_") %>%
      bind_rows(prism_bbs_data)
    }
    
    #Format the data a little and load into the sqlite database.
    prism_bbs_data$year <- as.numeric(prism_bbs_data$year)
    prism_bbs_data$month <- as.numeric(prism_bbs_data$month)
    
    #Spread out the climate variables ppt, tmean, etc into columns
    prism_bbs_data = prism_bbs_data %>%
      spread(clim_var, value)
    
    mydata <- copy_to(database, prism_bbs_data, temporary = FALSE,
                      indexes = list(c("cellID", "year", "month")))
    
    #Now return the data as asked for
    return(prism_bbs_data)
    
  }
}

###################################################################
#Helper function to calculate some of the bioclim variables,
#like "precip in coldest month"
###################################################################
max_min_combo=function(vec1,vec2,max=TRUE){
  #Return the value in vec1 in the position where
  #vec2 is either highest or lowest. But 1st check for na 
  #values. 
  if(any(is.na(vec1)) | any(is.na(vec2))){
    return(NA)
  } else  if(max){
    return(vec1[which.max(vec2)])
  } else {
    return(vec1[which.min(vec2)])
  }
}

#################################################################
#From raw prism monthly values calculate all the bioclim variables.
#You should not call this directly to load bioclim vars. Instead call
#get_bioclim_data(), which will 1st try to load the data from the sqlite
#db before processing it all from scratch. 
###################################################################
process_bioclim_data=function(){
  #Get the prism data. 
  prism_bbs_data=get_prism_data()
  

  
  #Process the quarter ones first.
  quarter_info=data.frame(month=1:12, quarter=c(1,1,1,2,2,2,3,3,3,4,4,4))
  bioclim_quarter_data= prism_bbs_data %>%
    left_join(quarter_info, by='month') %>%
    group_by(siteID, year, quarter) %>%
    summarize(precip=sum(ppt), temp=mean(tmean)) %>%
    ungroup() %>%
    group_by(siteID,year) %>%
    summarize(bio8=max_min_combo(temp, precip, max=TRUE),
              bio9=max_min_combo(temp, precip, max=FALSE),
              bio10=max(temp),
              bio11=min(temp),
              bio16=max(precip),
              bio17=min(precip),
              bio18=max_min_combo(precip, temp, max=TRUE),
              bio19=max_min_combo(precip, temp, max=FALSE)) %>%
    ungroup()
  
  #Next the yearly ones, joining the quartely ones  back in at the end. 
  bioclim_data=prism_bbs_data %>%
    group_by(siteID, year) %>%
    mutate(monthly_temp_diff=tmax-tmin) %>%
    summarize(bio1=mean(tmean),
              bio2=mean(monthly_temp_diff),
              bio4=sd(tmean)*100,
              bio5=max_min_combo(tmax,tmean,max=TRUE),
              bio6=max_min_combo(tmin,tmean,max=FALSE),
              bio12=sum(ppt),
              bio13=max(ppt),
              bio14=min(ppt),
              bio15=cv(ppt)) %>%
    ungroup() %>%
    mutate(bio7=bio5-bio6,
           bio3=(bio2/bio7)*100) %>%
    full_join(bioclim_quarter_data, by=c('siteID','year'))
  
  return(bioclim_data)
}


####################################################################
#Load the bioclim variables from sqlite. if they aren't available,
#load the prism data, process bioclim, load it in sqlite, then return 
#bioclim as requested.
#
#TODO: Possibly put a check here to make sure all years are accounted for?
####################################################################

get_bioclim_data=function(){
  #Query sqlite database for the bioclim_bbs_data table. If it exists, return it.
  #Otherwise create it from the raw prism data.
  if('bioclim_bbs_data' %in% src_tbls(database)){
    return(collect(tbl(database, sql('SELECT * from bioclim_bbs_data'))))
  } else { 
    print("bioclim data table not found, processing from scratch. ")
    bioclim_bbs_data=process_bioclim_data()
    
    copy_to(database, bioclim_bbs_data, temporary = FALSE, 
              indexes = list(c('cellID','year')))
     
    return(bioclim_bbs_data)
    
  } 
}