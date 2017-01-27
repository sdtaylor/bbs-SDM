# Extract bioclim info for all years and all spatial scales, and store in
# the database. This upscales the bioclim data as per bbs method2.
# Also assign all the bbs sites to grid cells at each spatial scale
# and store that info in the database

library(raster)
library(dplyr)
library(tidyr)
library(stringr)
library(rgeos)
library(magrittr)

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

#######################################################
#Takes a grid cell and site spatial df.
#Returns a df of colnames('spatial_cell_id','spatial_scale','siteID')
#######################################################
assign_sites_to_grid=function(g, sites, spatial_scale){
  g=as(g, 'SpatialPolygons')
  x=data.frame(spatial_cell_id=over(sites, g), siteID=sites@data$siteID) %>%
    filter(!is.na(spatial_cell_id))
  x$spatial_scale=spatial_scale
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

########################################################
#Return prism data for each site from the sqlite DB. If it's not loaded in
#the DB, then extract values from the raw prism rasters, store in DB, and return them.
########################################################
process_bioclim_data=function(){

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
    
    all_bioclim_files = list.files(bioclim_data_folder, full.names = TRUE)
    all_bioclim_files = all_bioclim_files[!grepl('xml', all_bioclim_files)]
    bioclim_stack = raster::stack(all_bioclim_files[1:5])
    
    locations = spTransform(locations, crs(bioclim_stack))
    
    #Build the different grids at each spatial scale
    #And list of which site is within each cell at each scale
    #Store in db for later use
    spatial_grids=list()
    spatial_grid_info=data.frame()
    for(this_spatial_scale in spatial_scales){
      this_grid=create_grid(this_spatial_scale, bioclim_stack)
      this_grid=subset_grid(this_grid, locations)
      spatial_grids[paste('size',this_spatial_scale, sep='-')]=this_grid
      
      spatial_grid_info = spatial_grid_info %>%
        bind_rows(assign_sites_to_grid(this_grid, locations, this_spatial_scale))
    }
    copy_to(database, spatial_grid_info, temporary = FALSE)
    

    #Extract the bioclim values for each cell size
    bioclim_bbs_data=data.frame()
    for(this_spatial_scale in spatial_scales){
      extracted <- extract_polygon(bioclim_stack, spatial_grids[[paste('size',this_spatial_scale,sep='-')]]) %>%
        group_by(id) %>%
        summarize_each(funs(mean(.,na.rm=T))) %>%
        rename(spatial_cell_id=id)

      extracted$spatial_scale=this_spatial_scale
    
      bioclim_bbs_data <- extracted %>%
        gather(bioclim_filename, value, -spatial_cell_id, -spatial_scale) %>%
        tidyr::extract(bioclim_filename, c("bioclim_var", "year"),
                       "([:alnum:]*)_([:digit:]{4})", convert=TRUE) %>%
        bind_rows(bioclim_bbs_data)
    }

    mydata <- copy_to(database, bioclim_bbs_data, temporary = FALSE,
                      indexes = list(c('spatial_scale','bioclim_var')))
    
}

#########################################

process_bioclim_data()

