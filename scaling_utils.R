
#Take a raster stack of yearly predictions and increase the grain size
#in both time and space.
scale_up_predictions = function(p, spatial_scale, temporal_scale){
  #Scaling up in space
  p = raster::stack(aggregate(p, fact=spatial_scale, fun=mean))
  
  #Scaling up in time
  num_years = length(p@layers)

  #the effective temporal sample size after increasing the grain. 
  #ie. 9 years of data with grain size of 3 = 3 sets
  num_sets=floor(num_years/temporal_scale)
  
  
  larger_temporal_scale = raster::stack()
  for(i in (1:num_sets)-1){
    layers_to_use = (1:temporal_scale) + (temporal_scale*i)
    set_scaled = p[[layers_to_use]]
    
    set_scaled = raster::overlay(p[[layers_to_use]], fun=mean)
    
    larger_temporal_scale = raster::stack(larger_temporal_scale, set_scaled)
  }

  return(larger_temporal_scale)  
}


#######################################################
#Create a continuous grid of a certain cell size across the extent of the template raster
#######################################################
templateRaster=raster('~/data/yearly_bioclim/bio1_1971.tif')
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
#Takes a grid cell and site spatial df.
#Returns a df of colnames('cellID','siteID')
#######################################################
assign_sites_to_grid=function(g, sites){
  g=as(g, 'SpatialPolygons')
  if(!identicalCRS(g, sites)){
    sites = spTransform(sites, CRSobj = crs(g))
  }
  x=data.frame(spatial_cell_id=over(sites, g), siteID=sites@data$siteID) %>%
    filter(!is.na(spatial_cell_id))
  return(x)
}

#Aggregate 
scale_up = function(p, spatial_scale, temporal_scale){
  spatial_grid = create_grid(spatial_scale, templateRaster)
  
  
}