
#Take a raster stack of yearly predictions and increase the grain size
#in both time and space.
scale_up_predictions(p, spatial_scale, temporal_scale){
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