library(tidyverse)
source('utils.R')

results = read_csv('./results/bbs_results_method2_verbose.csv')

spatial_scale_conversion = data.frame(spatial_scale_km=c(40, 80, 160, 320),
                                      spatial_scale=c(0.4, 0.8, 1.6, 3.2))

results = results %>%
  left_join(spatial_scale_conversion, by='spatial_scale')

x= results %>%
  filter(Aou==3320, set_id==16)


outcomes = data.frame(presence=c(0,1,1,0),
                      prediction=c(1,1,0,0),
                      type=c('fp','tp','fn','tn'))

#Gives cost per ha/yr in terms of prescribed treatment vs actual losses
calculate_cost = function(df, treatment_cost, loss_cost, threshold, forecast_type){
  spatial_scale_km = unique(df$spatial_scale_km)
  temporal_scale = unique(df$temporal_scale)
  
   if(forecast_type=='perfect'){
    df$prediction = df$presence
  } else if(forecast_type=='always') {
    df$prediction = 1
  } else if(forecast_type=='never') {
    df$prediction = 0
  } else if(forecast_type=='forecast'){
    df$prediction = (df$prediction>=threshold)*1
  } else {
    stop('No forecast type')
  }
   

  df = df %>%
    filter(data_type=='testing') %>%
    left_join(outcomes, by=c('presence','prediction'))
  
  total_treated_area_years = sum(df$type %in% c('fp','tp')) * spatial_scale_km^2 * temporal_scale
  total_fn_area_years      = sum(df$type == 'fn') * spatial_scale_km^2 * temporal_scale
  total_cost = (treatment_cost * total_treated_area_years) + (loss_cost * total_fn_area_years)

  total_sqkm = length(unique(df$spatial_cell_id)) * spatial_scale_km^2
  total_sqkm_yr = total_sqkm * temporal_scale
  
  return(total_cost/total_sqkm_yr)  
}

#x = read.csv('test_forecast.csv')

#possible_treatment_costs = seq(1,500,10)
#possible_loss_costs = seq(1,500,10)

treatment_cost = 10
possible_loss_costs = 10 / seq(0.01, 1, 0.01)

training_data = x %>%
  filter(data_type == 'training')

optimal_threshold = get_threshold(training_data$presence, training_data$prediction)

cost_results = data.frame()

for(loss_cost in possible_loss_costs){
  #for(treatment_cost in possible_treatment_costs){
    cost_results = cost_results %>%
      bind_rows(data.frame('a'= treatment_cost/ loss_cost,
                           'cost_forecast' = calculate_cost(x, treatment_cost=treatment_cost, loss_cost=loss_cost, threshold = optimal_threshold, forecast_type='forecast'),
                           'cost_perfect' = calculate_cost(x, treatment_cost=treatment_cost, loss_cost=loss_cost, threshold = optimal_threshold, forecast_type='perfect'),
                           'cost_never' = calculate_cost(x, treatment_cost=treatment_cost, loss_cost=loss_cost, threshold = optimal_threshold, forecast_type='never'),
                           'cost_always' = calculate_cost(x, treatment_cost=treatment_cost, loss_cost=loss_cost, threshold = optimal_threshold, forecast_type='always')))
  #}
}

#eq. 8.3
cost_results$cost_climate = with(cost_results, pmin(cost_always, cost_never))
#eq 8.5
cost_results$value = with(cost_results, (cost_climate - cost_forecast) / (cost_climate - cost_perfect))
