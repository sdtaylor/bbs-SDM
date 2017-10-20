library(tidyverse)
source('config.R')
source('utils.R')

site_level_predictions = read_csv('results/bbs_results_method1.csv')

#Disregard temporal scaling for the moment
#site_level_predictions = site_level_predictions %>%
#  filter(temporal_scale==1, temporal_cell_id==1)

#Convert from deg. lat/lon on a side to km
spatial_scale_conversion = data.frame(spatial_scale_km=c(40, 80, 160, 320),
                                      spatial_scale=c(0.4, 0.8, 1.6, 3.2))
site_level_predictions = site_level_predictions %>%
  left_join(spatial_scale_conversion, by='spatial_scale')

#forecast outcomes given binary predictions
outcomes = data.frame(presence=c(0,1,1,0),
                      prediction=c(1,1,0,0),
                      type=c('fp','tp','fn','tn'))

#Gives cost per km^2/yr in terms of prescribed treatment vs actual losses
calculate_cost = function(df, treatment_cost, loss_cost, threshold=0, expense_type){
  #spatial_scale_km = unique(df$spatial_scale_km)

  if(expense_type=='perfect'){
    df$prediction = df$presence
  } else if(expense_type=='always') {
    df$prediction = 1
  } else if(expense_type=='never') {
    df$prediction = 0
  } else if(expense_type=='forecast'){
    #Apply upscaling to the forecasts only
    df = df %>%
      group_by(spatial_cell_id) %>%
      mutate(prediction=max(prediction)) %>%
      ungroup()
  } else {
    stop('No forecast type')
  }
  

  
  
  df = df %>%
    left_join(outcomes, by=c('presence','prediction'))
  
  total_tp_fp_routes = sum(df$type %in% c('fp','tp')) 
  total_fn_routes      = sum(df$type == 'fn') 
  total_cost = (treatment_cost * total_tp_fp_routes) + (loss_cost * total_fn_routes)
  
  total_routes = length(unique(df$siteID))
  print(paste(total_routes))
  
  return(total_cost/total_routes)  
}

#######################################################################################

treatment_cost = 10
possible_loss_costs = treatment_cost / seq(0.01, 1, 0.01)

value_scores = data.frame()
model_stats = data.frame()

for(this_aou in unique(site_level_predictions$Aou)){
#for(this_aou in c(3870, 4090, 4650, 4740, 3200)){
  smallest_grain_subset = site_level_predictions %>%
    filter(Aou==this_aou, spatial_scale_km==min(spatial_scales_km))
  
  #Save some metrics about the model
  prevalence = round(mean(smallest_grain_subset$presence),2)
  spec = with(smallest_grain_subset, specificity(presence, prediction))
  sens = with(smallest_grain_subset, sensitivity(presence, prediction))
    model_stats = model_stats %>%
    bind_rows(data.frame('Aou' = this_aou, prevalence = prevalence, mss = spec+sens))
  
  #loss_cost and threshold here do not matter when expense_type is perfect or always
  smallest_grain_expense_perfect = calculate_cost(smallest_grain_subset, treatment_cost = treatment_cost,
                                                  loss_cost=0, threshold=0, expense_type='perfect')
  smallest_grain_expense_always = calculate_cost(smallest_grain_subset, treatment_cost = treatment_cost,
                                                  loss_cost=0, threshold=0, expense_type='always')

  for(this_spatial_scale in spatial_scales_km){
    data_subset = site_level_predictions %>%
      filter(Aou==this_aou, spatial_scale_km==this_spatial_scale)
    
    for(this_loss_cost in possible_loss_costs){
      smallest_grain_expense_never = calculate_cost(smallest_grain_subset, treatment_cost = treatment_cost,
                                                      loss_cost=this_loss_cost, threshold=0, expense_type='never')
      a = treatment_cost/this_loss_cost
      subset_expense = calculate_cost(data_subset, treatment_cost = treatment_cost, loss_cost = this_loss_cost,
                                     threshold=a, expense_type = 'forecast')
      
      value_scores = value_scores %>%
        bind_rows(data.frame('Aou'=this_aou, 'spatial_scale_km'=this_spatial_scale, 'a'=a, threshold = a,
                             'expense_forecast' = subset_expense, 'expense_perfect'=smallest_grain_expense_perfect, 
                             'expense_always'=smallest_grain_expense_always, 'expense_never'=smallest_grain_expense_never))
      
    }
  }
}


model_stats$prevalence = round(model_stats$prevalence, 2)
model_stats$mss = round(model_stats$mss,2)

value_scores$expense_max = with(value_scores, pmin(expense_always, expense_never))
value_scores$value = with(value_scores, (expense_max - expense_forecast)/(expense_max - expense_perfect))


write_csv(value_scores, 'results/bbs_value_scores.csv')

##########################################################################
# A bunch of stuff for debugging/exploring the value models
# 
# a_ratios_of_interest = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95)
# 
# 
# for(this_aou in unique(value_scores$Aou)){
#   prevalence = model_stats %>%
#     filter(Aou==this_aou) %>%
#     pull(prevalence)
#   mss = model_stats %>%
#     filter(Aou==this_aou) %>%
#     pull(mss)
#   value_plot_title = paste(' Aou: ',this_aou,", prevalence: ",prevalence,", mss: ",mss)
#   ratio_value_plot=ggplot(filter(value_scores, Aou==this_aou), aes(y=value, x=a, color=as.factor(spatial_scale_km), group=as.factor(spatial_scale_km))) + 
#     #geom_point() + 
#     geom_line(size=1.5) +
#     ylim(0,1) +
#     theme(plot.title = element_text(size = 20),
#           axis.text = element_text(size = 30),
#           axis.title = element_text(size = 25),
#           legend.text = element_text(size = 15), 
#           legend.title = element_text(size = 20),
#           strip.text.x=element_text(size=22),
#           strip.text.y=element_text(size=22),
#           legend.key.width = unit(15, units = 'mm')) +
#     labs(title = value_plot_title,
#          color = 'Spatial Grain') 
#   
#   # ethans_plot = value_scores %>%
#   #   filter(a %in% a_ratios_of_interest, Aou==this_aou) %>%
#   #   ggplot(aes(x=spatial_scale_km, y=value, group=as.factor(a), color=as.factor(a))) + 
#   #   geom_line(size=2) +
#   #   geom_hline(yintercept = 0, size=1.5) +
#   #   scale_color_brewer(type='qual')
#   
#   
#   print(ratio_value_plot)
#   #print(ethans_plot)
# }
# 
# a_ratios = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95)
# 
# percentage_plot = value_scores %>%
#   filter(a %in% a_ratios) %>%
#   #filter(value>0) %>%
#   group_by(a, Aou) %>%
#   filter(value == max(value)) %>%
#   ungroup() %>%
#   mutate(spatial_scale_km = ifelse(value<0, 'None', spatial_scale_km)) %>%
#   group_by(a, spatial_scale_km) %>%
#   tally()
# 
# grain_sizes=c(40,80,160,320)
# pretty_grain_sizes = c('40')
# 
# 
# ggplot(percentage_plot, aes(x=as.factor(a), y=n, fill=as.factor(spatial_scale_km))) +
#   geom_bar(position = 'fill', stat='identity') +
#   scale_fill_brewer(type='qual') + 
#   theme(plot.title = element_text(size = 30),
#         axis.text = element_text(size = 20),
#         axis.title = element_text(size = 22),
#         legend.text = element_text(size = 15), 
#         legend.title = element_text(size = 20),
#         strip.text.x=element_text(size=22),
#         strip.text.y=element_text(size=22),
#         legend.key.width = unit(15, units = 'mm')) +
#   labs(fill = 'Spatial Grain',
#        x='a') 
# 
