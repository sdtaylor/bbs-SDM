library(dplyr)
library(ggplot2)

#Graphs for gordon conference 2016. Uses the fraction skill score. 
##############################################################

#Spatial scale on the x axis, skill score on the y, different lines for temporal scale. 
results=read.csv('results/bbsSDMResults_with_timelags_all.csv') %>%
  filter(!is.na(Aou), !is.na(fss))


lags_to_keep=data.frame(temporal_scale=c(1,5,10), 
                        time_lag=c(35.5, 37.5, 35.0),
                        keep='yes')

results_summarized =results %>%
  group_by(temporal_scale, cellSize, windowID, modelName, Aou) %>%
  summarize(skill=mean(fss), skill_sd=sd(fss)) %>%
  ungroup() %>%
  mutate(time_lag=(temporal_scale * windowID) - (temporal_scale/2)) %>%
  left_join(lags_to_keep, by=c('temporal_scale','time_lag')) %>%
  filter(keep=='yes')

results_summarized$modelName = factor(results_summarized$modelName, levels=c('gbm','naive'), labels = c('Gradient Boosting Model','Naive Model'), ordered = TRUE)
results_summarized$cellSize = factor(results_summarized$cellSize, levels=c(0.1,1.0,2.0), labels = c('0.1°','1.0°','2.0°'), ordered = TRUE)



ggplot(results_summarized, aes(x=temporal_scale, y=skill, group=as.factor(cellSize), color=as.factor(cellSize))) +
  #geom_point()+
  geom_smooth(se=T, alpha=0.2, method='loess', size=3) +
  scale_x_continuous(breaks=c(1,5,10))+
  scale_colour_manual(values=c('#E69F00','#009E73','#CC79A7')) +
  theme_bw() +
  facet_grid(.~modelName) + 
  theme(legend.position = c(0.8, 0.20), 
        legend.direction = "horizontal",
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 20),
        strip.text.x=element_text(size=22),
        strip.text.y=element_text(size=22),
        panel.grid.major = element_line(colour = "gray64"),
        panel.grid.minor = element_line(colour = NA)
  ) +
  labs(colour = "Spatial Scale (deg. lat/long)", 
       y = "Fractions Skill Score", 
       x = "Temporal Scale (yrs)")  

##############################################################
#Single species graph over time and scales

this_sp=6840 #Hooded Warbler
results_this_sp =results %>%
  filter(Aou==this_sp,temporal_scale %in% c(5,10), cellSize %in% c(0.1, 1)) %>%
  #filter(Aou==this_sp) %>%
  group_by(temporal_scale, cellSize, windowID, modelName) %>%
  summarize(skill=mean(fss), skill_sd=sd(fss)) %>%
  ungroup() %>%
  mutate(time_lag=(temporal_scale * windowID) - (temporal_scale/2))


results_this_sp$temporal_scale = factor(results_this_sp$temporal_scale, levels=c(1,3,5,10), labels = c('Temporal Scale: 1 yrs','Temporal Scale: 3yrs','Temporal Scale: 5yrs','Temporal Scale: 10yrs'), ordered = TRUE)
results_this_sp$cellSize = factor(results_this_sp$cellSize, levels=c(0.1,1.0,2.0), labels = c('Spatial Scale: 0.1°','Spatial Scale: 1.0°','Spatial Scale: 2.0°'), ordered = TRUE)


ggplot(results_this_sp, aes(x=time_lag, y=skill, colour=modelName, group=modelName)) +
  geom_point(size=5) +
  geom_line(linetype='dashed', size=1.5) +
  geom_smooth(method='lm', se=FALSE, linetype='solid', size=1.8, aes(group=modelName)) +
  scale_colour_manual(values=c('#0072B2','#E69F00')) +
  #geom_hline(yintercept = 0.75) +
  theme_bw() +
  facet_grid(cellSize~temporal_scale, labeller = label_value) +
  theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 20), 
        plot.title = element_text(size = 21), 
        legend.text = element_text(size = 23), 
        legend.key = element_rect(fill = "white"), 
        legend.background = element_rect(fill = "white"), 
        strip.text.x=element_text(size=22),
        strip.text.y=element_text(size=22),
        legend.position = "none", legend.direction = "horizontal") +labs(x = "Years Into Future", 
                                                                           y = "Fractions Skill Score", colour = NULL)



#############################################################
library(tidyr)
#Histogram of slopes of the skill trends over time
trend_lines = results %>%
  group_by(Aou, modelName, temporal_scale, cellSize, windowID) %>%
  summarize(skill=mean(fss), skill_sd=sd(fss)) %>%
  ungroup() %>%
  mutate(time_lag=(temporal_scale * windowID) - (temporal_scale/2)) %>%
  group_by(Aou, modelName, temporal_scale, cellSize) %>%
  do(model=lm(skill~time_lag, data=.)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(slope=unlist(model)[2]$'coefficients.time_lag', intercept=unlist(model)[1]$`coefficients.(Intercept)`) %>%
  dplyr::select(-model) %>%
  ungroup() %>%
  filter(!is.na(slope))

#Extrapolatin of trend into future
future_skill = trend_lines %>%
  mutate(skill_at_80yrs=slope*150+ intercept) %>%
  mutate(skill_at_80yrs=ifelse(skill_at_80yrs<0, 0, ifelse(skill_at_80yrs>1,1, skill_at_80yrs)))


ggplot(future_skill, aes(x=cellSize, y=skill_at_80yrs, group=as.factor(temporal_scale), color=as.factor(temporal_scale))) +
  #geom_point()+
  geom_smooth(se=T, alpha=0.2, method='loess', size=3) +
  scale_colour_manual(values=c('#E69F00','#009E73','#CC79A7')) +
  theme_bw() +
  facet_grid(.~modelName) + 
  theme(legend.position = c(0.8, 0.20), 
        legend.direction = "horizontal",
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 19),
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 20),
        strip.text.x=element_text(size=20),
        strip.text.y=element_text(size=20),
        panel.grid.major = element_line(colour = "gray64"),
        panel.grid.minor = element_line(colour = NA)
  ) +
  labs(x = "Spatial Scale (deg. lat/long)", 
       y = "Fractions Skill Score", 
       colour = "Temporal Scale (yrs)")  




















##Histograms of slopes

slopes=trend_lines %>%
  dplyr::select(-intercept) %>%
  spread(modelName, slope) %>%
  mutate(gbm_slope_is_more=ifelse(gbm>naive, 'yes','no')) %>%
  group_by(temporal_scale, cellSize, gbm_slope_is_more) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  filter(!is.na(gbm_slope_is_more)) %>%
  spread(gbm_slope_is_more, n) %>%
  mutate(percent=yes/(yes+no))


slopes=trend_lines %>%
  filter( cellSize %in% c(0.1, 1.0), temporal_scale %in% c(5,10))

slopes$cellSize = factor(slopes$cellSize, levels=c(0.1,1.0,2.0), labels = c('Spatial Scale: 0.1°','Spatial Scale: 1.0°','Spatial Scale: 2.0°'), ordered = TRUE)
slopes$temporal_scale = factor(slopes$temporal_scale, levels=c(1,3,5,10), labels = c('Temporal Scale: 1 yrs','Temporal Scale: 3yrs','Temporal Scale: 5yrs','Temporal Scale: 10yrs'), ordered = TRUE)



ggplot(slopes, aes(slope, fill=modelName))+
  #stat_density( kernel='gaussian')+
  geom_density(alpha=0.6) +
  #geom_histogram(bins=50)+
  scale_fill_manual(values=c('#0072B2','#E69F00')) +
  xlim(-0.01, 0.01)+
  #geom_freqpoly()+
  theme_bw()  +
  geom_vline(xintercept = 0)+
  facet_grid(cellSize~temporal_scale) +
  theme(panel.grid.major = element_line(colour = NA),
        legend.title = element_text(colour = NA), 
        legend.title = element_text(size = 1),
        legend.position = 'none',
        strip.text.x=element_text(size=22),
        strip.text.y=element_text(size=22),
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 20)
        )+
  labs(x = "Trend in skill over time", 
       y = "Density Estimate", 
       colour = "Temporal Scale (yrs)",
       title=NULL,
       fill='.')   




################################################################
#Change in difference of skill between gbm and naive model over time.
library(tidyr)

skill_diff=results %>%
  group_by(Aou, modelName, temporal_scale, cellSize, windowID) %>%
  summarize(skill=mean(fss)) %>%
  ungroup() %>%
  mutate(time_lag=(temporal_scale * windowID) - (temporal_scale/2)) %>%
  spread(modelName, skill) %>%
  mutate(diff=gbm-naive) %>%
  group_by(temporal_scale, cellSize, time_lag) %>%
  summarize(diff_sd=sd(diff, na.rm=T), diff=mean(diff, na.rm=T)) %>%
  ungroup() %>%
  filter( cellSize %in% c(0.1, 1.0), temporal_scale %in% c(5,10))

indicator_labels=data.frame(x=c(9,9),
                            y=c(0.25, -0.25),
                            t=c('GBM Better', 'Naive Better'),
                            temporal_scale=10, cellSize=0.1)

indicator_lines=data.frame(x=c(9,9), xend=c(9,9), 
                           y=c(0.10,-0.10), yend=c(0.20, -0.2),
                           temporal_scale=10, cellSize=0.1)

skill_diff$cellSize = factor(skill_diff$cellSize, levels=c(0.1,1.0,2.0), labels = c('Spatial Scale: 0.1°','Spatial Scale: 1.0°','Spatial Scale: 2.0°'), ordered = TRUE)
skill_diff$temporal_scale = factor(skill_diff$temporal_scale, levels=c(1,3,5,10), labels = c('Temporal Scale: 1 yrs','Temporal Scale: 3yrs','Temporal Scale: 5yrs','Temporal Scale: 10yrs'), ordered = TRUE)
indicator_labels$cellSize = factor(indicator_labels$cellSize, levels=c(0.1,1.0,2.0), labels = c('Spatial Scale: 0.1°','Spatial Scale: 1.0°','Spatial Scale: 2.0°'), ordered = TRUE)
indicator_labels$temporal_scale = factor(indicator_labels$temporal_scale, levels=c(1,3,5,10), labels = c('Temporal Scale: 1 yrs','Temporal Scale: 3yrs','Temporal Scale: 5yrs','Temporal Scale: 10yrs'), ordered = TRUE)
indicator_lines$cellSize = factor(indicator_lines$cellSize, levels=c(0.1,1.0,2.0), labels = c('Spatial Scale: 0.1°','Spatial Scale: 1.0°','Spatial Scale: 2.0°'), ordered = TRUE)
indicator_lines$temporal_scale = factor(indicator_lines$temporal_scale, levels=c(1,3,5,10), labels = c('Temporal Scale: 1 yrs','Temporal Scale: 3yrs','Temporal Scale: 5yrs','Temporal Scale: 10yrs'), ordered = TRUE)


ggplot(skill_diff, aes(x=time_lag, y=diff)) +
  #geom_smooth(se=T, alpha=0.2, method='loess') +
  geom_point(size=5) +
  geom_line(size=2)+
  geom_errorbar(aes(ymax=diff+(diff_sd*2), ymin=diff-(diff_sd*2)), alpha=0.7, size=2) +
  geom_hline(yintercept = 0)+
  #ylim(-0.25, 0.1)+
  theme_bw() +
  geom_text(data=indicator_labels, aes(x=x, y=y, label=t), size=6, inherit.aes = F) +
  geom_segment(data=indicator_lines, aes(x=x, xend=xend, y=y, yend=yend), size=1, arrow=arrow(length=unit(0.25,'cm')), inherit.aes = F)+
  #annotate(geom='segment', x=c(11,11), xend=c(11,11), y=c(0.35,-0.35), yend=c(0.45, -0.45), arrow=arrow(length=unit(0.25,'cm')) )+
  
  facet_grid(cellSize~temporal_scale) + 
  #annotate('text', x=c(10,10), y=c(0.5, -0.5), label=c('GBM Better', 'Naive Better')) +
  
  theme(legend.position = c(0.8, 0.20), 
        legend.direction = "horizontal",
        axis.title = element_text(size = 30), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16),
        strip.text.x=element_text(size=18),
        strip.text.y=element_text(size=18),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA)
  ) +
  labs(x = "Time Lag", 
       y = "Difference between GBM and Naive models") 










