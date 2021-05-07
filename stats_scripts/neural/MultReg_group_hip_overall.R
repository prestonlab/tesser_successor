#Set working directory
#at work
#setwd("/Users/athula/Dropbox/Writing/Tesser Behavioral/analyses")

#on laptop
setwd("/Users/athulapudhiyidath/Dropbox/Experiments/TesserScan/analyses/neural")

library(dplyr)
#import within similarity
sim_all <- read.csv("LSMeans_scrambled_roi_WithinAcross_Subject_remixed.csv")
sim_rois <- select(sim_all, within_b_hip_ant, across_b_hip_ant, within_b_hip_tail, across_b_hip_tail)

#import behaviors 
group_all <- read.csv("TesserScan_grouping_perf.csv")

#gather all struct data
struct_all <- read.delim("Struct_acc_rt_dprime_OVERALL_scan.txt")
struct_behavior <- select(struct_all, overall_dprime)

group_sim <- cbind(group_all, struct_behavior,sim_rois)

#########  #########  #########  #########

library(psych)
library(tidyverse)
library(car)
library(jtools)
library(interactions)
library(emmeans)

#b_hip: ant + tail
#dist: within 
lin_model_5x <- lm(within_dist ~ within_b_hip_ant + within_b_hip_tail + overall_dprime, data = group_sim)
summary(lin_model_5x)
Anova(lin_model_5x)
vif(lin_model_5x)

#b_hip: ant + tail
#dist: across 
lin_model_6x <- lm(across_dist ~ within_b_hip_ant + within_b_hip_tail + overall_dprime, data = group_sim)
summary(lin_model_6x)
Anova(lin_model_6x)
vif(lin_model_6x)

##########################################################################################
#b_hip: ant + tail
#dist: within + across
lin_model1 <- lm(within_dist ~ within_b_hip_ant + within_b_hip_tail + across_b_hip_ant + across_b_hip_tail + overall_dprime, data = group_sim)
summary(lin_model1)
Anova(lin_model1)
vif(lin_model1)

#b_hip: ant + tail
#dist: within 
lin_model2 <- lm(within_dist ~ within_b_hip_ant + within_b_hip_tail + overall_dprime, data = group_sim)
summary(lin_model2)
Anova(lin_model2)
vif(lin_model2)

#b_hip: ant + tail
#dist: across 
lin_model2 <- lm(within_dist ~ across_b_hip_ant + across_b_hip_tail, data = group_sim)
summary(lin_model2)
Anova(lin_model2)
vif(lin_model2)
