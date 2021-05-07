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
parse_all <- read.delim("Parse_OverWalks_Results_Overall_scan.txt")
parse_behavior <- select(parse_all, SUBJECT, Bound_Prop, Other_Prop, Parse_Diff)

#gather all struct data
struct_all <- read.delim("Struct_acc_rt_dprime_OVERALL_scan.txt")
struct_behavior <- select(struct_all, overall_dprime)

#####
#putting together dist_roi
parse_sim <- cbind(parse_behavior, sim_rois, struct_behavior)
library(tidyr)
parse_sim_long <- gather(parse_sim, dist_type_roi, sim, within_b_hip_ant, across_b_hip_ant, within_b_hip_tail, across_b_hip_tail)
  
#separate by within/across + separate by bilateral anterior/posterior
library(reshape2)
parse_sim_split_long <- cbind(parse_sim_long, colsplit(parse_sim_long$dist_type_roi, "_", names = c("dist_type", "roi")))
parse_sim_split_long$dist_type <- as.factor(parse_sim_split_long$dist_type)
parse_sim_split_long$roi <- as.factor(parse_sim_split_long$roi)

#########  #########  #########  #########

library(psych)
library(tidyverse)
library(car)
library(jtools)
library(interactions)
library(emmeans)

#b_hip: ant + tail
#dist: within
lin_model_4x <- lm(Parse_Diff ~ within_b_hip_ant + within_b_hip_tail + overall_dprime, data = parse_sim)
summary(lin_model_4x)
Anova(lin_model_4x)
vif(lin_model_4x)

##########################################################################################
#b_hip: ant + tail
#dist: within + across
lin_model1 <- lm(Parse_Diff ~ within_b_hip_ant + within_b_hip_tail + across_b_hip_ant + across_b_hip_tail + overall_dprime, data = parse_sim)
summary(lin_model1)
Anova(lin_model1, type=3)
vif(lin_model1)

#b_hip: ant + tail
#dist: cross
lin_model3 <- lm(Parse_Diff ~ across_b_hip_ant + across_b_hip_tail + overall_dprime, data = parse_sim)
summary(lin_model3)
Anova(lin_model3)
vif(lin_model3)

