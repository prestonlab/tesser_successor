#Set working directory
#at work
#setwd("/Users/athula/Dropbox/Writing/Tesser Behavioral/analyses")

#on laptop
setwd("/Users/athulapudhiyidath/Dropbox/Experiments/TesserScan/analyses/neural")

library(dplyr)
#import within similarity
sim_all <- read.csv("LSMeans_scrambled_roi_WithinAcross_Subject_remixed.csv")
sim_bil_rois <- select(sim_all, within_b_hip_ant, across_b_hip_ant, within_b_hip_tail, across_b_hip_tail)
sim_right_rois <- select(sim_all, within_r_hip_ant, across_r_hip_ant, within_r_hip_tail, across_r_hip_tail)
sim_left_rois <- select(sim_all, within_l_hip_ant, across_l_hip_ant, within_l_hip_tail, across_l_hip_tail)

#import behaviors 
induct_all <- read.delim("Induct_Results_Bias_Overall_scan.txt")
struct_all <- read.delim("Struct_acc_rt_dprime_OVERALL_scan.txt")
struct_behavior <- select(struct_all, overall_dprime)

induct_sim_bil <- cbind(induct_all, struct_behavior, sim_bil_rois)
induct_sim_right <- cbind(induct_all, struct_behavior, sim_right_rois)
induct_sim_left <- cbind(induct_all, struct_behavior, sim_left_rois)
induct_sim_all <- cbind(induct_all, struct_behavior, sim_all)
#########  #########  #########  #########

library(psych)
library(tidyverse)
library(car)
library(jtools)
library(interactions)
library(emmeans)

#b_hip: ant + tail
#dist: within 
lin_model_1x <- lm(Overall_Bias ~ within_b_hip_ant + within_b_hip_tail + overall_dprime, data = induct_sim_bil)
summary(lin_model_1x)
Anova(lin_model_1x)
vif(lin_model_1x)

#within, left right, ant, post overall bias
lin_model_2x <- lm(Overall_Bias ~ within_l_hip_ant + within_l_hip_tail + within_r_hip_ant + within_r_hip_tail +  overall_dprime, data = induct_sim_all)
summary(lin_model_2x)
Anova(lin_model_2x)
vif(lin_model_2x)

##########################################################################################
#b_hip: ant + tail
#dist: within + across
lin_model1 <- lm(Overall_Bias ~ within_b_hip_ant + within_b_hip_tail + across_b_hip_ant + across_b_hip_tail + overall_dprime, data = induct_sim_bil)
summary(lin_model1)
Anova(lin_model1)
vif(lin_model1)

#b_hip: ant + tail
#dist: across 
lin_model3 <- lm(Overall_Bias ~ across_b_hip_ant + across_b_hip_tail + overall_dprime, data = induct_sim_bil)
summary(lin_model3)
Anova(lin_model3)
vif(lin_model3)

#r_hip: ant + tail
#dist: within + across
lin_model4 <- lm(Overall_Bias ~ within_r_hip_ant + within_r_hip_tail + across_r_hip_ant + across_r_hip_tail + overall_dprime, data = induct_sim_right)
summary(lin_model4)
Anova(lin_model4)
vif(lin_model4)

#r_hip: ant + tail
#dist: within
lin_model5 <- lm(Overall_Bias ~ within_r_hip_ant + within_r_hip_tail + overall_dprime, data = induct_sim_right)
summary(lin_model5)
Anova(lin_model5)
vif(lin_model5)

#l_hip: ant + tail
#dist: within + across
lin_model4 <- lm(Overall_Bias ~ within_l_hip_ant + within_l_hip_tail + across_l_hip_ant + across_l_hip_tail + overall_dprime, data = induct_sim_left)
summary(lin_model4)
Anova(lin_model4)
vif(lin_model4)
##########################################################################################
#within, left right, overall bias
lin_model5 <- lm(Overall_Bias ~ within_l_hip_ant + within_r_hip_ant + overall_dprime, data = induct_sim_all)
summary(lin_model5)
Anova(lin_model5)
vif(lin_model5)

lin_model5 <- lm(Prim_Bias ~ within_l_hip_ant + within_r_hip_ant + overall_dprime, data = induct_sim_all)
summary(lin_model5)
Anova(lin_model5)
vif(lin_model5)

lin_model5 <- lm(Bound1_Bias ~ within_l_hip_ant + within_r_hip_ant + overall_dprime, data = induct_sim_all)
summary(lin_model5)
Anova(lin_model5)
vif(lin_model5)

lin_model5 <- lm(Bound2_Bias ~ within_l_hip_ant + within_r_hip_ant + overall_dprime, data = induct_sim_all)
summary(lin_model5)
Anova(lin_model5)
vif(lin_model5)

#within, left right, Bound1 bias
lin_model6 <- lm(Bound1_Bias ~ within_l_hip_ant + within_l_hip_tail + within_r_hip_ant + within_r_hip_tail + overall_dprime, data = induct_sim_all)
summary(lin_model6)
Anova(lin_model6)
vif(lin_model6)

#within, left right, Bound2 bias
lin_model7 <- lm(Bound2_Bias ~ within_l_hip_ant + within_l_hip_tail + within_r_hip_ant + within_r_hip_tail + overall_dprime, data = induct_sim_all)
summary(lin_model7)
Anova(lin_model7)
vif(lin_model7)
