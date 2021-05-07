#Set working directory
#at work
#setwd("/Users/athula/Dropbox/Writing/Tesser Behavioral/analyses")

#on laptop
setwd("/Users/athulapudhiyidath/Dropbox/Writing/TesserScan/analyses")

#import within similarity
bound_sim_all <- read.delim("within_bound_mpfcsub_masks.txt")
prim_sim_all <- read.delim("within_prim_mpfcsub_masks.txt")

#import behaviors 
behavior_all <- read.delim("TesserScan_behavior_bias_summary.txt")

#extract inference behaviors
library(dplyr)
parse_behavior <- select(behavior_all, Parse_Bound, Parse_NonBound, Parse_Diff)

#####
#putting together similarity + inference
bound_sim_parse <- cbind(bound_sim_all, parse_behavior)
prim_sim_parse <- cbind(prim_sim_all, parse_behavior)

####### ALL TOGETHER CHECK  #########
bound_mpfc <- bound_sim_parse$post_mpfc
prim_mpfc <- prim_sim_parse$post_mpfc
parse <- parse_behavior$Parse_Diff

alltog <- cbind(bound_mpfc, prim_mpfc, parse)
alltog_df <- data.frame(alltog)

library(psych)
library(tidyverse)
library(car)
library(jtools)
library(interactions)
library(emmeans)

check_tog <- lm(parse ~ bound_mpfc + prim_mpfc, data = alltog_df)
summary(check_tog)
Anova(check_tog)
vif(check_tog)
#########  #########  #########  #########
