#Set working directory
#on Desktop
#setwd("/Users/athulapudhiyidath/Dropbox/Writing/Tesser Behavioral/analyses")

#on laptop
setwd("/Users/athulapudhiyidath/Dropbox/Experiments/TesserScan/analyses/behavior")

#ANOVA to compare the overall Acc on cover task with age group
#gather all struct data
struct_all <- read.delim("Struct_acc_rt_dprime_OVERALL_scan.txt")
struct_all$SUBJECT <- as.factor(struct_all$SUBJECT)
struct_all$SUBJCOUNT <- as.factor(struct_all$SUBJCOUNT)

#gather all inference data
inference_all <- read.delim("Induct_Results_Bias_Overall_scan.txt")
inference_all$SUBJECT <- as.factor(inference_all$SUBJECT)
inference_all$SUBJCOUNT <- as.factor(inference_all$SUB)

#importing full parse data
parse_all <- read.delim("Parse_OverWalks_Results_Overall_scan.txt")
parse_all$SUBJECT <- as.factor(parse_all$SUBJECT)

#importing full group data
group_all <- read.csv("TesserScan_grouping_perf.csv")
group_all$SUBJECT <- as.factor(group_all$SUBJECT)

#########################
#MERGE STRUCT + INFERENCE + PARSE + GROUP
#########################
struct_infer <- merge(struct_all, inference_all, by=c("SUBJECT"))
struct_infer_parse <- merge(struct_infer, parse_all, by=c("SUBJECT"))
struct_infer_parse$SUBJECT <- as.factor(struct_infer_parse$SUBJECT)

struct_infer_parse_group <- merge(struct_infer_parse, group_all, by=c("SUBJECT"))
struct_infer_parse_group$SUBJECT <- as.factor(struct_infer_parse_group$SUBJECT)

#########################
#MAKE DATA LONG FOR QUESTTYPE
#########################
library(tidyr)
struct_infer_parse_long <- gather(struct_infer_parse, quest_type, quest_bias, Prim_Bias, Bound1_Bias, Bound2_Bias)
struct_infer_parse_long$quest_type <- as.factor(struct_infer_parse_long$quest_type)


#########################
#MAKE DATA LONG FOR DISTTYPE
#########################
struct_infer_parse_group_long <- gather(struct_infer_parse_group, dist_type, dist, within_dist, across_dist)
struct_infer_parse_group_long$dist_type <- as.factor(struct_infer_parse_group_long$dist_type)


################################################################################################
#CONDUCTING MULTIPLE REGRESSION FOR TEMPORAL KNOWLEDGE TASK BIAS x PARSE PERF x GROUP PERF
library(psych)
library(tidyverse)
library(car)
library(jtools)
library(interactions)

############################################################
#running the multiple regression for inductive generalization
#Model 1: 
fit_mr1_overall <- lm(quest_bias ~ quest_type + Parse_Diff + (quest_type*Parse_Diff) + overall_dprime, data = struct_infer_parse_long)
summary(fit_mr1_overall)
Anova(fit_mr1_overall, type=3)

group_prim_bias <- lm(Prim_Bias ~ Parse_Diff, data = struct_infer_parse_group)
summary(group_prim_bias)
corr.test(select(struct_infer_parse_group, Prim_Bias, Parse_Diff))$r
corr.test(select(struct_infer_parse_group, Prim_Bias, Parse_Diff))$r^2
corr.test(select(struct_infer_parse_group, Prim_Bias, Parse_Diff))$t
corr.test(select(struct_infer_parse_group, Prim_Bias, Parse_Diff))$p

group_bound1_bias <- lm(Bound1_Bias ~ Parse_Diff, data = struct_infer_parse_group)
summary(group_bound1_bias)
corr.test(select(struct_infer_parse_group, Bound1_Bias, Parse_Diff))$r
corr.test(select(struct_infer_parse_group, Bound1_Bias, Parse_Diff))$r^2
corr.test(select(struct_infer_parse_group, Bound1_Bias, Parse_Diff))$t
corr.test(select(struct_infer_parse_group, Bound1_Bias, Parse_Diff))$p

group_bound2_bias <- lm(Bound2_Bias ~ Parse_Diff, data = struct_infer_parse_group)
summary(group_bound2_bias)
corr.test(select(struct_infer_parse_group, Bound2_Bias, Parse_Diff))$r
corr.test(select(struct_infer_parse_group, Bound2_Bias, Parse_Diff))$r^2
corr.test(select(struct_infer_parse_group, Bound2_Bias, Parse_Diff))$t
corr.test(select(struct_infer_parse_group, Bound2_Bias, Parse_Diff))$p

#Model 1: 
fit_mr2_overall <- lm(Overall_Bias ~ Parse_Diff + overall_dprime, data = struct_infer_parse)
summary(fit_mr2_overall)
Anova(fit_mr2_overall, type=3)

############################################################
#running the multiple regression for grouping task
#Model 1: 
fit_mr1_overall <- lm(dist ~ dist_type + Parse_Diff + (dist_type*Parse_Diff) + overall_dprime, data = struct_infer_parse_group_long)
summary(fit_mr1_overall)
Anova(fit_mr1_overall, type=3)

within_corr <- lm(within_dist ~ Parse_Diff, data = struct_infer_parse_group)
summary(within_corr)
corr.test(select(struct_infer_parse_group, within_dist, Parse_Diff))
corr.test(select(struct_infer_parse_group, within_dist, Parse_Diff))$r
corr.test(select(struct_infer_parse_group, within_dist, Parse_Diff))$r^2
corr.test(select(struct_infer_parse_group, within_dist, Parse_Diff))$t
corr.test(select(struct_infer_parse_group, within_dist, Parse_Diff))$p

across_corr <- lm(within_dist ~ Parse_Diff, data = struct_infer_parse_group)
summary(across_corr)
corr.test(select(struct_infer_parse_group, across_dist, Parse_Diff))$r
corr.test(select(struct_infer_parse_group, across_dist, Parse_Diff))$r^2
corr.test(select(struct_infer_parse_group, across_dist, Parse_Diff))$t
corr.test(select(struct_infer_parse_group, across_dist, Parse_Diff))$p
