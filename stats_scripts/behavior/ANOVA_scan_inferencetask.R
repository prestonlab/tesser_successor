#Set working directory
#Desktop
#setwd("/Users/athula/Dropbox/Experiments/TesserScan/analyses/behavior")

#Laptop
setwd("/Users/athulapudhiyidath/Dropbox/Experiments/TesserScan/analyses/behavior")

########
#IMPORTING STRUCT LEARN
########
#importing full structure learning Acc + RT 
struct_perf_full <- read.delim("Struct_acc_rt_dprime_OVERALL_scan.txt")

struct_perf_full$SUBJECT <- as.factor(struct_perf_full$SUBJECT)

########
#IMPORTING INFERENCE
########
#importing full structure learning Acc + RT 
inference_full <- read.delim("Induct_Results_Bias_Overall_scan.txt")

#changin relevant columns to factors
inference_full$SUBJECT <- as.factor(inference_full$SUBJECT)

########
#IMPORTING INFERENCE RT
########
#importing full structure learning Acc + RT 
inference_rt_full <- read.delim("Induct_Results_Mean_Overall_scan.txt")

#changin relevant columns to factors
inference_rt_full$SUBJECT <- as.factor(inference_rt_full$SUBJECT)

########
#IMPORTING INFERENCE CORR/INCORR RT
########
#importing full structure learning Acc + RT 
inference_rt_corr_full <- read.delim("Induct_Results_RT_Overall_scan.txt")

#changin relevant columns to factors
inference_rt_corr_full$SUBJECT <- as.factor(inference_rt_corr_full$SUBJECT)

#########################
#MERGE STRUCT + INFERENCE
#########################
struct_infer <- merge(struct_perf_full, inference_full, by=c("SUBJECT"))
struct_infer$SUBJECT <- as.factor(struct_infer$SUBJECT)

##################################################
#MERGE STRUCT + INFERENCE RT
##################################################
struct_infer_rt <- merge(struct_perf_full, inference_rt_full, by=c("SUBJECT"))
struct_infer_rt$SUBJECT <- as.factor(struct_infer_rt$SUBJECT)

##################################################
#MERGE STRUCT + INFERENCE CORR RT
##################################################
struct_infer_corr_rt <- merge(struct_perf_full, inference_rt_corr_full, by=c("SUBJECT"))
struct_infer_corr_rt$SUBJECT <- as.factor(struct_infer_corr_rt$SUBJECT)


#########################
#MAKE DATA LONG FOR QUESTTYPE
#########################
library(tidyr)
inference_bias_long <- gather(struct_infer, quest_type, quest_bias, Prim_Bias, Bound1_Bias, Bound2_Bias)
inference_bias_long$quest_type <- as.factor(inference_bias_long$quest_type)

##################################################
#MAKE DATA LONG FOR QUESTTYPE: RT
##################################################
library(tidyr)
inference_rt_long <- gather(struct_infer_rt, quest_type, quest_rt, Prim_RT, Bound1_RT, Bound2_RT)
inference_rt_long$quest_type <- as.factor(inference_rt_long$quest_type)

##################################################
#MAKE DATA LONG FOR QUESTTYPE: CORR RT
##################################################
library(tidyr)
inference_corr_rt_long <- gather(struct_infer_corr_rt, quest_type, quest_rt, Prim_Corr_RT, Prim_Incorr_RT, Bound1_Corr_RT, Bound1_Incorr_RT, Bound2_Corr_RT, Bound2_Incorr_RT)
library(reshape2)
inference_split_rt_long <- cbind(inference_corr_rt_long, colsplit(inference_corr_rt_long$quest_type, "_", names = c("quest", "accuracy")))
inference_split_rt_long$quest <- as.factor(inference_split_rt_long$quest)
inference_split_rt_long$accuracy <- as.factor(inference_split_rt_long$accuracy)

options(contrasts = c("contr.sum","contr.poly"))

library(lme4) 
library(afex)  
my_lmer = lmer(formula = quest_rt ~  quest + accuracy + quest*accuracy + (1|SUBJECT) + overall_dprime
               ,data = inference_split_rt_long)
anova(my_lmer) #for mixed model
summary(my_lmer) #for p-values of mixed model

#Getting R^2 values for mixed model? 
library("MuMIn")
r.squaredGLMM(my_lmer)

library("r2glmm")
r2beta(my_lmer, method = "nsj")

#Comparing RT across quest types
library("multcomp")
summary(glht(my_lmer, linfct=mcp(quest="Tukey")), test = adjusted(type = "bonferroni"))

#
t.test(inference_split_rt_long$quest_rt[inference_split_rt_long$quest=="Prim"], inference_split_rt_long$quest_rt[inference_split_rt_long$quest=="Bound1"], paired=TRUE)
t.test(inference_split_rt_long$quest_rt[inference_split_rt_long$quest=="Prim"], inference_split_rt_long$quest_rt[inference_split_rt_long$quest=="Bound2"], paired=TRUE)
t.test(inference_split_rt_long$quest_rt[inference_split_rt_long$quest=="Bound1"], inference_split_rt_long$quest_rt[inference_split_rt_long$quest=="Bound2"], paired=TRUE)

#Comparing RT across accuracy
library("multcomp")
summary(glht(my_lmer, linfct=mcp(accuracy="Tukey")), test = adjusted(type = "bonferroni"))
#
t.test(inference_split_rt_long$quest_rt[inference_split_rt_long$accuracy=="Corr_RT"], inference_split_rt_long$quest_rt[inference_split_rt_long$accuracy=="Incorr_RT"], paired=TRUE)

#Summarizing data
library(plyr)
inference_quest_acc_rt_summary <- ddply(inference_split_rt_long, c("quest", "accuracy"), summarise,
                                        N    = sum(!is.na(quest_rt)),
                                        mean = mean(quest_rt, na.rm=TRUE),
                                        sd   = sd(quest_rt, na.rm=TRUE),
                                        se   = sd / sqrt(N)
)
inference_quest_acc_rt_summary

inference_quest_rt_summary <- ddply(inference_split_rt_long, c("quest"), summarise,
                                    N    = sum(!is.na(quest_rt)),
                                    mean = mean(quest_rt, na.rm=TRUE),
                                    sd   = sd(quest_rt, na.rm=TRUE),
                                    se   = sd / sqrt(N)
)
inference_quest_rt_summary

inference_acc_rt_summary <- ddply(inference_split_rt_long, c("accuracy"), summarise,
                                  N    = sum(!is.na(quest_rt)),
                                  mean = mean(quest_rt, na.rm=TRUE),
                                  sd   = sd(quest_rt, na.rm=TRUE),
                                  se   = sd / sqrt(N)
)
inference_acc_rt_summary

#########################
#CONDUCTING ANOVA
#########################
library(ez)
inference_ancova_overall_ez = ezANOVA(data=inference_bias_long
                                      , dv=.(quest_bias)
                                      , wid=.(SUBJCOUNT)
                                      , within=.(quest_type)
                                      , between_covariates = overall_dprime
                                      , detailed=T
                                      , type=3
                                      , return_aov = TRUE)
print(inference_ancova_overall_ez)

#########################
#CONDUCTING t-tests
#########################
inference_bias_anova_means <- aov(quest_bias ~ quest_type, data=inference_bias_long) #repeated measures, replicating ezanova
summary(inference_bias_anova_means)

#Comparing bias across types
TukeyHSD(inference_bias_anova_means)
model.tables(inference_bias_anova_means, "means") 

#Comparing bias to 0
t.test(inference_full$Prim_Bias, mu=0, alternative="greater")
t.test(inference_full$Bound1_Bias, mu=0, alternative="greater")
t.test(inference_full$Bound2_Bias, mu=0, alternative="greater")

## [OVERALL STRUCT RT, means, SD, SE] ##
#can do means like this or above

library(doBy)
# broken down by Age Group AND Block
inf_bias_data <- summaryBy(quest_bias ~ quest_type, data=inference_bias_long, FUN=c(length,mean,sd))
inf_bias_data

# Rename column change.length to just N
names(inf_bias_data)[names(inf_bias_data)=="quest_bias.length"] <- "N"
inf_bias_data

# Calculate standard error of the mean
inf_bias_data$quest_bias.se <- inf_bias_data$quest_bias.sd / sqrt(inf_bias_data$N)
inf_bias_data

#########################
#CONDUCTING ANOVA: RT
#########################
library(ez)
inference_ancova_overall_ez = ezANOVA(data=inference_rt_long
                                      , dv=.(quest_rt)
                                      , wid=.(SUBJCOUNT)
                                      , within=.(quest_type)
                                      , between_covariates = overall_dprime
                                      , detailed=T
                                      , type=3
                                      , return_aov = TRUE)
print(inference_ancova_overall_ez)

