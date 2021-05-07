#Set working directory
#at work
#setwd("/Users/athula/Dropbox/Writing/Tesser Behavioral/analyses")

#on laptop
setwd("/Users/athulapudhiyidath/Dropbox/Writing/TesserScan/analyses")

#ANOVA to compare the overall Acc on cover task with age group
#gather all struct data
within_sim_all <- read.delim("within_overall_masks.txt")

#########################
#BILATERAL HIPPO
#########################
library(dplyr)
bil_within <- select(within_sim_all, SUBJECT, b_hip_ant, b_hip_tail)

library(tidyr)
bil_within_long <- gather(bil_within, region, sim, b_hip_ant, b_hip_tail)
within_type = rep("within", nrow(bil_within_long))
bil_within_long <- cbind(bil_within_long, within_type)
colnames(bil_within_long)[which(names(bil_within_long) == "within_type")] <- "sim_type"

bil_within_long$region <- as.factor(bil_within_long$region)
bil_within_long$sim_type <- as.factor(bil_within_long$sim_type)
bil_within_long$SUBJECT <- as.factor(bil_within_long$SUBJECT)

#####
#merge bilateral within and across df
library(ez)
bil_hip_anova = ezANOVA(data=bil_within_long
                                      , dv=.(sim)
                                      , wid=.(SUBJECT)
                                      , between=.(region)
                                      , detailed=T
                                      , type=3
                                      , return_aov = TRUE)
print(bil_hip_anova)

#conduct 1-sample t-test
t.test(within_sim_all$b_hip_ant, mu=0)
t.test(within_sim_all$b_hip_tail, mu=0)

#########################
#RIGHT HIPPO
#########################
library(dplyr)
right_within <- select(within_sim_all, SUBJECT, r_hip_ant, r_hip_tail)

library(tidyr)
right_within_long <- gather(right_within, region, sim, r_hip_ant, r_hip_tail)
within_type = rep("within", nrow(right_within_long))
right_within_long <- cbind(right_within_long, within_type)
colnames(right_within_long)[which(names(right_within_long) == "within_type")] <- "sim_type"

right_within_long$region <- as.factor(right_within_long$region)
right_within_long$sim_type <- as.factor(right_within_long$sim_type)
right_within_long$SUBJECT <- as.factor(right_within_long$SUBJECT)

#####
library(ez)
right_hip_anova = ezANOVA(data=right_within_long
                        , dv=.(sim)
                        , wid=.(SUBJECT)
                        , between=.(region)
                        , detailed=T
                        , type=3
                        , return_aov = TRUE)
print(right_hip_anova)

#conduct 1-sample t-test
t.test(within_sim_all$r_hip_ant, mu=0)
t.test(within_sim_all$r_hip_tail, mu=0)

#########################
#LEFT HIPPO
#########################
library(dplyr)
left_within <- select(within_sim_all, SUBJECT, l_hip_ant, l_hip_tail)

library(tidyr)
left_within_long <- gather(left_within, region, sim, l_hip_ant, l_hip_tail)
within_type = rep("within", nrow(left_within_long))
left_within_long <- cbind(left_within_long, within_type)
colnames(left_within_long)[which(names(left_within_long) == "within_type")] <- "sim_type"

left_within_long$region <- as.factor(left_within_long$region)
left_within_long$sim_type <- as.factor(left_within_long$sim_type)
left_within_long$SUBJECT <- as.factor(left_within_long$SUBJECT)

#####
#merge bilateral within and across df
library(ez)
left_hip_anova = ezANOVA(data=left_within_long
                          , dv=.(sim)
                          , wid=.(SUBJECT)
                          , between=.(region)
                          , detailed=T
                          , type=3
                          , return_aov = TRUE)
print(left_hip_anova)

#conduct 1-sample t-test
t.test(within_sim_all$l_hip_ant, mu=0)
t.test(within_sim_all$l_hip_tail, mu=0)
