#Set working directory
#at work
#setwd("/Users/athula/Dropbox/Writing/Tesser Behavioral/analyses")

#on laptop
setwd("/Users/athulapudhiyidath/Dropbox/Writing/TesserScan/analyses")

#ANOVA to compare the overall Acc on cover task with age group
#gather all struct data
diff_sim_all <- read.delim("diff_overall_masks.txt")

#########################
#BILATERAL HIPPO
#########################
library(dplyr)
bil_diff <- select(diff_sim_all, SUBJECT, b_hip_ant, b_hip_tail)

library(tidyr)
bil_diff_long <- gather(bil_diff, region, sim, b_hip_ant, b_hip_tail)
within_type = rep("within", nrow(bil_diff_long))
bil_diff_long <- cbind(bil_diff_long, within_type)
colnames(bil_diff_long)[which(names(bil_diff_long) == "within_type")] <- "sim_type"

bil_diff_long$region <- as.factor(bil_diff_long$region)
bil_diff_long$sim_type <- as.factor(bil_diff_long$sim_type)
bil_diff_long$SUBJECT <- as.factor(bil_diff_long$SUBJECT)

#####
#merge bilateral within and across df
library(ez)
bil_hip_anova = ezANOVA(data=bil_diff_long
                                      , dv=.(sim)
                                      , wid=.(SUBJECT)
                                      , between=.(region)
                                      , detailed=T
                                      , type=3
                                      , return_aov = TRUE)
print(bil_hip_anova)

#comparing ant/post to 0
t.test(diff_sim_all$b_hip_ant, mu=0)
t.test(diff_sim_all$b_hip_tail, mu=0)

#########################
#RIGHT HIPPO
#########################
library(dplyr)
right_diff <- select(diff_sim_all, SUBJECT, r_hip_ant, r_hip_tail)

library(tidyr)
right_diff_long <- gather(right_diff, region, sim, r_hip_ant, r_hip_tail)
within_type = rep("within", nrow(right_diff_long))
right_diff_long <- cbind(right_diff_long, within_type)
colnames(right_diff_long)[which(names(right_diff_long) == "within_type")] <- "sim_type"

right_diff_long$region <- as.factor(right_diff_long$region)
right_diff_long$sim_type <- as.factor(right_diff_long$sim_type)
right_diff_long$SUBJECT <- as.factor(right_diff_long$SUBJECT)

#####
library(ez)
right_hip_anova = ezANOVA(data=right_diff_long
                        , dv=.(sim)
                        , wid=.(SUBJECT)
                        , between=.(region)
                        , detailed=T
                        , type=3
                        , return_aov = TRUE)
print(right_hip_anova)

#conduct 2-sample t-test
t.test(diff_sim_all$r_hip_ant, diff_sim_all$r_hip_tail, alternative = "two.sided", var.equal = FALSE)

#conduct 1-sample t-test
t.test(diff_sim_all$r_hip_ant, mu=0)
t.test(diff_sim_all$r_hip_tail, mu=0)

#########################
#LEFT HIPPO
#########################
library(dplyr)
left_diff <- select(diff_sim_all, SUBJECT, l_hip_ant, l_hip_tail)

library(tidyr)
left_diff_long <- gather(left_diff, region, sim, l_hip_ant, l_hip_tail)
within_type = rep("within", nrow(left_diff_long))
left_diff_long <- cbind(left_diff_long, within_type)
colnames(left_diff_long)[which(names(left_diff_long) == "within_type")] <- "sim_type"

left_diff_long$region <- as.factor(left_diff_long$region)
left_diff_long$sim_type <- as.factor(left_diff_long$sim_type)
left_diff_long$SUBJECT <- as.factor(left_diff_long$SUBJECT)

#####
#merge bilateral within and across df
library(ez)
left_hip_anova = ezANOVA(data=left_diff_long
                          , dv=.(sim)
                          , wid=.(SUBJECT)
                          , between=.(region)
                          , detailed=T
                          , type=3
                          , return_aov = TRUE)
print(left_hip_anova)

#conduct 1-sample t-test
t.test(diff_sim_all$l_hip_ant, mu=0)
t.test(diff_sim_all$l_hip_tail, mu=0)