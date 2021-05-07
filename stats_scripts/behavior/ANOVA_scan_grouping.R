#Set working directory
#Desktop
setwd("/Users/athulapudhiyidath/Dropbox/Experiments/TesserScan/analyses/behavior")
#ANOVA to compare the overall RTs on cover task with age group


#ANOVA to compare the overall RTs on cover task with age group

########
#IMPORTING STRUCT LEARN
########
#importing all struct
struct_perf_full <- read.delim("Struct_acc_rt_dprime_OVERALL_scan.txt")

########
#IMPORTING GROUP
########
#importing all parse
group_full <- read.csv("TesserScan_grouping_perf.csv")

#changin relevant columns to factors
group_full$SUBJECT <- as.factor(group_full$SUBJECT)

#########################
#MERGE STRUCT + INFERENCE
#########################
struct_group <- merge(struct_perf_full, group_full, by=c("SUBJECT"))
struct_group$SUBJECT <- as.factor(struct_group$SUBJECT)

#########################
#MAKE DATA LONG FOR PARSETYPE
#########################
library(tidyr)
struct_group_long <- gather(struct_group, group_type, dist, within_dist, across_dist)
struct_group_long$group_type <- as.factor(struct_group_long$group_type)
struct_group_long$SUBJCOUNT <- as.factor(struct_group_long$SUBJCOUNT)


#########################
#CONDUCTING ANOVA
#########################
library(ez)
group_ancova_overall_ez = ezANOVA(data=struct_group_long
                                  , dv=.(dist)
                                  , wid=.(SUBJCOUNT)
                                  , within=.(group_type)
                                  , between_covariates = overall_dprime
                                  , detailed=T
                                  , type=3
                                  , return_aov = TRUE)

print(group_ancova_overall_ez)

#########################
#CONDUCTING t-tests
#########################
t.test(struct_group$within_dist, struct_group$across_dist, paired = TRUE, alternative = "two.sided")


## [OVERALL STRUCT RT, means, SD, SE] ##
#can do means like this or above

library(doBy)
# broken down by Age Group AND Block
group_data <- summaryBy(dist ~ group_type, data=struct_group_long, FUN=c(length,mean,sd))
group_data

# Rename column change.length to just N
names(group_data)[names(group_data)=="dist.length"] <- "N"
group_data

# Calculate standard error of the mean
group_data$dist.se <- group_data$dist.sd / sqrt(group_data$N)
group_data

