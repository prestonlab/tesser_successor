#Set working directory
#at work
setwd("/Users/athula/Dropbox/Experiments/TesserScan/Neural Results/rsa_analyses")

#on laptop
#setwd("/Users/athulapudhiyidath/Dropbox/Writing/TesserScan/analyses")

#ANOVA to compare the overall Acc on cover task with age group
#gather all struct data
scan_allbehave <- read.delim("TesserScan_behavior_bias_summary.txt")
scan_allbehave$SUBJECTID <- as.factor(scan_allbehave$SUBJECT)

#########################
#MAKE DATA LONG FOR INFERENCE
#########################
library(tidyr)
group_long <- gather(scan_allbehave, group_type, dist, Within, Across)
group_long$group_type <- as.factor(group_long$group_type)
group_long$SUBJECT <- as.factor(group_long$SUBJECT)

#PropParsed <- ParseType x [dprime -- TBD]
library(ez)
###########[INCLUDED IN MANUSCRIPT]###########
group_ancova_exps_ez = ezANOVA(data=group_long
                               , dv=.(dist)
                               , wid=.(SUBJECT)
                               , within=.(group_type)
                               , detailed=T
                               , type=3
                               , return_aov = TRUE)
print(group_ancova_exps_ez)

#conduct 2-sample t-test
t.test(scan_allbehave$Within, scan_allbehave$Across, paired = TRUE, alternative = "two.sided")

## [OVERALL PRIM x EXP, means, SD, SE] ##
library(doBy)
# broken down by Age Group AND Block
group_type_data <- summaryBy(dist ~ group_type, data=group_long, FUN=c(length,mean,sd))
group_type_data

# Rename column change.length to just N
names(group_type_data)[names(group_type_data)=="dist.length"] <- "N"
group_type_data

# Calculate standard error of the mean
group_type_data$dist.se <- group_type_data$dist.sd / sqrt(group_type_data$N)
group_type_data
