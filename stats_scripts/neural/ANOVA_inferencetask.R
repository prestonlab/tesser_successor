#Set working directory
#at work
#setwd("/Users/athula/Dropbox/Writing/Tesser Behavioral/analyses")

#on laptop
setwd("/Users/athulapudhiyidath/Dropbox/Writing/TesserScan/analyses")

#ANOVA to compare the overall Acc on cover task with age group
#gather all struct data
scan_allbehave <- read.delim("TesserScan_behavior_bias_summary.txt")
scan_allbehave$SUBJECTID <- as.factor(scan_allbehave$SUBJECT)

#########################
#MAKE DATA LONG FOR INFERENCE
#########################
library(tidyr)
infer_long <- gather(scan_allbehave, quest_type, quest_bias, PrimBias, Bound1Bias, Bound2Bias)
infer_long$quest_type <- as.factor(infer_long$quest_type)
infer_long$SUBJECT <- as.factor(infer_long$SUBJECT)

#Bias <- QuestType x [dprime -- TBD]
library(ez)
###########[INCLUDED IN MANUSCRIPT]###########
inf_ancova_exps_ez = ezANOVA(data=infer_long
                                        , dv=.(quest_bias)
                                        , wid=.(SUBJECT)
                                        , within=.(quest_type)
                                        , detailed=T
                                        , type=3
                                        , return_aov = TRUE)
print(inf_ancova_exps_ez)

#conduct 2-sample t-test
t.test(scan_allbehave$PrimBias, scan_allbehave$Bound1Bias, paired = TRUE, alternative = "two.sided")
t.test(scan_allbehave$PrimBias, scan_allbehave$Bound2Bias, paired = TRUE, alternative = "two.sided")
t.test(scan_allbehave$Bound1Bias, scan_allbehave$Bound2Bias, paired = TRUE, alternative = "two.sided")

## [OVERALL PRIM x EXP, means, SD, SE] ##
library(doBy)
# broken down by Age Group AND Block
bias_type_data <- summaryBy(quest_bias ~ quest_type, data=infer_long, FUN=c(length,mean,sd))
bias_type_data

# Rename column change.length to just N
names(bias_type_data)[names(bias_type_data)=="quest_bias.length"] <- "N"
bias_type_data

# Calculate standard error of the mean
bias_type_data$quest_bias.se <- bias_type_data$quest_bias.sd / sqrt(bias_type_data$N)
bias_type_data