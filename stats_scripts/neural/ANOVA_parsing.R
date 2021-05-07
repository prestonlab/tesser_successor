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
parse_long <- gather(scan_allbehave, parse_type, parse_prop, Parse_Bound, Parse_NonBound)
parse_long$parse_type <- as.factor(parse_long$parse_type)
parse_long$SUBJECT <- as.factor(parse_long$SUBJECT)

#PropParsed <- ParseType x [dprime -- TBD]
library(ez)
###########[INCLUDED IN MANUSCRIPT]###########
parse_ancova_exps_ez = ezANOVA(data=parse_long
                             , dv=.(parse_prop)
                             , wid=.(SUBJECT)
                             , within=.(parse_type)
                             , detailed=T
                             , type=3
                             , return_aov = TRUE)
print(parse_ancova_exps_ez)

#conduct 2-sample t-test
t.test(scan_allbehave$Parse_Bound, scan_allbehave$Parse_NonBound, paired = TRUE, alternative = "two.sided")

## [OVERALL PRIM x EXP, means, SD, SE] ##
library(doBy)
# broken down by Age Group AND Block
parse_type_data <- summaryBy(parse_prop ~ parse_type, data=parse_long, FUN=c(length,mean,sd))
parse_type_data

# Rename column change.length to just N
names(parse_type_data)[names(parse_type_data)=="parse_prop.length"] <- "N"
parse_type_data

# Calculate standard error of the mean
parse_type_data$quest_bias.se <- parse_type_data$parse_prop.sd / sqrt(parse_type_data$N)
parse_type_data
