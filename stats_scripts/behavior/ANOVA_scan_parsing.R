#Set working directory
#Desktop
setwd("/Users/athulapudhiyidath/Dropbox/Experiments/TesserScan/analyses/behavior")
#ANOVA to compare the overall RTs on cover task with age group

########
#IMPORTING STRUCT LEARN
########
#importing all struct
struct_perf_full <- read.delim("Struct_acc_rt_dprime_OVERALL_scan.txt")

########
#IMPORTING PARSE
########
#importing all parse
parse_full <- read.delim("Parse_OverWalks_Results_Overall_scan.txt")

#changin relevant columns to factors
parse_full$SUBJECT <- as.factor(parse_full$SUBJECT)

#########################
#MERGE STRUCT + INFERENCE
#########################
struct_parse <- merge(struct_perf_full, parse_full, by=c("SUBJECT"))
struct_parse$SUBJECT <- as.factor(struct_parse$SUBJECT)
struct_parse$SUBJCOUNT <- as.factor(struct_parse$SUBJCOUNT)

#########################
#MAKE DATA LONG FOR PARSETYPE
#########################
library(tidyr)
struct_parse_long <- gather(struct_parse, parse_type, parse_prop, Bound_Prop, Other_Prop)
struct_parse_long$parse_type <- as.factor(struct_parse_long$parse_type)
struct_parse_long$SUBJCOUNT <- as.factor(struct_parse_long$SUBJCOUNT)


#########################
#CONDUCTING ANOVA
#########################
library(ez)
parse_ancova_overall_ez = ezANOVA(data=struct_parse_long
                                  , dv=.(parse_prop)
                                  , wid=.(SUBJCOUNT)
                                  , within=.(parse_type)
                                  , between_covariates = overall_dprime
                                  , detailed=T
                                  , type=3
                                  , return_aov = TRUE)

print(parse_ancova_overall_ez)

#########################
#CONDUCTING t-tests
#########################
t.test(struct_parse$Bound_Prop, struct_parse$Other_Prop, paired = TRUE, alternative = "two.sided")


## [OVERALL STRUCT RT, means, SD, SE] ##
#can do means like this or above

library(doBy)
# broken down by Age Group AND Block
parse_data <- summaryBy(parse_prop ~ parse_type, data=struct_parse_long, FUN=c(length,mean,sd))
parse_data

# Rename column change.length to just N
names(parse_data)[names(parse_data)=="parse_prop.length"] <- "N"
parse_data

# Calculate standard error of the mean
parse_data$parse_prop.se <- parse_data$parse_prop.sd / sqrt(parse_data$N)
parse_data

