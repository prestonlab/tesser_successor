#Set working directory
#Laptop
setwd("/Users/athulapudhiyidath/Dropbox/Experiments/TesserScan/analyses/behavior")

##################################################################################################################################################################

###############################################################
#ANOVA to compare the overall accuracy and RT on strut learning
###############################################################
#importing all blocks of learning Acc + RT 
struct_acc_rt_block1 <- read.delim("Struct_acc_rt_dprime_Part1_BLOCK_1.txt")
struct_acc_rt_block2 <- read.delim("Struct_acc_rt_dprime_Part1_BLOCK_2.txt")
struct_acc_rt_block3 <- read.delim("Struct_acc_rt_dprime_Part1_BLOCK_3.txt")
struct_acc_rt_block4 <- read.delim("Struct_acc_rt_dprime_Part1_BLOCK_4.txt")
struct_acc_rt_block5 <- read.delim("Struct_acc_rt_dprime_Part1_BLOCK_5.txt")

#binding together all the blocks
struct_acc_dprime_rt_all_blocks <- rbind(struct_acc_rt_block1, struct_acc_rt_block2, struct_acc_rt_block3, struct_acc_rt_block4, struct_acc_rt_block5)

#make some columns factors
names_toconvert <- c('SUBJECT' ,'BLOCK', 'SUBJCOUNT')
struct_acc_dprime_rt_all_blocks[names_toconvert] <- lapply(struct_acc_dprime_rt_all_blocks[names_toconvert] , factor)

# < NOW THAT EVERYTHING IS CONVERTED, CAN PERFORM SUBSEQUENT ANALYSES > #

########################################################
########### [OVERALL STRUCTURE LEARNING ACC] ########### 
#now doing a repeated measures 1-way within subjects ANOVA for accuracy
library(ez)
struct_blocks_acc_anova_v1 <- ezANOVA(data=struct_acc_dprime_rt_all_blocks
                                      , dv=overall_acc
                                      , wid=SUBJCOUNT
                                      , within=BLOCK
                                      , detailed=T
                                      , type=3)
print(struct_blocks_acc_anova_v1)

#trying to do a mixed model the aov way to ultimately save out the means and differences between conditions:
# which would translate in our case as:
struct_blocks_acc_anova_v2 <- aov(overall_acc ~ BLOCK + Error(SUBJCOUNT/BLOCK), data=struct_acc_dprime_rt_all_blocks) #repeated measures, replicating ezanova
summary(struct_blocks_acc_anova_v2)

struct_blocks_acc_anova_means <- aov(overall_acc ~ BLOCK, data=struct_acc_dprime_rt_all_blocks) #repeated measures, replicating ezanova
TukeyHSD(struct_blocks_acc_anova_means)

## [OVERALL STRUCT ACC MEANS, SD, SE] ##
model.tables(struct_blocks_acc_anova_means, "means") #can do means like this or above

library(doBy)
# broken down by Age Group AND Block
acc_blocks_data <- summaryBy(overall_acc ~ BLOCK + BLOCK, data=struct_acc_dprime_rt_all_blocks, FUN=c(length,mean,sd))
acc_blocks_data

# Rename column change.length to just N
names(acc_blocks_data)[names(acc_blocks_data)=="overall_acc.length"] <- "N"
acc_blocks_data

# Calculate standard error of the mean
acc_blocks_data$overall_rt.se <- acc_blocks_data$overall_acc.sd / sqrt(acc_blocks_data$N)
acc_blocks_data

########################################################
########### [OVERALL STRUCTURE LEARNING D-PRIME] ########### 
library(ez)
struct_blocks_dprime_anova_v1 <- ezANOVA(data=struct_acc_dprime_rt_all_blocks
                                         , dv=overall_dprime
                                         , wid=SUBJCOUNT
                                         , within=BLOCK
                                         , detailed=T
                                         , type=3)
print(struct_blocks_dprime_anova_v1)

#trying to do a mixed model the aov way to ultimately save out the means and differences between conditions:
# which would translate in our case as:
struct_blocks_dprime_anova_v2 <- aov(overall_dprime ~ BLOCK + Error(SUBJCOUNT/BLOCK), data=struct_acc_dprime_rt_all_blocks) #repeated measures, replicating ezanova
summary(struct_blocks_dprime_anova_v2)

struct_blocks_dprime_anova_means <- aov(overall_dprime ~ BLOCK, data=struct_acc_dprime_rt_all_blocks) #repeated measures, replicating ezanova
TukeyHSD(struct_blocks_dprime_anova_means)
pairwise.t.test(struct_acc_dprime_rt_all_blocks$overall_dprime, struct_acc_dprime_rt_all_blocks$BLOCK, p.adj="bonferroni")

## [OVERALL STRUCT ACC MEANS, SD, SE] ##
model.tables(struct_blocks_dprime_anova_means, "means") #can do means like this or above

library(doBy)
# broken down by Age Group AND Block
dprime_blocks_data <- summaryBy(overall_dprime ~ BLOCK + BLOCK, data=struct_acc_dprime_rt_all_blocks, FUN=c(length,mean,sd))
dprime_blocks_data

# Rename column change.length to just N
names(dprime_blocks_data)[names(dprime_blocks_data)=="overall_dprime.length"] <- "N"
dprime_blocks_data

# Calculate standard error of the mean
dprime_blocks_data$overall_dprime.se <- dprime_blocks_data$overall_dprime.sd / sqrt(dprime_blocks_data$N)
dprime_blocks_data

#######################################################
########### [OVERALL STRUCTURE LEARNING RT] ########### 
#now doing a repeated measures 1-way within subjects ANOVA for RT
library(ez)
struct_blocks_rt_anova_v1 <- ezANOVA(data=struct_acc_dprime_rt_all_blocks
                                     , dv=overall_rt
                                     , wid=SUBJCOUNT
                                     , within=BLOCK
                                     , detailed=T
                                     , type=3)
print(struct_blocks_rt_anova_v1)

#trying to do a mixed model the aov way to ultimately save out the means and differences between conditions:
# which would translate in our case as:
struct_blocks_rt_anova_v2 <- aov(overall_rt ~ BLOCK + Error(SUBJCOUNT/BLOCK), data=struct_acc_dprime_rt_all_blocks) #repeated measures, replicating ezanova
summary(struct_blocks_rt_anova_v2)

struct_blocks_rt_anova_means <- aov(overall_rt ~ BLOCK, data=struct_acc_dprime_rt_all_blocks) #repeated measures, replicating ezanova
TukeyHSD(struct_blocks_rt_anova_means)

## [OVERALL STRUCT RT, means, SD, SE] ##
model.tables(struct_blocks_rt_anova_means, "means") #can do means like this or above

library(doBy)
# broken down by Age Group AND Block
rt_blocks_data <- summaryBy(overall_rt ~ BLOCK + BLOCK, data=struct_acc_dprime_rt_all_blocks, FUN=c(length,mean,sd))
rt_blocks_data

# Rename column change.length to just N
names(rt_blocks_data)[names(rt_blocks_data)=="overall_rt.length"] <- "N"
rt_blocks_data

# Calculate standard error of the mean
rt_blocks_data$overall_rt.se <- rt_blocks_data$overall_rt.sd / sqrt(rt_blocks_data$N)
rt_blocks_data

