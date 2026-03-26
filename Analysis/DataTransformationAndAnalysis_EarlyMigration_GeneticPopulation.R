### Early Migration Genetic Population Data Transformation and AICc Analysis Code of
### Host et al. Escapement quality of sockeye salmon (Oncorhynchus nerka) returning to a glacially dominated watershed: relationships among body size, somatic energy, reproductive investment, and migratory difficulty
### Published with Canadian Journal of Fisheries and Aquatic Sciences
# March 26, 2026

#Important notes on the genetic population data: It only has data from 2020 and 2021. Raw data file is called "EarlyMigration_GeneticPopulation_RawData.csv"

# ================= Body Size (PC1 Score) Data Transformation =============
#Clear Environment, load important packages
rm(list = ls())
library(AICcmodavg)
library(lubridate)

#download original data
OriginalData<-read.csv("EarlyMigration_RunTimingGroup_BodyCondition.csv")
  #if you wish to skip this step, you can use the transformed data "EarlyMigration_RunTimingGroup_BodySize.csv" and jump straight to the AICc analysis

### Removing NAs that will mess up PCA code
any(is.na(OriginalData$Group_Assignment)) #this is the genetic population assignment variable
any(is.na(OriginalData$Year))
any(is.na(OriginalData$Sex)) #if any NAs, remove them!


### Now, we can create the actual PC score for the body size that will be used in the AICc analysis ####
#OriginalData$Fish_Leng_1_mm 
#OriginalData$Fish_Ht_mm 
#OriginalData$Fish_Grth_mm 

bc_pca2 <- princomp(OriginalData[,c(75,76,77)], cor=T, scores=T, covmat = NULL) #should work now without any NAs, identical scores to that of the one above.  units do not matter!
summary(bc_pca2, loadings=T, cutoff=0.0001) #summary of PC analysis
bc_pca2$scores
screeplot(bc_pca2, type=c('lines'))
bodysize_pc_g_mm <- bc_pca2$scores[,1] #creates a variable that is the body size index scores

#combine dataframe with column for body size PC scores
OriginalData <- cbind(OriginalData, bodysize_pc_g_mm) 

#Look at Genetic Population variable and determine if trimming needs to happen
OriginalData$Group_Assignment <- as.factor(OriginalData$Group_Assignment)
levels(OriginalData$Group_Assignment) #should be 9 levels

#Now, we have a dataset that can be saved as fully transformed and ready to be analyzed for body size AICc
write.csv(OriginalData, "EarlyMigration_GeneticPopulation_BodyCondition.csv")



# ================= Energy Density Data Transformation =============
rm(list = ls()) #clear environment
OriginalData <-read.csv("EarlyMigration_GeneticPopulation_BodyCondition.csv")

## Again, we need to check for the presence of NAs that will mess up levels / model weight
any(is.na(OriginalData$Group_Assignment))
any(is.na(OriginalData$Year))
any(is.na(OriginalData$Sex)) #if any NAs, remove!

#Look at Genetic Population variable and determine if trimming needs to happen
OriginalData$Group_Assignment <- as.factor(OriginalData$Group_Assignment)
levels(OriginalData$Group_Assignment) #should be 9 levels

#transform energy density to correct units (mJ/g)
OriginalData$EnergyPDry_1_mJ_g <- OriginalData$EnergyPDry_1 / 1000

#save data as "transformed" 
write.csv(OriginalData, "EarlyMigration_GeneticPopulation_BodyCondition.csv")



# ================= Total Energy Data Transformation =============
rm(list = ls()) #clear environment
OriginalData <-read.csv("EarlyMigration_GeneticPopulation_BodyCondition.csv")

## Again, we need to check for the presence of NAs that will mess up levels / model weight
any(is.na(OriginalData$Collection_Location))
any(is.na(OriginalData$Group_Assignment))
any(is.na(OriginalData$Year))
any(is.na(OriginalData$Sex)) #remove any NAs


#Look at Genetic Population variable and determine if trimming needs to happen
OriginalData$Group_Assignment <- as.factor(OriginalData$Group_Assignment)
levels(OriginalData$Group_Assignment) #should have 9 levels

#Calculate total energy variable (mJ)
OriginalData$EnergyPDry_1_mJ_g <- OriginalData$EnergyPDry_1 / 1000
OriginalData$Carcass.guts_Wt_g <- OriginalData$Carcass.guts_Wt * 1000
OriginalData$TotalEnergy_carcass_mJ <- OriginalData$EnergyPDry_1_mJ_g * OriginalData$Carcass.guts_Wt_g
str(OriginalData)

write.csv(OriginalData, "EarlyMigration_GeneticPopulation_BodyCondition.csv")



# ================= Fecundity Data Transformation =============
rm(list = ls()) #clear environment
##for fecundity, we will use previously transformed body size data, because the PC score will be an explanatory variable in the resulting AICc analysis
OriginalData <-read.csv("EarlyMigration_GeneticPopulation_BodySize.csv")

## Removing NAs that could disrupt model analysis, check each explanatory variable
any(is.na(OriginalData$Group_Assignment))
any(is.na(OriginalData$Year))
any(is.na(OriginalData$Sex)) #should all be False, as we previously checked for this and this is the transformed data

#subset data for females only
OriginalDataF <- subset(OriginalData, Sex == "F")

#calculate a fecundity metric 
# fecundity = (# of eggs / weight of that sample in g) * (total gonad weight in kg * 1000)
OriginalDataF$Fecun_num_by_wt1 <- (OriginalDataF$Fem_Fecun_1_num/OriginalDataF$Fem_Fecun_1_wt)*(OriginalDataF$Gonad_Wt * 1000)

#save data as "transformed" 
write.csv(OriginalData, "EarlyMigration_GeneticPopulation_Fecundity.csv")




# ================= Fineness Ratio Data Transformation =============
rm(list = ls()) #clear environment
#for fineness, we will use previously transformed body size data
OriginalData <-read.csv("EarlyMigration_GeneticPopulation_BodyCondition.csv")

#All data is already transformed, so there should be no NAs, but still check
## Removing NAs that could disrupt model analysis, check each explanatory variable
any(is.na(OriginalData$Group_Assignment))
any(is.na(OriginalData$Year))
any(is.na(OriginalData$Sex)) #should all be False, as we previously checked for this and this is the transformed data

### Fineness Ratio = A fish's "fineness ratio" refers to the ratio of its body length to its maximum cross-sectional diameter
#Need to calculate fineness according to sex
library(dplyr)
# Fineness Ratio 2: fish length / fish width, dependent on sex
OriginalData <- OriginalData %>%
  group_by(Sex) %>%
  mutate(fineness = Fish_Leng_1_mm / Fish_Wdth_mm)

str(OriginalData$fineness) #fineness is numeric

#save data as "transformed" 
write.csv(OriginalData, "EarlyMigration_GeneticPopulation_BodyCondition.csv")




# ================= AICc Analysis and Model Averaging Example for Response Variable ~ Genetic Population Analysis =============
# I will provide one full walk through example of how we completed our AICc analysis and resulting model average approach for this manuscript. This code can be reflected across all response variables. I will also show the model set for fecundity, as it is slightly different 
### Data Analysis for energy density and total energy follow the same exact analysis as above, just replacing the transformed data sets with respective data sets and replacing body size with the response variable of interest
rm(list = ls()) #clear environment
#for fineness, we will use previously transformed body size data
DataSet <-read.csv("EarlyMigration_GeneticPopulation_BodyCondition.csv")

DataSet$Group_Assignment <- as.factor(DataSet$Group_Assignment)
DataSet$Year<- as.factor(DataSet$Year)
DataSet$Sex<- as.factor(DataSet$Sex)
str(DataSet)

model1 <- lm(bodysize_pc_g_mm~1, data=DataSet)
model2 <- lm(bodysize_pc_g_mm~Collection_RiverMile_m + Group_Assignment, data = DataSet)
model3 <- lm(bodysize_pc_g_mm~Collection_RiverMile_m + Sex, data = DataSet)
model4 <- lm(bodysize_pc_g_mm~Collection_RiverMile_m + Group_Assignment + Sex, data = DataSet)


#AICc model selection
library(AICcmodavg)
AIC(model1, model2, model3, model4)

#table of AIC results
mynames1 <- paste("model", as.character(1:4), sep = "")
models <- list(model1, model2, model3, model4)

# Generate AIC table
myaicc1 <- aictab(models, modnames = mynames1)
print(myaicc1)
# Convert AIC table to a data frame for easier manipulation
aic_df <- as.data.frame(myaicc1)

colnames(aic_df)[colnames(aic_df) == "Modnames"] <- "Model #"

#Now, we can model average these results, Model averaging based on AICc, conditional averaging (default is full=TRUE)
library(MuMIn)
model_avg <- model.avg(models, rank = "AICc", full = TRUE)

# Summary with coefficients, SE, CIs
summary(model_avg)

# Full averaged coefficients
avg_coefs_full <- coef(model_avg, full = TRUE)

# Full averaged confidence intervals
avg_confint_full <- confint(model_avg, full = TRUE)

# Full averaged variance-covariance matrix and SE
avg_vcov_full <- vcov(model_avg, full = TRUE)
avg_se_full <- sqrt(diag(avg_vcov_full))

# Create tidy data frame with full averages
coefs_df_full <- as.data.frame(avg_coefs_full)
colnames(coefs_df_full) <- "Avg_Coef"
coefs_df_full$CI_Lower <- avg_confint_full[, 1]
coefs_df_full$CI_Upper <- avg_confint_full[, 2]
coefs_df_full$SE <- avg_se_full

# Calculate parameter importance (parameter likelihoods)
param_likelihoods <- sw(model_avg)
as.data.frame(param_likelihoods)

# Manually specify parameter likelihoods from sw(model_avg)
manual_likelihoods <- c(
  "(Intercept)" = 1.0,
  "Collection_RiverMile_m" = 0.9979668,
  "Group_AssignmentChitina" = 0.2343056,
  "Group_AssignmentGulkana" = 0.2343056,
  "Group_AssignmentGulkanaHatchery" = 0.2343056,
  "Group_AssignmentKlutinaLake" = 0.2343056,
  "Group_AssignmentKlutinaTonsinaOutlet" = 0.2343056,
  "Group_AssignmentSlana" = 0.2343056,
  "Group_AssignmentTanadaCopperLakes" = 0.2343056,
  "Group_AssignmentTazlina" = 0.2343056,
  "SexM" = 0.9502672 #these numbers will all change according to the analysis, but this is the code used to create the table.
)


# Add ParamLikelihoods by matching rownames
coefs_df_full$ParamLikelihood <- manual_likelihoods[rownames(coefs_df_full)]

# View final table
print(coefs_df_full) #this dataframe can be saved as a .csv for future reference and further analysis



# ================= AICc Analysis and Model Averaging Example for Fecundity ~ Genetic Population Analysis =============
# All of the model analysis is the same, but the model set is slightly larger and looks like this, so replace with respective model set and be sure to adjust the AICc analysis and model averaging accordingly
DataSetF <- read.csv("EarlyMigration_GeneticPopulation_Fecundity.csv")

DataSet$Group_Assignment <- as.factor(DataSet$Group_Assignment)
DataSet$Year<- as.factor(DataSet$Year)
DataSet$Sex<- as.factor(DataSet$Sex)
str(DataSet)

#8 models for fecundity analysis
model1 <- lm(Fecun_num_by_wt1~1, data=DataSetF)
model2 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + bodysize_pc_g_mm, data=DataSetF)
model3 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + EnergyPDry_1_mJ_g, data=DataSetF)
model4 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + Group_Assignment, data=DataSetF)
model5 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + bodysize_pc_g_mm + EnergyPDry_1_mJ_g, data=DataSetF)
model6 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + bodysize_pc_g_mm + Group_Assignment, data=DataSetF)
model7 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + EnergyPDry_1_mJ_g + Group_Assignment, data=DataSetF)
model8 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + bodysize_pc_g_mm + EnergyPDry_1_mJ_g + Group_Assignment, data=DataSetF)

#so there are now 8 models, rather than 4. This needs to be reflected in the AICc and model lists that are made to produce tables/results from AICc model selection and model averaging
#For the manual likelihood calculation, there will also be likelihoods for bodysize_pc_g_mm and EnergyPDry_1_mJ_g


# ================= AICc Analysis and Model Averaging Example for Fineness ~ Genetic Population Analysis =============
# All of the model analysis is the same, but the model set is slightly smaller and looks like this, so replace with respective model set and be sure to adjust the AICc analysis and model averaging accordingly
DataSet <- read.csv("EarlyMigration_GeneticPopulation_BodyCondition.csv")

DataSet$Group_Assignment <- as.factor(DataSet$Group_Assignment)
DataSet$Year<- as.factor(DataSet$Year)
DataSet$Sex<- as.factor(DataSet$Sex)
str(DataSet)

model1 <- lm(fineness~1, data=DataSet)
model2 <- lm(fineness~Collection_RiverMile_m + Sex, data = DataSet)
model3 <- lm(fineness~Collection_RiverMile_m + Group_Assignment, data = DataSet)
model4 <- lm(fineness~Collection_RiverMile_m + Sex + Group_Assignment, data = DataSet)
#so there are 4 models. Will follow the same structure as body condition analysis, adjust AICc model averaging and model lists as needed.


# ================= Final Analysis Notes =============
# Each AICc analysis for a response variable ~ Genetic Population will result in an AICc model selection results table as well as a table of model averaged results. These were used to interpet results, make conclusions, and produce figures for the manuscript


