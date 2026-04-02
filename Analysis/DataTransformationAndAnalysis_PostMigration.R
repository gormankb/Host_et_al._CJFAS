### Post Migration Spawning Population Data Transformation and AICc Analysis Code of
### Host et al. Escapement quality of sockeye salmon (Oncorhynchus nerka) returning to a glacially dominated watershed: relationships among body size, somatic energy, reproductive investment, and migratory difficulty
### Published with Canadian Journal of Fisheries and Aquatic Sciences
# April 2, 2026

#Important notes on the genetic population data: It only has data from 2019 and 2021. 

# ================= Body Size (PC1 Score) Data Transformation =============
#Clear Environment, load important packages
rm(list = ls())
library(AICcmodavg)
library(lubridate)

#download original data
OriginalData<-read.csv("PostMigration_BodyCondition.csv")

### Removing NAs that will mess up PCA code
any(is.na(OriginalData$Collection_Location)) #this is the genetic population assignment variable
any(is.na(OriginalData$Year))
any(is.na(OriginalData$Sex)) #if any NAs, remove them!


### Now, we can create the actual PC score for the body size that will be used in the AICc analysis ####
#OriginalData$Fish_Leng_1_mm 
#OriginalData$Fish_Ht_mm 
#OriginalData$Fish_Grth_mm 

bc_pca2 <- princomp(OriginalData[,c(69,71,72)], cor=T, scores=T, covmat = NULL) #should work now without any NAs, identical scores to that of the one above.
summary(bc_pca2, loadings=T, cutoff=0.0001) #summary of PC analysis
bc_pca2$scores
screeplot(bc_pca2, type=c('lines'))
bodysize_pc_g_mm <- bc_pca2$scores[,1] #creates a variable that is the body size index scores

#combine dataframe with column for body size PC scores
OriginalData <- cbind(OriginalData, bodysize_pc_g_mm) 


#Now, we have a dataset that can be saved as fully transformed and ready to be analyzed for body size AICc
write.csv(OriginalData, "PostMigration_BodyCondition.csv")



# ================= Energy Density Data Transformation =============
rm(list = ls()) #clear environment
OriginalData <-read.csv("PostMigration_BodyCondition.csv")

## Again, we need to check for the presence of NAs that will mess up levels / model weight
any(is.na(OriginalData$Collection_Location))
any(is.na(OriginalData$Year))
any(is.na(OriginalData$Sex)) #if any NAs, remove!

#transform energy density to correct units (mJ/g)
OriginalData$EnergyPDry_1_mJ_g <- OriginalData$EnergyPDry_1 / 1000

#save data as "transformed" 
write.csv(OriginalData, "PostMigration_BodyCondition.csv")



# ================= Total Energy Data Transformation =============
rm(list = ls()) #clear environment
OriginalData <-read.csv("PostMigration_BodyCondition.csv")

## Again, we need to check for the presence of NAs that will mess up levels / model weight
any(is.na(OriginalData$Collection_Location))
any(is.na(OriginalData$Year))
any(is.na(OriginalData$Sex)) #remove any NAs


#Calculate total energy variable (mJ)
OriginalData$EnergyPDry_1_mJ_g <- OriginalData$EnergyPDry_1 / 1000
OriginalData$Carcass.guts_Wt_g <- OriginalData$Carcass.guts_Wt * 1000
OriginalData$TotalEnergy_carcass_mJ <- OriginalData$EnergyPDry_1_mJ_g * OriginalData$Carcass.guts_Wt_g
str(OriginalData)

write.csv(OriginalData, "PostMigration_BodyCondition.csv")



# ================= Egg Size Data Transformation =============
rm(list = ls()) #clear environment
##for egg size, we will use previously transformed body size data, because the PC score will be an explanatory variable in the resulting AICc analysis
OriginalData <-read.csv("PostMigration_BodyCondition.csv")

## Removing NAs that could disrupt model analysis, check each explanatory variable
any(is.na(OriginalData$Collection_Location))
any(is.na(OriginalData$Year))
any(is.na(OriginalData$Sex)) #should all be False, as we previously checked for this and this is the transformed data

#subset data for females only
OriginalDataF <- subset(OriginalData, Sex == "F")

# convert gonad_wt to grams from kilograms
OriginalDataF$Gonad_Wt_g <- OriginalDataF$Gonad_Wt * 1000

OriginalDataF$Fem_Egg_Sz_1 #egg size metric

#Subset data to make sure all egg size values are above 0 and there is no NAs present in gonad weight or egg size
OriginalDataF <- subset(OriginalDataF, Fem_Egg_Sz_1 > 0 & !is.na(Fem_Egg_Sz_1) & !is.na(Gonad_Wt_g))


#save data as "transformed" 
write.csv(OriginalData, "PostMigration_EggSize.csv")




# ================= Fineness Ratio Data Transformation =============
rm(list = ls()) #clear environment
#for fineness, we will use previously transformed body size data
OriginalData <-read.csv("PostMigration_BodyCondition.csv")

#All data is already transformed, so there should be no NAs, but still check
## Removing NAs that could disrupt model analysis, check each explanatory variable
any(is.na(OriginalData$Collection_Location))
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
write.csv(OriginalData, "PostMigration_BodyCondition.csv")




# ================= AICc Analysis and Model Averaging Example for Response Variable ~ Spawning Population Analysis =============
# I will provide one full walk through example of how we completed our AICc analysis and resulting model average approach for this manuscript. This code can be reflected across all response variables. I will also show the model set for egg size, as it is slightly different 
### Data Analysis for energy density and total energy follow the same exact analysis as above, just replacing the transformed data sets with respective data sets and replacing body size with the response variable of interest
rm(list = ls()) #clear environment
#for fineness, we will use previously transformed body size data
DataSet <-read.csv("PostMigration_BodyCondition.csv")

DataSet$Collection_Location <- as.factor(DataSet$Collection_Location) #metric for spawning population
DataSet$Year<- as.factor(DataSet$Year)
DataSet$Sex<- as.factor(DataSet$Sex)
str(DataSet)

model1 <- lm(bodysize_pc_g_mm~1, data=DataSet)
model2 <- lm(bodysize_pc_g_mm~ Gonad_Wt_g + Collection_Location, data = DataSet)
model3 <- lm(bodysize_pc_g_mm~ Gonad_Wt_g + Sex, data = DataSet)
model4 <- lm(bodysize_pc_g_mm~ Gonad_Wt_g + Year, data = DataSet)
model5 <- lm(bodysize_pc_g_mm~ Gonad_Wt_g + Collection_Location + Sex, data = DataSet)
model6 <- lm(bodysize_pc_g_mm~ Gonad_Wt_g + Collection_Location + Year, data = DataSet)
model7 <- lm(bodysize_pc_g_mm~ Gonad_Wt_g + Sex + Year, data = DataSet)
model8 <- lm(bodysize_pc_g_mm~ Gonad_Wt_g + Collection_Location + Sex + Collection_Location:Sex, data = DataSet)
model9 <- lm(bodysize_pc_g_mm~Gonad_Wt_g + Sex + Year + Sex:Year, data = DataSet)
model10 <- lm(bodysize_pc_g_mm~ Gonad_Wt_g + Collection_Location + Sex + Year, data = DataSet)


#AICc model selection
library(AICcmodavg)
AIC(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10)

#table of AIC results
mynames1 <- paste("model", as.character(1:10), sep = "")
models <- list(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10)

# Generate AIC table
myaicc1 <- aictab(models, modnames = mynames1)
print(myaicc1)
# Convert AIC table to a data frame for easier manipulation
aic_df <- as.data.frame(myaicc1)

colnames(aic_df)[colnames(aic_df) == "Modnames"] <- "Model #" #can save this dataframe as a .csv now to look back on AICc results

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
  "Gonad_Wt_g" = 1.000000,
  "Collection_LocationGulkana Hatchery" = 0.41159041, 
  "Collection_LocationLong Lake" = 0.41159041,
  "Collection_LocationMahlo" = 0.41159041,
  "Collection_LocationMentasta" = 0.41159041,
  "Collection_LocationSt Anne" = 0.41159041,
  "Collection_LocationTanada" = 0.41159041,
  "Year2020" = 0.96892340,
  "Year2021" = 0.96892340, 
  "SexM" = 1.00000000,
  "SexM:Year2020" = 0.25346865,
  "SexM:Year2021" = 0.25346865,
  "Collection_LocationGulkana Hatchery:SexM" = 0.00140456,
  "Collection_LocationLong Lake:SexM" = 0.00140456,
  "Collection_LocationMahlo:SexM" = 0.00140456,
  "Collection_LocationMentasta:SexM" = 0.00140456,
  "Collection_LocationSt Anne:SexM" = 0.00140456,
  "Collection_LocationTanada:SexM" = 0.00140456
) #these numbers will all change according to the analysis, but this is the code used to create the table.



# Add ParamLikelihoods by matching rownames
coefs_df_full$ParamLikelihood <- manual_likelihoods[rownames(coefs_df_full)]

# View final table
print(coefs_df_full) #this dataframe can be saved as a .csv for future reference and further analysis

#reminder, total energy and energy density AICc analysis are the exact same as above, with the exception of replacing the response variable



# ================= AICc Analysis and Model Averaging Example for Egg Size ~ Spawning Population Analysis =============
# All of the model analysis is the same, but the model set is slightly larger and looks like this, so replace with respective model set and be sure to adjust the AICc analysis and model averaging accordingly
DataSetF <- read.csv("PostMigration_Fecundity.csv")

DataSet$Collection_Location <- as.factor(DataSet$Collection_Location)
DataSet$Year<- as.factor(DataSet$Year)
DataSet$Sex<- as.factor(DataSet$Sex)
str(DataSet)

#models for egg size
model1 <- lm(Fem_Egg_Sz_1~1, data=DataSetF)
model2 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g + bodysize_pc_g_mm, data=DataSetF)
model3 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g + EnergyPDry_1_mJ_g, data=DataSetF)
model4 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g + Collection_Location, data=DataSetF)
model5 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g + Year, data=DataSetF)
model6 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g + bodysize_pc_g_mm + EnergyPDry_1_mJ_g, data=DataSetF)
model7 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g + bodysize_pc_g_mm + Collection_Location, data=DataSetF)
model8 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g + bodysize_pc_g_mm + Year, data=DataSetF)
model9 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g + EnergyPDry_1_mJ_g + Collection_Location, data=DataSetF)
model10 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g + EnergyPDry_1_mJ_g*Collection_Location, data=DataSetF)
model11 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g + EnergyPDry_1_mJ_g + Year, data=DataSetF)
model12 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g + Collection_Location + Year, data=DataSetF)
model13 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g + bodysize_pc_g_mm + EnergyPDry_1_mJ_g + Collection_Location, data=DataSetF)
model14 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g + bodysize_pc_g_mm + EnergyPDry_1_mJ_g + Year, data=DataSetF)
model15 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g + bodysize_pc_g_mm + Collection_Location + Year, data=DataSetF)
model16 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g +  EnergyPDry_1_mJ_g + Collection_Location + Year, data=DataSetF)
model17 <- lm(Fem_Egg_Sz_1~Gonad_Wt_g + bodysize_pc_g_mm + EnergyPDry_1_mJ_g + Collection_Location + Year, data=DataSetF)

#so there are now 17 models, rather than 10. This needs to be reflected in the AICc and model lists that are made to produce tables/results from AICc model selection and model averaging
#For the manual likelihood calculation, there will also be likelihoods for bodysize_pc_g_mm and EnergyPDry_1_mJ_g


# ================= AICc Analysis and Model Averaging Example for Fineness ~ Spawning Population Analysis =============
# All of the model analysis is the same, but the model set is slightly smaller and looks like this, so replace with respective model set and be sure to adjust the AICc analysis and model averaging accordingly
DataSet <- read.csv("PostMigration_BodyCondition.csv")

DataSet$Collection_Location <- as.factor(DataSet$Collection_Location)
DataSet$Year<- as.factor(DataSet$Year)
DataSet$Sex<- as.factor(DataSet$Sex)
str(DataSet)

model1 <- lm(fineness~1, data=DataSet)
model2 <- lm(fineness~ Gonad_Wt_g + Collection_Location, data = DataSet)
model3 <- lm(fineness~ Gonad_Wt_g + Sex, data = DataSet)
model4 <- lm(fineness~ Gonad_Wt_g + Collection_Location + Sex, data = DataSet)
#so there are 4 models. AICc analysis and Model Averaging will follow the same structure as body condition analysis, adjust AICc model averaging and model lists as needed.


# ================= Final Analysis Notes =============
# Each AICc analysis for a response variable ~ Spawning Population will result in an AICc model selection results table as well as a table of model averaged results. These were used to interpet results, make conclusions, and produce figures for the manuscript.
