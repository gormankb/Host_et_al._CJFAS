### Early Migration Seasonal Run Timing Data Transformation and AICc Analysis Code
### Host et al. Escapement quality of sockeye salmon (Oncorhynchus nerka) returning to a glacially dominated watershed: relationships among body size, somatic energy, reproductive investment, and migratory difficulty
### Published with Canadian Journal of Fisheries and Aquatic Sciences
# March 23, 2026

#Important note on Seasonal Run Timing Raw Data: contains data from 2019-2021, raw data file is named "EarlyMigration_RunTimingGroup_RawData.csv"

# ================= Body Size (PC1 Score) Data Transformation =============
#Clear Environment, load important packages
rm(list = ls())
library(AICcmodavg)
library(lubridate)

#download original data
OriginalData<-read.csv("EarlyMigration_RunTimingGroup_BodyCondition.csv")
  #if you wish to skip this step, you can use the transformed data and jump straight to the AICc analysis

### Removing NAs that could disrupt model analysis, check each explanatory variable
any(is.na(OriginalData$Collection_Location)) 
any(is.na(OriginalData$RunTimingGroup))
any(is.na(OriginalData$Year))
any(is.na(OriginalData$Sex)) #remove any NAs present in these explanatory variables


### To create the actual PC score for the body size that will be used in the AICc analysis ####
## Convert leng_1, Ht, and Grth -- for bodysize_pc_g_mm
#OriginalData$Fish_Leng_1_mm 
#OriginalData$Fish_Ht_mm 
#OriginalData$Fish_Grth_mm 

bc_pca2 <- princomp(OriginalData[,c(70,71,72)], cor=T, scores=T, covmat = NULL) #PC score uses using leng_1, height, girth (double check column numbers, if need be)
summary(bc_pca2, loadings=T, cutoff=0.0001) #summary of PC analysis
bc_pca2$scores
screeplot(bc_pca2, type=c('lines'))
bodysize_pc_g_mm <- bc_pca2$scores[,1]

#combine dataframe with column for body size PC scores
OriginalData <- cbind(OriginalData, bodysize_pc_g_mm) 

#Now, we have a dataset that can be saved as fully transformed and ready to be analyzed for body size AICc
write.csv(OriginalData, "EarlyMigration_RunTimingGroup_BodyCondition.csv")




# ================= Energy Density Data Transformation =============
rm(list = ls()) #clear environment
OriginalData <-read.csv("EarlyMigration_RunTimingGroup_BodyCondition.csv")
#again, important to check for NAs
### Removing NAs that could disrupt model analysis, check each explanatory variable
any(is.na(OriginalData$Collection_Location)) 
any(is.na(OriginalData$RunTimingGroup))
any(is.na(OriginalData$Year))
any(is.na(OriginalData$Sex)) #remove any NAs if needed

#convert energy density to mJ/g
OriginalData$EnergyPDry_1_mJ_g <- OriginalData$EnergyPDry_1 / 1000

#save data as "transformed" 
write.csv(OriginalData, "EarlyMigration_RunTimingGroup_BodyCondition.csv")





# ================= Total Energy Data Transformation =============
rm(list = ls()) #clear environment
OriginalData <-read.csv("EarlyMigration_RunTimingGroup_BodyCondition.csv")
#again, important to check for NAs
### Removing NAs that could disrupt model analysis, check each explanatory variable
any(is.na(OriginalData$Collection_Location)) 
any(is.na(OriginalData$RunTimingGroup))
any(is.na(OriginalData$Year))
any(is.na(OriginalData$Sex)) #remove NAs present 

#need to add total energy (mJ) as a variable to this data set
# Total energy = PDry_1 * Fish_Wt == mJ energy
OriginalData$EnergyPDry_1_mJ_g <- OriginalData$EnergyPDry_1 / 1000
OriginalData$Carcass.guts_Wt_g <- OriginalData$Carcass.guts_Wt * 1000
OriginalData$TotalEnergy_carcass_g <- OriginalData$EnergyPDry_1_mJ_g * OriginalData$Carcass.guts_Wt_g

#save data as "transformed" 
write.csv(OriginalData, "EarlyMigration_RunTimingGroup_BodyCondition.csv")





# ================= Fecundity Data Transformation =============
rm(list = ls()) #clear environment
##for fecundity, we will use previously transformed body size data, because the PC score will be an explanatory variable in the resulting AICc analysis
OriginalData <-read.csv("EarlyMigration_RunTimingGroup_Fecundity.csv")
#again, important to check for NAs
### Removing NAs that could disrupt model analysis, check each explanatory variable
any(is.na(OriginalData$Collection_Location)) 
any(is.na(OriginalData$RunTimingGroup))
any(is.na(OriginalData$Year))
any(is.na(OriginalData$Sex)) #should all be False, as we previously checked for this

#subset data for females only
OriginalData <- subset(OriginalData, Sex == "F")

#calculate a fecundity metric 
  # fecundity = (# of eggs / weight of that sample in g) * (total gonad weight in kg * 1000)
OriginalData$Fecun_num_by_wt1 <- (OriginalData$Fem_Fecun_1_num/OriginalData$Fem_Fecun_1_wt)*(OriginalData$Gonad_Wt * 1000)

#save data as "transformed" 
write.csv(OriginalData, "EarlyMigration_RunTimingGroup_Fecundity.csv")




# ================= Fineness Ratio Data Transformation =============
rm(list = ls()) #clear environment
#for fineness, we will use previously transformed body size data
OriginalData <-read.csv("EarlyMigration_RunTimingGroup_BodyCondition.csv")
#OriginalData$Fish_Wdth_mm used for fineness

# Check NAs
any(is.na(OriginalData$Fish_Leng_1_mm)) #FALSE
any(is.na(OriginalData$Fish_Wdth_mm)) #FALSE

#Now, we can calculate fineness ratio
DataSet <- DataSet %>%
  group_by(Sex) %>%
  mutate(fineness = Fish_Leng_1_mm / Fish_Wdth_mm)

#save data as "transformed" 
write.csv(OriginalData, "EarlyMigration_RunTimingGroup_BodyCondition.csv")






# ================= AICc Analysis and Model Averaging Example for Response Variable ~ RunTimingGroup Analysis =============
# I will provide one full walk through example of how we completed our AICc analysis and resulting model average approach for this manuscript. This code can be reflected across all response variables. I will also show the model set for fecundity as well as fineness ratio, as they are slightly different 
### Data Analysis for energy density and total energy follow the same exact analysis as above, just replacing the transformed data sets with respective data sets and replacing body size with the response variable of interest
rm(list = ls()) #clear environment
DataSet <- read.csv("EarlyMigration_RuntimingGroup_BodyCondition.csv")

#### CHECK DIRECTION OF PC SCORE ###
library(ggplot2)
ggplot(DataSet, aes(x=Fish_Leng_1, y=bodysize_pc_g_mm)) +
  geom_point()
#so direction of PC score is positive, can check with other variables

# View DataSet.
head(DataSet)
summary(DataSet)
nrow(DataSet)
str(DataSet)

#make response variables factors
DataSet$Collection_Location <- as.factor(DataSet$Collection_Location)
DataSet$RunTimingGroup <- as.factor(DataSet$RunTimingGroup)
DataSet$Year<- as.factor(DataSet$Year)
DataSet$Sex<- as.factor(DataSet$Sex)
str(DataSet) #collection_rivermile_m is numeric


##Making Models for AICc Analysis, using supplementary material model sets
model1 <- lm(bodysize_pc_g_mm~1, data=DataSet)
model2 <- lm(bodysize_pc_g_mm~Collection_RiverMile_m + RunTimingGroup, data = DataSet)
model3 <- lm(bodysize_pc_g_mm~Collection_RiverMile_m + Sex, data = DataSet)
model4 <- lm(bodysize_pc_g_mm~Collection_RiverMile_m + Year, data = DataSet)
model5 <- lm(bodysize_pc_g_mm~Collection_RiverMile_m + RunTimingGroup + Sex, data = DataSet)
model6 <- lm(bodysize_pc_g_mm~Collection_RiverMile_m + RunTimingGroup + Year, data = DataSet)
model7 <- lm(bodysize_pc_g_mm~Collection_RiverMile_m + Sex + Year, data = DataSet)
model8 <- lm(bodysize_pc_g_mm~Collection_RiverMile_m + RunTimingGroup + Sex + RunTimingGroup:Sex, data = DataSet)
model9 <- lm(bodysize_pc_g_mm~Collection_RiverMile_m + Sex + Year + Sex:Year, data = DataSet)
model10 <- lm(bodysize_pc_g_mm~Collection_RiverMile_m + RunTimingGroup + Sex + Year, data = DataSet)


#AICc model selection
library(AICcmodavg)
AIC(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10)

#table of AIC results
mynames1 <- paste("model", as.character(1:10), sep = "")
models <- list(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10)

# Generate AIC table
myaicc1 <- aictab(models, modnames = mynames1)
print(myaicc1) #can save this as a .csv if you would like, to refer back to later
aic_df <- as.data.frame(myaicc1)

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

# Create tidy dataframe with full averages
coefs_df_full <- as.data.frame(avg_coefs_full)
colnames(coefs_df_full) <- "Avg_Coef"
coefs_df_full$CI_Lower <- avg_confint_full[, 1]
coefs_df_full$CI_Upper <- avg_confint_full[, 2]
coefs_df_full$SE <- avg_se_full

# Calculate parameter importance (parameter likelihoods)
param_likelihoods <- sw(model_avg)

# Manually specify parameter likelihoods from sw(model_avg)
manual_likelihoods <- c(
  "(Intercept)" = 1.0,
  "Collection_RiverMile_m" = 0.9999989367,
  "RunTimingGroupLate" = 0.9943143436,
  "RunTimingGroupMiddle" = 0.9943143436,
  "SexM" = 0.9938102633,
  "Year2020" = 0.6518253418,
  "Year2021" = 0.6518253418,
  "RunTimingGroupLate:SexM" = 0.0745580097,
  "RunTimingGroupMiddle:SexM" = 0.0745580097,
  "SexM:Year2020" = 0.0005110073,
  "SexM:Year2021" = 0.0005110073 #these numbers will all change according to the analysis, but this is the code used to create the table.
)

# Add ParamLikelihoods by matching rownames
coefs_df_full$ParamLikelihood <- manual_likelihoods[rownames(coefs_df_full)]

# View final table
print(coefs_df_full) #this data frame can be saved as a .csv for future reference and further analysis



# ================= AICc Analysis and Model Averaging Example for Fecundity ~ RunTimingGroup Analysis =============
# All of the model analysis is the same, but the model set is slightly larger and looks like this, so replace with respective model set and be sure to adjust the AICc analysis and model averaging accordingly
DataSetF <- read.csv("EarlyMigration_RunTimingGroup_Fecundity.csv")

# View DataSet.
head(DataSetF)
summary(DataSetF)
nrow(DataSetF)
str(DataSetF)

#make response variables factors
DataSetF$RunTimingGroup <- as.factor(DataSetF$RunTimingGroup)
DataSetF$Year<- as.factor(DataSetF$Year)

#model set for Fecundity AICc
model1 <- lm(Fecun_num_by_wt1~1, data=DataSetF)
model2 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + bodysize_pc_g_mm, data=DataSetF)
model3 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + EnergyPDry_1_mJ_g, data=DataSetF)
model4 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + RunTimingGroup, data=DataSetF)
model5 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + Year, data=DataSetF)
model6 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + bodysize_pc_g_mm + EnergyPDry_1_mJ_g, data=DataSetF)
model7 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + bodysize_pc_g_mm + RunTimingGroup, data=DataSetF)
model8 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + bodysize_pc_g_mm + Year, data=DataSetF)
model9 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + EnergyPDry_1_mJ_g + RunTimingGroup, data=DataSetF)
model10 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + EnergyPDry_1_mJ_g*RunTimingGroup, data=DataSetF)
model11 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + EnergyPDry_1_mJ_g + Year, data=DataSetF)
model12 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + RunTimingGroup + Year, data=DataSetF)
model13 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + bodysize_pc_g_mm + EnergyPDry_1_mJ_g + RunTimingGroup, data=DataSetF)
model14 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + bodysize_pc_g_mm + EnergyPDry_1_mJ_g + Year, data=DataSetF)
model15 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + bodysize_pc_g_mm + RunTimingGroup + Year, data=DataSetF)
model16 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m +  EnergyPDry_1_mJ_g + RunTimingGroup + Year, data=DataSetF)
model17 <- lm(Fecun_num_by_wt1~Collection_RiverMile_m + bodysize_pc_g_mm + EnergyPDry_1_mJ_g + RunTimingGroup + Year, data=DataSetF)

#so there are now 17 models, rather than 10. This needs to be reflected in the AICc and model lists that are made to produce tables/results from AICc model selection and model averaging



# ================= AICc Analysis and Model Averaging Example for Fineness ~ RunTimingGroup Analysis =============
# All of the model analysis is the same, but the model set is slightly smaller and looks like this, so replace with respective model set and be sure to adjust the AICc analysis and model averaging accordingly
DataSet <- read.csv("EarlyMigration_RunTimingGroup_BodyCondition.csv")

DataSet$Collection_Location <- as.factor(DataSet$Collection_Location)
DataSet$RunTimingGroup <- as.factor(DataSet$RunTimingGroup)
DataSet$Sex<- as.factor(DataSet$Sex)

model1 <- lm(fineness~1, data=DataSet)
model2 <- lm(fineness~Collection_RiverMile_m + RunTimingGroup, data = DataSet)
model3 <- lm(fineness~Collection_RiverMile_m + Sex, data = DataSet)
model4 <- lm(fineness~Collection_RiverMile_m + RunTimingGroup + Sex, data = DataSet)

#so there is 4 models now rather than 10. This should be reflected in the AICc and model lists that are made to produce tables/results from AICc model selection and model averaging



# ================= Final Analysis Notes =============
# Each AICc analysis for a response variable ~ RunTimingGroup will result in an AICc model selection results table as well as a table of model averaged results. These were used to interpet results, make conclusions, and produce figures for the manuscript

