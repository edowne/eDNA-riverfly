# BIOS0034 Research Project
# Student number: HVVY6
# 23/08/2024

# Data Analysis Part 3: Exploration of relationships and predictive model ####


# Set-up ####

# Clearing environment
rm(list=ls())

# Setting working directory: 
setwd("C:/Users/elisa/OneDrive - University College London/BIOS0034 Research Project/MY DATA/Analysis")

getwd()

# Loading packages 
library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(stringr)
library(drc)
library(scales)
library(lattice)
library(cowplot)
library(GGally)
library(lme4)

#Load packages
library(arm)
library(car)
library(ggbeeswarm)
library(ggplot2)
library(lattice)
library(lawstat)
library(outliers)

# Reading in dataset:
data_0 <- read.csv("COMP_Proj_Data.csv")

# Checking dataset: 
str(data_0)

# Preparing dataset to include only variables of interest: 
data_env <- data_0 %>% dplyr::select(Sample_Ref, Date_Collected, Water_Temp_DegC, pH, ORP_mV,
                                   Conductivity, Phos_AV, Nitr_AV, DepthMean_cm, Met_Temp_DegC,
                                   Met_WindSpd_mph, Met_WindDir, Met_UV, Met_Pres_mb, Met_Hum_.,
                                   Sample_Type, Date_Extracted, PCR_Plate, PCR_Code, Adj_conc, KS2_HOG_Count)


# Filtering out samples that weren't analysed: 
data_env <- data_env %>% dplyr::filter(PCR_Code!="")

# Checking dataset structure: 
str(data_env)

# Make sure categorical variables are set as factors: 
data_env$Sample_Type <- as.factor(data_env$Sample_Type) # sample type
data_env$Met_WindDir <- as.factor(data_env$Met_WindDir) # Wind direction
data_env$Met_UV <- as.factor(data_env$Met_UV) # UV level     
data_env$PCR_Plate <- as.factor(data_env$PCR_Plate) # PCR plate
data_env$PCR_Code <- as.factor(data_env$PCR_Code) # PCR sample code

# Convert dates to date format: 
data_env <- data_env %>% mutate(Date_Collected = as.Date(Date_Collected, format = "%d/%m/%Y")) %>% 
  mutate(Date_Extracted = as.Date(Date_Extracted, format = "%d/%m/%Y"))
# Check: 
str(data_env)

# Convert "NA" concentrations to zero, to indicate no target DNA detected in sample: 
data_env <- data_env %>% mutate(Adj_conc = ifelse(is.na(Adj_conc), 0, Adj_conc))

# Filter out controls - I am only looking at the water samples in this analysis: 
data_env <- data_env %>% filter(Sample_Type %in% c("Edge_Depth", "Edge_Surface", 
                                               "Transect_Depth"))


# DATA EXPLORATION ####

# Steps were followed from Smith C, Uzal A & Warren M (2020) "Statistics in R for 
# Biodiversity Conservation". Amazon, UK, with data exploration steps based on 
# recommended protocol from Zuur AF, Ieno EN & Elphick CS (2010) A protocol for data
# exploration to avoid common statistical problems. "Methods in Ecology and Evolution"
# 1, 3-14.

# Initially a Gaussian GLM model was selected because GLMs can be used to predict 
# a continuous and normally distributed response variable (hoglouse counts) from one
# or more independent variables, which can be continuous or categorical. 

# Aims of the analysis: 
  # The aims of this analysis were to test whether hoglouse abundance as indicated 
  # by qPCR DNA concentration was: 
    # (1) Positively associated with hoglouse abundance from kick-sampling count data
    # (2) Associated with any of the environmental variables measured 
    # (3) Associated with any technical variables (extraction date and PCR plate)
  # with the dependent variable being disaggregated by sampling method. 

# 1. DATA EXPLORATION ####

# 1(a) Check for missing data ####

# Screening for NAs in the dataset:
colSums(is.na(data_env))
# 5 pH measurements are missing (due to probes not available in September, and pH
      # probe not working in January)
# 2 ORP and conductivity measurements are missing (due to probes not available in 
      # September)
# 2 phosphates and nitrates measurements are missing (due to tests not being carried
      # out in September)


# 1(b) Check for outliers ####

# Use Cleveland plots to look for outliers in continuous variables:
Names <- c("Water_Temp_DegC", "pH", "ORP_mV", "Conductivity", "Phos_AV", "Nitr_AV",
           "DepthMean_cm", "Met_Temp_DegC", "Met_WindSpd_mph", "Met_Pres_mb", 
           "Met_Hum_.", "Adj_conc", "KS2_HOG_Count")
dotplot(as.matrix(as.matrix(data_env[,Names])),
        groups=FALSE,
        strip = strip.custom(bg = 'white',
                             par.strip.text = list(cex = 1)),
        scales = list(x = list(relation = "free", draw = TRUE),
                      y = list(relation = "free", draw = FALSE)),
        col = "#228886", cex  = 1, pch = 20,
        xlab = list(label = "Data range (value of the variable)", cex = 1.2),
        ylab = list(label = "Order of the data", 
                    cex = 1.2))

# Cleveland dotplots don't indicate any clear outliers in the continuous variables. 

# Use boxplots to check for outliers in categorical variables: 
par(mfrow = c(1,2), mar=c(5,5,2,2), cex.lab=1.5)
# Drop unused levels from the Sample_Type factor
data_env$Sample_Type <- droplevels(data_env$Sample_Type)
boxplot(Adj_conc~Sample_Type, 
        ylab="DNA concentration", 
        xlab="Sampling method",
        pch=16, cex=1, col="#33CCCC", 
        data=data_env)
# Possible outlier in Edge-Surface - already knew about this - but not sure it makes 
# sense to remove it for this analysis. 

boxplot(Adj_conc~Met_WindDir, 
        ylab="DNA concentration", 
        xlab="Wind direction",
        pch=16, cex=1, col="#33CCCC", 
        data=data_env)

boxplot(Adj_conc~Met_UV, 
        ylab="DNA concentration", 
        xlab="UV level",
        pch=16, cex=1, col="#33CCCC", 
        data=data_env)

boxplot(Adj_conc~PCR_Plate, 
        ylab="DNA concentration", 
        xlab="PCR plate number",
        pch=16, cex=1, col="#33CCCC", 
        data=data_env)

# Potential outlier in PCR plate number, but doesn't make sense to remove it. 


# 1(c) Check for data balance in categorical variables ####

table(data_env$Sample_Type) # Sample Types are well balanced (7, 6, 7)
table(data_env$Met_WindDir) # Wind Direction is fairly well balanced (6, 5, 3, 6)
table(data_env$Met_UV) # UV levels are not well balanced (3 high, 11 medium, 6 low)
                      # but do reflect normal distribution and expected for an 
                      # environmental variable 
table(data_env$PCR_Plate) # Data is balanced over the three PCR plates (8,6,6)
table(data_env$Date_Collected) # Data is balanced over sample collection date
table(data_env$Date_Extracted) # Data is fairly balanced over extraction date, 
                      # although with higher number (6) extracted on two dates


# There are a lot of variables here, so care must be taken fitting the data to 
# such a complex model. 


# 1(d) Check for lots of zeros in the response variable ####

# Calculate % of zeros in qPCR DNA conc: 
sum(data_env$Adj_conc == 0)*100/nrow(data_env)

# 55% of samples did not detect the target DNA in the qPCR experiment. This figure
# seems high and might cause problems with the analysis. This can be addressed through
# additional statistical processes. 


# BUT qPCR data will not be the response variable - count data will be: 
# Calculate % of zeros in kick-sampling count data: 
sum(data_env$KS2_HOG_Count == 0)*100/nrow(data_env)

# No zeros, so there should not be a problem with zero inflation in model. 


# 1(e) Check for multicollinearity among environmental covariates ####

# Visualise multicollinearity using a correlation matrix with corresponding pairplots: 
Coll <- c("Adj_conc", "Water_Temp_DegC", "pH", "ORP_mV", "Conductivity", "Phos_AV",
          "Nitr_AV", "DepthMean_cm", "Met_Temp_DegC", "Met_WindSpd_mph", "Met_WindDir", 
          "Met_UV", "Met_Pres_mb", "Met_Hum_.", "Sample_Type")

panel.cor <- function(x, y, digits=1, prefix="", cex.cor = 6)
{usr <- par("usr"); on.exit(par(usr))
par(usr = c(0, 1, 0, 1))
r1 = cor(x,y,use="pairwise.complete.obs")
r <- abs(cor(x, y,use="pairwise.complete.obs"))
txt <- format(c(r1, 0.1), digits=digits)[1]
txt <- paste(prefix, txt, sep="")
if(missing(cex.cor)) { cex <- 0.6/strwidth(txt) } else {
  cex = cex.cor}
text(0.5, 0.5, txt, cex = cex * r)}
pairs(data_env[, Coll], lower.panel = panel.cor, cex.labels = 0.5)

# Pairplots indicate strong collinearity between several variables, most notably: 
  # Air temperature & humidity (-1.0)
  # Water temperature & Air temperature (0.9)
  # Water temperature & humidity (-0.8)
  # ORP and air pressure (-0.8) 
  # pH and nitrates (0.8)

# I will remove Met_Temp_DegC (air temp), humidity and pressure as strongly correlated 
# with other variables

# Remove the covariates with missing values (pH, ORP, Conductivity, phosphates, nitrates)

# Check variance inflation for the remaining covariates using the 'vif' function 
# from the 'car' package: 
# options(digits=3)
# car::vif(glm(Adj_conc ~ Water_Temp_DegC + #pH + ORP_mV + Conductivity + Phos_AV + Nitr_AV +
  #      DepthMean_cm + Met_WindSpd_mph + Met_WindDir +
  #      #Met_Pres_mb + Met_Hum_. + Met_Temp_DegC +
  #      Sample_Type,
  #      family = quasipoisson,
  #      data = data_env))

# Output of VIF test shows significant problems with multicollinearity, indicating 
# the model might not be appropriate. 


# 1(f) Check for normality of the response variable #### 

# For Adj_conc: 

# Histogram:
data_env %>% ggplot(aes(x=Adj_conc)) +
  geom_histogram(fill="#33CCCC", colour="white")+
  theme_bw()

# Shapiro-Wilk test for normality: 
shapiro.test(data_env$Adj_conc) # p<0.01 Data is significantly skewed, i.e. not normally distributed. 

# QQplot: 
qqnorm(data_env$Adj_conc)
qqline(data_env$Adj_conc) # Broadly follows the line

# Conclusion is that the Adj_conc data is not normally distributed as it is skewed 
# (and has already been log transformed). 


# For hoglouse counts: 
# Histogram: 
data_env %>% ggplot(aes(x=KS2_HOG_Count))+
  geom_histogram(fill="#F4B183", colour="white")+
  theme_bw()

# Shapiro-Wilk test for normality: 
shapiro.test(data_env$KS2_HOG_Count) # p>0.05, Data is not significantly different to 
# normal distribution.  

# QQplot: 
qqnorm(data_env$KS2_HOG_Count)
qqline(data_env$KS2_HOG_Count) # Far from perfect normal distribution but does seem to 
# broadly follow the normal line 


# 1(g) Check for homogeneity of variance ####

# For Adj_conc: 

# Boxplots with each grouping variable:
# Sample Type: 
data_env %>% ggplot(aes(x=Sample_Type, y=Adj_conc, fill=Sample_Type))+
  geom_boxplot()+
  theme_bw(base_size = 16)+
  scale_fill_viridis(discrete=TRUE, option="viridis", alpha=0.6, name="Sample type")+
  geom_point(data = data_env, aes(x = Sample_Type, y = Adj_conc), 
             pch = 23, colour = "black", fill="red", size = 3)+
  xlab("Sample type")+
  ylab("Relative concentration of target DNA detected by qPCR")

kruskal.test(Adj_conc~Sample_Type, data_env) # p>0.05, no sig differences between sample types


# Wind direction: 
data_env %>% ggplot(aes(x=Met_WindDir, y=Adj_conc, fill=Met_WindDir))+
  geom_boxplot()+
  theme_bw(base_size = 16)+
  scale_fill_viridis(discrete=TRUE, option="viridis", alpha=0.6, name="Wind direction")+
  geom_point(data = data_env, aes(x = Met_WindDir, y = Adj_conc), 
             pch = 23, colour = "black", fill="red", size = 3)+
  xlab("Wind direction")+
  ylab("Relative concentration of target DNA detected by qPCR")

kruskal.test(Adj_conc~Met_WindDir, data_env) # p>0.05, no sig differences between wind direction

# UV levels:
data_env %>% ggplot(aes(x=Met_UV, y=Adj_conc, fill=Met_UV))+
  geom_boxplot()+
  theme_bw(base_size = 16)+
  scale_fill_viridis(discrete=TRUE, option="viridis", alpha=0.6, name="UV level")+
  geom_point(data = data_env, aes(x = Met_UV, y = Adj_conc), 
             pch = 23, colour = "black", fill="red", size = 3)+
  xlab("UV level")+
  ylab("Relative concentration of target DNA detected by qPCR")

kruskal.test(Adj_conc~Met_UV, data_env) # Not significant

# PCR plate: 
data_env %>% ggplot(aes(x=PCR_Plate, y=Adj_conc, fill=PCR_Plate))+
  geom_boxplot()+
  theme_bw(base_size = 16)+
  scale_fill_viridis(discrete=TRUE, option="viridis", alpha=0.6, name="PCR plate no.")+
  geom_point(data = data_env, aes(x = PCR_Plate, y = Adj_conc), 
             pch = 23, colour = "black", fill="red", size = 3)+
  xlab("PCR plate no.")+
  ylab("Relative concentration of target DNA detected by qPCR")

kruskal.test(Adj_conc~PCR_Plate, data_env) # Not significant. 


# Scatterplots with continuous variables: 

# Calculating means for each continuous variable: 
Mean_Adj_conc <- mean(data_env$Adj_conc, na.rm=TRUE) # Response variable
Mean_WaterTemp <- mean(data_env$Water_Temp_DegC, na.rm=TRUE)
Mean_pH <- mean(data_env$pH, na.rm=TRUE)
Mean_ORP <- mean(data_env$ORP_mV, na.rm=TRUE)
Mean_cond <- mean(data_env$Conductivity, na.rm=TRUE)
Mean_phos <- mean(data_env$Phos_AV, na.rm=TRUE)
Mean_nitr <- mean(data_env$Nitr_AV, na.rm=TRUE)
Mean_depth <- mean(data_env$DepthMean_cm, na.rm=TRUE)
Mean_Airtemp <- mean(data_env$Met_Temp_DegC, na.rm=TRUE)
Mean_Windspd <- mean(data_env$Met_WindSpd_mph, na.rm=TRUE)
Mean_Pres <- mean(data_env$Met_Pres_mb, na.rm=TRUE)
Mean_Hum <- mean(data_env$Met_Hum_., na.rm=TRUE)
Mean_HogCount <- mean(data_env$KS2_HOG_Count, na.rm=TRUE)

scatter_watertemp <- data_env %>% ggplot(aes(x=Water_Temp_DegC, y=Adj_conc))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_WaterTemp, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_Adj_conc, linetype="dashed", col="#C45A12")+
  xlab("Water temperature (degrees Celsius")+
  ylab("Relative concentration of target DNA detected by qPCR")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_watertemp

scatter_pH <- data_env %>% ggplot(aes(x=pH, y=Adj_conc))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_pH, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_Adj_conc, linetype="dashed", col="#C45A12")+
  xlab("Water pH")+
  ylab("Relative concentration of target DNA detected by qPCR")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_pH

scatter_ORP <- data_env %>% ggplot(aes(x=ORP_mV, y=Adj_conc))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_ORP, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_Adj_conc, linetype="dashed", col="#C45A12")+
  xlab("Water oxidative reductive potential (mV)")+
  ylab("Relative concentration of target DNA detected by qPCR")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_ORP

scatter_cond <- data_env %>% ggplot(aes(x=Conductivity, y=Adj_conc))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_cond, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_Adj_conc, linetype="dashed", col="#C45A12")+
  xlab("Water conductivity (µS/cm)")+
  ylab("Relative concentration of target DNA detected by qPCR")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_cond

scatter_phos <- data_env %>% ggplot(aes(x=Phos_AV, y=Adj_conc))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_phos, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_Adj_conc, linetype="dashed", col="#C45A12")+
  xlab("Phosphate levels (mg/L PO[4]^3-)") +
  ylab("Relative concentration of target DNA detected by qPCR")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_phos

scatter_nitr <- data_env %>% ggplot(aes(x=Nitr_AV, y=Adj_conc))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_nitr, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_Adj_conc, linetype="dashed", col="#C45A12")+
  xlab("Nitrate levels (mg/L NO[2]-[-N])") +
  ylab("Relative concentration of target DNA detected by qPCR")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_nitr

scatter_depth <- data_env %>% ggplot(aes(x=DepthMean_cm, y=Adj_conc))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_depth, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_Adj_conc, linetype="dashed", col="#C45A12")+
  xlab("Water depth (cm)") +
  ylab("Relative concentration of target DNA detected by qPCR")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_depth

scatter_airtemp <- data_env %>% ggplot(aes(x=Met_Temp_DegC, y=Adj_conc))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_Airtemp, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_Adj_conc, linetype="dashed", col="#C45A12")+
  xlab("Air temperature (degrees Celsius)") +
  ylab("Relative concentration of target DNA detected by qPCR")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_airtemp

scatter_windspd <- data_env %>% ggplot(aes(x=Met_WindSpd_mph, y=Adj_conc))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_Windspd, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_Adj_conc, linetype="dashed", col="#C45A12")+
  xlab("Wind speed (mph)") +
  ylab("Relative concentration of target DNA detected by qPCR")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_windspd

scatter_pres <- data_env %>% ggplot(aes(x=Met_Pres_mb, y=Adj_conc))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_Pres, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_Adj_conc, linetype="dashed", col="#C45A12")+
  xlab("Pressure (mb)") +
  ylab("Relative concentration of target DNA detected by qPCR")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_pres

scatter_hum <- data_env %>% ggplot(aes(x=Met_Hum_., y=Adj_conc))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_Hum, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_Adj_conc, linetype="dashed", col="#C45A12")+
  xlab("Humidity (%)") +
  ylab("Relative concentration of target DNA detected by qPCR")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_hum

scatter_HogCount <- data_env %>% ggplot(aes(x=KS2_HOG_Count, y=Adj_conc))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_HogCount, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_Adj_conc, linetype="dashed", col="#C45A12")+
  xlab("Hoglouse count from kick-sampling") +
  ylab("Relative concentration of target DNA detected by qPCR")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_HogCount

#combining the scatters onto one plot:
scatter_combo <- plot_grid(scatter_airtemp, scatter_cond, scatter_depth, scatter_hum,
                           scatter_nitr, scatter_ORP, scatter_pH, scatter_phos,
                           scatter_pres, scatter_watertemp, scatter_windspd, scatter_HogCount, 
                           align = 'v', labels = 'AUT0')
scatter_combo  

# Findings: Quite a bit of variance around the means, but very few data points so 
# hard to draw any conclusions. 

# Repeating for kick-sampling count data as response variable: 


scatter_watertemp2 <- data_env %>% ggplot(aes(x=Water_Temp_DegC, y=KS2_HOG_Count))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_WaterTemp, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_HogCount, linetype="dashed", col="#C45A12")+
  xlab("Water temperature (degrees Celsius")+
  ylab("Hoglouse counts from kick-sampling")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_watertemp2

scatter_pH2 <- data_env %>% ggplot(aes(x=pH, y=KS2_HOG_Count))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_pH, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_HogCount, linetype="dashed", col="#C45A12")+
  xlab("Water pH")+
  ylab("Hoglouse counts from kick-sampling")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_pH2

scatter_ORP2 <- data_env %>% ggplot(aes(x=ORP_mV, y=KS2_HOG_Count))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_ORP, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_HogCount, linetype="dashed", col="#C45A12")+
  xlab("Water oxidative reductive potential (mV)")+
  ylab("Hoglouse counts from kick-sampling")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_ORP2

scatter_cond2 <- data_env %>% ggplot(aes(x=Conductivity, y=KS2_HOG_Count))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_cond, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_HogCount, linetype="dashed", col="#C45A12")+
  xlab("Water conductivity (µS/cm)")+
  ylab("Hoglouse counts from kick-sampling")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_cond2

scatter_phos2 <- data_env %>% ggplot(aes(x=Phos_AV, y=KS2_HOG_Count))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_phos, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_HogCount, linetype="dashed", col="#C45A12")+
  xlab("Phosphate levels (mg/L PO[4]^3-)") +
  ylab("Hoglouse counts from kick-sampling")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_phos2

scatter_nitr2 <- data_env %>% ggplot(aes(x=Nitr_AV, y=KS2_HOG_Count))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_nitr, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_HogCount, linetype="dashed", col="#C45A12")+
  xlab("Nitrate levels (mg/L NO[2]-[-N])") +
  ylab("Hoglouse counts from kick-sampling")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_nitr2

scatter_depth2 <- data_env %>% ggplot(aes(x=DepthMean_cm, y=KS2_HOG_Count))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_depth, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_HogCount, linetype="dashed", col="#C45A12")+
  xlab("Water depth (cm)") +
  ylab("Hoglouse counts from kick-sampling")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_depth2

scatter_airtemp2 <- data_env %>% ggplot(aes(x=Met_Temp_DegC, y=KS2_HOG_Count))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_Airtemp, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_HogCount, linetype="dashed", col="#C45A12")+
  xlab("Air temperature (degrees Celsius)") +
  ylab("Hoglouse counts from kick-sampling")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_airtemp2

scatter_windspd2 <- data_env %>% ggplot(aes(x=Met_WindSpd_mph, y=KS2_HOG_Count))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_Windspd, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_HogCount, linetype="dashed", col="#C45A12")+
  xlab("Wind speed (mph)") +
  ylab("Hoglouse counts from kick-sampling")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_windspd2

scatter_pres2 <- data_env %>% ggplot(aes(x=Met_Pres_mb, y=KS2_HOG_Count))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_Pres, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_HogCount, linetype="dashed", col="#C45A12")+
  xlab("Pressure (mb)") +
  ylab("Hoglouse counts from kick-sampling")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_pres2

scatter_hum2 <- data_env %>% ggplot(aes(x=Met_Hum_., y=KS2_HOG_Count))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_Hum, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_HogCount, linetype="dashed", col="#C45A12")+
  xlab("Humidity (%)") +
  ylab("Hoglouse counts from kick-sampling")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_hum2

scatter_qPCR <- data_env %>% ggplot(aes(x=Adj_conc, y=KS2_HOG_Count))+
  geom_point(col="#228886")+
  theme_bw()+
  geom_vline(xintercept=Mean_Adj_conc, linetype="dashed", col="#C45A12")+
  geom_hline(yintercept=Mean_HogCount, linetype="dashed", col="#C45A12")+
  ylab("Hoglouse count from kick-sampling") +
  xlab("Relative concentration of target DNA detected by qPCR")+
  theme(plot.title = element_text(hjust = 0.5))
scatter_qPCR

#combining the scatters onto one plot:
scatter_combo <- plot_grid(scatter_airtemp2, scatter_cond2, scatter_depth2, scatter_hum2,
                           scatter_nitr2, scatter_ORP2, scatter_pH2, scatter_phos2,
                           scatter_pres2, scatter_watertemp2, scatter_windspd2, scatter_qPCR, 
                           align = 'v', labels = 'AUT0')
scatter_combo  

# Findings: Quite a bit of variance around the means, but very few data points so 
# hard to draw any conclusions. Broadly even. 


# 1(f) Check for relationships among dependent and independent variables ####

# Since such a high degree of multicollinearity in previous step, I'm only going to 
# include water temperature (as a broad indicator of all weather conditions), pH (
# as a general indicator of water quality/chemical conditions), and depth as my 
# environmental variables

# Looking at Adj_conc as dependent variable: 

par(mfrow=c(2,2), mar=c(5,5,1,1), cex.lab = 1)
plot(Adj_conc ~ Water_Temp_DegC,  data = data_env,
     xlab = "Water temp (deg C)", 
     ylab = "Hoglouse eDNA concentration (relative qty, log10)",
     pch = 16, cex = 1.3)
plot(Adj_conc ~ pH,  data = data_env,
     xlab = "Water pH", 
     ylab = "Hoglouse eDNA concentration (relative qty, log10)",
     pch = 16, cex = 1.3)
plot(Adj_conc ~ ORP_mV,  data = data_env,
     xlab = "Water ORP", 
     ylab = "Hoglouse eDNA concentration (relative qty, log10)",
     pch = 16, cex = 1.3)
plot(Adj_conc ~ Conductivity,  data = data_env,
     xlab = "Water conductivity", 
     ylab = "Hoglouse eDNA concentration (relative qty, log10)",
     pch = 16, cex = 1.3)
plot(Adj_conc ~ Phos_AV,  data = data_env,
     xlab = "Phosphate level", 
     ylab = "Hoglouse eDNA concentration (relative qty, log10)",
     pch = 16, cex = 1.3)
plot(Adj_conc ~ Nitr_AV,  data = data_env,
     xlab = "Nitrate level", 
     ylab = "Hoglouse eDNA concentration (relative qty, log10)",
     pch = 16, cex = 1.3)
plot(Adj_conc ~ DepthMean_cm,  data = data_env,
     xlab = "Mean depth (cm)", 
     ylab = "Hoglouse eDNA concentration (relative qty, log10)",
     pch = 16, cex = 1.3)
plot(Adj_conc ~ KS2_HOG_Count,  data = data_env,
     xlab = "Hoglouse counts from kick-sampling", 
     ylab = "Hoglouse eDNA concentration (relative qty, log10)",
     pch = 16, cex = 1.3)

# Drop unused levels from the Sample_Type factor
data_env$Sample_Type <- droplevels(data_env$Sample_Type)
boxplot(Adj_conc ~ Sample_Type, data = data_env, 
        xlab = "Sampling method", 
        ylab = "Hoglouse eDNA concentration (relative qty, log10)",
        pch = 16, cex = 1.3, col = "#33CCCC")
boxplot(Adj_conc ~ Date_Extracted, data = data_env, 
        xlab = "DNA extraction date", 
        ylab = "Hoglouse eDNA concentration (relative qty, log10)",
        pch = 16, cex = 1.3, col = "#F4B183")
boxplot(Adj_conc ~ PCR_Plate, data = data_env, 
        xlab = "PCR plate number", 
        ylab = "Hoglouse eDNA concentration (relative qty, log10)",
        pch = 16, cex = 1.3, col = "#5B9BD5")

# No clear relationships emerging; however, if the zeros are ignored, then 
# some relationships might be present. 


# Looking at hoglouse counts as dependent variable: 

par(mfrow=c(2,2), mar=c(5,5,1,1), cex.lab = 1)
plot(KS2_HOG_Count ~ Water_Temp_DegC,  data = data_env,
     xlab = "Water temp (deg C)", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, cex = 1.3)
plot(KS2_HOG_Count ~ pH,  data = data_env,
     xlab = "Water pH", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, cex = 1.3)
plot(KS2_HOG_Count ~ ORP_mV,  data = data_env,
     xlab = "Water ORP", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, cex = 1.3)
plot(KS2_HOG_Count ~ Conductivity,  data = data_env,
     xlab = "Water conductivity", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, cex = 1.3)
plot(KS2_HOG_Count ~ Phos_AV,  data = data_env,
     xlab = "Phosphate level", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, cex = 1.3)
plot(KS2_HOG_Count ~ Nitr_AV,  data = data_env,
     xlab = "Nitrate level", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, cex = 1.3)
plot(KS2_HOG_Count ~ DepthMean_cm,  data = data_env,
     xlab = "Mean depth (cm)", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, cex = 1.3)
plot(KS2_HOG_Count ~ Adj_conc,  data = data_env,
     xlab = "Relative DNA concentration from qPCR", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, cex = 1.3)

# Drop unused levels from the Sample_Type factor
data_env$Sample_Type <- droplevels(data_env$Sample_Type)
boxplot(KS2_HOG_Count ~ Sample_Type, data = data_env, 
        xlab = "Sampling method", 
        ylab = "Hoglouse counts from kick-sampling",
        pch = 16, cex = 1.3, col = "#33CCCC")
boxplot(KS2_HOG_Count ~ Date_Extracted, data = data_env, 
        xlab = "DNA extraction date", 
        ylab = "Hoglouse counts from kick-sampling",
        pch = 16, cex = 1.3, col = "#F4B183")
boxplot(KS2_HOG_Count ~ PCR_Plate, data = data_env, 
        xlab = "PCR plate number", 
        ylab = "Hoglouse counts from kick-sampling",
        pch = 16, cex = 1.3, col = "#5B9BD5")

# Indication of a negative relationship between the hoglouse counts and water 
# temperature. No other clear relationships emerging; however, if the zero DNA
# detections are ignored, then there may be a positive relationship there. 


# Re-running with no absence data ####

# If we re-run the plots after removing Adj_conc == 0, what relationships emerge? 

# Convert "0" concentrations back to "NA", to ignore DNA absences: 
data_env_1 <- data_env %>% mutate(Adj_conc = ifelse(Adj_conc == 0, NA, Adj_conc))

# 1(a) Check for missing data ####

# Screening for NAs in the dataset:
colSums(is.na(data_env_1))
# 5 pH measurements are missing (due to probes not available in September, and pH
# probe not working in January)
# 2 ORP and conductivity measurements are missing (due to probes not available in 
# September)
# 2 phosphates and nitrates measurements are missing (due to tests not being carried
# out in September)
# In addition this time, we have 11 missing values for DNA detections.


# 1(b) Check for outliers ####

# Use Cleveland plots to look for outliers in continuous variables:
Names <- c("Water_Temp_DegC", "pH", "ORP_mV", "Conductivity", "Phos_AV", "Nitr_AV",
           "DepthMean_cm", "Met_Temp_DegC", "Met_WindSpd_mph", "Met_Pres_mb", 
           "Met_Hum_.", "Adj_conc", "KS2_HOG_Count")
dotplot(as.matrix(as.matrix(data_env_1[,Names])),
        groups=FALSE,
        strip = strip.custom(bg = 'white',
                             par.strip.text = list(cex = 1)),
        scales = list(x = list(relation = "free", draw = TRUE),
                      y = list(relation = "free", draw = FALSE)),
        col = "#228886", cex  = 1, pch = 20,
        xlab = list(label = "Data range (value of the variable)", cex = 1.2),
        ylab = list(label = "Order of the data", 
                    cex = 1.2))

# Cleveland dotplots indicate one outlier in the DNA concentration at relative conc of ~1. 
# This is also a lower concentration than found in the negative extraction control 
# sample that amplified - so might be indicative of this level of concentration 
# being indistinguishable from contamination.  

# Remove the outlier in question: 
data_env_1_clean <- data_env_1 %>% filter(is.na(Adj_conc) | Adj_conc <= 1.0 | Adj_conc >= 1.2)


# 1(c) Check for data balance in categorical variables ####

table(data_env_1_clean$Sample_Type) # Sample Types are well balanced (6,7,6)
table(data_env_1_clean$Met_WindDir) # Wind Direction is fairly well balanced (6, 5, 3, 5)
table(data_env_1_clean$Met_UV) # UV levels are not well balanced (3 high, 10 medium, 6 low)
# but do reflect normal distribution and expected for an environmental variable 
table(data_env_1_clean$PCR_Plate) # Data is fairly balanced over the three PCR plates (8,6,5)
table(data_env_1_clean$Date_Collected) # Data is balanced over sample collection date
table(data_env_1_clean$Date_Extracted) # Data is fairly balanced over extraction date


# There are a lot of variables here, so care must be taken fitting the data to 
# such a complex model. 


# 1(d) Check for lots of zeros in the response variable ####

# Calculate % of zeros in Adj_conc: 
sum(data_env_1_clean$Adj_conc == 0)*100/nrow(data_env)

# Zeros in the response variable have been removed, so this is no longer relevant
# as we are now only looking at presence data and ignoring absence data. 

# There were already no zeros in the kick-sampling data. 


# 1(e) Check for multicollinearity among environmental covariates ####

# Visualise multicollinearity using a correlation matrix with corresponding pairplots: 
Coll <- c("Water_Temp_DegC", "pH", "ORP_mV", "Conductivity", "Phos_AV",
          "Nitr_AV", "DepthMean_cm", "Met_Temp_DegC", "Met_WindSpd_mph", "Met_WindDir", 
          "Met_UV", "Met_Pres_mb", "Met_Hum_.", "Sample_Type")

panel.cor <- function(x, y, digits=1, prefix="", cex.cor = 6)
{usr <- par("usr"); on.exit(par(usr))
par(usr = c(0, 1, 0, 1))
r1 = cor(x,y,use="pairwise.complete.obs")
r <- abs(cor(x, y,use="pairwise.complete.obs"))
txt <- format(c(r1, 0.1), digits=digits)[1]
txt <- paste(prefix, txt, sep="")
if(missing(cex.cor)) { cex <- 0.6/strwidth(txt) } else {
  cex = cex.cor}
text(0.5, 0.5, txt, cex = cex * r)}
pairs(data_env_1_clean[, Coll], lower.panel = panel.cor, cex.labels = 0.5)

# Pairplots indicate strong collinearity among variables again 


# 1(f) Check for relationships among dependent and independent variables ####

par(mfrow=c(2,2), mar=c(5,5,1,1), cex.lab = 1)

plot(KS2_HOG_Count ~ Adj_conc,  data = data_env_1_clean,
     xlab = "Hoglouse eDNA concentration (relative qty, log10)", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, col = "#228886", cex = 1.3)
abline(lm(KS2_HOG_Count ~ Adj_conc, data = data_env_1_clean), col = "#C45A12", lwd = 1)
cor(data_env_1_clean$KS2_HOG_Count, data_env_1_clean$Adj_conc, use = "complete.obs")

# Strong correlation, 0.877
cor.test(data_env_1_clean$KS2_HOG_Count, data_env_1_clean$Adj_conc, use = "complete.obs")

plot(KS2_HOG_Count ~ Water_Temp_DegC,  data = data_env_1_clean,
     xlab = "Water temp (deg C)", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, col = "#228886", cex = 1.3)
abline(lm(KS2_HOG_Count ~ Water_Temp_DegC, data = data_env_1_clean), col = "#C45A12", lwd = 1)
cor(data_env_1_clean$KS2_HOG_Count, data_env_1_clean$Water_Temp_DegC, use = "complete.obs")

# Negative correlation, r = -0.702
cor.test(data_env_1_clean$KS2_HOG_Count, data_env_1_clean$Water_Temp_DegC, use = "complete.obs")

plot(KS2_HOG_Count  ~ pH,  data = data_env_1_clean,
     xlab = "Water pH", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, col = "#228886", cex = 1.3)
abline(lm(KS2_HOG_Count ~ pH, data = data_env_1_clean), col = "#C45A12", lwd = 1)
cor(data_env_1_clean$KS2_HOG_Count, data_env_1_clean$pH, use = "complete.obs")

# Weak correlation

plot(KS2_HOG_Count ~ ORP_mV,  data = data_env_1_clean,
     xlab = "Water ORP", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, col = "#228886", cex = 1.3)
abline(lm(KS2_HOG_Count ~ ORP_mV, data = data_env_1_clean), col = "#C45A12", lwd = 1)
cor(data_env_1_clean$KS2_HOG_Count, data_env_1_clean$ORP_mV, use = "complete.obs")

# No correlation

plot(KS2_HOG_Count  ~ Conductivity,  data = data_env_1_clean,
     xlab = "Water conductivity", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, col = "#228886", cex = 1.3)
abline(lm(KS2_HOG_Count ~ Conductivity, data = data_env_1_clean), col = "#C45A12", lwd = 1)
cor(data_env_1_clean$KS2_HOG_Count, data_env_1_clean$Conductivity, use = "complete.obs")

# No/very weak correlation

plot(KS2_HOG_Count  ~ Phos_AV,  data = data_env_1_clean,
     xlab = "Phosphate level", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, col = "#228886", cex = 1.3)
abline(lm(KS2_HOG_Count ~ Phos_AV, data = data_env_1_clean), col = "#C45A12", lwd = 1)
cor(data_env_1_clean$KS2_HOG_Count, data_env_1_clean$Phos_AV, use = "complete.obs")

# Weak negative correlation

plot(KS2_HOG_Count ~ Nitr_AV,  data = data_env_1_clean,
     xlab = "Nitrate level", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, col = "#228886", cex = 1.3)
abline(lm(KS2_HOG_Count ~ Nitr_AV, data = data_env_1_clean), col = "#C45A12", lwd = 1)
cor(data_env_1_clean$KS2_HOG_Count, data_env_1_clean$Nitr_AV, use = "complete.obs")

# Weak negative correlation

plot(KS2_HOG_Count ~ DepthMean_cm,  data = data_env_1_clean,
     xlab = "Mean depth (cm)", 
     ylab = "Hoglouse counts from kick-sampling",
     pch = 16, col = "#228886", cex = 1.3)
abline(lm(KS2_HOG_Count ~ DepthMean_cm, data = data_env_1_clean), col = "#C45A12", lwd = 1)
cor(data_env_1_clean$KS2_HOG_Count, data_env_1_clean$DepthMean_cm, use = "complete.obs")

# Weak negative correlation


# The plots indicate potential relationship between count data from kick-sampling 
# and qPCR DNA detection levels (strong positive relationship) and a strong negative 
# relationship with water temperature. Other indications are weak relationships. 


# 2. MODEL FITTING #### 

# The data exploration showed: 
  # 1. No outliers in the response variable (=count data from kick-sampling). 
  # 2. A broadly normal response variable and broadly conforming to homogeneity.
  # 3. No zeros in the response variable. 
  # 4. Significant collinearity in independent variables. 
  # 5. No interactions were tested. 
  # 6. Independence of response variable was assumed but not tested.

# Fitting a Gaussian General Linear Model (GLM)

# First we need to ignore the missing values in the dataset. 

data_model <- data_env_1_clean %>% filter(!is.na(Adj_conc))

# Additive model with water temperature: 
glm1 <- lm(KS2_HOG_Count ~ Adj_conc + Water_Temp_DegC, data=data_model)
summary(glm1)

# p>0.05 for water temp; p<0.05 for Adj_conc. 

# Backward stepwise selection: 

# Assessing AICs to see if covariates should be dropped from the model: 
drop1(glm1)

# AIC is lower if water temperature is removed. 


# Re-fit model with only one predictor variable:  

glm2 <- lm(KS2_HOG_Count ~ Adj_conc, data=data_model)
summary(glm2)


# MODEL VALIDATION ####

# Checking homogeneity of residual variance

Fitted <- fitted(glm2)
Resid <- resid(glm2)
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = Fitted, y=Resid, xlab="Fitted values", ylab="Pearson residuals", 
     pch=16, cex=1.5)
abline(h=0, lty=2)

# Distribution of residuals around 0 is broadly consistent along the fitted values


# Checking normality of residuals

# Adding residuals to dataframe: 
p <- ggplot()
p <- p + ylab("Frequency")
p <- p + xlab("Pearson residuals")
p <- p + theme(panel.background = element_blank())
p <- p + theme(panel.border = element_rect(fill = NA, 
                                           colour = "black", size = 1))
p <- p + theme(strip.background = element_rect(fill = "white", 
                                               colour = "white", size = 1))
p <- p + theme(text = element_text(size=15))
p <- p + geom_histogram(colour = "black", fill = "white", 
                        data = data_model, aes(x=Resid), bins = 10)
p

# Doesn't look normal, but little data to go on. 

#Normality of residuals
shapiro.test(Resid)

# p>0.05 therefore assumption of normality of the model residuals in satisfied. 


# Checking for influential observations

par(mfrow = c(1,1))
plot(cooks.distance(glm2),
     xlab="Observation", ylab="Cook's distance", type="h", 
     ylim=c(0,1.1), cex.lab=1.5)
abline(h=1, lty=2)

# There are no clear influential observations where Cook's distance exceeds 1. 


# MODEL PRESENTATION ####

# Summarise and visualise the model
options(digits=3)
summary(glm2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)   
#(Intercept)   -117.5       34.6   -3.40   0.0145 * 
#  Adj_conc        53.2       11.9    4.48   0.0042 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 8.56 on 6 degrees of freedom
#Multiple R-squared:  0.77,	Adjusted R-squared:  0.731 
#F-statistic:   20 on 1 and 6 DF,  p-value: 0.00421


# Plotting the model: 

#Fit linear regression to data
ggplot(data_model, aes(x = Adj_conc, y = KS2_HOG_Count)) +
  geom_point(shape = 16, size = 3, alpha = 0.7) +
  geom_smooth(method = 'lm', colour = 'red', se = FALSE, size = 1) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, size = 1)) +
  theme(strip.background = element_rect(fill = "white", 
                                        colour = "white", size = 1)) +
  theme(text = element_text(size=13)) +
  xlab("Log-10 relative DNA concentration from qPCR") + 
  ylab("Hoglouse counts from kick-sampling")

#Obtain model summary
summary(glm2)

# Visualise the model as a figure

range(data_model$Adj_conc)
# 2.46 to 3.26

MyData <- expand.grid(
  Adj_conc = seq(2.46, 3.26, 
               length = 25))

X <- model.matrix(~ Adj_conc, data = MyData)

#Calculate predicted kick-sampling count values
MyData$Pred <- X %*% coef(glm2)

#Calculate standard errors
MyData$SE <- sqrt(  diag(X %*% vcov(glm2) %*% t(X))  )

#And using the predicted values and standard errors, 
#calculate 95% confidence intervals for predicted fecundities
MyData$SeUp <- MyData$Pred + 1.96 * MyData$SE
MyData$SeLo <- MyData$Pred - 1.96 * MyData$SE

# PLOT: 
p <- ggplot()
p <- p + geom_point(data = data_model, aes(y = KS2_HOG_Count, x = Adj_conc), 
                    shape = 16, size = 2.5,  alpha = 0.7, colour="#228886")
p <- p + xlab("Log-10 relative DNA concentration from qPCR")
p <- p + ylab("Hoglouse counts from kick-sampling")
p <- p + theme(text = element_text(size=13)) 
p <- p + theme(panel.background = element_blank())
p <- p + theme(panel.border = element_rect(fill = NA, 
                                           colour = "black", size = 1))
p <- p + theme(strip.background = element_rect(fill = "white", 
                                               color = "white", size = 1))
p <- p + geom_line(data = MyData,aes(x = Adj_conc, y = Pred), 
                   colour = "#C45A12", size = 1.2)
p <- p + geom_ribbon(data = MyData, aes(x = Adj_conc, ymax = SeUp, 
                                        ymin = SeLo ), fill = "#F4B183", alpha=0.2)
p

# CONCLUSIONS #### 

# An additive Gaussian GLM was fitted to data to predict hoglouse counts from 
# hoglouse DNA detections through qPCR and water temperature. Backwards stepwise model
# selection was carried out and water temperature dropped from the model. This 
# resulted in a linear regression model being selected as the best-fitting model
# which took the form: Hoglouse counts = -117.5 + 53.2 x log10 DNA concentration
# (r^2 = 0.77, df = 6, p < 0.001)

