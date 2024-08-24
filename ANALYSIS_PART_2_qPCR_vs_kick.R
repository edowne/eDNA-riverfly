# BIOS0034 Research Project
# Student number: HVVY6
# 23/08/2024

# Data Analysis Part 2: qPCR/eDNA data vs kick-sampling data, and methods comparison ####


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

# Reading in dataset:
data_0 <- read.csv("COMP_Proj_Data.csv")

# Checking dataset: 
str(data_0)

# Reducing dataset to columns of interest: 
data_1 <- data_0 %>% dplyr::select(Sample_Ref, Date_Collected, Sample_Type, PCR_Code, 
                            Adj_conc, KS2_HOG_Count, KS2_HOG_Score, KS2_TotalScore, 
                            KS2_SpRich)

# Filtering out samples that weren't analysed: 
data_1 <- data_1 %>% dplyr::filter(PCR_Code!="")


# Make sure "Sample type" is a grouping factor: 
data_1$Sample_Type <- as.factor(data_1$Sample_Type)
# Check: 
str(data_1)

# Convert dates to date format: 
data_1 <- data_1 %>% mutate(Date_Collected = as.Date(Date_Collected, format = "%d/%m/%Y"))
# Check: 
str(data_1)

# Convert "NA" concentrations to zero, to indicate no target DNA detected in sample: 
data_1 <- data_1 %>% mutate(Adj_conc = ifelse(is.na(Adj_conc), 0, Adj_conc))

# Filter out controls - I am only looking at the water samples in this analysis: 
data_1 <- data_1 %>% filter(Sample_Type %in% c("Edge_Depth", "Edge_Surface", 
                                               "Transect_Depth"))


# Comparison of qPCR data with kick-sampling data ####

# Plotting the qPCR data against the kick-sampling data: 

plot_A <- ggplot(data_1, aes(x = Date_Collected)) +
  geom_line(data=data_1, aes(y = Adj_conc, colour = Sample_Type)) + # line for sample type
  geom_line(data=data_1, aes(y = KS2_HOG_Count * (max(Adj_conc) / max(KS2_HOG_Count)), 
                colour = "Hoglouse counts from kick-sampling"), linetype = "dashed") +# line for kick-sampling counts
  # Scale the secondary y-axis to match the original Concentration scale
  scale_y_continuous(sec.axis = sec_axis(~ . * (max(data_1$KS2_HOG_Count) / max(data_1$Adj_conc)), 
                                         name = "Hoglouse counts from kick-sampling")) +
    # Format x-axis to display months and set limits
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y", 
               limits = as.Date(c("2023-09-01", "2024-08-01"))) +
  labs(y = "Relative concentration of hoglouse DNA", x = "Date", 
       colour = "Legend") +
  theme_bw() +
  theme(axis.title.y.right = element_text(colour = "black"),
        axis.title.y.left = element_text(colour = "black"),
        axis.text.x = element_text(angle=45, hjust = 1),
        legend.position = "bottom")

#theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

plot_A

# Create scatterplots for each comparison: 

# Overall data: 
scatter_1 <- ggplot(data=data_1, aes(x=Adj_conc, y=KS2_HOG_Count))+
  geom_point()+
  theme_bw()
scatter_1

cor.test(~Adj_conc + KS2_HOG_Count, data=data_1, method="spearman")
allpots.lm <- lm(KS2_HOG_Count ~ Adj_conc, data=data_1)
coef(allpots.lm)
plot_allpots_lm <- ggplot(data=data_1, aes(x=Adj_conc, y=KS2_HOG_Count))+
  geom_point()+
  theme_bw()+
  geom_smooth(method="lm")+
  ylab("Kick-sampling counts")+
  xlab("Relative concentration of target DNA detected by qPCR")
#annotate("text",label="y = ")
plot_allpots_lm



# Filter for Pot 1(Edge surface samples): 
data_Pot1 <- data_1 %>% filter(Sample_Type=="Edge_Surface")

# Scatterplot for Pot 1 samples: 
scatter_pot1 <- ggplot(data=data_Pot1, aes(x=Adj_conc, y=KS2_HOG_Count)) +
                         geom_point() +
  geom_smooth(method=lm)+
                         theme_bw()
scatter_pot1                       

cor.test(~Adj_conc + KS2_HOG_Count, data=data_Pot1, method="spearman")
Pot1.lm <- lm(KS2_HOG_Count ~ Adj_conc, data=data_Pot1)
coef(Pot1.lm)
plot_Pot1_lm <- ggplot(data=data_Pot1, aes(x=Adj_conc, y=KS2_HOG_Count))+
  geom_point()+
  theme_bw()+
  geom_smooth(method="lm")+
  ylab("Kick-sampling counts")+
  xlab("Relative concentration of target DNA detected by qPCR")
  #annotate("text",label="y = ")
plot_Pot1_lm


# Filter for Pot 2 (Edge Depth samples):
data_Pot2 <- data_1 %>% filter(Sample_Type=="Edge_Depth")

# Scatterplot for Pot 2 samples: 
scatter_pot2 <- ggplot(data=data_Pot2, aes(x=Adj_conc, y=KS2_HOG_Count)) +
  geom_point()+
  theme_bw()
scatter_pot2                 

cor.test(~Adj_conc + KS2_HOG_Count, data=data_Pot2, method="spearman")
Pot2.lm <- lm(KS2_HOG_Count ~ Adj_conc, data=data_Pot2)
coef(Pot2.lm)
plot_Pot2_lm <- ggplot(data=data_Pot2, aes(x=Adj_conc, y=KS2_HOG_Count))+
  geom_point()+
  theme_bw()+
  geom_smooth(method="lm")+
  ylab("Kick-sampling counts")+
  xlab("Relative concentration of target DNA detected by qPCR")
#annotate("text",label="y = ")
plot_Pot2_lm


# Filter for Pot 3 (Transect Depth samples):
data_Pot3 <- data_1 %>% filter(Sample_Type=="Transect_Depth")

# Scatterplot for Pot 3 samples: 
scatter_pot3 <- ggplot(data=data_Pot3, aes(x=Adj_conc, y=KS2_HOG_Count)) +
  geom_point()+
  theme_bw()
scatter_pot3    

cor.test(~Adj_conc + KS2_HOG_Count, data=data_Pot3, method="spearman")
Pot3.lm <- lm(KS2_HOG_Count ~ Adj_conc, data=data_Pot3)
coef(Pot3.lm)
plot_Pot3_lm <- ggplot(data=data_Pot3, aes(x=Adj_conc, y=KS2_HOG_Count))+
  geom_point()+
  theme_bw()+
  geom_smooth(method="lm")+
  ylab("Kick-sampling counts")+
  xlab("Relative concentration of target DNA detected by qPCR")
#annotate("text",label="y = ")
plot_Pot3_lm


# Correlational analysis = not significant for any of the three sample types, 
# nor across all samples.  



# Data distribution checks ####

# Checking distribution of the data: 
hist(data_1$Adj_conc) # Not normally distributed due to high number of zeros 
hist(data_1$KS2_HOG_Count) # Hard to tell if normally distributed with data
  # available, but doesn't seem to follow Poisson distr as is common with other
  # count data. Appears skewed. Try some transformations: 

# Trying a log (base 10) transformation: 
hist(log(data_1$Adj_conc,10))
hist(log(data_1$KS2_HOG_Count, 10))

# Trying a log (base 10) + 1 transformation: 
hist(log(data_1$Adj_conc+1, 10))
hist(log(data_1$KS2_HOG_Count+1, 10))

# Trying a square root transformation: 
hist(sqrt(data_1$Adj_conc+0.5))
hist(sqrt(data_1$KS2_HOG_Count+0.5))

# Trying de'logged data for the conc: 
hist(10^data_1$Adj_conc)

# Logging this again with +1: 
hist(log((10^data_1$Adj_conc)+1, 10))


# Checking for zero inflation: 
# Zeros in Adj_conc: 
sum(data_1$Adj_conc == 0) # 11 zeros, out of 20 observations, more than 50% zeros

# Zeros in kick-sampling data: 
sum(data_1$KS2_HOG_Count == 0) # No zeros in the KS count data


# Sample type comparisons #### 

# Boxplots for each sample type: 


boxplot_pots <- ggplot(data_1, aes(x=Sample_Type, y=Adj_conc))+
  geom_boxplot(fill="#33CCCC")+
  geom_jitter(shape = 4, size = 2, stroke=1.1, width = 0.2, colour="#C45A12")+
  xlab("Sampling method")+
  ylab("Relative concentration of hoglouse DNA found in sample")+
  scale_x_discrete(limits = c("Edge_Surface", "Edge_Depth", "Transect_Depth"),
                   labels = c("Edge_Surface" = "Edge-Surface (Pot 1)", 
                              "Edge_Depth" = "Edge-Depth (Pot 2)", 
                              "Transect_Depth" = "Transect-Depth (Pot 3)"))+
  theme_bw()

boxplot_pots

# Potential outlier in Edge-surface (Pot 1) which, if removed, would reduce the
# detections in this group to zero. 

# ANOVA to test for differences between the sample types, with outlier kept in: 

# Ignoring lack of normality for the moment: 
anova_sample_types <- aov(Adj_conc ~ Sample_Type, data_1)
summary(anova_sample_types)

# Post-hoc Tukey HSD test: 
TukeyHSD(anova_sample_types)


# Trying non-parametric tests of difference: 
kruskal.test(Adj_conc ~ Sample_Type, data_1)
# No significant different (p>0.05)


# Now removing outlier and re-running analysis ####

# Removing the outlier: 
data_1_clean <- data_1 %>% filter(is.na(Adj_conc) | Adj_conc <= 3.2)

# Boxplots for each sample type: 

boxplot(Adj_conc ~ Sample_Type, data = data_1_clean)

boxplot_pots_2 <- ggplot(data_1_clean, aes(x=Sample_Type, y=Adj_conc))+
  geom_boxplot(fill="#33CCCC")+
  geom_jitter(shape = 4, size = 2, stroke=1.1, width = 0.2, colour="#C45A12")+
  xlab("Sample type")+
  ylab("Relative concentration of hoglouse DNA found in sample")+
  scale_x_discrete(limits = c("Edge_Surface", "Edge_Depth", "Transect_Depth"),
                   labels = c("Edge_Surface" = "Edge-Surface (Pot 1)", 
                              "Edge_Depth" = "Edge-Depth (Pot 2)", 
                              "Transect_Depth" = "Transect-Depth (Pot 3)"))+
  theme_bw()

boxplot_pots_2

# ANOVA to test for differences between the sample types, now without outlier: 

# Ignoring lack of normality for the moment: 
anova_sample_types_2 <- aov(Adj_conc ~ Sample_Type, data_1_clean)
summary(anova_sample_types_2)

# Difference is now significant, p<0.05

# Post-hoc Tukey HSD test: 
TukeyHSD(anova_sample_types_2)

# sig (p<0.05) difference between edge-surface and edge-depth, but not sig
# between other comparisons 


# Trying non-parametric tests of difference: 
kruskal.test(Adj_conc ~ Sample_Type, data_1_clean)
# No significant different (p>0.05)

# CONCLUSIONS: Some indication that edge-surface sampling is less likely to detect 
# target DNA, but not between other pots. Since data did not conform to parametric assumption, 
# non-parametric test ran and results not significant. 
# Insufficient datapoints here to demonstrate conclusive differences, although some 
# indication taht there is something here to explore futher and a depth sampling device might 
# be warranted when sampling benthic macroinvertebrates. 

