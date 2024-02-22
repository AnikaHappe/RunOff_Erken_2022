#_______________________________________________________________________________
#
# R SCRIPT FOR ANALYSIS OF THE RUN-OFF MESOCOSM EXPERIMENT AT LAKE ERKEN 2022
# Author: Anika Happe
# Date: 22nd February 2024
#
#_______________________________________________________________________________

rm(list=ls()) #Empty the environment
setwd("/Users/Anika/Library/Mobile Documents/com~apple~CloudDocs/PhD/Experiment in Schweden")

#Load packages
library(readxl)
library(tidyverse)
library(PKNCA) #For AUC function
library(readr)
library(tidyr)
library(dplyr)
library(drc)
library(DescTools)
library(bayestestR)
library(writexl)
library(lubridate)

#Important Papers for Methodology: 
#Hillebrand et al. 2018 (DOI: 10.1111/ele.12867), 
#White et al. 2020 (DOI: 10.1038/s41559-020-01315-w),
#Urrutia-Cordero et al. 2021 (DOI: 10.1111/1365-2745.13804)

#_______________________________________________________________________________

#### FIG1: Nutrient Additions ####

Table_Additions <- read_excel("Data/Table_Additions.xlsx")
Table_Additions <- pivot_longer(Table_Additions, cols = starts_with("Treatment_"), names_to = "treatment", values_to = "value")

plot1 <- ggplot(Table_Additions, aes(x = Day, y = value, fill = treatment)) +
  annotate("rect", xmin = 21, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey87")+
  geom_bar(stat = "identity", position = "dodge", width = 1.0) +
  labs(x = "Day", y = "Value", fill = "Treatment") +
  scale_fill_manual(values = c("Treatment_E" = "indianred", "Treatment_I" = "skyblue3", "Treatment_D" = "palegreen4"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))
plot1

ggsave("Fig1_Additions.eps", plot1 , unit = "cm", device = "eps", width = 35, height = 13, dpi = 300)

#_______________________________________________________________________________

#### FIG2: Time Series of DO and C ####

#The data table for the time series (except DO) can be found on PANGAEA!

#Read data and preparation of main data set
data <- read_excel("R Scripts/Tidy Graphs/Data_PANGAEA_SITES2022.xlsx")
data <- dplyr::select(data, Level, Meso = Mesocosm, Rep = Replicate, Treat = `Treatment: run-off`, S_No = `Sampling Number`, Day, 
                      C_µmol_L = `Carbon, intracellular [µmol/L]`, N_µmol_L = `Nitrogen, intracellular [µmol/L]`, P_µmol_L = `Phosphorus, intracellular [µmol/L]`, Si_µmol_L = `Silicate, intracellular [µmol/L]`,
                      CN = `C:N ratio, intracellular`, CP = `C:P ratio, intracellular`, CSi = `C:Si ratio, intracellular`,
                      NP = `N:P ratio, intracellular`, SiN = `Si:N ratio, intracellular`, SiP = `Si:P ratio, intracellular`)

str(data) #Check the data structure and transform is necessary
data$Meso <- as.numeric(data$Meso) 
data$Level <- as.factor(data$Level)
data$Treat <- as.factor(data$Treat)
data$CN <- as.numeric(data$CN)
data$CP <- as.numeric(data$CP)
data$CSi <- as.numeric(data$CSi)
data$NP <- as.numeric(data$NP)
data$SiN <- as.numeric(data$SiN)
data$SiP <- as.numeric(data$SiP)
data$P_µmol_L <- as.numeric(data$P_µmol_L)
data$Si_µmol_L <- as.numeric(data$Si_µmol_L)

data_phy1 <- filter(data, Level == "Smaller105µm") #We only need this level in this step
data_phy1 <- filter(data_phy1, Treat != "L") #We exclude the lake samples
data_phy1 <- filter(data_phy1, Meso != 2) #Mesocosm 2 is excluded

#Data from mesocosm #2 is excluded since it received one wrong nutrient addition
#early on and behaved quite differently leading to large standard deviations.

data_phy1_P <- filter(data_phy1, P_µmol_L > 0) #Excludes the NA value only
data_phy1_P <- filter(data_phy1_P, CP < 1000) #ID 295 and 296 (see below)

#Mesocosms/Timepoints are excluded that have NAs for P because it makes trouble
#when calculating means, standard deviations and ratios. And an extreme outlier
#with the ID 295 and 296 where the other three values are very close together and this
#one is extremly low, it not seemed to have worked in the filter processing.

data_phy_CN <- filter(data_phy1, C_µmol_L < 400) #Excludes one extreme data point
data_phy_CN <- mutate(data_phy_CN, SiC = Si_µmol_L/C_µmol_L)
data_phy_CN <- group_by(data_phy_CN, Day, Treat)
data_phy_CN <- mutate(data_phy_CN, Mean_N_µmol_L = mean(N_µmol_L), Sd_N_µmol_L = sd(N_µmol_L), Mean_C_µmol_L = mean(C_µmol_L), Sd_C_µmol_L = sd(C_µmol_L), Mean_CN = mean(CN), Sd_CN = sd(CN),
                      Mean_CSi = mean(CSi), Sd_CSi = sd(CSi), Mean_SiN = mean(SiN), sd_SiN = sd(SiN),
                      Mean_SiC = mean(SiC), Sd_SiC = sd(SiC))

data_phy1_P <- group_by(data_phy1_P, Day, Treat)
data_phy1_P <- mutate(data_phy1_P, Mean_P_µmol_L = mean(P_µmol_L), Sd_P_µmol_L = sd(P_µmol_L), Mean_CP = mean(CP), Sd_CP = sd(CP), Mean_NP = mean(NP), Sd_NP = sd(NP), Mean_SiP = mean(SiP), Sd_SiP = sd(SiP))

#CREATION OF PLOTS FOR THE RATIOS

data_phy_CN$Treat <- factor(data_phy_CN$Treat,levels = c("C","D","I","E"))
plot_C <- ggplot(data_phy_CN, aes(Day-1, Mean_C_µmol_L, color = Treat))+
  annotate("rect", xmin = 21, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey87")+
  geom_line(size=1.2)+
  geom_point(size=1.7)+
  geom_errorbar(aes(ymin = Mean_C_µmol_L - Sd_C_µmol_L, ymax = Mean_C_µmol_L + Sd_C_µmol_L), width = 0.4, size=0.9)+
  scale_colour_manual(values = c("grey60", "palegreen4", "skyblue3", "indianred"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
plot_C

ggsave("Fig2_C.eps", plot_C , unit = "cm", device = "eps", width = 15, height = 10, dpi = 300)

plot_CN <- ggplot(data_phy_CN, aes(Day-1, Mean_CN, color = Treat))+
  annotate("rect", xmin = 21, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey87")+
  geom_line(size=1.2)+
  geom_point(size=1.7)+
  geom_errorbar(aes(ymin = Mean_CN - Sd_CN, ymax = Mean_CN + Sd_CN), width = 0.4, size=0.9)+
  scale_colour_manual(values = c("grey60", "palegreen4", "skyblue3", "indianred"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
plot_CN

ggsave("Fig2_CN.eps", plot_CN , unit = "cm", device = "eps", width = 15, height = 10, dpi = 300)

plot_SiC <- ggplot(data_phy_CN, aes(Day-1, Mean_SiC, color = Treat))+
  annotate("rect", xmin = 21, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey87")+
  geom_line(size=1.2)+
  geom_point(size=1.7)+
  geom_errorbar(aes(ymin = Mean_SiC - Sd_SiC, ymax = Mean_SiC + Sd_SiC), width = 0.4, size=0.9)+
  scale_colour_manual(values = c("grey60", "palegreen4", "skyblue3", "indianred"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
plot_SiC

ggsave("Fig2_SiC.eps", plot_SiC , unit = "cm", device = "eps", width = 15, height = 10, dpi = 300)

plot_CP <- ggplot(data_phy1_P, aes(Day-1, Mean_CP, color = Treat))+
  annotate("rect", xmin = 21, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey87")+
  geom_line(size=1.2)+
  geom_point(size=1.7)+
  geom_errorbar(aes(ymin = Mean_CP - Sd_CP, ymax = Mean_CP + Sd_CP), width = 0.4, size=0.9)+
  scale_colour_manual(values = c("grey60", "palegreen4", "skyblue3", "indianred"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
plot_CP

ggsave("Fig2_CP.eps", plot_CP , unit = "cm", device = "eps", width = 15, height = 10, dpi = 300)

plot_NP <- ggplot(data_phy1_P, aes(Day-1, Mean_NP, color = Treat))+
  annotate("rect", xmin = 21, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey87")+
  geom_line(size=1.2)+
  geom_point(size=1.7)+
  geom_errorbar(aes(ymin = Mean_NP - Sd_NP, ymax = Mean_NP + Sd_NP), width = 0.4, size=0.9)+
  scale_colour_manual(values = c("grey60", "palegreen4", "skyblue3", "indianred"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
plot_NP

ggsave("Fig2_NP.eps", plot_NP , unit = "cm", device = "eps", width = 15, height = 10, dpi = 300)

plot_SiP <- ggplot(data_phy1_P, aes(Day-1, Mean_SiP, color = Treat))+
  annotate("rect", xmin = 21, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey87")+
  geom_line(size=1.2)+
  geom_point(size=1.7)+
  geom_errorbar(aes(ymin = Mean_SiP - Sd_SiP, ymax = Mean_SiP + Sd_SiP), width = 0.4, size=0.9)+
  scale_colour_manual(values = c("grey60", "palegreen4", "skyblue3", "indianred"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
plot_SiP

ggsave("Fig2_SiP.eps", plot_SiP , unit = "cm", device = "eps", width = 15, height = 10, dpi = 300)

plot_SiN <- ggplot(data_phy_CN, aes(Day-1, Mean_SiN, color = Treat))+
  annotate("rect", xmin = 21, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey87")+
  geom_line(size=1.2)+
  geom_point(size=1.7)+
  geom_errorbar(aes(ymin = Mean_SiN - sd_SiN, ymax = Mean_SiN + sd_SiN), width = 0.4, size=0.9)+
  scale_colour_manual(values = c("grey60", "palegreen4", "skyblue3", "indianred"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
plot_SiN

ggsave("Fig2_SiN.eps", plot_SiN , unit = "cm", device = "eps", width = 15, height = 10, dpi = 300)

#DISSOLVED OXYGEN DATA

#This data is available from the SITES AQUANET network
data_DO <- read_csv("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Experiment in Schweden/Data/6) Sensor Data SITES Erken/Cleaned data + R scripts/Erken_Daily_avg_final_new_clean.csv")

data.1 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(1)", "DOconc_optode(1)")
data.1 <- mutate(data.1, Meso = 1, Treat = "C", Rep = 1)
data.1 <- rename(data.1, "DOsat" = "DOsat_optode(1)", "DOconc" = "DOconc_optode(1)")

data.2 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(2)", "DOconc_optode(2)")
data.2 <- mutate(data.2, Meso = 2, Treat = "D", Rep = 1)
data.2 <- rename(data.2, "DOsat" = "DOsat_optode(2)", "DOconc" = "DOconc_optode(2)")

data.3 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(3)", "DOconc_optode(3)")
data.3 <- mutate(data.3, Meso = 3, Treat = "I", Rep = 1)
data.3 <- rename(data.3, "DOsat" = "DOsat_optode(3)", "DOconc" = "DOconc_optode(3)")

data.4 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(4)", "DOconc_optode(4)")
data.4 <- mutate(data.4, Meso = 4, Treat = "E", Rep = 1)
data.4 <- rename(data.4, "DOsat" = "DOsat_optode(4)", "DOconc" = "DOconc_optode(4)")

data.5 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(5)", "DOconc_optode(5)")
data.5 <- mutate(data.5, Meso = 5, Treat = "I", Rep = 2)
data.5 <- rename(data.5, "DOsat" = "DOsat_optode(5)", "DOconc" = "DOconc_optode(5)")

data.6 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(6)", "DOconc_optode(6)")
data.6 <- mutate(data.6, Meso = 6, Treat = "E", Rep = 2)
data.6 <- rename(data.6, "DOsat" = "DOsat_optode(6)", "DOconc" = "DOconc_optode(6)")

data.7 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(7)", "DOconc_optode(7)")
data.7 <- mutate(data.7, Meso = 7, Treat = "C", Rep = 2)
data.7 <- rename(data.7, "DOsat" = "DOsat_optode(7)", "DOconc" = "DOconc_optode(7)")

data.8 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(8)", "DOconc_optode(8)")
data.8 <- mutate(data.8, Meso = 8, Treat = "D", Rep = 2)
data.8 <- rename(data.8, "DOsat" = "DOsat_optode(8)", "DOconc" = "DOconc_optode(8)")

data.9 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(9)", "DOconc_optode(9)")
data.9 <- mutate(data.9, Meso = 9, Treat = "E", Rep = 3)
data.9 <- rename(data.9, "DOsat" = "DOsat_optode(9)", "DOconc" = "DOconc_optode(9)")

data.10 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(10)", "DOconc_optode(10)")
data.10 <- mutate(data.10, Meso = 10, Treat = "C", Rep = 3)
data.10 <- rename(data.10, "DOsat" = "DOsat_optode(10)", "DOconc" = "DOconc_optode(10)")

data.11 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(11)", "DOconc_optode(11)")
data.11 <- mutate(data.11, Meso = 11, Treat = "D", Rep = 3)
data.11 <- rename(data.11, "DOsat" = "DOsat_optode(11)", "DOconc" = "DOconc_optode(11)")

data.12 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(12)", "DOconc_optode(12)")
data.12 <- mutate(data.12, Meso = 12, Treat = "I", Rep = 3)
data.12 <- rename(data.12, "DOsat" = "DOsat_optode(12)", "DOconc" = "DOconc_optode(12)")

data.13 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(13)", "DOconc_optode(13)")
data.13 <- mutate(data.13, Meso = 13, Treat = "D", Rep = 4)
data.13 <- rename(data.13, "DOsat" = "DOsat_optode(13)", "DOconc" = "DOconc_optode(13)")

data.14 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(14)", "DOconc_optode(14)")
data.14 <- mutate(data.14, Meso = 14, Treat = "I", Rep = 4)
data.14 <- rename(data.14, "DOsat" = "DOsat_optode(14)", "DOconc" = "DOconc_optode(14)")

data.15 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(15)", "DOconc_optode(15)")
data.15 <- mutate(data.15, Meso = 15, Treat = "E", Rep = 4)
data.15 <- rename(data.15, "DOsat" = "DOsat_optode(15)", "DOconc" = "DOconc_optode(15)")

data.16 <- dplyr::select(data_DO, TIMESTAMP, "DOsat_optode(16)", "DOconc_optode(16)")
data.16 <- mutate(data.16, Meso = 16, Treat = "C", Rep = 4)
data.16 <- rename(data.16, "DOsat" = "DOsat_optode(16)", "DOconc" = "DOconc_optode(16)")

data_full <- rbind(data.1, data.2, data.3, data.4, data.5, data.6, data.7, data.8, data.9, data.10, data.11, data.12, data.13, data.14, data.15, data.16)
#write_xlsx(data_full,"DO_Data_Tidy_Erken2022.xlsx")

#Calculate mean and variance

data_full$TIMESTAMP <- as.Date(data_full$TIMESTAMP)
data_full$Day <- as.integer(data_full$TIMESTAMP - min(data_full$TIMESTAMP))

data_full$Treat <- as.factor(data_full$Treat)
data_full <- filter(data_full, Meso != 2) #Exclude mesocosm 2 again
data_full <- group_by(data_full, TIMESTAMP, Treat)
data_full <- mutate(data_full, Mean_DOsat = mean(DOsat), SD_DOsat = sd(DOsat), Mean_DOconc = mean(DOconc), SD_DOconc = sd(DOconc))
data_full <- filter(data_full, Rep == 4) #All replicates contain the same data now

#Plotting of DO data

data_full$Treat <- factor(data_full$Treat,levels = c("C","D","I","E"))
Fig2_DO <- ggplot(data_full, aes(Day, Mean_DOsat, color = Treat))+
  annotate("rect", xmin = 21, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey87")+
  geom_line(size=1.2)+
  geom_point(size=1.7)+
  geom_errorbar(aes(ymin = Mean_DOsat - SD_DOsat, ymax = Mean_DOsat + SD_DOsat), width = 0.4, size=0.9)+
  scale_colour_manual(values = c("grey60", "palegreen4", "skyblue3", "indianred"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
Fig2_DO

ggsave("Fig2_DO.eps", Fig2_DO , unit = "cm", device = "eps", width = 15, height = 10, dpi = 300)

#PHOTOSYNTHETIC ACTIVE RADIATION

Data_C <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Experiment in Schweden/Data/6) Sensor Data SITES Erken/Light_Sensor/Erken_hourly_avg_control.csv")
Data_D <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Experiment in Schweden/Data/6) Sensor Data SITES Erken/Light_Sensor/Erken_hourly_avg_daily.csv")
Data_I <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Experiment in Schweden/Data/6) Sensor Data SITES Erken/Light_Sensor/Erken_hourly_avg_intermediate.csv")
Data_E <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Experiment in Schweden/Data/6) Sensor Data SITES Erken/Light_Sensor/Erken_hourly_avg_extreme.csv")

Data_C <- dplyr::select(Data_C, TIMESTAMP, PAR_C = PAR_Apogee_avg)
Data_D <- dplyr::select(Data_D, TIMESTAMP, PAR_D = PAR_Apogee_avg)
Data_I <- dplyr::select(Data_I, TIMESTAMP, PAR_I = PAR_Apogee_avg)
Data_E <- dplyr::select(Data_E, TIMESTAMP, PAR_E = PAR_Apogee_avg)
Data_Full <- Data_C %>% full_join(Data_D, by = "TIMESTAMP") %>% full_join(Data_I, by = "TIMESTAMP") %>% full_join(Data_E, by = "TIMESTAMP")
Data_Tidy <- Data_Full %>% pivot_longer(cols = starts_with("PAR_"), names_to = "Treat", values_to = "Value")
Data_Tidy <- Data_Tidy %>% mutate(Julian_Day = yday(TIMESTAMP))

Data_Mean_Day <- group_by(Data_Tidy, Julian_Day, Treat)
Data_Mean_Day <- mutate(Data_Mean_Day, Mean_Day = mean(Value)) #No standard deviation available here
Data_Mean_Day <- filter(Data_Mean_Day, Julian_Day != 188) #Measures started at 19h
Data_Mean_Day <- dplyr::select(Data_Mean_Day, Treat, Julian_Day, Mean_Day)
Data_Mean_Day <- distinct(Data_Mean_Day)
Data_Mean_Day <- mutate(Data_Mean_Day, Incub_Day = Julian_Day - 189) #Translate Julian Day into Incubation Day (07.07.2022 = Day 0)

Data_Mean_Day$Treat <- factor(Data_Mean_Day$Treat,levels = c("PAR_C","PAR_D","PAR_I","PAR_E"))
Fig2_PAR <- ggplot(Data_Mean_Day, aes(Incub_Day, Mean_Day, color = Treat))+
  annotate("rect", xmin = 21, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey87")+
  geom_line(size=1.2)+
  geom_point(size=1.7)+
  #geom_errorbar(aes(ymin = Mean_Day - Sd_Day, ymax = Mean_Day + Sd_Day), width = 0.4, size=0.9)+
  scale_colour_manual(values = c("grey60", "palegreen4", "skyblue3", "indianred"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
Fig2_PAR

ggsave("Fig2_PAR.eps", Fig2_PAR , unit = "cm", device = "eps", width = 15, height = 10, dpi = 300)

#_______________________________________________________________________________

#### FIG 3: LRRs, OEV and recovery of phytoplankton ####

#### RECOVERY ####

data_phy_CN_recovery <- data_phy_CN
data_phy_CN_recovery <- filter(data_phy_CN_recovery, S_No == 10)

#To include the variance, we calculate the LRR for each replicate. So,
#ln(disturbed from replicate 1 / mean control) and so on. With this
#we can calculate a LRR in the end.

#Paste the value of the mean undisturbed (control) treatment in a new column
data_phy_CN_recovery <- mutate(data_phy_CN_recovery, Control_C = mean(data_phy_CN_recovery$C_µmol_L[data_phy_CN_recovery$Treat == "C"]))
data_phy_CN_recovery <- mutate(data_phy_CN_recovery, Control_CN = mean(data_phy_CN_recovery$CN[data_phy_CN_recovery$Treat == "C"]))
data_phy_CN_recovery <- mutate(data_phy_CN_recovery, Control_SiN = mean(data_phy_CN_recovery$SiN[data_phy_CN_recovery$Treat == "C"]))
data_phy_CN_recovery <- mutate(data_phy_CN_recovery, Control_SiC = mean(data_phy_CN_recovery$CSi[data_phy_CN_recovery$Treat == "C"]))
data_phy_CN_recovery <- filter(data_phy_CN_recovery, Treat != "C") #Not necessary anymore!

#Calculate the LRR for Recovery = ln(Disturbed/Control) for C
data_phy_CN_recovery_C <- mutate(data_phy_CN_recovery, Recovery = log(C_µmol_L/Control_C), Stab_Para = "C")
data_phy_CN_recovery_C <- group_by(data_phy_CN_recovery_C, Treat)
data_phy_CN_recovery_C_mean <- mutate(data_phy_CN_recovery_C, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
data_phy_CN_recovery_C_mean <- filter(data_phy_CN_recovery_C_mean, Rep == 4) #They should all be the same now!

#Calculate the LRR for Recovery = ln(Disturbed/Control) for C:N Ratio
data_phy_CN_recovery_CN <- mutate(data_phy_CN_recovery, Recovery = log(CN/Control_CN), Stab_Para = "CN")
data_phy_CN_recovery_CN <- group_by(data_phy_CN_recovery_CN, Treat)
data_phy_CN_recovery_CN_mean <- mutate(data_phy_CN_recovery_CN, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
data_phy_CN_recovery_CN_mean <- filter(data_phy_CN_recovery_CN_mean, Rep == 4) #They should all be the same now!

#Calculate the LRR for Recovery = ln(Disturbed/Control) for Si:N Ratio
data_phy_CN_recovery_SiN <- mutate(data_phy_CN_recovery, Recovery = log(SiN/Control_SiN), Stab_Para = "SiN")
data_phy_CN_recovery_SiN <- group_by(data_phy_CN_recovery_SiN, Treat)
data_phy_CN_recovery_SiN_mean <- mutate(data_phy_CN_recovery_SiN, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
data_phy_CN_recovery_SiN_mean <- filter(data_phy_CN_recovery_SiN_mean, Rep == 4) #They should all be the same now!

#Calculate the LRR for Recovery = ln(Disturbed/Control) for Si:N Ratio
data_phy_CN_recovery_SiC <- mutate(data_phy_CN_recovery, Recovery = log(SiC/Control_SiC), Stab_Para = "SiC")
data_phy_CN_recovery_SiC <- group_by(data_phy_CN_recovery_SiC, Treat)
data_phy_CN_recovery_SiC_mean <- mutate(data_phy_CN_recovery_SiC, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
data_phy_CN_recovery_SiC_mean <- filter(data_phy_CN_recovery_SiC_mean, Rep == 4) #They should all be the same now!

#Now: The same for P ratios with the P-data set

data_phy_P_recovery <- data_phy1_P
data_phy_P_recovery <- filter(data_phy_P_recovery, S_No == 10)

#Paste the value of the mean undisturbed (control) treatment in a new column
data_phy_P_recovery <- mutate(data_phy_P_recovery, Control_NP = mean(data_phy_P_recovery$NP[data_phy_P_recovery$Treat == "C"]))
data_phy_P_recovery <- mutate(data_phy_P_recovery, Control_CP = mean(data_phy_P_recovery$CP[data_phy_P_recovery$Treat == "C"]))
data_phy_P_recovery <- mutate(data_phy_P_recovery, Control_SiP = mean(data_phy_P_recovery$SiP[data_phy_P_recovery$Treat == "C"]))
data_phy_P_recovery <- filter(data_phy_P_recovery, Treat != "C") #Not necessary anymore!

#Calculate the LRR for Recovery = ln(Disturbed/Control) for CP
data_phy_P_recovery_CP <- mutate(data_phy_P_recovery, Recovery = log(CP/Control_CP), Stab_Para = "CP")
data_phy_P_recovery_CP <- group_by(data_phy_P_recovery_CP, Treat)
data_phy_P_recovery_CP_mean <- mutate(data_phy_P_recovery_CP, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
data_phy_P_recovery_CP_mean <- filter(data_phy_P_recovery_CP_mean, Rep == 4) #They should all be the same now!

#Calculate the LRR for Recovery = ln(Disturbed/Control) for NP
data_phy_P_recovery_NP <- mutate(data_phy_P_recovery, Recovery = log(NP/Control_NP), Stab_Para = "NP")
data_phy_P_recovery_NP <- group_by(data_phy_P_recovery_NP, Treat)
data_phy_P_recovery_NP_mean <- mutate(data_phy_P_recovery_NP, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
data_phy_P_recovery_NP_mean <- filter(data_phy_P_recovery_NP_mean, Rep == 4) #They should all be the same now!

#Calculate the LRR for Recovery = ln(Disturbed/Control) for SiP
data_phy_P_recovery_SiP <- mutate(data_phy_P_recovery, Recovery = log(SiP/Control_SiP), Stab_Para = "SiP")
data_phy_P_recovery_SiP <- group_by(data_phy_P_recovery_SiP, Treat)
data_phy_P_recovery_SiP_mean <- mutate(data_phy_P_recovery_SiP, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
data_phy_P_recovery_SiP_mean <- filter(data_phy_P_recovery_SiP_mean, Rep == 4) #They should all be the same now!

#Put all the data sets together
data_recovery_mean <- rbind(data_phy_CN_recovery_C_mean, data_phy_CN_recovery_CN_mean, data_phy_CN_recovery_SiC_mean, data_phy_CN_recovery_SiN_mean, data_phy_P_recovery_CP_mean, data_phy_P_recovery_NP_mean, data_phy_P_recovery_SiP_mean)
data_recovery_mean <- dplyr::select(data_recovery_mean, Mesocosm = Meso, Treat, Stab_Para, mean_Recovery, sd_Recovery)
data_recovery <- rbind(data_phy_CN_recovery_C, data_phy_CN_recovery_CN, data_phy_CN_recovery_SiC, data_phy_CN_recovery_SiN, data_phy_P_recovery_CP, data_phy_P_recovery_NP, data_phy_P_recovery_SiP)
data_recovery <- dplyr::select(data_recovery, Mesocosm = Meso, Treat, Rep, Stab_Para, Recovery)

#write_xlsx(data_recovery, "SITES_2022_SmallPhyto_Recovery.xlsx")

#### CALCULATE LOG-RESPONSE RATIOS (LRR) ####

data_phy_CN_resilience <- data_phy_CN
data_phy_CN_resilience$S_No <- as.numeric(data_phy_CN_resilience$S_No)
data_phy_P_resilience <- data_phy1_P
data_phy_P_resilience$S_No <- as.numeric(data_phy_P_resilience$S_No)

#Put the mean value of the control for each sampling day in a new column
data_phy_CN_resilience_DIE <- mutate(data_phy_CN_resilience, Mean_C_Control = NA, Mean_CN_Control = NA, Mean_SiC_Control = NA, Mean_SiN_Control = NA)
data_phy_P_resilience_DIE <- mutate(data_phy_P_resilience, Mean_CP_Control = NA, Mean_NP_Control = NA, Mean_SiP_Control = NA)

for (i in unique(data_phy_CN_resilience_DIE$S_No)) {
  data_phy_CN_resilience_DIE$Mean_C_Control[data_phy_CN_resilience_DIE$S_No == i] <-
    data_phy_CN_resilience_DIE$Mean_C_µmol_L[data_phy_CN_resilience_DIE$S_No == i & data_phy_CN_resilience_DIE$Treat == "C" & data_phy_CN_resilience_DIE$Rep == 4]
}
for (i in unique(data_phy_CN_resilience_DIE$S_No)) {
  data_phy_CN_resilience_DIE$Mean_CN_Control[data_phy_CN_resilience_DIE$S_No == i] <-
    data_phy_CN_resilience_DIE$Mean_CN[data_phy_CN_resilience_DIE$S_No == i & data_phy_CN_resilience_DIE$Treat == "C" & data_phy_CN_resilience_DIE$Rep == 4]
}
for (i in unique(data_phy_CN_resilience_DIE$S_No)) {
  data_phy_CN_resilience_DIE$Mean_SiC_Control[data_phy_CN_resilience_DIE$S_No == i] <-
    data_phy_CN_resilience_DIE$Mean_SiC[data_phy_CN_resilience_DIE$S_No == i & data_phy_CN_resilience_DIE$Treat == "C" & data_phy_CN_resilience_DIE$Rep == 4]
}
for (i in unique(data_phy_CN_resilience_DIE$S_No)) {
  data_phy_CN_resilience_DIE$Mean_SiN_Control[data_phy_CN_resilience_DIE$S_No == i] <-
    data_phy_CN_resilience_DIE$Mean_SiN[data_phy_CN_resilience_DIE$S_No == i & data_phy_CN_resilience_DIE$Treat == "C" & data_phy_CN_resilience_DIE$Rep == 4]
}
for (i in unique(data_phy_P_resilience_DIE$S_No)) {
  data_phy_P_resilience_DIE$Mean_CP_Control[data_phy_P_resilience_DIE$S_No == i] <-
    data_phy_P_resilience_DIE$Mean_CP[data_phy_P_resilience_DIE$S_No == i & data_phy_P_resilience_DIE$Treat == "C" & data_phy_P_resilience_DIE$Rep == 1]
}
for (i in unique(data_phy_P_resilience_DIE$S_No)) {
  data_phy_P_resilience_DIE$Mean_NP_Control[data_phy_P_resilience_DIE$S_No == i] <-
    data_phy_P_resilience_DIE$Mean_NP[data_phy_P_resilience_DIE$S_No == i & data_phy_P_resilience_DIE$Treat == "C" & data_phy_P_resilience_DIE$Rep == 1]
}
for (i in unique(data_phy_P_resilience_DIE$S_No)) {
  data_phy_P_resilience_DIE$Mean_SiP_Control[data_phy_P_resilience_DIE$S_No == i] <-
    data_phy_P_resilience_DIE$Mean_SiP[data_phy_P_resilience_DIE$S_No == i & data_phy_P_resilience_DIE$Treat == "C" & data_phy_P_resilience_DIE$Rep == 1]
}

data_phy_CN_resilience_DIE <- mutate(data_phy_CN_resilience_DIE, LRR_C = log(C_µmol_L/Mean_C_Control), LRR_CN = log(CN/Mean_CN_Control), LRR_SiC = log(SiC/Mean_SiC_Control), LRR_SiN = log(SiN/Mean_SiN_Control))
data_phy_CN_resilience_DIE_mean <- group_by(data_phy_CN_resilience_DIE, Treat, S_No)
data_phy_CN_resilience_DIE_mean <- mutate(data_phy_CN_resilience_DIE_mean, mean_LRR_C = mean(LRR_C), sd_LRR_C = sd(LRR_C), mean_LRR_CN = mean(LRR_CN), sd_LRR_CN = sd(LRR_CN), mean_LRR_SiC = mean(LRR_SiC), sd_LRR_SiC = sd(LRR_SiC), mean_LRR_SiN = mean(LRR_SiN), sd_LRR_SiN = sd(LRR_SiN))
data_phy_CN_resilience_DIE_mean <- filter(data_phy_CN_resilience_DIE_mean, Rep == 4)
data_phy_CN_resilience_DIE_mean <- filter(data_phy_CN_resilience_DIE_mean, Treat != "C")
data_phy_CN_resilience_DIE_mean <- na.omit(data_phy_CN_resilience_DIE_mean)

data_phy_P_resilience_DIE <- mutate(data_phy_P_resilience_DIE, LRR_CP = log(CP/Mean_CP_Control), LRR_NP = log(NP/Mean_NP_Control), LRR_SiP = log(SiP/Mean_SiP_Control))
data_phy_P_resilience_DIE_mean <- group_by(data_phy_P_resilience_DIE, Treat, S_No)
data_phy_P_resilience_DIE_mean <- mutate(data_phy_P_resilience_DIE_mean, mean_LRR_CP = mean(LRR_CP), sd_LRR_CP = sd(LRR_CP), mean_LRR_NP = mean(LRR_NP), sd_LRR_NP = sd(LRR_NP), mean_LRR_SiP = mean(LRR_SiP), sd_LRR_SiP = sd(LRR_SiP))
data_phy_P_resilience_DIE_mean <- filter(data_phy_P_resilience_DIE_mean, Rep == 2)
data_phy_P_resilience_DIE_mean <- filter(data_phy_P_resilience_DIE_mean, Treat != "C")
data_phy_P_resilience_DIE_mean <- na.omit(data_phy_P_resilience_DIE_mean)

#At this point it would be possible to create time series plots of the LRRs

#### OVERALL ECOLOGICAL VULNERABILITY (OEV) ####

#Take the LRR of the parameters that I want to look at (over entire time series)
AUC_CN <- data_phy_CN_resilience_DIE #From the resilience calculations
AUC_P <- data_phy_P_resilience_DIE #From the resilience calculations

#Fill the gaps of missing data, otherwise the function does not work!
#I manually took the mean between the point before and after.
M15_Rep4_TreatE_SNo3 <- data.frame(Meso = 15, Rep = 4, Treat = "E", S_No = 3, LRR_CP = -0.2790574, LRR_NP = -0.3579183, LRR_PSi = (((-0.614737865)+(0.200229520))/2))
M16_Rep4_TreatC_SNo3 <- data.frame(Meso = 16, Rep = 4, Treat = "C", S_No = 3, LRR_CP = -0.1189456, LRR_NP = -0.1337062, LRR_PSi = (((-0.111085221)+(-0.061122558))/2))
M12_Rep3_TreatI_SNo8 <- data.frame(Meso = 12, Rep = 3, Treat = "I", S_No = 8, LRR_CP = -0.4139685, LRR_NP = -0.2880625, LRR_PSi = (((-0.985547248)+(-1.311735031))/2))
AUC_P$Rep <- as.numeric(AUC_P$Rep)
AUC_P <- rbind(AUC_P, M15_Rep4_TreatE_SNo3, M16_Rep4_TreatC_SNo3, M12_Rep3_TreatI_SNo8)

auc_CN <- AUC_CN %>% mutate(mesocosmID = paste(Meso, Treat, Rep, sep = "_"))
auc_CN <- auc_CN[with(auc_CN, order(mesocosmID, S_No)),]
auc_P <- AUC_P %>% mutate(mesocosmID = paste(Meso, Treat, Rep, sep = "_"))
auc_P <- auc_P[with(auc_P, order(mesocosmID, S_No)),]

#The first sampling needs to be 0 in order to use the AUC function
auc_CN$S_No <- auc_CN$S_No-1
auc_P$S_No <- auc_P$S_No-1

#Create an USI (identifier for the experimental unit)
USI_CN <- unique(auc_CN$mesocosmID)
USI_P <- unique(auc_P$mesocosmID)

#Create empty data frame for AUC
area_C <- data.frame()
area_CN <- data.frame()
area_SiC <- data.frame()
area_SiN <- data.frame()
area_CP <- data.frame()
area_NP <- data.frame()
area_SiP <- data.frame()

#Conduct a loop for calculating AUC
for(i in 1:length(USI_CN)){
  temp <- auc_CN[auc_CN$mesocosmID==USI_CN[i], ] #temp is a temporary data frame
  area_value <- pk.calc.auc(abs(temp$LRR_C),temp$S_No)
  area_C <- rbind(area_C,data.frame(temp, area_value))
} #For C
for(i in 1:length(USI_CN)){
  temp <- auc_CN[auc_CN$mesocosmID==USI_CN[i], ] #temp is a temporary data frame
  area_value <- pk.calc.auc(abs(temp$LRR_CN),temp$S_No)
  area_CN <- rbind(area_CN,data.frame(temp, area_value))
} #For CN
for(i in 1:length(USI_CN)){
  temp <- auc_CN[auc_CN$mesocosmID==USI_CN[i], ] #temp is a temporary data frame
  area_value <- pk.calc.auc(abs(temp$LRR_SiC),temp$S_No)
  area_SiC <- rbind(area_SiC,data.frame(temp, area_value))
} #For SiC
for(i in 1:length(USI_CN)){
  temp <- auc_CN[auc_CN$mesocosmID==USI_CN[i], ] #temp is a temporary data frame
  area_value <- pk.calc.auc(abs(temp$LRR_SiN),temp$S_No)
  area_SiN <- rbind(area_SiN,data.frame(temp, area_value))
} #For SiN
for(i in 1:length(USI_P)){
  temp <- auc_P[auc_P$mesocosmID==USI_P[i], ] #temp is a temporary data frame
  area_value <- pk.calc.auc(abs(temp$LRR_CP),temp$S_No)
  area_CP <- rbind(area_CP,data.frame(temp, area_value))
} #For CP
for(i in 1:length(USI_P)){
  temp <- auc_P[auc_P$mesocosmID==USI_P[i], ] #temp is a temporary data frame
  area_value <- pk.calc.auc(abs(temp$LRR_NP),temp$S_No)
  area_NP <- rbind(area_NP,data.frame(temp, area_value))
} #For NP
for(i in 1:length(USI_P)){
  temp <- auc_P[auc_P$mesocosmID==USI_P[i], ] #temp is a temporary data frame
  area_value <- pk.calc.auc(abs(temp$LRR_SiP),temp$S_No)
  area_SiP <- rbind(area_SiP,data.frame(temp, area_value))
} #For SiP

#Remove duplicated rows
area_C <- area_C %>% distinct(mesocosmID, .keep_all = TRUE) %>% mutate(Stab_Para = "C") %>% dplyr::select(Mesocosm = Meso, Treat, Rep, Stab_Para, AUC = area_value)
area_CN <- area_CN %>% distinct(mesocosmID, .keep_all = TRUE) %>% mutate(Stab_Para = "CN") %>% dplyr::select(Mesocosm = Meso, Treat, Rep, Stab_Para, AUC = area_value)
area_SiC <- area_SiC %>% distinct(mesocosmID, .keep_all = TRUE) %>% mutate(Stab_Para = "SiC") %>% dplyr::select(Mesocosm = Meso, Treat, Rep, Stab_Para, AUC = area_value)
area_SiN <- area_SiN %>% distinct(mesocosmID, .keep_all = TRUE) %>% mutate(Stab_Para = "SiN") %>% dplyr::select(Mesocosm = Meso, Treat, Rep, Stab_Para, AUC = area_value)
area_CP <- area_CP %>% distinct(mesocosmID, .keep_all = TRUE) %>% mutate(Stab_Para = "CP") %>% dplyr::select(Mesocosm = Meso, Treat, Rep, Stab_Para, AUC = area_value)
area_NP <- area_NP %>% distinct(mesocosmID, .keep_all = TRUE) %>% mutate(Stab_Para = "NP") %>% dplyr::select(Mesocosm = Meso, Treat, Rep, Stab_Para, AUC = area_value)
area_SiP <- area_SiP %>% distinct(mesocosmID, .keep_all = TRUE) %>% mutate(Stab_Para = "SiP") %>% dplyr::select(Mesocosm = Meso, Treat, Rep, Stab_Para, AUC = area_value)

#Check if it worked with one example

area_CP$Treat <- factor(area_CN$Treat,levels = c("C","D","I","E"))
ggplot(area_CP, aes(x = Treat, y = AUC, fill = Treat)) +
theme_bw()+ 
  geom_hline(yintercept=0, linetype = "dashed")+
  geom_boxplot()

#Put it into form to merge with final stability table
AUC_full <- rbind(area_C, area_CN, area_CP, area_NP, area_SiN, area_SiP, area_SiC)
AUC_full <- filter(AUC_full, Treat != "C")
final_stability <- left_join(AUC_full, data_recovery, by = c("Mesocosm", "Treat", "Rep", "Stab_Para"))

#write_xlsx(final_stability, "SITES_2022_SmallPhyto_Stability_Feb2024.xlsx")

#PREPARE FOR

final_stability$Stab_Para <- as.factor(final_stability$Stab_Para)
final_stability$Rep <- as.numeric(final_stability$Rep)
final_stability_mean1 <- group_by(final_stability, Treat, Stab_Para)
final_stability_mean1 <- mutate(final_stability_mean1, Mean_AUC = mean(AUC), Sd_AUC = sd(AUC), Mean_Recovery = mean(Recovery), Sd_Recovery = sd(Recovery))
final_stability_mean1 <- filter(final_stability_mean1, Rep == 4)
final_stability_mean1 <- dplyr::select(final_stability_mean1, -AUC, -Recovery)
final_stability_mean1 <- pivot_longer(final_stability_mean1, cols = starts_with(c("Mean_", "Sd_")), names_to = c(".value", "variable"), names_pattern = "^(Mean|Sd)_([A-Za-z_]+)")

final_stability_mean1$Treat <- factor(final_stability_mean1$Treat, levels = c("D", "I", "E"))
final_stability_mean1$Stab_Para <- factor(final_stability_mean1$Stab_Para, levels = c("SiP", "SiN", "NP", "SiC", "CP", "CN", "C"))
Fig3 <- ggplot(dataxy, aes(x = Mean, y = Stab_Para, color = Treat)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_point(position = position_dodge(width = -0.5), size = 5) +
  geom_errorbar(aes(xmin = Mean - Sd, xmax = Mean + Sd), width = 0.3, position = position_dodge(width = -0.5))+
  scale_colour_manual(values = c("palegreen4", "skyblue3", "indianred"))+
  facet_grid(.~Stab_Para) +
  xlab("") +
  ylab("") +
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))+
  facet_grid(.~variable, scales = "free")
Fig3

ggsave("Fig3.eps", Fig3 , unit = "cm", device = "eps", width = 20, height = 6, dpi = 300)

#_______________________________________________________________________________

#FIG 4: Recovery across levels and within zooplankton

#### PICKED ZOOPLANKTON & CALCULATION OF ZOOPLANKTON FILTERS ####

#DATA HANDLING

data_zoo_picked <- read_excel("Data/3) Gepicktes Zooplankton/Zoo_Pick_SITES_2022_CN.xlsx")
data_zoo_picked$Individuals <- as.numeric(data_zoo_picked$Individuals)
data_zoo_picked <- filter(data_zoo_picked, Individuals > 0) #To exclude empty tin capsules!
data_zoo_picked <- dplyr::select(data_zoo_picked, Sampling_ID, Mesocosm_No, Treatment, Replicate, Species, Individuals, Parameter, DryWeight_mg, N_mg, C_mg, N_µg, C_µg)
data_zoo_picked$Parameter <- as.factor(data_zoo_picked$Parameter)
data_zoo_picked$Species <- as.factor(data_zoo_picked$Species)
data_zoo_picked$Replicate <- as.numeric(data_zoo_picked$Replicate)
data_zoo_picked$Treatment <- as.factor(data_zoo_picked$Treatment)
data_zoo_picked$Mesocosm_No <- as.numeric(data_zoo_picked$Mesocosm_No)
data_zoo_picked$Sampling_ID <- as.factor(data_zoo_picked$Sampling_ID)
data_zoo_picked_CN <- filter(data_zoo_picked, Parameter == "CN")
data_zoo_picked_P <- filter(data_zoo_picked, Parameter == "POP")

data_zoo_picked_CN$N_µg[data_zoo_picked_CN$N_µg == 0] <- 0.09
data_zoo_picked_CN <- mutate(data_zoo_picked_CN, Weight_mg_per_Individual = DryWeight_mg/Individuals)
data_zoo_picked_CN <- mutate(data_zoo_picked_CN, C_µg_per_Individual = C_µg/Individuals)
data_zoo_picked_CN <- mutate(data_zoo_picked_CN, N_µg_per_Individual = N_µg/Individuals)
data_zoo_picked_CN <- filter(data_zoo_picked_CN, N_µg < 5)

#Calculate mean and standard deviations
data_zoo_picked_CN <- ungroup(data_zoo_picked_CN)
data_zoo_picked_CN_mean <- dplyr::group_by(data_zoo_picked_CN, Sampling_ID, Treatment, Species)
data_zoo_picked_CN_mean <- dplyr::mutate(data_zoo_picked_CN_mean, Mean_N_µg_Ind = mean(N_µg_per_Individual), Sd_N_µg_Ind = sd(N_µg_per_Individual))
data_zoo_picked_CN_mean <- dplyr::mutate(data_zoo_picked_CN_mean, Mean_C_µg_Ind = mean(C_µg_per_Individual), Sd_C_µg_Ind = sd(C_µg_per_Individual))
data_zoo_picked_CN_mean$Sampling_ID <- factor(data_zoo_picked_CN_mean$Sampling_ID, levels = c("Start", "Mid", "End"))

#1) How much C, N, P was there per individual zooplankton? (calculated above)

#2) How many copepods and cladocerans were there per liter (Beni)?

data_zoo_beni <- read_excel("Data/8) Zooplankton Beni/Zooplankton_Beni.xlsx")
data_zoo_beni <- dplyr::select(data_zoo_beni, Lake, Mesocosm, Treatment, Pool, Sample, Group, Individual_L)
data_zoo_beni <- filter(data_zoo_beni, Lake == "Erken")
data_zoo_beni <- filter(data_zoo_beni, Group != "Rotifer")
data_zoo_beni <- filter(data_zoo_beni, Group != "cop_nauplii")
data_zoo_beni <- group_by(data_zoo_beni, Mesocosm, Treatment, Pool, Group, Sample)
data_zoo_beni_sum <- mutate(data_zoo_beni, sum_Ind_L = sum(Individual_L)) #Sum up the cladoceran groups (different species!)
data_zoo_beni_sum <- mutate(data_zoo_beni_sum, UniqueID = paste(Lake, Mesocosm, Treatment, Pool, Sample, Group, sep = "_"))
duplicated_rows <- duplicated(data_zoo_beni_sum$UniqueID)
data_zoo_beni_sum <- data_zoo_beni_sum[!duplicated_rows, ]
data_zoo_beni_end <- filter(data_zoo_beni_sum, Sample == "end")

data_zoo_beni_end$Lake <- as.factor(data_zoo_beni_end$Lake)
data_zoo_beni_end$Treatment <- as.factor(data_zoo_beni_end$Treatment)
data_zoo_beni_end$Pool <- as.numeric(data_zoo_beni_end$Pool)
data_zoo_beni_end$Sample <- as.factor(data_zoo_beni_end$Sample)
data_zoo_beni_end$Group <- as.factor(data_zoo_beni_end$Group)

#3) Calculate C, N, P in µmol/L for cladocerans/copepods!

data_zoo_1_end <- filter(data_zoo_picked_CN_mean, Sampling_ID == "End") 
data_zoo_1_end <- dplyr::select(data_zoo_1_end, Mesocosm = Mesocosm_No, Treatment, Pool = Replicate, Sample = Sampling_ID, Group = Species, C_µg_per_Individual, N_µg_per_Individual) 
data_zoo_1_end$Group <- gsub("Cladocerans", "Cladocera", data_zoo_1_end$Group) #Change spelling
data_zoo_1_end$Group <- gsub("Copepods", "Copepoda", data_zoo_1_end$Group) #Change spelling

data_zoo_beni_end <- filter(data_zoo_beni_end, Treatment != "LE")
data_zoo_beni_end <- data_zoo_beni_end[!(data_zoo_beni_end$Treatment == "D" & data_zoo_beni_end$Pool == 3 & data_zoo_beni_end$Group == "Copepoda"), ]
data_zoo_calc <- merge(data_zoo_1_end, data_zoo_beni_end, by = c("Mesocosm", "Treatment", "Pool", "Group"))

#Now: Calculate C, N, P in µg/L and µmol/L!
data_zoo_calc <- mutate(data_zoo_calc, C_µg_L = C_µg_per_Individual*sum_Ind_L, N_µg_L = N_µg_per_Individual*sum_Ind_L)
data_zoo_calc <- mutate(data_zoo_calc, C_µmol_L = C_µg_L/12.0107, N_µmol_L = N_µg_L/14.0067)
data_zoo_calc_mean <- group_by(data_zoo_calc, Treatment, Group)
data_zoo_calc_mean <- mutate(data_zoo_calc_mean, Mean_C_µmol_L = mean(C_µmol_L), Sd_C_µmol_L = sd(C_µmol_L), Mean_N_µmol_L = mean(N_µmol_L), Sd_N_µmol_L = sd(N_µmol_L))

#4) Calculate C, N, P for total zooplankton (needed for later calculations)!

data_zoo_calc_total <- group_by(data_zoo_calc, Mesocosm)
data_zoo_calc_total <- mutate(data_zoo_calc_total, C_µmol_L_Sum = sum(C_µmol_L), N_µmol_L_Sum = sum(N_µmol_L))
data_zoo_calc_total <- filter(data_zoo_calc_total, Group == "Cladocera") #Should be the same now!
data_zoo_calc_total_LRR <- data_zoo_calc_total

#Paste the value of the mean undisturbed (control) treatment in a new column
data_zoo_calc_total_LRR <- mutate(data_zoo_calc_total_LRR, Control_C = mean(data_zoo_calc_total_LRR$C_µmol_L_Sum[data_zoo_calc_total_LRR$Treatment == "C"]))
data_zoo_calc_total_LRR <- mutate(data_zoo_calc_total_LRR, CN = (C_µmol_L_Sum)/(N_µmol_L_Sum))
data_zoo_calc_total_LRR <- mutate(data_zoo_calc_total_LRR, Control_CN = mean(data_zoo_calc_total_LRR$CN[data_zoo_calc_total_LRR$Treatment == "C"]))
data_zoo_calc_total_LRR <- filter(data_zoo_calc_total_LRR, Treatment != "C") #Not necessary anymore!

#Calculate the LRR for Recovery = ln(Disturbed/Control) for C
data_zoo_calc_total_LRR_C <- mutate(data_zoo_calc_total_LRR, Recovery = log(C_µmol_L_Sum/Control_C), Stab_Para = "C")
data_zoo_calc_total_LRR_CN <- mutate(data_zoo_calc_total_LRR, Recovery = log(CN/Control_CN), Stab_Para = "CN")
data_zoo_total_LRR <- rbind(data_zoo_calc_total_LRR_C, data_zoo_calc_total_LRR_CN)

#5) Calculate LRR (for plot) and recovery for C and CN for both zooplankton groups!

#Paste the value of the mean undisturbed (control) treatment in a new column
data_zoo_calc_LRR <- mutate(data_zoo_calc, Control_C = mean(data_zoo_calc$C_µmol_L[data_zoo_calc$Treatment == "C"]))
data_zoo_calc_LRR <- mutate(data_zoo_calc_LRR, CN = (C_µmol_L)/(N_µmol_L))
data_zoo_calc_LRR <- mutate(data_zoo_calc_LRR, Control_CN = mean(data_zoo_calc_LRR$CN[data_zoo_calc_LRR$Treatment == "C"]))
data_zoo_calc_LRR <- filter(data_zoo_calc_LRR, Treatment != "C") #Not necessary anymore!

#Calculate the LRR for Recovery = ln(Disturbed/Control) for C
data_zoo_calc_LRR_C <- mutate(data_zoo_calc_LRR, Recovery = log(C_µmol_L/Control_C), Stab_Para = "C")
data_zoo_calc_LRR_CN <- mutate(data_zoo_calc_LRR, Recovery = log(CN/Control_CN), Stab_Para = "CN")
zoo_recovery <- rbind(data_zoo_calc_LRR_C, data_zoo_calc_LRR_CN)

#Calculate the mean
zoo_recovery <- filter(zoo_recovery, Mesocosm != 2)
#write_xlsx(zoo_recovery, "SITES_2022_Zoo_Recovery.xlsx")
zoo_recovery_mean <- group_by(zoo_recovery, Treatment, Group, Stab_Para)
zoo_recovery_mean <- mutate(zoo_recovery_mean, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
zoo_recovery_mean <- filter(zoo_recovery_mean, Pool == 4) #They should all be the same now!
zoo_recovery_mean <- dplyr::select(zoo_recovery_mean, Treatment, Pool, Group, Sample.x, Stab_Para, mean_Recovery, sd_Recovery)
zoo_recovery_mean1 <- pivot_longer(zoo_recovery_mean, cols = starts_with(c("mean_", "sd_")), names_to = c(".value", "variable"), names_pattern = "^(mean|sd)_([A-Za-z_]+)")

zoo_recovery_mean1$Treatment <- factor(zoo_recovery_mean1$Treatment, levels = c("D", "I", "E"))
Fig4_Part2 <- ggplot(zoo_recovery_mean1, aes(x = mean, y = variable, color = Treatment, shape = Group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_point(position = position_dodge(width = -0.5), size = 5) +
  scale_shape_manual(values=c(15, 17))+
  geom_errorbar(aes(xmin = mean - sd, xmax = mean + sd), width = 0.3, position = position_dodge(width = -0.5))+
  scale_colour_manual(values = c("palegreen4", "skyblue3", "indianred"))+
  facet_grid(.~Stab_Para, scales = "free") +
  xlab("") +
  ylab("") +
  #annotate("text", x = -0.45, y = 2.5, label = "Destabilising", color = "black", size = 3, family = "Times")+
  #annotate("text", x = 0.45, y = 2.5, label = "Stabilising", color = "black", size = 3, family = "Times")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))
Fig4_Part2

ggsave("Fig4_Part2.eps", Fig4_Part2, unit = "cm", device = "eps", width = 25, height = 5, dpi = 300)

#6) How much of the C, N, P on the large filter-fraction is NOT zooplankton?
#C, N, P on zooplankton filter (µmol/L) - C, N, P upscaled from picked zooplankton!

#Read data and preparation of main data set
data <- read_excel("R Scripts/Packed_CN_Phy_Erken_Uebertragen.xlsx")
data <- filter(data, Meso != "Add")
data <- filter(data, Meso != "Lake")
data <- filter(data, Meso != "Blank")
data <- filter(data, Level == "Zoo")
data <- filter(data, S_No == 10)
data <- mutate(data, Rep = 1)
data$Rep[data$Meso == 5 | data$Meso == 6 | data$Meso == 7 | data$Meso == 8] <- 2
data$Rep[data$Meso == 9 | data$Meso == 10 | data$Meso == 11 | data$Meso == 12] <- 3
data$Rep[data$Meso == 13 | data$Meso == 14 | data$Meso == 15 | data$Meso == 16] <- 4

data <- dplyr::select(data, Meso, Rep, Treat, S_No, N_µmol_L_Filter = N_µmol_L, C_µmol_L_Filter = C_µmol_L)
data$Meso <- as.numeric(data$Meso)
data$Treat <- as.factor(data$Treat)
data$N_µmol_L_Filter <- as.numeric(data$N_µmol_L_Filter)
data$C_µmol_L_Filter <- as.numeric(data$C_µmol_L_Filter)

data_zoo_calc_total <- dplyr::select(data_zoo_calc_total, Meso = Mesocosm, Rep = Pool, Treat = Treatment, C_µmol_L_Sum, N_µmol_L_Sum)
large_phyto <- merge(data, data_zoo_calc_total, by = c("Meso", "Rep", "Treat"))

#Substract the calculated zooplankton fraction from the zooplankton filters to get the
#rest on the filter (Gloeotrichia + large phytoplankton)

large_phyto <- mutate(large_phyto, Rest_C_µmol_L = C_µmol_L_Filter-C_µmol_L_Sum, Rest_N_µmol_L = N_µmol_L_Filter-N_µmol_L_Sum)
large_phyto$Rest_C_µmol_L[large_phyto$Rest_C_µmol_L < 0] <- 0 #Set everything smaller than 0 to 0 (1 case)

#Paste the value of the mean undisturbed (control) treatment in a new column
large_phyto_LRR <- mutate(large_phyto, Control_C = mean(large_phyto$Rest_C_µmol_L[large_phyto$Treat == "C"]))
large_phyto_LRR <- mutate(large_phyto_LRR, CN = (Rest_C_µmol_L)/(Rest_N_µmol_L))
large_phyto_LRR <- mutate(large_phyto_LRR, Control_CN = mean(large_phyto_LRR$CN[large_phyto_LRR$Treat == "C"]))
large_phyto_LRR <- filter(large_phyto_LRR, Treat != "C") #Not necessary anymore!

#Calculate the LRR for Recovery = ln(Disturbed/Control) for C
large_phyto_LRR_C <- mutate(large_phyto_LRR, Recovery = log(Rest_C_µmol_L/Control_C), Stab_Para = "C")
large_phyto_LRR_CN <- mutate(large_phyto_LRR, Recovery = log(CN/Control_CN), Stab_Para = "CN")
large_phyto_recovery <- rbind(large_phyto_LRR_C, large_phyto_LRR_CN)

#write_xlsx(large_phyto_recovery, "SITES_2022_LargePhyto_Recovery.xlsx")

#Calculate the mean
large_phyto_mean <- group_by(large_phyto_recovery, Treat, Stab_Para)
large_phyto_mean <- mutate(large_phyto_mean, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
large_phyto_mean <- filter(large_phyto_mean, Rep == 4) #They should all be the same now!
large_phyto_mean <- dplyr::select(large_phyto_mean, Treat, Rep, Stab_Para, mean_Recovery, sd_Recovery)
large_phyto_mean1 <- pivot_longer(large_phyto_mean, cols = starts_with(c("mean_", "sd_")), names_to = c(".value", "variable"), names_pattern = "^(mean|sd)_([A-Za-z_]+)")

#Put data sets together
AllGroups_Recovery <- read_excel("R Scripts/AllGroups_Recovery.xlsx")
AllGroups_Recovery$Stab_Para <- as.factor(AllGroups_Recovery$Stab_Para)
AllGroups_Recovery$Size_Fraction <- as.factor(AllGroups_Recovery$Size_Fraction)
AllGroups_Recovery$Treatment <- as.factor(AllGroups_Recovery$Treatment)
AllGroups_Recovery$Mesocosm <- as.factor(AllGroups_Recovery$Mesocosm)
AllGroups_Recovery$Pool <- as.factor(AllGroups_Recovery$Pool)
AllGroups_Recovery$Recovery <- as.numeric(AllGroups_Recovery$Recovery)

#Standardisation: Calculate absolute values (see Hillebrand et al. 2019)
AllGroups_Recovery <- mutate(AllGroups_Recovery, Recovery_Abs = abs(Recovery))
AllGroups_Recovery_C <- filter(AllGroups_Recovery, Stab_Para == "C")
AllGroups_Recovery_CN <- filter(AllGroups_Recovery, Stab_Para == "CN")

#Calculate mean
AllGroups_Recovery <- filter(AllGroups_Recovery, Mesocosm != 2)
AllGroups_Recovery_mean <- mutate(AllGroups_Recovery, ID = paste(Size_Fraction, Mesocosm, Stab_Para, sep = "_"))
AllGroups_Recovery_mean <- filter(AllGroups_Recovery_mean, ID != "Zoo_2_CN")
AllGroups_Recovery_mean <- group_by(AllGroups_Recovery_mean, Size_Fraction, Treatment, Stab_Para) 
AllGroups_Recovery_mean <- mutate(AllGroups_Recovery_mean, Recovery_Mean = mean(Recovery), Recovery_Sd = sd(Recovery))
AllGroups_Recovery_mean <- filter(AllGroups_Recovery_mean, Pool == 4)

#Plot recovery of all groups
AllGroups_Recovery_mean$Treatment <- factor(AllGroups_Recovery_mean$Treatment, levels = c("D", "I", "E"))
AllGroups_Recovery_mean$Size_Fraction <- factor(AllGroups_Recovery_mean$Size_Fraction, levels = c("Zoo", "LargePhyto", "SmallPhyto"))
Fig4_Part1 <- ggplot(AllGroups_Recovery_mean, aes(x = Recovery_Mean, y = Size_Fraction, color = Treatment)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_point(position = position_dodge(width = -0.5), size = 9) +
  scale_shape_manual(values=c(15, 17))+
  geom_errorbar(aes(xmin = Recovery_Mean - Recovery_Sd, xmax = Recovery_Mean + Recovery_Sd), width = 0.4, position = position_dodge(width = -0.5))+
  scale_colour_manual(values = c("palegreen4", "skyblue3", "indianred"))+
  facet_grid(.~Stab_Para, scales = "free") +
  xlab("") +
  ylab("") +
  #annotate("text", x = -0.45, y = 2.5, label = "Destabilising", color = "black", size = 3, family = "Times")+
  #annotate("text", x = 0.45, y = 2.5, label = "Stabilising", color = "black", size = 3, family = "Times")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18))
Fig4_Part1

ggsave("Fig4_Part1.eps", Fig4_Part1, unit = "cm", device = "eps", width = 25, height = 12, dpi = 300)

#_______________________________________________________________________________

