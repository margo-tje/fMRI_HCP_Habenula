rm(list=ls()) #remove all variables in workspace
source("summarySE.R") #load the user defined function summarySE 

install.packages("reshape2")
install.packages("ggpubr")
install.packages("effectsize")
install.packages("rstatix")

#load all the libraries you need (after installing them)
library(reshape2)
library(lme4)
library(lmerTest)
library(readxl)
library(writexl)
library(xlsx)
library(ggplot2)  
library(coda)
library(plyr)
library(gridExtra)
library(car)
library(tidyverse)
library(knitr)
library(effects)
library(psych)
library(candisc)
library(ggpubr)
library(effectsize)
library(rstatix)

cbPalette <- c(  "#CC79A7","#b3b3b3","#0072B2", "#D55E00","#56B4E9","#E69F00",  "#009E73", "#F0E442")
cbPalette2 <- c( "#56B4E9","#E69F00",  "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000")


#BMI Dataset (n=300)
Copes <- read_excel("L:/basic/divg/MSlomp_fMRI/tfMRI_MixedModel/GamblingTask_Copes1-6_BMI.xlsx")
Copes$FD_After <- Copes$FDAfter
Copes$Hb_SegAnatom_Bi <- Copes$Hb_SegAnatom_L + Copes$Hb_SegAnatom_R #Calculate bilateral volume

#OR HbA1c dataset (n=72subs)
Copes <- read_excel("L:/basic/divg/MSlomp_fMRI/tfMRI_MixedModel/GamblingTask_Copes1-6_HbA1C.xlsx")
Copes$Hb_SegAnatom_Bi <- Copes$Hb_volL + Copes$Hb_volR #Calculate bilateral volume

#Make sure all variables have the right type
Copes$SubID <- as.character(Copes$SubID)
Copes$Smoking = factor(Copes$SSAGA_TB_Still_Smoking, levels = c(0,1), labels = c("No","Yes"))
Copes$HbA1C <- as.numeric(Copes$HbA1C) #Only a few values in set; is it okay to treat as numeric?
Copes$Hematocrit_1 <- as.numeric(Copes$Hematocrit_1)
Copes$Hematocrit_2 <- as.numeric(Copes$Hematocrit_2)
Copes$HematocritAvg <- (Copes$Hematocrit_1 + Copes$Hematocrit_2 ) / 2 #Average Hematocrit over the two days
Copes$Gender <- factor(Copes$Gender)
Copes$HbCat <- ifelse(Copes$HbA1C >= 5.7, "1","0")
Copes$HbCat = factor(Copes$HbCat, levels = c(0,1), labels = c("LowHb","HighHb"))
Copes$BMICat <- ifelse(Copes$BMI >= 30, "2", ifelse(Copes$BMI < 25, "0","1"))
Copes$BMICat = factor(Copes$BMICat, levels = c(0,1,2), labels = c("Lean","Overweight", "Obese"))

Copes$Gambling_Task_Punish_Median_RT_Smaller <- as.numeric(Copes$Gambling_Task_Punish_Median_RT_Smaller)
Copes$Gambling_Task_Punish_Median_RT_Larger <- as.numeric(Copes$Gambling_Task_Punish_Median_RT_Larger)
Copes$Gambling_Task_Reward_Median_RT_Smaller <- as.numeric(Copes$Gambling_Task_Reward_Median_RT_Smaller)
Copes$Gambling_Task_Reward_Median_RT_Larger <- as.numeric(Copes$Gambling_Task_Reward_Median_RT_Larger)

Copes$Gambling_Task_Punish_Perc_Smaller <- as.numeric(Copes$Gambling_Task_Punish_Perc_Smaller)
Copes$Gambling_Task_Punish_Perc_Larger <-  as.numeric(Copes$Gambling_Task_Punish_Perc_Larger)
Copes$Gambling_Task_Reward_Perc_Smaller <- as.numeric(Copes$Gambling_Task_Reward_Perc_Smaller)
Copes$Gambling_Task_Reward_Perc_Larger <-  as.numeric(Copes$Gambling_Task_Reward_Perc_Larger)

## From wide to long format
data_long <- melt(Copes,id.vars=c("SubID","Age_in_Yrs","Gender","Handedness","Race","BMI","BMICat","HematocritAvg","BPSystolic","BPDiastolic","ThyroidHormone","HbA1C","HbCat","SES","THC","Smoking","Sadness_Unadj","FD_After","Hb_SegAnatom_Bi","Gambling_Task_Punish_Perc_Larger","Gambling_Task_Punish_Perc_Smaller","Gambling_Task_Punish_Median_RT_Smaller","Gambling_Task_Punish_Median_RT_Larger"),measure.vars = c("Cope1_HbL","Cope1_HbR","Cope1_HbBi","Cope2_HbL","Cope2_HbR","Cope2_HbBi","Cope3_HbL","Cope3_HbR","Cope3_HbBi"))
data_long$value <- as.numeric(data_long$value)

#Cope 1 = Punishmnt, Cope 2 = Reward, Cope 3= Punishment vs Reward
Cope1_Bi = subset(data_long, variable=="Cope1_HbBi" )
Cope2_Bi = subset(data_long, variable=="Cope2_HbBi" )
Cope3_Bi = subset(data_long, variable=="Cope3_HbBi" )
Cope123_Bi = rbind(Cope1_Bi,Cope2_Bi,Cope3_Bi)
Cope12_Bi = rbind(Cope1_Bi,Cope2_Bi)

print(t.test(Cope1_Bi$value,Cope2_Bi$value, paired=TRUE)) #Paired t-test : Activity during punishment different than during reward?

cohens_d(Cope1_Bi$value,Cope2_Bi$value,paired=TRUE)

#Test if punishment and reward are significantly different.

##################################
# Create bar plot with mean +- SD for punishment and reward
tBMI <- summarySE(Cope12_Bi, measurevar="value", groupvars=c("variable"),na.rm=TRUE)  #Summary for the differnt groups
tBMI$variable <- factor(tBMI$variable,labels = c("Punishment","Reward"))

ggplot(tBMI, aes(y=value, x=variable, fill = variable)) + geom_bar(position=position_dodge(), stat="identity", show.legend= FALSE) + geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2, position=position_dodge(.9)) + xlab("") + ylab("Z-statistic") +  theme_bw() + theme(text = element_text(size=30) )   + scale_fill_manual(values=cbPalette,labels = labels)


##   300 pps, BMI set
#Save
ggsave("L:/basic/divg/MSlomp_fMRI/tfMRI_MixedModel/NoOutlierRm_Gambling_300pps_PunRew-Bi_Barplot.png",dpi=300,width=5,height=6)

#Run model
model1 <- lm(value ~ BMI + HbA1C  + Age_in_Yrs    + Gender   ,data=Cope3_Bi) 
print(summary(model1))


#Same, but for 72 pps HbA1c set
ggsave("L:/basic/divg/MSlomp_fMRI/tfMRI_MixedModel/NoOutlierRm_Gambling_72pps_PunRew-Bi_Barplot.png",dpi=300,width=5,height=6)

model1 <- lm(value ~ BMI * HbA1C  + Age_in_Yrs    + Gender   ,data=Cope3_Bi) 
print(summary(model1))

#Optional:
model2 <- lm(value ~ BMI * HbA1C  + Age_in_Yrs    + Gender  + FD_After ,data=Cope3_Bi) 
print(summary(model2))

anova(model1,model2) #Compared the model with and without FD

