rm(list=ls()) #remove all variables in workspace

#load all the libraries you need (after installing them)
install.packages("lme4")
install.packages("lmerTest")
install.packages("readxl")
install.packages("writexl")
install.packages("xls")
install.packages("ggplot2")
install.packages("coda")
install.packages("plyr")
install.packages("gridExtra")
install.packages("car")
install.packages("tidyverse")
install.packages("knitr")
install.packages("effects")
install.packages("psych")
install.packages("candisc")
install.packages("Rmisc")
install.packages("reshape2")
install.packages("sjPlot")
install.packages("sjmisc")
install.packages("ggplot2")
install.packages("MASS")
install.packages("faraway")
install.packages("plyr")
install.packages("MHTdiscrete")
install.packages("car")

library(MHTdiscrete)
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
library(reshape2)
library(Rmisc)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(MASS)
library(faraway)
library(car)

cbPalette <- c( "#56B4E9","#E69F00",  "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000")
data <- read_excel({file}) # Load the data

data$HbCat <- ifelse(data$HbA1C >= 5.7, "1","0")
data$HbCat = factor(data$HbCat, levels = c(0,1), labels = c("LowHb","HighHb"))
data$BMICat <- ifelse(data$BMI >= 30, "2", ifelse(data$BMI < 25, "0","1"))
data$BMICat = factor(data$BMICat, levels = c(0,1,2), labels = c("Lean","Overweight", "Obese"))

#Make sure all variables have the right type
data$SubID <- as.character(data$SubID)
data$HbA1C <- as.numeric(data$HbA1C) #Only a few values in set; is it okay to treat as numeric?
data$Thyroid <- as.numeric(data$Thyroid) #Only a few values in set; is it okay to treat as numeric?
data$Hematocrit_Avg <- as.numeric(data$Hematocrit_Avg)
data$Gender <- factor(data$Gender)
data$BPSys <- as.numeric(data$BPSys) 
data$BPDia <- as.numeric(data$BPDia) 
data$Source <- factor(data$Source)
data$HbVol_Bi <- data$HabenulaVol_L + data$HabenulaVol_R
data$Smoking <- replace(data$Smoking, data$Smoking == NaN, 0) # REPLACE THE SMOKING NAN WITH 0
data$Smoking <- as.factor(data$Smoking)

data$DepresCat <- ifelse(data$Depres >= 45, "1","0")
data$DepresCat <- factor(data$DepresCat, levels = c(0,1), labels = c("LowDepres","HighDepres"))

## From wide to long format#
data_long <- melt(data,id.vars=c("SubID","Age","Gender","Handedness","Race","BMI","BMICat","Hematocrit_Avg","BPSys","THC","BPDia","Thyroid","HbA1C","HbCat","SES","Smoking","Depres","DepresCat","MeanFD","HbVol_Bi","Source"),measure.vars = c("SeedHb","Ctr_DM","Ctr_CM","pgACC","Insula","NAc","AMY","HPC","VTA","HYP"))
data_long$value <- as.numeric(data_long$value) 

#The sources and connections to loop through
connections <- c("pgACC", "Insula", "NAc","AMY","VTA","HYP" )   ##HPC excluded
connectionNames  <- c("Pregenual Anterior Cingulate Cortex", "Insula", "Nucleus Accumbens (NAc)","Amygdala","Ventral Tegmental Area (VTA)","Hypothalamus" )   ##HPC excluded
sources <- c("SeedHb","Ctr_CM","Ctr_DM")

#Prepare data frames to save values out of the loop
pvals <- data.frame(matrix(nrow = length(connections), ncol = 5))
fstat <- data.frame(matrix(nrow = length(connections), ncol = 4))
names(fstat) <- c("Connection","FStat","df1","df2")
overall_pval <- data.frame(matrix(nrow = length(connections), ncol = 2))
overall_rsquare <- data.frame(matrix(nrow = length(connections), ncol = 3))
names(overall_rsquare) <- c("Connection","R2","AdjustedR2")
names(pvals) <- c("Intercept","BMI","HbA1C","GenderM","Age")
tvals <- pvals
sderror <- pvals
estimates <- pvals
cor_pvals <- pvals #p-values that will be corrected
cor <- data.frame(matrix(nrow = length(connections), ncol = 8))
onesample <- data.frame(matrix(nrow = length(connections), ncol = 3))
names(onesample) <- c("Connection","t-value","p-value")
names(cor) <- c("","BMI","HbA1c","Lowhb_BMI","HighHb_BMI","Lean_HbA1c","Overw_HbA1c","Obese_HbA1c")

s=1 #1=Hb, 2=CM, 3=DM
SourceData = subset(data_long, Source==sources[s])
SourceData <- na.omit(SourceData)

for (j in 1:length(connections)){
  connectionData = subset(SourceData, variable==connections[j])
  connectionData <- na.omit(connectionData) #remove rows with missing values (Hb-Hb, CM-CM and DM-DM)

  #Run the linear model
  print(connections[j])
  model.full1 <- lm(value ~ BMI + HbA1C  + Gender +  Age   ,data=connectionData)
  print(summary(model.full1))

  # Plot data including regression line using intercept + slope
  coeff<-coefficients(model.full1)     #get intercept and slope value
  intercept<-coeff[1]
  slope<- coeff[2]

  # Plot
  p5 <- ggplot(connectionData, aes(x=BMI,  y=value))  + geom_point(size=2) +  ylab('Correlation Strength') + scale_colour_manual(values=cbPalette) + theme_bw(base_size=25) + theme(legend.position = c(0.8,0.8)) +theme(text = element_text(family = "sans")) + labs(title =connectionNames[j]) + xlab("BMI (kg/mm2)") + ylim(-0.06, 0.06) +
    geom_abline(intercept = intercept, slope = slope, color="red", linetype="dashed", size=1.5)
   plot(p5)
  
  #Save the plot
   ggsave(paste({outputfile},connections[j],"Scatter.png",sep="_"),dpi=300,width=8,height=6)
   
  #Save p values to correct for multiple comparisons
  pvals[j,] <- coef(summary(model.full1))[,4] #p-value
  tvals[j,] <- coef(summary(model.full1))[,3] #t-value
  sderror[j,] <- coef(summary(model.full1))[,2] #standard error
  estimates[j,] <- coef(summary(model.full1))[,1] #estimates
  
  # Save F-statistic, df etc
  fstat[j,2:4] <- summary(model.full1)$fstatistic #F-statistic, df's
  overall_pval[j,2] <- pf(fstat[j,2], fstat[j,3], fstat[j,4], lower.tail=FALSE) #pvalue
  overall_rsquare[j,2] <- summary(model.full1)$r.squared
  overall_rsquare[j,3] <- summary(model.full1)$adj.r.squared
 
  # Add names
  sderror[j,1] <- connections[j]
  estimates[j,1] <- connections[j]
  tvals[j,1] <- connections[j]
  pvals[j,1] <- connections[j]
  cor_pvals[j,1] <- connections[j]
  fstat[j,1]  <- connections[j]
  overall_pval[j,1] <- connections[j]
  overall_rsquare[j,1] <- connections[j]
  
}
  
#Bonferoni corrected
cor_pvals[1:6,2] <- p.adjust(pvals[1:6,2],method="fdr") 
cor_pvals[1:6,3] <- p.adjust(pvals[1:6,3],method="fdr") 
cor_pvals[1:6,4] <- p.adjust(pvals[1:6,4],method="fdr") 
cor_pvals[1:6,5] <- p.adjust(pvals[1:6,5],method="fdr")  

#Save all output to 1 Excel file (append sheets)
write.xlsx(estimates,{outputfile},sheetName="Estimates")
write.xlsx(tvals,{outputfile},sheetName="Tvalues",append=TRUE)
write.xlsx(sderror,{outputfile},sheetName="SdError",append=TRUE)
write.xlsx(pvals,{outputfile},sheetName="UncorrectedP",append=TRUE)
write.xlsx(cor_pvals,{outputfile},sheetName="BonferoniCorrectedP",append=TRUE)
write.xlsx(fstat,{outputfile},sheetName="FStats",append=TRUE)
write.xlsx(overall_pval,{outputfile},sheetName="ModelP",append=TRUE)
write.xlsx(overall_rsquare,{outputfile},sheetName="R2",append=TRUE)