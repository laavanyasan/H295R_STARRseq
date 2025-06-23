###qPCR for H295R-dCas9-p300
##2024-11

#clear workspace
rm(list=ls())

#set directory
setwd("")

#load packages
library(gdata)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(viridis)
library(rstatix)
library(DescTools)
library(devtools)
library(readxl)

##Note sheets in the source data file are : "Figure 5c". "Figure 5d", "Supplementary Figure 28" and "H295R Samples" respectively
#for loading the following tables. 

data <- read_excel("H295Rsamples.xlsx", col_names = TRUE, sheet="Psamples")
SampleNames <- read_excel("Sample_Names.xlsx", col_names=c("Row", "Rep", "Sample", "Treatment"), sheet="Psamples")

data$SampleName <- NA
data$SampleTrt <- NA
data$SampleRep <- NA
for (i in 1:nrow(data)){
  for (j in 1:nrow(SampleNames)){
    if (data$Sample.Name[i]==toString(SampleNames$Row[j])){
      data$SampleName[i] <- toString(SampleNames$Sample[j])
      data$SampleTrt[i] <- toString(SampleNames$Treatment[j])
      data$SampleRep[i] <- toString(SampleNames$Rep[j])
      #data$Quantity[i] <- as.numeric(SampleNames$Dilution[j])
    }
  }
}


##MAKE ALL ct'S NUMERICAL
##REMOVING SAMPLES THAT ARE PROBLEMS

df.filter <- data %>% filter (Sample.Name != "P1")
df.filter$CT.STRBP <- as.numeric(df.filter$CT.STRBP)
df.filter$CT.LHX2 <- as.numeric(df.filter$CT.LHX2)
df.filter$CT.GAPDH <- as.numeric(df.filter$CT.GAPDH)
df.filter$CT.DENND1A <- as.numeric(df.filter$CT.DENND1A)
df.filter$CT.CYP17A1 <- as.numeric(df.filter$CT.CYP17A1)


##Convert CT's to avergae CT's per tech replicate. Then do the first internal normalization by subtracting GAPDH for each sample

df.1 <- df.filter %>% group_by(SampleRep,SampleTrt) %>% summarise(aCT.STRBP = mean(CT.STRBP, na.rm=TRUE), aCT.LHX2 = mean(CT.LHX2, na.rm=TRUE), aCT.GAPDH = mean(CT.GAPDH,na.rm=TRUE), aCT.DENND1A = mean(CT.DENND1A,na.rm=TRUE),aCT.CYP17A1 = mean(CT.CYP17A1,na.rm=TRUE),Sample.Rep = SampleRep, Sample.Name = SampleName)
df.1 <- unique(df.1)
df.1 <- df.1 %>% mutate(delCT.STRBP = aCT.STRBP - aCT.GAPDH) %>% mutate(delCT.LHX2 = aCT.LHX2 - aCT.GAPDH) %>% mutate(delCT.DENND1A = aCT.DENND1A - aCT.GAPDH) %>% mutate(delCT.CYP17A1 = aCT.CYP17A1 - aCT.GAPDH)


##Get the second normalizer (refCT) for each neg control sample


refCT1 <- mean((df.1 %>% filter(SampleTrt=="DMSO") %>% filter(Sample.Name == "SCR"))$delCT.STRBP, na.rm = TRUE)
refCT2 <- mean((df.1 %>% filter(SampleTrt=="DMSO") %>% filter(Sample.Name == "SCR"))$delCT.LHX2, na.rm = TRUE)
refCT3 <- mean((df.1 %>% filter(SampleTrt=="DMSO") %>% filter(Sample.Name == "SCR"))$delCT.DENND1A, na.rm = TRUE)
refCT4 <- mean((df.1 %>% filter(SampleTrt=="DMSO") %>% filter(Sample.Name == "SCR"))$delCT.CYP17A1, na.rm = TRUE)

##Get 2^ fold change for each gene in each samples

d1 <- df.1

d1 <- d1 %>% mutate(raisedFC.STRBP = 2^(-(delCT.STRBP - refCT1))) %>% mutate(raisedFC.CYP17A1 = 2^(-(delCT.CYP17A1 - refCT4))) %>% mutate(raisedFC.DENND1A = 2^(-(delCT.DENND1A - refCT3))) %>% mutate(raisedFC.LHX2 = 2^(-(delCT.LHX2 - refCT2)))

##Gather df to then make boxplot

df.6 <- d1 %>% gather("Primer", "FC", 14:17)
df.6 <- df.6 %>% mutate(l2fc = log2(FC))

##Make boxplots per gene
boxplot <- ggplot(data = (df.6 %>% filter(Primer %in% c("raisedFC.STRBP", "raisedFC.LHX2")))) +
  geom_boxplot(aes(y = log2(FC), x=Sample.Name, colour = SampleTrt)) + theme_classic()  + facet_wrap(~Primer) +
  geom_point(aes(y = log2(FC), x=Sample.Name, colour = SampleTrt), position = position_dodge(width = 0.8)) + xlab("H295R-dCas9-KRAB") + ylab("log2 FC") +
  theme_classic() +   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(-4,4) +
  scale_colour_brewer(palette = "Dark2")

boxplot

########################
########################
###qPCR for H295R-dCas9-KRAB
##2024-11

#clear workspace
rm(list=ls())

#set directory

setwd("")

#load packages
library(gdata)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(viridis)
library(rstatix)
library(DescTools)
library(devtools)
library(readxl)


##Note sheets in the source data file are : "Supplementary FIgure 28" and 
#"Supplementary FIgure 28 Samples" respectively
#for loading the following tables. 

data <- read_excel("H295Rsamples.xlsx", col_names = TRUE, sheet="Ksamples")

SampleNames <- read_excel("Sample_Names.xlsx", col_names=c("Row", "Rep", "Sample", "Treatment"), sheet="Ksamples")

data$SampleName <- NA
data$SampleTrt <- NA
data$SampleRep <- NA
for (i in 1:nrow(data)){
  for (j in 1:nrow(SampleNames)){
    if (data$Sample.Name[i]==toString(SampleNames$Row[j])){
      data$SampleName[i] <- toString(SampleNames$Sample[j])
      data$SampleTrt[i] <- toString(SampleNames$Treatment[j])
      data$SampleRep[i] <- toString(SampleNames$Rep[j])
      #data$Quantity[i] <- as.numeric(SampleNames$Dilution[j])
    }
  }
}


df.filter <- data
df.filter$CT.STRBP <- as.numeric(df.filter$CT.STRBP)
df.filter$CT.LHX2 <- as.numeric(df.filter$CT.LHX2)
df.filter$CT.GAPDH <- as.numeric(df.filter$CT.GAPDH)
df.filter$CT.DENND1A <- as.numeric(df.filter$CT.DENND1A)
df.filter$CT.CYP17A1 <- as.numeric(df.filter$CT.CYP17A1)
#df.filter <- na.omit(df.filter)


df.1 <- df.filter %>% group_by(SampleRep,SampleTrt) %>% summarise(aCT.STRBP = mean(CT.STRBP, na.rm=TRUE), aCT.LHX2 = mean(CT.LHX2, na.rm=TRUE), aCT.GAPDH = mean(CT.GAPDH,na.rm=TRUE), aCT.DENND1A = mean(CT.DENND1A,na.rm=TRUE),aCT.CYP17A1 = mean(CT.CYP17A1,na.rm=TRUE),Sample.Rep = SampleRep, Sample.Name = SampleName)
df.1 <- unique(df.1)
df.1 <- df.1 %>% mutate(delCT.STRBP = aCT.STRBP - aCT.GAPDH) %>% mutate(delCT.LHX2 = aCT.LHX2 - aCT.GAPDH) %>% mutate(delCT.DENND1A = aCT.DENND1A - aCT.GAPDH) %>% mutate(delCT.CYP17A1 = aCT.CYP17A1 - aCT.GAPDH)

refCT1 <- mean((df.1 %>% filter(SampleTrt=="DMSO") %>% filter(Sample.Name == "SCR"))$delCT.STRBP, na.rm = TRUE)
refCT2 <- mean((df.1 %>% filter(SampleTrt=="DMSO") %>% filter(Sample.Name == "SCR"))$delCT.LHX2, na.rm = TRUE)
refCT3 <- mean((df.1 %>% filter(SampleTrt=="DMSO") %>% filter(Sample.Name == "SCR"))$delCT.DENND1A, na.rm = TRUE)
refCT4 <- mean((df.1 %>% filter(SampleTrt=="DMSO") %>% filter(Sample.Name == "SCR"))$delCT.CYP17A1, na.rm = TRUE)

d1 <- df.1

d1 <- d1 %>% mutate(raisedFC.STRBP = 2^(-(delCT.STRBP - refCT1))) %>% mutate(raisedFC.CYP17A1 = 2^(-(delCT.CYP17A1 - refCT4))) %>% mutate(raisedFC.DENND1A = 2^(-(delCT.DENND1A - refCT3))) %>% mutate(raisedFC.LHX2 = 2^(-(delCT.LHX2 - refCT2)))

##Remove SS1-2 LHX value (37!) for qPCR
d1$aCT.LHX2 <- replace(d1$aCT.LHX2, d1$aCT.LHX2 > 36, NA)
(d1 %>% filter(SampleRep =="SS1-2") %>% filter(SampleTrt == "Forskolin"))$aCT.LHX2


##Gather df to then make boxplot

df.6 <- d1 %>% gather("Primer", "FC", 14:17)
df.6 <- df.6 %>% mutate(l2fc = log2(FC))


boxplot <- ggplot(data = (df.6 %>% filter(Primer=="raisedFC.DENND1A"))) +
  geom_boxplot(aes(y = log2(FC), x=Sample.Name, colour = SampleTrt)) + theme_classic()  + facet_wrap(~Primer) +
  geom_point(aes(y = log2(FC), x=Sample.Name, colour = SampleTrt), position = position_dodge(width = 0.8)) + xlab("H295R-dCas9-KRAB") + ylab("log2 FC") +
  theme_classic() +   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(-3,3) +
  scale_colour_brewer(palette = "Dark2")

boxplot

