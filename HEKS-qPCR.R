###qPCR for HEK293T-dCas9-p300
##2023-08

#clear workspace
rm(list=ls())

#set directory
setwd("~/OneDrive - Duke University/ReddyLab/Experiments/20230814-qPCRs")

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

data <- read_excel("HEK-qpcrs.xlsx", col_names = TRUE)
SampleNames <- read_excel("Sample_Names.xlsx", col_names=c("Row","Sample", "Rep"), sheet="Sample2")

data$SampleName <- NA
data$SampleRep <- NA
for (i in 1:nrow(data)){
  for (j in 1:nrow(SampleNames)){
    if (data$Sample.Name[i]==toString(SampleNames$Row[j])){
      data$SampleName[i] <- toString(SampleNames$Sample[j])
      data$SampleRep[i] <- toString(SampleNames$Rep[j])
      #data$Quantity[i] <- as.numeric(SampleNames$Dilution[j])
    }
  }
}


df.filter <- data
df.filter$CT.SRBP <- as.numeric(df.filter$CT.SRBP)
df.filter$CT.LHX <- as.numeric(df.filter$CT.LHX)
df.filter$CT.GAPDH <- as.numeric(df.filter$CT.GAPDH)
df.filter$CT.DENND1A <- as.numeric(df.filter$CT.DENND1A)


df.1 <- df.filter %>% group_by(SampleName) %>% summarise(aCT.SRBP = mean(CT.SRBP, na.rm=TRUE), aCT.LHX = mean(CT.LHX, na.rm=TRUE), aCT.GAPDH = mean(CT.GAPDH,na.rm=TRUE), aCT.DENND1A = mean(CT.DENND1A,na.rm=TRUE), SampleRep = SampleRep)
df.2 <- unique(df.1)
df3 <- df.filter %>% group_by(SampleRep) %>% summarise(aCT.SRBP = mean(CT.SRBP, na.rm=TRUE), aCT.LHX = mean(CT.LHX, na.rm=TRUE), aCT.GAPDH = mean(CT.GAPDH,na.rm=TRUE), aCT.DENND1A = mean(CT.DENND1A,na.rm=TRUE), SampleRep = SampleRep)
df.3 <- unique(df3)
df.3 <- df.3 %>% mutate(delCT.SRBP = aCT.SRBP - aCT.GAPDH) %>% mutate(delCT.LHX = aCT.LHX - aCT.GAPDH) %>% mutate(delCT.DENND1A = aCT.DENND1A - aCT.GAPDH)

refCT1 <- 6.897949
refCT2 <- 14.04702
refCT3 <- 8.292393

df.3 <- df.3 %>% mutate(raisedFC.SRBP = 2^(-(delCT.SRBP - refCT1))) %>% mutate(raisedFC.LHX = 2^(-(delCT.LHX - refCT2))) %>% mutate(raisedFC.DENND1A = 2^(-(delCT.DENND1A - refCT3)))

df.4 <- df.3 %>% gather("Primer", "FC", 9:11) %>% filter(Primer != "raisedFC.SRBP")


df5 <- df.filter %>% group_by(SampleName) %>% summarise(aCT.SRBP = mean(CT.SRBP, na.rm=TRUE), aCT.LHX = mean(CT.LHX, na.rm=TRUE), aCT.GAPDH = mean(CT.GAPDH,na.rm=TRUE), aCT.DENND1A = mean(CT.DENND1A,na.rm=TRUE), SampleRep = SampleRep)
df.5 <- unique(df5)
df.5 <- df.5 %>% mutate(delCT.SRBP = aCT.SRBP - aCT.GAPDH) %>% mutate(delCT.LHX = aCT.LHX - aCT.GAPDH) %>% mutate(delCT.DENND1A = aCT.DENND1A - aCT.GAPDH)
df.5 <- df.5 %>% mutate(raisedFC.SRBP = 2^(-(delCT.SRBP - refCT1))) %>% mutate(raisedFC.LHX = 2^(-(delCT.LHX - refCT2))) %>% mutate(raisedFC.DENND1A = 2^(-(delCT.DENND1A - refCT3)))
df.5 <- df.5 %>% gather("Primer", "FC", 10:12)
CT.df.wide <- df.5  %>% group_by(SampleRep, Primer) %>% summarise(avg.FC = mean((FC), na.rm=TRUE), sd = sd((FC), na.rm=TRUE), sem = (sd/sqrt(length(df.5))))
CT.df.wide.2 <- CT.df.wide %>% mutate(li = avg.FC - sd, ui = avg.FC + sd)

boxplot <- ggplot(data = df.5) +
  geom_boxplot(aes(y = log2(FC), x=SampleRep, colour = SampleRep)) + theme_classic()  + facet_wrap(~Primer) +
  geom_point(aes(y = log2(FC), x=SampleRep)) + xlab("HEK") + ylab("log2 FC") +
  theme_classic() +   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_brewer(palette = "Dark2")

boxplot