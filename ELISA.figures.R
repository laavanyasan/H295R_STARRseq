#2024-11
#Revision ELISAs
#Boxplots
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
library(NatParksPalettes)
library(ggsignif)
library(readxl)

#Testosterone - H295R-dCas9-p300
elisa.data <- read_excel("p300-T.xls",sheet="forR",col_names = TRUE)
elisa.data <- elisa.data %>% filter(Sample !="NONE")


boxplot <- ggplot(data = (elisa.data)) +
  geom_boxplot(aes(y = (TestosteroneConc/1000), x=Sample, colour = Drug)) + theme_classic() + facet_wrap(~Day) +
  geom_point(aes(y = (TestosteroneConc/1000), x=Sample, colour = Drug, shape = Day), position = position_dodge(width = 0.8)) +
  ylab("Testosterone concentration (ng/ml)") 

boxplot

##Change in T between day 2 and day 4 H295R-dCas9-p300
elisa.df.wide <- elisa.data2 %>% group_by(Sample, Day, Drug) %>% summarise(avg.Test = mean(TestosteroneConc), sd = sd(TestosteroneConc), sem = (sd/sqrt(length(elisa.data))))
elisa.df.wide <- elisa.df.wide %>% mutate(li = avg.Test - sem, ui = avg.Test + sem)

p.ELISA.3 <- ggplot(data = elisa.df.wide, aes(x = Day, y = avg.Test)) +
  geom_point(aes(x = Day, y = avg.Test, colour = Sample), size = 4) +
  geom_line(aes(group = Sample, colour = Sample)) +
  xlab("H295R-dcas9p300 with gRNA agianst DENND1A promoter") + ylab("Testosterone conc(pg/ml or ng/l)") +
  theme_classic() +   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_grid(~Drug)

p.ELISA.3

########
#Testosterone - H295R-dCas9-KRAB
elisa.data <- read_excel("KRAB-T.xls",sheet="forR",col_names = TRUE)


boxplot <- ggplot(data = (elisa.data)) +
  geom_boxplot(aes(y = (TestosteroneConc/1000), x=Sample, colour = Drug)) + theme_classic() + facet_wrap(~Day) +
  geom_point(aes(y = (TestosteroneConc/1000), x=Sample, colour = Drug, shape = Day), position = position_dodge(width = 0.8)) +
  ylab("Testosterone concentration (ng/ml)")

boxplot
#######

#Estradiol - H295R-dCas9-p300

elisa.data <- read_excel("p300-E.xls",sheet="forR",col_names = TRUE)


elisa.data$TestosteroneConc <- as.numeric(elisa.data$TestosteroneConc)


boxplot <- ggplot(data = (elisa.data)) +
  geom_boxplot(aes(y = (TestosteroneConc/1000), x=Sample, colour = Drug)) + theme_classic() + facet_wrap(~Day+Drug) +
  geom_point(aes(y = (TestosteroneConc/1000), x=Sample, colour = Drug, shape = Day), position = position_dodge(width = 0.8)) +
  ylab("Estradiol concentration pg/ml")

boxplot


#########

#Estradiol - H295R-dCas9-KRAB

elisa.data <- read_excel("KRAB-E.xls",sheet="forR",col_names = TRUE)


elisa.data$TestosteroneConc <- as.numeric(elisa.data$TestosteroneConc)


boxplot <- ggplot(data = (elisa.data)) +
  geom_boxplot(aes(y = (TestosteroneConc/1000), x=Sample, colour = Drug)) + theme_classic() + facet_wrap(~Day+Drug) +
  geom_point(aes(y = (TestosteroneConc/1000), x=Sample, colour = Drug, shape = Day), position = position_dodge(width = 0.8)) +
  ylab("Estradiol concentration pg/ml")

boxplot

#####
##T test generate
filt.data <- elisa.data %>% filter(Sample %in% c("S5", "S6")) %>% filter(Drug=="Forskolin")
t.test(filt.data$TestosteroneConc ~ filt.data$Sample) 
###


