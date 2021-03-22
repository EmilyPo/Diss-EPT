# script to generate visualization of ct and ept trends from SKCPH quarterly reports 
# 2012 Q1 to 2018 Q1 
# EPT percentatges based on graphs in report - not raw data 

library(readxl)
library(ggplot2)
library(tidyverse)

dat <- read_xlsx("~/Documents/Dissertation/Data/chlamydia_stats/chlamydia_ept_trends_skphc.xlsx")
colnames(dat)[2:5] <- c("Women", "Men", "EPT Offered", "Treated >=1 Partner")
dat$time <- 1:nrow(dat)

diag <- dat %>% 
  select(-c("EPT Offered", "Treated >=1 Partner")) %>% 
  gather(key="Sex", value = "Diagnoses", -c(`Year/Quarter`, time))

ept <- dat %>% 
  select(-c(Women, Men)) %>% 
  gather(key="EPT Provision/Usage", value = "Percent", -c(`Year/Quarter`, time))

# diagnoses
diag %>% ggplot(aes(x = time, y = Diagnoses, color=Sex)) +
  geom_point() +
  geom_smooth(method="loess", formula=y~x) +
  theme(axis.text.x = element_text(angle=45)) +
  scale_x_continuous(breaks=1:25, labels = dat$`Year/Quarter`) + 
  scale_y_continuous(name="Reported Chlamydia Diagnoses", limits = c(0,1500))

ept %>% ggplot(aes(x=time, y=Percent, color=`EPT Provision/Usage`)) +
          geom_point() +
          geom_smooth(method="loess", formula=y~x) +
          theme(axis.text.x = element_text(angle=45)) +
          scale_x_continuous(breaks=1:25, labels = dat$`Year/Quarter`) + 
          scale_y_continuous(name="Percent of Patients", limits=c(0,100))


