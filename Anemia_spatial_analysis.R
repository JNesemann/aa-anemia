# ALTO AMAZONAS ANEMIA SPATIAL ANALYSIS
# JOHN NESEMANN

options(scipen = 999)

#### packages ####
library(tidyverse)
library(haven)
library(foreign)

#### data ####
setwd("~/Desktop/DATA/aa-anemia")

# anemia data
data <- read_dta("anemia_alto_amazonas_abril.dta") %>%
  select(unique_id:hh_numb, pplinhh:gps_alt, dbs_hg, dbs_hg_na, anemia, anemia_level_sa, anemia_orig) %>%
  rename(anemia.f=anemia_orig)
data

# population data from 2017 alto amazonas health network census
# this will be used to calculate age and sex adjusted prevalence estimates 
AAcensus <- read_csv("2018AA Census.csv") %>%
  # filtering out yurimaguas and total since we did not include yurimaguas in our sampling frame
  filter(MICRORED!="YURIMAGUAS" & MICRORED!="TOTAL") %>%
  # selecting relevant columns 
  select(MICRORED, `1a`,`2a`,`3a`,`4a`,`5a`,`6a`,`7a`,`8a`,`9a`) %>% select(-MICRORED) %>%
  mutate(`1`=sum(`1a`),
         `2`=sum(`2a`),
         `3`=sum(`3a`),
         `4`=sum(`4a`),
         `5`=sum(`5a`),
         `6`=sum(`6a`),
         `7`=sum(`7a`),
         `8`=sum(`8a`),
         `9`=sum(`9a`),
         total=sum(`1a`)+sum(`2a`)+sum(`3a`)+sum(`4a`)+sum(`5a`)+sum(`6a`)+sum(`7a`)+sum(`8a`)+sum(`9a`)) %>% select(-(`1a`:`9a`)) %>% unique() %>%
  pivot_longer(cols = `1`:total, names_to ="age.group19") %>%
  mutate(weight=value/11599)
# creating the weights
weights19 <- AAcensus %>% filter(age.group19!="total") %>% mutate(age.group19=as.numeric(age.group19))

# 










