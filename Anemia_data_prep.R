# ALTO AMAZONAS ANEMIA SPATIAL ANALYSIS
# JOHN NESEMANN

#### packages ####
library(tidyverse)
library(foreign)
library(sf)
library(gridExtra)
library(here)

#### configurations and objects of constant value ####
options(scipen = 999)
crs <- 4326 # this is the coordinate reference system for the project, 
# these numbers correspond to WGS84

#### data ####
# creating a uncleaned/untransformed data object
ts <- haven::read_dta(here("data","anemia_alto_amazonas_abril.dta")) %>%
  group_by(hh_numb) %>%
  fill(gps_lat, .direction = c("updown")) %>%
  fill(gps_long, .direction = c("updown"))

# making sure variables were correctly coded
xtabs(data=ts, ~dbs_hg_na+anemia_orig, addNA=T)
xtabs(data=ts, ~dbs_hg_na+anemia_level_sa, addNA=T)
# looks about right, some overlap due to age specific anemia cut offs

# anemia data
# village level
data <- haven::read_dta(here("data","anemia_alto_amazonas_abril.dta")) %>%
  select(unique_id:hh_numb, pplinhh:gps_alt, dbs_hg, dbs_hg_na, anemia, anemia_level_sa, 
         anemia_orig) %>%
  dplyr::rename(anemia.f=anemia_orig) %>%
  # the code JK send requires that the data be aggregated at a certain level, 
  # will do village level here
  group_by(community) %>%
  summarise(n_tested=sum(!is.na(age)),
            n_any_anemia=sum(anemia.f %in% c("mild","moderate","severe")),
            # n_mild=sum(anemia.f=="mild"),
            # grouping moderate and severe together due to the low number of severe cases
            n_modsev=sum(anemia.f %in% c("moderate","severe")),
            p_any=n_any_anemia/n_tested,
            # p_mild=n_mild/n_tested,
            p_modsev=n_modsev/n_tested,
            mean_hg=mean(dbs_hg_na, na.rm=T),
            median_hg=median(dbs_hg_na, na.rm=T),
            lat=median(gps_lat, na.rm = T),
            lon=median(gps_long, na.rm = T)) %>%
  # ggplot(data=data,aes(x=p_any)) + geom_histogram() # does not appear normal
  # ggplot(data=data,aes(x=p_modsev)) + geom_histogram() # does not appear normal
  # transforming to the empirical logit scale
  mutate(any.logit=log((n_any_anemia + 0.5) / (n_tested - n_any_anemia + 0.5)),
         modsev.logit=log((n_modsev + 0.5) / (n_tested - n_modsev + 0.5)))

# checking normality of transformed data
# ggplot(data, aes(x=p_any)) + geom_histogram(bins = 10)
# ggplot(data, aes(x=any.logit)) + geom_histogram(bins = 10) # maybe it is more normal??
# ggplot(data, aes(x=modsev.logit)) + geom_histogram(bins = 10)
# ggplot(data, aes(x=p_modsev)) + geom_histogram(bins = 10) 
# this is definitely more normal

# household level
data.hhs <- haven::read_dta(here("data","anemia_alto_amazonas_abril.dta")) %>%
  select(unique_id:hh_numb, pplinhh:gps_alt, dbs_hg, dbs_hg_na, anemia, 
         anemia_level_sa, anemia_orig, unique_id, age, sex) %>%
  select(-age_cat, -anemia, -anemia_level_sa) %>%
  dplyr::rename(anemia.f=anemia_orig) %>% ungroup() %>%
  # grouping the few households with the same GPS coordinates together
  # mutate(hh_numb=case_when(hh_numb %in% c("H-2512384","H-2512385","H-2512387")~"H-2512384",
  #                          TRUE ~ hh_numb)) %>%
  group_by(hh_numb, community) %>%
  fill(gps_lat, .direction = c("updown")) %>%
  fill(gps_long, .direction = c("updown"))
data.hhs

# setting seed
set.seed(1072021)

# selecting a random hg value per household
data.hhs <- data.hhs %>% 
  group_by(hh_numb, community) %>%
  sample_n(., 1) %>% 
  # creating a random hg variable
  mutate(rand_hg=dbs_hg_na) %>%
  full_join(., data.hhs) %>%
  mutate(atleast1_mild=if_else(anemia.f == "normal", 0, 1), 
         # xtabs(data=data.hhs, ~atleast1_mild+anemia.f, addNA=T)
         atleast1_mod=if_else(anemia.f %in% c("moderate","severe"), 1, 0)) %>% 
  # xtabs(data=data.hhs, ~atleast1_mod+anemia.f, addNA=T)
  summarise(# village=community,
    n_ppl=sum(!is.na(unique_id)),
    mean_age=mean(age, na.rm=T),
    n_female=sum(sex==0),
    n_tested=sum(!is.na(age)),
    n_any_anemia=sum(anemia.f %in% c("mild","moderate","severe")),
    n_mild=sum(anemia.f=="mild"),
    # grouping moderate and severe together due to the low number of severe cases
    n_modsev=sum(anemia.f %in% c("moderate","severe")),
    p_any=n_any_anemia/n_tested,
    p_mild=n_mild/n_tested,
    p_modsev=n_modsev/n_tested,
    lat=median(gps_lat, na.rm = T),
    lon=median(gps_long, na.rm = T),
    mean_hg=mean(dbs_hg_na, na.rm=T), 
    median_hg=median(dbs_hg_na, na.rm=T), 
    atleast1_mild=sum(atleast1_mild), # xtabs(data=data.hhs, ~atleast1_mild, addNA=T)
    atleast1_mod=sum(atleast1_mod), # xtabs(data=data.hhs, ~atleast1_mod, addNA=T)
    rand_hg=sum(rand_hg, na.rm=T)) %>% # xtabs(data=data.hhs, ~rand_hg, addNA=T)
  # transforming to the empirical logit scale
  mutate(any.logit=log((n_any_anemia + 0.5) / (n_tested - n_any_anemia + 0.5)),
         modsev.logit=log((n_modsev + 0.5) / (n_tested - n_modsev + 0.5)),
         # scaling atleast1 variables to binary
         atleast1_mild=if_else(atleast1_mild>0,1,0),
         # xtabs(data=data.hhs, ~atleast1_mild, addNA=T)
         atleast1_mild.f=if_else(atleast1_mild==1,"case","control"),
         # xtabs(data=data.hhs, ~atleast1_mild + atleast1_mild.f, addNA=T)
         atleast1_mod=if_else(atleast1_mod>0,1,0),
         # xtabs(data=data.hhs, ~atleast1_mod, addNA=T)
         atleast1_mod.f=if_else(atleast1_mod==1,"case","control")) # %>%
  # xtabs(data=data.hhs, ~atleast1_mod + atleast1_mod.f, addNA=T)
  # filtering out missing gps values
  # filter(!is.na(lat) & !is.na(lon)) %>% ungroup() # %>%
  # filtering out communities with small numbers of observations
  # filter(community != "Nuevo Barranquita" & community != "Corazon de Jesus")

data.hhs

# checking normality
# ggplot(data.hhs, aes(x=p_any)) + geom_histogram(bins = 10)
# ggplot(data.hhs, aes(x=any.logit)) + geom_histogram(bins = 10) # more normal
# ggplot(data.hhs, aes(x=p_modsev)) + geom_histogram(bins = 10)
# ggplot(data.hhs, aes(x=modsev.logit)) + geom_histogram(bins = 10)
# ggplot(data.hhs, aes(x=mean_hg)) + geom_histogram() # looks normal to me
# ggplot(data.hhs, aes(x=median_hg)) + geom_histogram() # looks normal to me
# ggplot(data.hhs, aes(x=rand_hg)) + geom_histogram()

# checking for duplicates
data.hhs %>% group_by(hh_numb) %>% mutate(dups=n()) %>% filter(dups>1) # no duplicates

#### summary stats ####

# summary stats of # of villages
ts %>% 
  ungroup() %>%
  summarise(n=sum(!is.na(unique_id)),
                 n_hh=length(unique(hh_numb)),
                 na_gps=sum(is.na(gps_lat)),
                 v_villages=length(unique(community)))

ts %>% group_by(hh_numb) %>% summarize(na_lat=sum(is.na(gps_lat))) %>% filter(na_lat>0) # 18 households missing gps data
# but totaling 37 individuals

data.hhs %>% 
  ungroup() %>%
  summarise(n_hh=length(unique(hh_numb)),
                   na_gps=sum(is.na(lat)),
            p_nagps=na_gps/n_hh * 100,
                   v_villages=length(unique(community)))

# characteristics of the households missing GPS
data.hhs %>%
  ungroup() %>%
  mutate(na_gps=if_else(is.na(lat) | is.na(lon), 1, 0)) %>%
  group_by(na_gps) %>%
  summarise(n=sum(n_ppl),
            n_hhs=sum(!is.na(hh_numb)),
            mean_n=mean(n_ppl),
            sd_n=sd(n_ppl),
            median_n=median(n_ppl),
            p25=quantile(n_ppl, 0.25),
            p75=quantile(n_ppl, 0.75),
            mean.age=mean(mean_age),
            sd_age=sd(mean_age),
            n_fem=sum(n_female),
            p_fem=n_fem/n*100,
            n_mild=sum(atleast1_mild),
            p_mild=n_mild/n_hhs*100,
            n_mod=sum(atleast1_mod),
            p_mod=n_mod/n_hhs*100)

# do I need a statistical test to compare the two?
# can do chi-square test later if JK requires
ts.xt <- data.hhs %>%
  ungroup() %>%
  mutate(na_gps=if_else(is.na(lat) | is.na(lon), 1, 0))

xtabs(data=ts.xt, ~na_gps+atleast1_mod, addNA=T)
chisq.test(xtabs(data=ts.xt, ~na_gps+atleast1_mod, addNA=T))



# distribution of hg values

# overall
supfig1 <- ts %>%
  mutate(vline=case_when(age_group2 == "1-4" ~ 11,
                         age_group2 == "5-9" ~ 11.5),
         age_group2=case_when(age_group2 == "1-4" ~ "1-4 years old",
                              age_group2 == "5-9" ~ "5-9 years old")) %>% 
  ggplot(data=., aes(dbs_hg_na)) +
  geom_histogram(bins=40) +
  facet_grid(.~age_group2) +
  geom_vline(aes(xintercept = vline, color = "red"), linetype = 2) +
  theme_bw() + 
  labs(x="Hemoglobin", y = "Count") +
  theme(legend.position = "none") 
supfig1

# ggsave(here("figures","supfig1.eps"), plot = supfig1)

# village level
supfig2 <- ts %>%
  mutate(vline=case_when(age_group2 == "1-4" ~ 11,
                         age_group2 == "5-9" ~ 11.5),
         age_group2=case_when(age_group2 == "1-4" ~ "1-4 years old",
                              age_group2 == "5-9" ~ "5-9 years old")) %>% 
  ggplot(data=., aes(dbs_hg_na)) +
  geom_histogram(bins=20) +
  facet_grid(age_group2~community) +
  geom_vline(aes(xintercept = vline, color = "red"), linetype = 2) +
  theme_bw() + 
  labs(x="Hemoglobin", y = "Count") +
  theme(legend.position = "none")

# OR
supfig2 <- ts %>%
  mutate(vline=case_when(age_group2 == "1-4" ~ 11,
                         age_group2 == "5-9" ~ 11.5),
         age_group2=case_when(age_group2 == "1-4" ~ "1-4 years old",
                              age_group2 == "5-9" ~ "5-9 years old")) %>% 
  ggplot(data=., aes(dbs_hg_na)) +
  geom_histogram(bins=20) +
  facet_grid(.~community) +
  # geom_vline(aes(xintercept = vline, color = "red"), linetype = 2) +
  theme_bw() + 
  labs(x="Hemoglobin", y = "Count") +
  theme(legend.position = "none")
supfig2
# ggsave(here("figures","supfig2.eps"), plot = supfig2)

ts %>% group_by(community) %>%
  summarise(n=sum(!is.na(unique_id)),
            n_hg=sum(!is.na(dbs_hg_na))) %>%
  arrange(n_hg)

# remove nuevo barranquito and corazon de jesus due to sparse data
# data <- data %>%
#   filter(community != "Nuevo Barranquita" & community != "Corazon de Jesus")
# 
# data.hhs <- data.hhs %>%
#   filter(community != "Nuevo Barranquita" & community != "Corazon de Jesus")
