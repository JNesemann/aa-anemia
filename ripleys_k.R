#### Ripley's K analysis ####

library(spatstat)
library(maptools)

# using database generated in Anemia_data_prep.R but filtering out missing GPS values
data.hhs <- data.hhs %>% 
  filter(!is.na(lat) & !is.na(lon))

# when using spatstat we need to create a planar point pattern object (ppp)

# first we need an observation window corresponding to our study area

# importing the admin boundaries from https://data.humdata.org/dataset/limites-de-peru
setwd("~/Desktop/data/aa-anemia")
aa <- st_read("per_adm_ign_20200714_shp/per_admbnda_adm2_ign_20200714.shp") %>%
  filter(ADM2_ES == "Alto Amazonas") %>% st_as_sf(crs = crs) %>%
  # adding a buffer to fit hhs that lie slightly outside of the limit
  st_buffer(., 0.8)
ggplot(aa) + geom_sf()
st_is_longlat(aa) # TRUE

# creating an observation window
# note that spatstat need projected coordinates (i.e., planar ones)
# thus I cannot use longitude and latitude
# EPSG 3857 is the web mercator projection - consider changing to one local to peru (5387)
# trying another way
aa_crs <- st_transform(aa, 3857)
st_is_longlat(aa_crs) # FALSE
aa_crs_sp <- as(aa_crs, "Spatial")
aa_win <- as.owin(aa_crs_sp)
aa_win
# success! visualizing it
plot(aa_win)

# converting hhs data to new crs
hhs.merc <- data.hhs %>% 
  # first converting into an sf object
  st_as_sf(coords = c("lon", "lat"), crs = crs) %>%
  # and now transforming to a mercator projection
  st_transform(crs = 3857)
st_is_longlat(hhs.merc) # FALSE

# extracting the coordinates for each point from spatial object 
# and add them to original dataset
data.hhs[c("merc_X","merc_Y")] <- st_coordinates(hhs.merc)

# creating a (marked) planar point pattern (ppp) object
mild <- ppp(x = data.hhs$merc_X, y = data.hhs$merc_Y,
            window = aa_win,
            marks = as.factor(data.hhs$atleast1_mild.f))
mild.km <- rescale(mild, 1000, "km")

# visualizing
plot(mild.km, cols = c("red","navy"))
summary(mild.km) # this calculates the spatial intensity of cases and controls

# first step is to split the households into cases and controls
mildcases <- split(mild.km)$case
mildcases
mildcont <- split(mild.km)$control
mildcont 

# determining diameters of the villages
# per https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4580214/
# the max scan distance should be half the diameter of the village
# ts <- geosphere::distm(cbind(filter(data.hhs, community=="Vista Allegre")$lat, 
#                              filter(data.hhs, community=="Vista Allegre")$lon))
# max(ts/1000/2) # distance in meters so /1000 for Km and /2 for half the diameter
# # San antonio 0.50
# # 06 de Julio 0.30
# # Angamos 0.05
# # Bellavista (Balsapuerto) 0.28
# # Bellavista (Jeberos) 0.96
# # Bethel 0.45
# # Centro America 0.30
# # Corazon de Jesus 0.03
# # Dos de Mayo 0.37
# # Huancayo 0.16
# # Huatapi 0.17
# # Loma Linda 0.31
# # Maranatha 0.45
# # Nuevo Arica 0.19
# # Nuevo Barranquita 0.22
# # Nuevo Iquitos 0.08
# # Nuevo Papaplaya 0.23
# # Panam 0.86
# # San Antonio 0.50
# # Tamarate 0.36
# # Union Campesino 0.12
# # Vista Allegre 0.23
rmax <- ((0.5+0.3+0.28+0.96+0.45+0.3+0.37+0.16+0.17+0.31+0.45+0.19+0.22+0.23+0.86+0.5+0.36+0.12+0.23)/22)
# # 0.32 km or 320 meters

# compute the k-function for cases
Kmildcases <- Kest(mildcases, 
                   correction = "isotropic",
                   # r = rsplit,
                   rmax = rmax) # this is the maximum distance for the scan

# plot results
plot(Kmildcases)
# black line lies above the red line, suggesting evidence of clustering

# generating confidence bands for the H0 of CRS
plot(envelope(mildcases,
              nsim = 99)) # change to 9999 prn

# now testing to see if clustering exists after controling for distribution of controls
# this is done by subtracting the Kcontronl from Kcase (Kcase - Kcontrol)
# if this is close to 0 then the observed clustering is due to the underlying pop @ risk
# Kcontrol function
Kmildcont <- Kest(mildcont, 
                  correction = "isotropic",
                  # r = rsplit,
                  rmax = rmax)
plot(Kmildcont)

# computing difference 
Dmild <- Kmildcases$iso - Kmildcont$iso
rmild <- Kmildcases$r
plot(rmild, Dmild, type = "l") # not very pretty

# the confidence bands can be determined by randomly relabeling cases/controls
source(here("envelopeD_rmax.R"))
envelopeD_rmax(cases = mildcases, contr = mildcont, 
                       nsims = 99, # can increase this PRN
                       # rsplit = rsplit,
                       Rmax = rmax) 

#### repeating for moderate anemia ####
# generating ppp
mod <- ppp(x = data.hhs$merc_X, y = data.hhs$merc_Y,
           window = aa_win,
           marks = as.factor(data.hhs$atleast1_mod.f))
# rescaling
mod.km <- rescale(mod, 1000, "km")

# splitting
modcases <- split(mod.km)$case
modcases
modcont <- split(mod.km)$control
modcont

# computing k for cases and controls
Kmodcases <- Kest(modcases, 
                  correction = "isotropic",
                  # r = rsplit,
                  rmax = rmax)
plot(Kmodcases)
Kmodcont <- Kest(modcont, 
                 correction = "isotropic",
                 # r = rsplit,
                 rmax = rmax)
plot(Kmodcont)
plot(Kmodcases, add = TRUE)

# computing difference 
Dmod <- Kmodcases$iso - Kmodcont$iso
rmod <- Kmodcases$r
plot(rmod, Dmod, type = "l") # not very pretty

# repeating it with confidence intervals
envelopeD_rmax(cases = modcases, contr = modcont,
                      nsims = 99, # can increase this PRN
                      # rsplit = rsplit,
                      Rmax = rmax)

#### plotting with discrete comparison distances ####
# essentially this will just smooth out the graphs described above

# defining distances at which to calculate K
rsplit <- seq(0, rmax, by = 0.05)

#### mild-or-worse plot ####

# defining all my objects used in the function
win <- mildcases$window
win

x <- c(mildcases$x, mildcont$x)
x

y <- c(mildcases$y, mildcont$y)
y

cc <- c(rep("case", mildcases$n), rep("control", mildcont$n))
cc

# monte carlo function
mc <- function(nsims, Rmax) {
  Dsim <- c() # we'll store our result in this object
  for (i in 1:nsims) {
    cc <- sample(cc)
    simcases <- ppp(x = x[cc == "case"], 
                    y = y[cc =="case"], 
                    window = win)
    simcontr <- ppp(x = x[cc == "control"], 
                    y = y[cc == "control"], 
                    window = win)
    
    Kcases <- Kest(simcases, 
                   correction = "isotropic", 
                   r = rsplit,
                   rmax = Rmax)$iso 
    Kcontrols <- Kest(simcontr, 
                      correction = "isotropic", 
                      r = rsplit,
                      rmax = Rmax)$iso
    dsim <- Kcases - Kcontrols 
    Dsim <- cbind(Dsim, dsim)
  }
  Dsim
}

# testing it out  
Dsim <- mc(nsims = 999, Rmax = rmax)
Dsim  

# getting centiles 
qts <- apply(Dsim, 1, quantile, probs = c(0.025, 0.975))
qts

Kcases <- Kest(mildcases, 
               correction = "isotropic", 
               r = rsplit,
               rmax = rmax)

Kcontrols <- Kest(mildcont, 
                  correction = "isotropic", 
                  r = rsplit,
                  rmax = rmax)

D <- Kcases$iso - Kcontrols$iso
D

r <- Kcases$r
r

# creating a dataframe
mild.plot <- rbind(r, D, qts) %>% as.data.frame() %>%
  rownames_to_column() %>%
  pivot_longer(., V1:V7, 
               names_to = "name",
               values_to = "value") %>%
  pivot_wider(names_from = rowname,
              values_from = value) %>%
  rename(low=`2.5%`, high=`97.5%`) %>%
  select(-name) %>%
  mutate(abline.y = 0) %>%
  # plotting
  ggplot(data = ., aes(x=r, y=D)) + 
  geom_point(shape = 1) +
  geom_line() +
  # monte carlo limits
  geom_line(aes(x=r, y=high), linetype = 2) + 
  geom_line(aes(x=r, y=low), linetype = 2) +
  # reference line of 0
  geom_line(aes(x=r, y=abline.y), color = "red", linetype = 3) + 
  theme_bw() + 
  labs(y = "Difference in K-functions (case - control)", x = "Distance (km)",
       title = "A) Mild-or-worse anemia")
mild.plot

# saving
ggsave(here("figures", "kplots", "kmild.eps"), mild.plot,
       width = 6, height = 3)


#### moderate or worse plot ####

# redefining my objects
win <- modcases$window
win

x <- c(modcases$x, modcont$x)
x

y <- c(modcases$y, modcont$y)
y

cc <- c(rep("case", modcases$n), rep("control", modcont$n))
cc

# monte carlo simulation 
Dsim <- mc(nsims = 999, Rmax = rmax)
Dsim  

# getting centiles 
qts <- apply(Dsim, 1, quantile, probs = c(0.025, 0.975))
qts

Kcases <- Kest(modcases, 
               correction = "isotropic", 
               r = rsplit,
               rmax = rmax)

Kcontrols <- Kest(modcont, 
                  correction = "isotropic", 
                  r = rsplit,
                  rmax = rmax)

D <- Kcases$iso - Kcontrols$iso
D

r <- Kcases$r
r

# creating a dataframe
mod.plot <- rbind(r, D, qts) %>% as.data.frame() %>%
  rownames_to_column() %>%
  pivot_longer(., V1:V7, 
               names_to = "name",
               values_to = "value") %>%
  pivot_wider(names_from = rowname,
              values_from = value) %>%
  rename(low=`2.5%`, high=`97.5%`) %>%
  select(-name) %>%
  mutate(abline.y = 0) %>%
  # plotting
  ggplot(data = ., aes(x=r, y=D)) + 
  geom_point(shape = 1) +
  geom_line() +
  # monte carlo limits
  geom_line(aes(x=r, y=high), linetype = 2) + 
  geom_line(aes(x=r, y=low), linetype = 2) +
  # reference line of 0
  geom_line(aes(x=r, y=abline.y), color = "red", linetype = 3) + 
  theme_bw() + 
  labs(y = "Difference in K-functions (case - control)", x = "Distance (km)",
       title = "B) Moderate-or-worse anemia")
mod.plot

# saving
ggsave(here("figures", "kplots", "kmod.eps"), mod.plot,
       width = 6, height = 3)






# using R script designed by LSHTM with some modifications
source(here("envelopeD_JMN.R")) 
# source(here("envelopeD_rmax.R"))

# mild or worse anemia
envelopeD_jmn(cases = mildcases, contr = mildcont,
                      nsims = 99,
                      rsplit = rsplit,
                      Rmax = rmax)
# # comparing
# envelopeD_rmax(cases = mildcases, contr = mildcont, 
#                        nsims = 99,
#                        Rmax = rmax) 
# # great, these look essentially the same

# moderate or worse
envelopeD_jmn(cases = modcases, contr = modcont,
              nsims = 99,
              rsplit = rsplit,
              Rmax = rmax)
# # comparing
# envelopeD_rmax(cases = modcases, contr = modcont,
#                       nsims = 99,
#                       Rmax = 320) 
# # these also look essentially the same

# creating plots
envelopeD_jmn(cases = mildcases, contr = mildcont,
              nsims = 999,
              rsplit = rsplit,
              Rmax = 320)
envelopeD_jmn(cases = modcases, contr = modcont,
              nsims = 999,
              rsplit = rsplit,
              Rmax = 320)

# to do:
# fix axes 
# relabel axes
# change units to Km so the estimated D(r) is less extreme
# confirm interpretation of line below the confidence envelope


#### troubleshooting creating plots































