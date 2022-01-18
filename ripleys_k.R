#### Ripley's K analysis ####

library(spatstat)
library(maptools)

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
# not that spatstat need projected coordinates (i.e., planar ones)
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
glimpse(data.hhs)

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
ts <- geosphere::distm(cbind(filter(data.hhs, community=="Vista Allegre")$lat, 
                             filter(data.hhs, community=="Vista Allegre")$lon))
max(ts/1000/2) # distance in meters so /1000 for Km and /2 for half the diameter
# San antonio 0.50
# 06 de Julio 0.30
# Angamos 0.05
# Bellavista (Balsapuerto) 0.28
# Bellavista (Jeberos) 0.96
# Bethel 0.45
# Centro America 0.30
# Corazon de Jesus 0.03
# Dos de Mayo 0.37
# Huancayo 0.16
# Huatapi 0.17
# Loma Linda 0.31
# Maranatha 0.45
# Nuevo Arica 0.19
# Nuevo Barranquita 0.22
# Nuevo Iquitos 0.08
# Nuevo Papaplaya 0.23
# Panam 0.86
# San Antonio 0.50
# Tamarate 0.36
# Union Campesino 0.12
# Vista Allegre 0.23
(0.5+0.3+0.28+0.96+0.45+0.3+0.37+0.16+0.17+0.31+0.45+0.19+0.22+0.23+0.86+0.5+0.36+0.12+0.23)/22
# 0.32

# defining distances at which to calculate K
rsplit <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.32)

# compute the k-function for cases
Kmildcases <- Kest(mildcases, 
                   correction = "isotropic",
                   r = rsplit,
                   rmax = 0.32) # this is the maximum distance for the scan

# plot results
plot(Kmildcases)
# black line lies above the red line, suggesting evidence of clustering

# # generating confidence bands for the H0 of CRS
# plot(envelope(mildcases, 
#               nsim = 999)) # change to 9999 prn

# now testing to see if clustering exists after controling for distribution of controls
# this is done by subtracting the Kcontronl from Kcase (Kcase - Kcontrol)
# if this is close to 0 then the observed clustering is due to the underlying pop @ risk
# Kcontrol function
Kmildcont <- Kest(mildcont, 
                  correction = "isotropic",
                  r = rsplit,
                  rmax = 0.32)
plot(Kmildcont)

# computing difference 
Dmild <- Kmildcases$iso - Kmildcont$iso
rmild <- Kmildcases$r
plot(rmild, Dmild, type = "l") # not very pretty
# best so far: 0.96 or 0.2 or 0.8

# the confidence bands can be determined by randomly relabeling cases/controls
source(here("envelopeD_JMN.R")) # using R script designed by LSHTM with some modifications
mildKplot <- envelopeD(cases = mildcases, contr = mildcont, 
                       nsims = 999, # can increase this PRN
                       rsplit = rsplit,
                       Rmax = 0.32) 
mildKplot

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
                  r = rsplit,
                  rmax = 0.5)
Kmodcont <- Kest(modcont, 
                 correction = "isotropic",
                 r = rsplit,
                 rmax = 0.5)
# computing difference 
Dmod <- Kmodcases$iso - Kmodcont$iso
rmod <- Kmodcases$r
plot(rmod, Dmod, type = "l") # not very pretty

# repeating it with confidence intervals
modKplot <- envelopeD(cases = modcases, contr = modcont,
                      nsims = 999, # can increase this PRN
                      rsplit = rsplit,
                      Rmax = 0.32)
modKplot
