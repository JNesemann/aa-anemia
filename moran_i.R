# Moran's I Analysis

# use data created in Anemia_data_prep.R

#### Morans I at the village level ####
# all code below is based off that provided by JK

# converting my data to SF
data_spatial <- data %>%
  # filter(community != "Nuevo Barranquita" & community != "Corazon de Jesus") %>%
  st_as_sf(coords = c("lon","lat"), crs = crs) %>%
  as("Spatial")

# create symmetric matrix of distances between every point
dist_matrix <- sapply(1:nrow(data_spatial),
                      # using the distGeo method accounts for globe curvature when calculating distance
                      function(x) geosphere::distGeo(p1 = data_spatial, p2 = data_spatial[x,]))

# take inverse of distance 
dist_matrix_inv <- 1/dist_matrix
# replace the Inf values in the diagonals with 0
diag(dist_matrix_inv) <- 0

# estimate moran's I and statistical significance for each level of anemia
moran_pmild <- ape::Moran.I(data$p_any, dist_matrix_inv) %>% as_tibble(.) %>%
  mutate(measure = "Mild-or-worse prevalence")
moran_pmild.logit <- ape::Moran.I(data$any.logit, dist_matrix_inv) %>% as_tibble(.) %>%
  mutate(measure = "Mild-or-worse logit prevalence")
moran_pmod <- ape::Moran.I(data$p_modsev, dist_matrix_inv) %>% as_tibble(.) %>%
  mutate(measure = "Moderate-or-worse prevalence")
moran_pmod.logit <- ape::Moran.I(data$modsev.logit, dist_matrix_inv) %>% as_tibble(.) %>%
  mutate(measure = "Moderate-or-worse logit prevalence")
moran_meanhg <- ape::Moran.I(data$mean_hg, dist_matrix_inv) %>% as_tibble(.) %>%
  mutate(measure = "Mean Hg")
moran_medianhg <- ape::Moran.I(data$median_hg, dist_matrix_inv) %>% as_tibble(.) %>%
  mutate(measure = "Median Hg")

moran.village <- rbind(moran_pmild, moran_pmild.logit, moran_pmod, moran_pmod.logit, moran_meanhg, moran_medianhg) %>%
  select(measure, everything())
moran.village

write_csv(moran.village, here("tables","village_morans.csv"))

#### Moran's I at household level ####
# not sure if this will work as the prevalence estimates in the houses will be crude (ie 0, 0.5, or 1 given the low n in each household)

# converting my data to SF
data_spatial.hhs <- data.hhs %>% 
  # filter(community != "Nuevo Barranquita" & community != "Corazon de Jesus") %>%
  # just applying this to one village
  #filter(community=="Panam") %>%
  st_as_sf(coords = c("lon","lat"), crs = crs) %>%
  as("Spatial")

# create symmetric matrix of distances between every point
dist_matrix.hhs <- sapply(1:nrow(data_spatial.hhs),
                          # using the distGeo method accounts for globe curvature when calculating distance
                          function(x) geosphere::distGeo(p1 = data_spatial.hhs, p2 = data_spatial.hhs[x,]))

# take inverse of distance 
dist_matrix_inv.hhs <- 1/dist_matrix.hhs # due to duplicates with 0 distance between them this is leading to 1/0 = Inf
# replace the Inf values in the diagonals with 0
diag(dist_matrix_inv.hhs) <- 0

# estimate moran's I and statistical significance for each level of anemia using hh level
# hhmoran_any <- ape::Moran.I(data.hhs$p_any, dist_matrix_inv.hhs) %>% as_tibble(.) %>%
#   mutate(measure = "Mild-or-worse prevalence")
# hhmoran_anylogit <- ape::Moran.I(data.hhs$any.logit, dist_matrix_inv.hhs) %>% as_tibble(.) %>%
#   mutate(measure = "Mild-or-worse logit prevalence")
# hhmoran_mod <- ape::Moran.I(data.hhs$p_modsev, dist_matrix_inv.hhs) %>% as_tibble(.) %>%
#   mutate(measure = "Mod-or-worse prevalence")
# hhmoran_mod.logit <- ape::Moran.I(data.hhs$modsev.logit, dist_matrix_inv.hhs) %>% as_tibble(.) %>%
#   mutate(measure = "Mod-or-worse prevalence")
hhmoran_meanhg <- ape::Moran.I(data.hhs$mean_hg, dist_matrix_inv.hhs) %>% as_tibble(.) %>%
  mutate(measure = "Mean Hg")
hhmoran_medianhg <- ape::Moran.I(data.hhs$median_hg, dist_matrix_inv.hhs) %>% as_tibble(.) %>%
  mutate(measure = "Median Hg")
hhmoran_randhg <- ape::Moran.I(data.hhs$rand_hg, dist_matrix_inv.hhs) %>% as_tibble(.) %>%
  mutate(measure = "Random Hg")

moran.hh.overall <- rbind(hhmoran_meanhg, hhmoran_medianhg, hhmoran_randhg) %>% 
  select(measure, everything())

# hhmoran_any, hhmoran_anylogit, hhmoran_mod, hhmoran_mod.logit, 

moran.hh.overall
ggplot(data = data.hhs, aes(x=p_any)) + geom_histogram() #+ facet_grid(.~community)
ggplot(data = data.hhs, aes(x=any.logit)) + geom_histogram()
# the above plot shows that the prevalence distributions are not normal so using  
# prevalence outcomes at the household level are not valid

write_csv(moran.hh.overall, here("tables", "moran_hh_overall.csv"))


#### Morans I for all households grouped by village ####

data.hhs %>% group_by(community) %>%
  summarise(n=sum(!is.na(hh_numb))) %>% view(.)

# creating a function to apply the Moran's I to the households within each village
village.moran <- function(outcome_var, village) {
  
  # 1. create a df for just the village
  temp_df <- data.hhs %>%
    filter(community == village)
  
  # 2. convert to spatial format
  temp_df_spatial <- temp_df %>% 
    st_as_sf(coords = c("lon","lat"), crs = crs) %>%
    as("Spatial")
  
  # 2. create symmetric matrix of distances between every point
  dist_matrix <- sapply(1:nrow(temp_df_spatial),
                        function(x) geosphere::distGeo(p1 = temp_df_spatial, p2 = temp_df_spatial[x,]))
  
  # 3. take inverse of distance
  dist_matrix_inv <- 1/dist_matrix
  
  # 4. replace the Inf values with 0 
  diag(dist_matrix_inv) <- 0
  
  # 5. estimate Moran's I and p-value for each level of anemia
  ret <- ape::Moran.I(temp_df %>% pull(outcome_var), dist_matrix_inv) %>% 
    as_tibble() %>% 
    mutate(outcome_var = outcome_var, village = village) %>%
    select(village, outcome_var, observed, expected, p.value)
  
  # 6. bind all of these together
  # ret <- rbind(any, mild, mod) %>% select(village, level, observed, expected, p.value)
  
  # 6. return the results
  return(ret)
  
}

# now to iterate it over every village
nuevo_arica <- mapply(village.moran,
                      outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                      village = "Nuevo Arica",
                      SIMPLIFY = F) %>%
  bind_rows()

dos_mayo <- mapply(village.moran,
                   outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                   village = "Dos de Mayo",
                   SIMPLIFY = F) %>%
  bind_rows()

bellavista_b <- mapply(village.moran,
                       outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                       village = "Bellavista (Balsapuerto)",
                       SIMPLIFY = F) %>%
  bind_rows()

nuevo_papaplaya <- mapply(village.moran,
                          outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                          village = "Nuevo Papaplaya",
                          SIMPLIFY = F) %>%
  bind_rows()

bethel <- mapply(village.moran,
                 outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                 village = "Bethel",
                 SIMPLIFY = F) %>%
  bind_rows()

bellavista_j <- mapply(village.moran,
                       outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                       village = "Bellavista (Jeberos)",
                       SIMPLIFY = F) %>%
  bind_rows()

vista_allegre <- mapply(village.moran,
                        outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                        village = "Vista Allegre",
                        SIMPLIFY = F) %>%
  bind_rows()

huancayo <- mapply(village.moran,
                   outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                   village = "Huancayo",
                   SIMPLIFY = F) %>%
  bind_rows()

nuevo_arica <- mapply(village.moran,
                      outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                      village = "Nuevo Arica",
                      SIMPLIFY = F) %>%
  bind_rows()

sies_julio <- mapply(village.moran,
                     outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                     village = "06 de Julio",
                     SIMPLIFY = F) %>%
  bind_rows()

tamarate <- mapply(village.moran,
                   outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                   village = "Tamarate",
                   SIMPLIFY = F) %>%
  bind_rows()

huatapi <- mapply(village.moran,
                  outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                  village = "Huatapi",
                  SIMPLIFY = F) %>%
  bind_rows()

union_campesino <- mapply(village.moran,
                          outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                          village = "Union Campesino",
                          SIMPLIFY = F) %>%
  bind_rows()

nuevo_iquitos <- mapply(village.moran,
                        outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                        village = "Nuevo Iquitos",
                        SIMPLIFY = F) %>%
  bind_rows()

angamos <- mapply(village.moran,
                  outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                  village = "Angamos",
                  SIMPLIFY = F) %>%
  bind_rows()

san_antonio <- mapply(village.moran,
                      outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                      village = "San Antonio",
                      SIMPLIFY = F) %>%
  bind_rows()

panam <- mapply(village.moran,
                outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                village = "Panam",
                SIMPLIFY = F) %>%
  bind_rows()

maranatha <- mapply(village.moran,
                    outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                    village = "Maranatha",
                    SIMPLIFY = F) %>%
  bind_rows()

loma_linda <- mapply(village.moran,
                     outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                     village = "Loma Linda",
                     SIMPLIFY = F) %>%
  bind_rows()

centro_america <- mapply(village.moran,
                         outcome_var = c("mean_hg", "median_hg", "rand_hg"), # "p_any", "any.logit", "p_modsev", "modsev.logit", 
                         village = "Centro America",
                         SIMPLIFY = F) %>%
  bind_rows()

village.morans <- rbind(centro_america, loma_linda, maranatha, panam, san_antonio, angamos, 
                        nuevo_iquitos, union_campesino, huatapi, tamarate, sies_julio, nuevo_arica, 
                        huancayo, vista_allegre, bellavista_j, bethel, nuevo_papaplaya, bellavista_b, 
                        dos_mayo, nuevo_arica) %>%
  left_join((ts %>%
               group_by(community) %>%
               summarise(n_hh=length(unique(hh_numb))) %>% rename(village=community))) %>%
  mutate(sig.p=if_else(p.value>0.05,0,1))

village.morans

write_csv(village.morans, here("tables", "morans_per_village.csv"))

# ggplot(village.morans, aes(x=expected, y=observed, color = sig.p)) + 
#   geom_point() +
#   facet_grid(~outcome_var)

pervillage <- village.morans %>%
  group_by(outcome_var) %>%
  summarise(mean=mean(observed),
            sd=sd(observed),
            median=median(observed),
            p25=quantile(observed, 0.25),
            p75=quantile(observed, 0.75),
            min=min(observed),
            max=max(observed)) %>%
  mutate(p25=round(p25, digits=2),
         p75=round(p75, digits=2)) %>%
  unite(iqr, p25, p75, sep = " to ")
pervillage

write_csv(pervillage, here("tables", "morans_per_village_summary.csv"))

