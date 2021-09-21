# ALTO AMAZONAS ANEMIA SPATIAL ANALYSIS
# JOHN NESEMANN

#### packages ####
library(tidyverse)
library(haven)
library(foreign)
library(sf)
library(gridExtra)
library(here)

#### configurations and obects of constant value ####
options(scipen = 999)
crs <- 4326 # this is the coordinate reference system for the project, these numbers correspond to WGS84

#### data ####
setwd("~/Desktop/DATA/aa-anemia")

# anemia data
# village level
data <- read_dta("anemia_alto_amazonas_abril.dta") %>%
  select(unique_id:hh_numb, pplinhh:gps_alt, dbs_hg, dbs_hg_na, anemia, anemia_level_sa, anemia_orig) %>%
  dplyr::rename(anemia.f=anemia_orig) %>%
  # the code JK send requires that the data be aggregated at a certain level, will do village level here
  # might be able to do household level later to evaluate spatial distribution within each community
  group_by(community) %>%
  summarise(n_tested=sum(!is.na(age)),
            n_any_anemia=sum(anemia.f %in% c("mild","moderate","severe")),
            n_mild=sum(anemia.f=="mild"),
            # grouping moderate and severe together due to the low number of severe cases
            n_modsev=sum(anemia.f %in% c("moderate","severe")),
            p_any=n_any_anemia/n_tested,
            p_mild=n_mild/n_tested,
            p_modsev=n_modsev/n_tested,
            # going for medians to get village midpoint as often there were a few houses located very far away
            lat=median(gps_lat, na.rm = T),
            lon=median(gps_long, na.rm = T))
data

# household level
# not sure if this will work statistically given each household has a very small N but will give it a shot later if I have the time
data.hhs <- read_dta("anemia_alto_amazonas_abril.dta") %>%
  select(unique_id:hh_numb, pplinhh:gps_alt, dbs_hg, dbs_hg_na, anemia, anemia_level_sa, anemia_orig, unique_id) %>%
  select(-age_cat, -anemia, -anemia_level_sa) %>%
  dplyr::rename(anemia.f=anemia_orig) %>% ungroup() %>%
  # grouping the few households with the same GPS coordinates together
  mutate(hh_numb=case_when(hh_numb %in% c("H-2512384","H-2512385","H-2512387")~"H-2512384",
                           TRUE ~ hh_numb)) %>%
  group_by(hh_numb, community) %>%
  summarise(# village=community,
            n_ppl=sum(!is.na(unique_id)),
            n_tested=sum(!is.na(age)),
            n_any_anemia=sum(anemia.f %in% c("mild","moderate","severe")),
            n_mild=sum(anemia.f=="mild"),
            # grouping moderate and severe together due to the low number of severe cases
            n_modsev=sum(anemia.f %in% c("moderate","severe")),
            p_any=n_any_anemia/n_tested,
            p_mild=n_mild/n_tested,
            p_modsev=n_modsev/n_tested,
            # going for medians to get village midpoint as often there were a few houses located very far away
            lat=median(gps_lat, na.rm = T),
            lon=median(gps_long, na.rm = T)) %>%
  # filtering out missing gps values
  filter(!is.na(lat) & !is.na(lon)) %>%
  ungroup()
data.hhs

# checking for duplicates
data.hhs %>% mutate(dups=n()) %>% filter(dups>1) # no duplicates

#### checking with other data ####
ts <- read_csv("loretodata_swabs.csv", guess_max = 3000) %>%
  select(unique.id, community, age, sex, hh_numb, gps_lat:gps_long, dbs.int.key, dbs.filtercode, dbs.hg, trach.int.key, trach.refused, census,
         dbs.notes) %>%
  mutate(census=if_else(census==1&is.na(trach.int.key),1,0), # xtabs(data=loreto.data,~census+is.na(trach.int.key),addNA=T)
         refused.exam=if_else(is.na(trach.int.key) | trach.refused==1,1,0)) %>% # xtabs(data=loreto.data,~refused.exam+is.na(trach.int.key),addNA=T)
         filter(census==0) %>% 
  filter(hh_numb!="control") %>%
  mutate(dbs.hg.na=ifelse(dbs.hg == -999999999, NA_real_, dbs.hg),
         hemocue=case_when(grepl("hemocue",dbs.notes)~1, # we had hemocue malfunctions in these communities, sometimes we got it to work and other times we did not
                             is.na(dbs.hg.na) & !is.na(dbs.int.key) & community %in% c("Corazon de Jesus","Nuevo Barranquita","Loma Linda","Maranatha")~1,
                             TRUE ~ 0), # addmargins(xtabs(data=anemia,~hemocue+!is.na(dbs.hg.na),addNA=T))
         coag=if_else(grepl("coag",dbs.notes),1,0), # xtabs(data=anemia,~coag,addNA=T)
         refused=if_else(is.na(dbs.int.key),1,0), # addmargins(xtabs(data=anemia,~refused+!is.na(dbs.hg.na),addNA=T))
         missing=if_else(hemocue==0&coag==0&refused==0&is.na(dbs.hg.na),1,0), # addmargins(xtabs(data=anemia,~missing+!is.na(dbs.hg.na),addNA=T))
         # creating a variable with all individuals identified above
         exclude=if_else(hemocue==1 | coag==1 | refused==1 | missing==1, 1, 0)) %>% # addmargins(xtabs(data=anemia,~exclude+!is.na(dbs.hg.na),addNA=T))
  # now I am filtering out points to get the current anemia dataset
  filter(age>=1 & age<10) %>%
  filter(exclude==0)

ts %>% group_by(community, hh_numb) %>% 
  summarise(n_ppl=sum(!is.na(unique.id)),
            n_tested=sum(!is.na(age)),
            #n_any_anemia=sum(anemia.f %in% c("mild","moderate","severe")),
            #n_mild=sum(anemia.f=="mild"),
            # grouping moderate and severe together due to the low number of severe cases
            #n_modsev=sum(anemia.f %in% c("moderate","severe")),
            #p_any=n_any_anemia/n_tested,
            #p_mild=n_mild/n_tested,
            #p_modsev=n_modsev/n_tested,
            # going for medians to get village midpoint as often there were a few houses located very far away
            lat=median(gps_lat, na.rm = T),
            lon=median(gps_long, na.rm = T))

#### descriptive code ####
ggplot(data, aes(x=p_any)) + geom_histogram(bins=15)
ggplot(data, aes(x=p_mild)) + geom_histogram(bins=15)
ggplot(data, aes(x=p_modsev)) + geom_histogram(bins=15) # none of these look normally distributed

#### global morans I ####
# all code below is based off that provided by JK

# converting my data to SF
data_spatial <- data %>%
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
ape::Moran.I(data$p_any, dist_matrix_inv) %>% as_tibble(.)
ape::Moran.I(data$p_mild, dist_matrix_inv) %>% as_tibble(.)
ape::Moran.I(data$p_modsev, dist_matrix_inv) %>% as_tibble(.) # none are significant likely due to the small sample size

#### semi-variograms ####

# create a spatial points dataframe to allow for variogram building
hhs_spatial <- data.hhs %>%
  st_as_sf(coords = c("lon","lat"), crs = crs) %>%
  as("Spatial")

# rule of thumb: limit max distance to hald of max interpoint distance
# estimate distance from each point to all other points
clu_distances <- sapply(1:nrow(hhs_spatial), function(x) geosphere::distGeo(p1 = hhs_spatial, p2 = hhs_spatial[x,]))
# divide by 1000 to turn meters into KM
max_dist <- max(clu_distances) / 2 / 1000
  
# getvariogram function estimates variograms, fitted models, and monte carlo envelopes
# input_df: input dataframe
# input_formula: formula for variogram estimation
# input_survey: survey time point to concentrate on -- does not apply to me so will not use
# input_model: model to fit to variogram
# return line, fit (sserr), effective range, and empirical values (exp_var) for variograms
get_variogram <- function(input_df, input_formula, input_model = "Exp"){
  
  # create temporary spatial dataframe based on inputs
  temp_spatial <- input_df %>%
    st_as_sf(coords = c("lon", "lat"), remove = FALSE, crs = crs) %>%
    as("Spatial")
  
  ## get variogram model fits
  # note that automap creates fewer, but larger bins than gstat (gstat results were much more unstable)
  # note: `dist` in results is the average distance of all pairs in that bin (see `gstat` documentation)
  temp_fit <- automap::autofitVariogram(formula = as.formula(input_formula),
                                        input_data = temp_spatial,
                                        model = input_model, # ability to fit multiple models at once is deprecated?
                                        kappa = c(1.5, 2.5, 3.5, 5, 10, 15, 20), # only used for Matern
                                        # miscFitOptions = list(merge.small.bins = TRUE) # alternative way to create larger bins
                                        miscFitOptions = list(min.np.bin = 10)) # automap has this feature, gstat does not 
  temp_var_model <- temp_fit$var_model
  
  # save values for fitted model
  temp_line <- gstat::variogramLine(temp_var_model, maxdist = max_dist)
  variog_line <- temp_line %>% mutate(formula = input_formula, model = input_model)
  
  # estimate effective range using nugget + 95% of sill
  temp_sill <- temp_var_model[which(temp_var_model$model == input_model), "psill"]
  temp_nugget <- temp_var_model[which(temp_var_model$model == "Nug"), "psill"]
  temp_effrange_y <- temp_sill*0.95 + temp_nugget
  temp_effrange <- ifelse(max(temp_line$gamma) < temp_effrange_y,
                          NA,
                          temp_line$dist[min(which(temp_line$gamma >= temp_effrange_y))])
  variog_stats <- data.frame(formula = input_formula, model = input_model,
                             sserr = temp_fit$sserr, effrange = temp_effrange)
  
  # save values for empirical variogram
  variog_expvar <- temp_fit$exp_var %>% dplyr::select(np, dist, gamma) %>%
    mutate(formula = input_formula, model = input_model)
  
  ret <- list(variog_line = variog_line,
              variog_stats = variog_stats,
              variog_expvar = variog_expvar)
  
  return(ret)
}
  
## `geo_perm`: helper function to conduct permutations
# input_df: input dataframe
# input_seed: for replication of permutations
# return permuted df - fixed lat/lon values, permuted outcomes
geo_perm <- function(input_df, input_seed){
  fix_locations <- input_df %>% 
    dplyr::select(lon, lat)
  
  set.seed(input_seed)
  temp_permute <- input_df %>% 
    dplyr::select(-c(lon, lat)) %>% 
    sample_n(size = nrow(.), replace = FALSE)
  
  return(bind_cols(fix_locations, temp_permute))
}

## `get_variog_perms`: function to conduct Monte Carlo variogram permutations
# input_df: input dataframe
# input_formula: formula for variogram estimation
# input_survey: survey time point to concentrate on -- JN: does not apply to this situation so will remove from code
# n_permute: number of permutations, must be >0
# return empirical values (exp_var) and 95% interval for permutations
get_variog_perms <- function(input_df, input_formula, n_permute){
  
  # keep track of seeds; return list for replication
  rand_seed_list <- sample(1:10000, size = n_permute, replace = TRUE)
  
  permute_results <- lapply(
    c(1:n_permute),
    function(x){
      rand_seed <- rand_seed_list[x]
      
      temp_permute_df <- geo_perm(input_df = input_df,
                                  input_seed = rand_seed)
      
      temp_permute_variog <- get_variogram(input_df = temp_permute_df,
                                           input_formula = input_formula)
      
      temp_permute_variog$variog_expvar %>% rename(!!paste0("gamma_", x) := gamma)
    }) %>%
    reduce(left_join, by = c("np", "dist", "formula", "model"))
  
  # add observed results
  variog_results <- get_variogram(input_df = input_df,
                                  input_formula = input_formula)
  
  permute_results <- permute_results %>% 
    left_join(variog_results$variog_expvar %>% rename(gamma_0 = gamma),
              by = c("np", "dist", "formula", "model"))
  
  # get upper and lower bounds for 95% envelope
  permute_bounds <- apply(permute_results %>% dplyr::select(starts_with("gamma")), 1, quantile, c(0.025, 0.975)) %>%
    t() %>% 
    as.data.frame() %>% 
    rename(q0.025 = "2.5%", q0.975 = "97.5%") %>% 
    mutate(fill_flag = "1") # used for creating ggplot legend
  
  ret <- list(permute_results = permute_results %>%
                bind_cols(permute_bounds) %>%
                dplyr::select(-c(model)), # note that model is not important for permutations, uses only empirical
              rand_seed_list = rand_seed_list)
  
  return(ret)
}  
  
## iterate through formulas, surveys, and models
# 200 permutations takes ~1.1 minutes per formula
# 1000 permutations takes ~6.1 minutes per formula
n_permute_variog <- 200
variog_formula_list <- c("p_any~1", "p_mild~1", "p_modsev~1", # JN: I think this is the null hypothesis, i.e., prev unrelated to location
                         "p_any~lat+lon", "p_mild~lat+lon", "p_modsev~lat+lon") # and this is the alternative
variog_results <- mapply(get_variogram,
                         input_formula = rep(variog_formula_list, each = 8),
                         # input_survey = rep(rep(survey_list, 2), 6), 
                         input_model = rep(rep(c("Mat", "Exp"), each = 4), 6),
                         MoreArgs = list(input_df = data),
                         SIMPLIFY = TRUE)

variog_perm_results <- mapply(get_variog_perms,
                              input_formula = rep(variog_formula_list, each = 4),
                              # input_survey = rep(survey_list, 6),
                              n_permute = n_permute_variog,
                              MoreArgs = list(input_df = data),
                              SIMPLIFY = TRUE)

variog_line_df <- variog_results[1,] %>% bind_rows() # class(variog_line_df)
variog_stats_df <- variog_results[2,] %>% bind_rows() # class(variog_stats_df)
variog_expvar_df <- variog_results[3,] %>% bind_rows() # class(variog_expvar_df)
variog_perm_results_df <- variog_perm_results[1,] %>% bind_rows() # class(variog_perm_results_df)

## plot variograms
# setting colors
variog_colors <- c("Exp" = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[4],
                   "Mat" = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[5])
# any anemia
plot_any <- ggplot() +
  geom_ribbon(data = variog_perm_results_df %>% filter(formula == "p_any~1"), # originally temp_permute
              aes(x = dist, ymin = q0.025, ymax = q0.975, fill = fill_flag),
              alpha = 0.3) +
  geom_point(data = variog_expvar_df %>% filter(formula == "p_any~1"), aes(x = dist, y = gamma)) + # replaced temp_expvar with variog_expvar_df
  geom_text(data = variog_expvar_df %>% filter(formula == "p_any~1"), aes(x = dist, y = 0.06*0.95, label = np),
            size = 2) + # empirical variog, does not depend on model
  geom_line(data = variog_line_df %>% filter(formula == "p_any~1"), aes(x = dist, y = gamma, color = model)) +
  geom_vline(data = variog_stats_df %>% filter(formula == "p_any~1"), aes(xintercept = effrange, color = model),
             lty = "dashed", show.legend = FALSE) +
  scale_color_manual(values = variog_colors,
                     labels = c("Exp" = "Exponential", "Mat" = "Matern")) +
  scale_fill_manual(values = c("1" = "grey"), labels = c("1" = "Monte Carlo\nenvelope")) +
  labs(x = "Distance (km)",
       y = "Semivariance",
       color = "Fitted model",
       fill = "") +
  coord_cartesian(x = c(0, max((variog_expvar_df%>% filter(formula == "p_any~1")) %>% pull(dist))),
                  y = c(0, 0.062)) +
  theme_classic() +
  theme(title = element_text(size = 10),
        panel.grid.major = element_line(color = "grey95"),
        ) + # legend.position = "none"
  guides(fill  = guide_legend(order = 2),
         color = guide_legend(order = 1))
plot_any

# mild anemia
plot_mild <- ggplot() +
  geom_ribbon(data = variog_perm_results_df %>% filter(formula == "p_mild~1"), # originally temp_permute
              aes(x = dist, ymin = q0.025, ymax = q0.975, fill = fill_flag),
              alpha = 0.3) +
  geom_point(data = variog_expvar_df %>% filter(formula == "p_mild~1"), aes(x = dist, y = gamma)) + # replaced temp_expvar with variog_expvar_df
  geom_text(data = variog_expvar_df %>% filter(formula == "p_mild~1"), aes(x = dist, y = 0.06*0.95, label = np),
            size = 2) + # empirical variog, does not depend on model
  geom_line(data = variog_line_df %>% filter(formula == "p_mild~1"), aes(x = dist, y = gamma, color = model)) +
  geom_vline(data = variog_stats_df %>% filter(formula == "p_mild~1"), aes(xintercept = effrange, color = model),
             lty = "dashed", show.legend = FALSE) +
  scale_color_manual(values = variog_colors,
                     labels = c("Exp" = "Exponential", "Mat" = "Matern")) +
  scale_fill_manual(values = c("1" = "grey"), labels = c("1" = "Monte Carlo\nenvelope")) +
  labs(x = "Distance (km)",
       y = "Semivariance",
       color = "Fitted model",
       fill = "") +
  coord_cartesian(x = c(0, max((variog_expvar_df%>% filter(formula == "p_mild~1")) %>% pull(dist))),
                  y = c(0, 0.062)) +
  theme_classic() +
  theme(title = element_text(size = 10),
        panel.grid.major = element_line(color = "grey95"),
        ) + # legend.position = "none"
  guides(fill  = guide_legend(order = 2),
         color = guide_legend(order = 1))
plot_mild 

# moderate to severe anemia
plot_modsev <- ggplot() +
  geom_ribbon(data = variog_perm_results_df %>% filter(formula == "p_modsev~1"), # originally temp_permute
              aes(x = dist, ymin = q0.025, ymax = q0.975, fill = fill_flag),
              alpha = 0.3) +
  geom_point(data = variog_expvar_df %>% filter(formula == "p_modsev~1"), aes(x = dist, y = gamma)) + # replaced temp_expvar with variog_expvar_df
  geom_text(data = variog_expvar_df %>% filter(formula == "p_modsev~1"), aes(x = dist, y = 0.06*0.95, label = np),
            size = 2) + # empirical variog, does not depend on model
  geom_line(data = variog_line_df %>% filter(formula == "p_modsev~1"), aes(x = dist, y = gamma, color = model)) +
  geom_vline(data = variog_stats_df %>% filter(formula == "p_modsev~1"), aes(xintercept = effrange, color = model),
             lty = "dashed", show.legend = FALSE) +
  scale_color_manual(values = variog_colors,
                     labels = c("Exp" = "Exponential", "Mat" = "Matern")) +
  scale_fill_manual(values = c("1" = "grey"), labels = c("1" = "Monte Carlo\nenvelope")) +
  labs(x = "Distance (km)",
       y = "Semivariance",
       color = "Fitted model",
       fill = "") +
  coord_cartesian(x = c(0, max((variog_expvar_df%>% filter(formula == "p_modsev~1")) %>% pull(dist))),
                  y = c(0, 0.062)) +
  theme_classic() +
  theme(title = element_text(size = 10),
        panel.grid.major = element_line(color = "grey95"),
        ) + # legend.position = "none"
  guides(fill  = guide_legend(order = 2),
         color = guide_legend(order = 1))
plot_modsev

# putting them all together
semivar_plots <- ggpubr::ggarrange(plot_any, plot_mild, plot_modsev,
                                   nrow = 1, common.legend = T, legend = "right")
semivar_plots

# saving
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/R projects/aa-anemia")
ggsave("figures/village_semivariogram.pdf", semivar_plots)

#### repeating semivariogram at the hh level ####
n_permute_variog <- 200

variog_results.hhs <- mapply(get_variogram,
                         input_formula = rep(variog_formula_list, each = 8),
                         # input_survey = rep(rep(survey_list, 2), 6), 
                         input_model = rep(rep(c("Mat", "Exp"), each = 4), 6),
                         MoreArgs = list(input_df = data.hhs),
                         SIMPLIFY = TRUE)

variog_perm_results.hhs <- mapply(get_variog_perms,
                              input_formula = rep(variog_formula_list, each = 4),
                              # input_survey = rep(survey_list, 6),
                              n_permute = n_permute_variog,
                              MoreArgs = list(input_df = data.hhs),
                              SIMPLIFY = TRUE)

variog_linehhs_df <- variog_results.hhs[1,] %>% bind_rows() # class(variog_line_df)
variog_statshhs_df <- variog_results.hhs[2,] %>% bind_rows() # class(variog_stats_df)
variog_expvarhhs_df <- variog_results.hhs[3,] %>% bind_rows() # class(variog_expvar_df)
variog_perm_resultshhs_df <- variog_perm_results.hhs[1,] %>% bind_rows() # class(variog_perm_results_df)

## plot variograms

# any anemia
plot_any.hh <- ggplot() +
  geom_ribbon(data = variog_perm_resultshhs_df %>% filter(formula == "p_any~1"), # originally temp_permute
              aes(x = dist, ymin = q0.025, ymax = q0.975, fill = fill_flag),
              alpha = 0.3) +
  geom_point(data = variog_expvarhhs_df %>% filter(formula == "p_any~1"), aes(x = dist, y = gamma)) + # replaced temp_expvar with variog_expvar_df
  geom_text(data = variog_expvarhhs_df %>% filter(formula == "p_any~1"), aes(x = dist, y = 0.06*0.95, label = np),
            size = 2) + # empirical variog, does not depend on model
  geom_line(data = variog_linehhs_df %>% filter(formula == "p_any~1"), aes(x = dist, y = gamma, color = model)) +
  geom_vline(data = variog_statshhs_df %>% filter(formula == "p_any~1"), aes(xintercept = effrange, color = model),
             lty = "dashed", show.legend = FALSE) +
  scale_color_manual(values = variog_colors,
                     labels = c("Exp" = "Exponential", "Mat" = "Matern", "Sph" = "Spherical")) +
  scale_fill_manual(values = c("1" = "grey"), labels = c("1" = "Monte Carlo\nenvelope")) +
  labs(x = "Distance (km)",
       y = "Semivariance",
       color = "Fitted model",
       fill = "") +
  coord_cartesian(x = c(0, max((variog_expvarhhs_df%>% filter(formula == "p_any~1")) %>% pull(dist))),
                  y = c(0, 0.25)) +
  theme_classic() +
  theme(title = element_text(size = 10),
        panel.grid.major = element_line(color = "grey95"),
  ) + # legend.position = "none"
  guides(fill  = guide_legend(order = 2),
         color = guide_legend(order = 1))
plot_any.hh

# mild anemia
plot_mild.hh <- ggplot() +
  geom_ribbon(data = variog_perm_resultshhs_df %>% filter(formula == "p_mild~1"), # originally temp_permute
              aes(x = dist, ymin = q0.025, ymax = q0.975, fill = fill_flag),
              alpha = 0.3) +
  geom_point(data = variog_expvarhhs_df %>% filter(formula == "p_mild~1"), aes(x = dist, y = gamma)) + # replaced temp_expvar with variog_expvar_df
  geom_text(data = variog_expvarhhs_df %>% filter(formula == "p_mild~1"), aes(x = dist, y = 0.06*0.95, label = np),
            size = 2) + # empirical variog, does not depend on model
  geom_line(data = variog_linehhs_df %>% filter(formula == "p_mild~1"), aes(x = dist, y = gamma, color = model)) +
  geom_vline(data = variog_statshhs_df %>% filter(formula == "p_mild~1"), aes(xintercept = effrange, color = model),
             lty = "dashed", show.legend = FALSE) +
  scale_color_manual(values = variog_colors,
                     labels = c("Exp" = "Exponential", "Mat" = "Matern", "Sph" = "Spherical")) +
  scale_fill_manual(values = c("1" = "grey"), labels = c("1" = "Monte Carlo\nenvelope")) +
  labs(x = "Distance (km)",
       y = "Semivariance",
       color = "Fitted model",
       fill = "") +
  coord_cartesian(x = c(0, max((variog_expvarhhs_df%>% filter(formula == "p_mild~1")) %>% pull(dist))),
                  y = c(0, 0.25)) +
  theme_classic() +
  theme(title = element_text(size = 10),
        panel.grid.major = element_line(color = "grey95"),
  ) + # legend.position = "none"
  guides(fill  = guide_legend(order = 2),
         color = guide_legend(order = 1))
plot_mild.hh

# moderate to severe anemia
plot_modsev.hh <- ggplot() +
  geom_ribbon(data = variog_perm_resultshhs_df %>% filter(formula == "p_modsev~1"), # originally temp_permute
              aes(x = dist, ymin = q0.025, ymax = q0.975, fill = fill_flag),
              alpha = 0.3) +
  geom_point(data = variog_expvarhhs_df %>% filter(formula == "p_modsev~1"), aes(x = dist, y = gamma)) + # replaced temp_expvar with variog_expvar_df
  geom_text(data = variog_expvarhhs_df %>% filter(formula == "p_modsev~1"), aes(x = dist, y = 0.06*0.95, label = np),
            size = 2) + # empirical variog, does not depend on model
  geom_line(data = variog_linehhs_df %>% filter(formula == "p_modsev~1"), aes(x = dist, y = gamma, color = model)) +
  geom_vline(data = variog_statshhs_df %>% filter(formula == "p_modsev~1"), aes(xintercept = effrange, color = model),
             lty = "dashed", show.legend = FALSE) +
  scale_color_manual(values = variog_colors,
                     labels = c("Exp" = "Exponential", "Mat" = "Matern", "Sph" = "Spherical")) +
  scale_fill_manual(values = c("1" = "grey"), labels = c("1" = "Monte Carlo\nenvelope")) +
  labs(x = "Distance (km)",
       y = "Semivariance",
       color = "Fitted model",
       fill = "") +
  coord_cartesian(x = c(0, max((variog_expvarhhs_df%>% filter(formula == "p_modsev~1")) %>% pull(dist))),
                  y = c(0, 0.25)) +
  theme_classic() +
  theme(title = element_text(size = 10),
        panel.grid.major = element_line(color = "grey95"),
  ) + # legend.position = "none"
  guides(fill  = guide_legend(order = 2),
         color = guide_legend(order = 1))
plot_modsev.hh

# putting them all together
semivar_plots.hh <- ggpubr::ggarrange(plot_any.hh, plot_mild.hh, plot_modsev.hh,
                                   nrow = 1, common.legend = T, legend = "right")
semivar_plots.hh

# saving it
ggsave("figures/hh_variogram.pdf", semivar_plots.hh)

#### repeating Moran's I at village level ####
# not sure if this will work as the prevalence estimates in the houses will be crude (ie 0, 0.5, or 1 given the low n in each household)

# converting my data to SF
data_spatial.hhs <- data.hhs %>% 
  # just applying this to one village
  filter(community=="Panam") %>%
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

# estimate moran's I and statistical significance for each level of anemia
ape::Moran.I(filter(data.hhs, community=="Panam")$p_any, dist_matrix_inv.hhs) %>% as_tibble(.)
ape::Moran.I(filter(data.hhs, community=="Panam")$p_mild, dist_matrix_inv.hhs) %>% as_tibble(.)
ape::Moran.I(filter(data.hhs, community=="Panam")$p_modsev, dist_matrix_inv.hhs) %>% as_tibble(.) # none are significant likely due to the small sample size

# it works for one village, now I have to create a function that will apply it to all 21
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

village.moran("p_any", "Panam") # Works!

village_list <- data.hhs %>% 
  # filtering out these two since there are only two households with anemia/gps data and this is not enough for morans I
  filter(!(community %in% c("Nuevo Barranquita","Corazon de Jesus"))) %>%
  ungroup() %>%
  # getting a unique list of the remaining village
  distinct(.$community) %>% as.vector()
village_list


#### now to iterate it over every village ####
moran.village <- mapply(village.moran,
                        outcome_var = rep(c("p_any", "p_mild", "p_modsev"), each = 19),
                        village = village_list,
                        SIMPLIFY = F)

nuevo_arica <- mapply(village.moran,
       outcome_var = c("p_any", "p_mild", "p_modsev"),
       village = "Nuevo Arica",
       SIMPLIFY = F) %>%
  bind_rows()

dos_mayo <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Dos de Mayo",
                      SIMPLIFY = F) %>%
  bind_rows()

bellavista_b <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Bellavista (Balsapuerto)",
                      SIMPLIFY = F) %>%
  bind_rows()

nuevo_papaplaya <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Nuevo Papaplaya",
                      SIMPLIFY = F) %>%
  bind_rows()

nuevo_barranquita <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Nuevo Barranquita",
                      SIMPLIFY = F) %>%
  bind_rows() # duplicate error, simply two households so not enough data for Moran's I

corazon_de_jesus <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village ="Corazon de Jesus",
                      SIMPLIFY = F) %>%
  bind_rows() # duplicate error

bethel <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Bethel",
                      SIMPLIFY = F) %>%
  bind_rows()

bellavista_j <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Bellavista (Jeberos)",
                      SIMPLIFY = F) %>%
  bind_rows()

vista_allegre <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Vista Allegre",
                      SIMPLIFY = F) %>%
  bind_rows()

huancayo <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Huancayo",
                      SIMPLIFY = F) %>%
  bind_rows()

nuevo_arica <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Nuevo Arica",
                      SIMPLIFY = F) %>%
  bind_rows()

sies_julio <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "06 de Julio",
                      SIMPLIFY = F) %>%
  bind_rows()

tamarate <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Tamarate",
                      SIMPLIFY = F) %>%
  bind_rows()

huatapi <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Huatapi",
                      SIMPLIFY = F) %>%
  bind_rows()

union_campesino <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Union Campesino",
                      SIMPLIFY = F) %>%
  bind_rows()

nuevo_iquitos <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Nuevo Iquitos",
                      SIMPLIFY = F) %>%
  bind_rows()

angamos <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Angamos",
                      SIMPLIFY = F) %>%
  bind_rows()

san_antonio <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "San Antonio",
                      SIMPLIFY = F) %>%
  bind_rows()

panam <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Panam",
                      SIMPLIFY = F) %>%
  bind_rows()

maranatha <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Maranatha",
                      SIMPLIFY = F) %>%
  bind_rows()

loma_linda <- mapply(village.moran,
                      outcome_var = c("p_any", "p_mild", "p_modsev"),
                      village = "Loma Linda",
                      SIMPLIFY = F) %>%
  bind_rows()

centro_america <- mapply(village.moran,
                     outcome_var = c("p_any", "p_mild", "p_modsev"),
                     village = "Centro America",
                     SIMPLIFY = F) %>%
  bind_rows() # duplicate error?

village.morans <- rbind(centro_america, loma_linda, maranatha, panam, san_antonio, angamos, 
                        nuevo_iquitos, union_campesino, huatapi, tamarate, sies_julio, nuevo_arica, 
                        huancayo, vista_allegre, bellavista_j, bethel, nuevo_papaplaya, bellavista_b, 
                        dos_mayo, nuevo_arica)
village.morans

ggplot(village.morans, aes(x=observed)) + 
  geom_histogram(bins = 20) + 
  facet_grid(~outcome_var)

village.morans %>%
  group_by(outcome_var) %>%
  summarise(mean=mean(observed),
            sd=sd(observed))
  

#### question: to repeat semivariograms at village level ####
# https://link.springer.com/chapter/10.1007/978-94-011-1739-5_14
  # according to this you need about 50 data points to make the semivariogram worthwhile
  # and my class notes recommend >50 points in total and each bin (group of points of similar distances apart) should have ≥30 pairwise comparisons
  # lastly the data needs to be normally distributed, otherwise it needs to be transformed.

#### old code ####
# # population data from 2017 alto amazonas health network census
# # this will be used to calculate age and sex adjusted prevalence estimates 
# AAcensus <- read_csv("2018AA Census.csv") %>%
#   # filtering out yurimaguas and total since we did not include yurimaguas in our sampling frame
#   filter(MICRORED!="YURIMAGUAS" & MICRORED!="TOTAL") %>%
#   # selecting relevant columns 
#   select(MICRORED, `1a`,`2a`,`3a`,`4a`,`5a`,`6a`,`7a`,`8a`,`9a`) %>% select(-MICRORED) %>%
#   mutate(`1`=sum(`1a`),
#          `2`=sum(`2a`),
#          `3`=sum(`3a`),
#          `4`=sum(`4a`),
#          `5`=sum(`5a`),
#          `6`=sum(`6a`),
#          `7`=sum(`7a`),
#          `8`=sum(`8a`),
#          `9`=sum(`9a`),
#          total=sum(`1a`)+sum(`2a`)+sum(`3a`)+sum(`4a`)+sum(`5a`)+sum(`6a`)+sum(`7a`)+sum(`8a`)+sum(`9a`)) %>% select(-(`1a`:`9a`)) %>% unique() %>%
#   pivot_longer(cols = `1`:total, names_to ="age.group19") %>%
#   mutate(weight=value/11599)
# # creating the weights
# weights19 <- AAcensus %>% filter(age.group19!="total") %>% mutate(age.group19=as.numeric(age.group19))
# # trying to open a different dataset
# data.hhs2 <- read_csv("loretodata_swabs.csv", guess_max = 3000) %>% 
#   select(district:unique.id, dbs.hg, gps_lat:gps_long, census, trach.int.key, trach.refused, dbs.notes,
#          dbs.int.key) %>%
#   # applying the same selection criteria for the data that I sent to Jorge
#   mutate(census=if_else(census==1&is.na(trach.int.key),1,0), # xtabs(data=loreto.data,~census+is.na(trach.int.key),addNA=T)
#          refused.exam=if_else(is.na(trach.int.key) | trach.refused==1,1,0),
#          dbs.hg.na=if_else(dbs.hg == -999999999, NA_real_, dbs.hg)) %>% # xtabs(data=loreto.data,~refused.exam+is.na(trach.int.key),addNA=T)
#   # view(loreto.data %>% filter(refused.exam==1&!is.na(trach.int.key)))
#   mutate(hemocue=case_when(grepl("hemocue",dbs.notes)~1, # we had hemocue malfunctions in these communities, sometimes we got it to work and other times we did not
#                            is.na(dbs.hg.na) & !is.na(dbs.int.key) & community %in% c("Corazon de Jesus","Nuevo Barranquita","Loma Linda","Maranatha")~1,
#                            TRUE ~ 0), # addmargins(xtabs(data=anemia,~hemocue+!is.na(dbs.hg.na),addNA=T))
#          coag=if_else(grepl("coag",dbs.notes),1,0), # xtabs(data=anemia,~coag,addNA=T)
#          refused=if_else(is.na(dbs.int.key),1,0), # addmargins(xtabs(data=anemia,~refused+!is.na(dbs.hg.na),addNA=T))
#          missing=if_else(hemocue==0&coag==0&refused==0&is.na(dbs.hg.na),1,0), # addmargins(xtabs(data=anemia,~missing+!is.na(dbs.hg.na),addNA=T))
#          # creating a variable with all individuals identified above
#          exclude=if_else(hemocue==1 | coag==1 | refused==1 | missing==1, 1, 0))
# 
# 
# data.hhs2 %>% group_by(hh_numb, community) %>% 
#   # applying same selection criteria as those for the original data I sent to jorge
#   filter(census==0) %>% filter(hh_numb!="control") %>%
#   filter(exclude==0) %>% filter(age<10) %>% filter(age>=1) %>% 
#   summarise( # village=community,
#     n_ppl=sum(!is.na(unique.id)),
#     n_tested=sum(!is.na(age)),
#     # n_any_anemia=sum(anemia.f %in% c("mild","moderate","severe")),
#     # n_mild=sum(anemia.f=="mild"),
#     # grouping moderate and severe together due to the low number of severe cases
#     # n_modsev=sum(anemia.f %in% c("moderate","severe")),
#     # p_any=n_any_anemia/n_tested,
#     # p_mild=n_mild/n_tested,
#     # p_modsev=n_modsev/n_tested,
#     # going for medians to get village midpoint as often there were a few houses located very far away
#     lat=median(gps_lat, na.rm = T),
#     lon=median(gps_long, na.rm = T)) %>%
#   # filtering out missing gps values
#   filter(!is.na(lat) & !is.na(lon))
# 
# 
# 
# summarise(n=sum(!is.na(unique.id))) %>%
#   group_by(hh_numb) %>%
#   mutate(dups=n()) %>% filter(dups>1) # no dups when I do it this way
# 
# data.hhs2 %>% group_by(hh_numb) %>%
#   mutate(dups=n()) %>% filter(dups>1) # plenty of duplicates here....

