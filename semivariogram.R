# ALTO AMAZONAS ANEMIA SPATIAL ANALYSIS @ HH LEVEL
# JOHN NESEMANN

#### semivariogram code ####
# create a spatial points dataframe to allow for variogram building
hhs_spatial <- data.hhs %>%
  st_as_sf(coords = c("lon","lat"), crs = crs) %>%
  as("Spatial")

# trying to see what model the autofitVariogram function chooses
ts <- automap::autofitVariogram(formula = mean_hg ~ 1,
                          input_data = data.hhs %>% 
                            st_as_sf(coords = c("lon","lat"),
                                     remove = F,
                                     crs = crs) %>% as("Spatial"),
                          miscFitOptions = list(merge.small.bins = T,
                                                min.np.bin = 50))
ts
ts %>% plot(.)

ts2 <- automap::autofitVariogram(formula = median_hg ~ 1,
                          input_data = data.hhs %>% 
                            st_as_sf(coords = c("lon","lat"),
                                     remove = F,
                                     crs = crs) %>% as("Spatial"),
                          miscFitOptions = list(merge.small.bins = T,
                                                min.np.bin = 50)) 
ts2
ts2 %>% plot(.)

ts3 <- automap::autofitVariogram(formula = rand_hg ~ 1,
                          input_data = data.hhs %>% 
                            st_as_sf(coords = c("lon","lat"),
                                     remove = F,
                                     crs = crs) %>% as("Spatial"),
                          miscFitOptions = list(merge.small.bins = T,
                                                min.np.bin = 50)) 
ts3
ts3 %>% plot(.)

# rule of thumb: limit max distance to half of max interpoint distance
# estimate distance from each point to all other points
clu_distances <- sapply(1:nrow(hhs_spatial), 
                        function(x) geosphere::distGeo(p1 = hhs_spatial, 
                                                       p2 = hhs_spatial[x,]))
# divide by 1000 to turn meters into KM
max_dist <- max(clu_distances) / 2 / 1000

# minimum lag distance shoudl be > smallest inter-point distance
min(clu_distances)

# bins should have >30 pairwise comparisons, will set it to 150 for this analysis 
# given the large number of comparisons

# getvariogram function estimates variograms, fitted models, and monte carlo envelopes
# input_df: input dataframe
# input_formula: formula for variogram estimation
# input_survey: survey time point to concentrate on -- does not apply to me so will not use
# input_model: model to fit to variogram
# return line, fit (sserr), effective range, and empirical values (exp_var) for variograms
get_variogram <- function(input_df, input_formula, input_model = "Sph"){
  
  # create temporary spatial dataframe based on inputs
  temp_spatial <- input_df %>%
    st_as_sf(coords = c("lon", "lat"), remove = FALSE, crs = crs) %>%
    as("Spatial")
  
  ## get variogram model fits
  # note that automap creates fewer, but larger bins than gstat 
  # (gstat results were much more unstable)
  # note: `dist` in results is the average distance of all pairs in that bin 
  # (see `gstat` documentation)
  temp_fit <- automap::autofitVariogram(formula = as.formula(input_formula),
                                        input_data = temp_spatial,
                                        model = input_model, 
                                        # ability to fit multiple models is deprecated?
                                        kappa = c(1.5, 2.5, 3.5, 5, 10, 15, 20), 
                                        # only used for Matern
                                        # automap has this feature, gstat does not 
                                        miscFitOptions = list(min.np.bin = 150,
                                                              merge.small.bins = T)) 
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
# input_survey: survey time point to concentrate on J
# N: does not apply to this situation so will remove from code
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
  permute_bounds <- apply(permute_results %>% dplyr::select(starts_with("gamma")), 
                          1, quantile, c(0.025, 0.975)) %>%
    t() %>% 
    as.data.frame() %>% 
    rename(q0.025 = "2.5%", q0.975 = "97.5%") %>% 
    mutate(fill_flag = "1") # used for creating ggplot legend
  
  ret <- list(permute_results = permute_results %>%
                bind_cols(permute_bounds) %>%
                dplyr::select(-c(model)), 
              # note that model is not important for permutations, uses only empirical
              rand_seed_list = rand_seed_list)
  
  return(ret)
} 

#### semivariogram at the hh level ####
n_permute_variog <- 200
variog_formula_list <- c("mean_hg~1", "median_hg~1", "rand_hg~1",
                         # JN: I think this is the null hypothesis
                         "mean_hg~lat+lon", "median_hg~lat+lon","rand_hg~lat+lon") # and this is the alternative

variog_results.hhs <- mapply(get_variogram,
                             input_formula = rep(variog_formula_list, each = 8),
                             # input_survey = rep(rep(survey_list, 2), 6), 
                             input_model = rep(rep(c("Sph", "Ste"), each = 4), 6),
                             MoreArgs = list(input_df = data.hhs),
                             SIMPLIFY = TRUE)

variog_perm_results.hhs <- mapply(get_variog_perms,
                                  input_formula = rep(variog_formula_list, each = 8),
                                  # input_survey = rep(survey_list, 6),
                                  n_permute = n_permute_variog,
                                  MoreArgs = list(input_df = data.hhs),
                                  SIMPLIFY = TRUE) 

variog_linehhs_df <- variog_results.hhs[1,] %>% bind_rows() 
# class(variog_line_df)
variog_statshhs_df <- variog_results.hhs[2,] %>% bind_rows() 
# class(variog_stats_df)
variog_expvarhhs_df <- variog_results.hhs[3,] %>% bind_rows() 
# class(variog_expvar_df)
variog_perm_resultshhs_df <- variog_perm_results.hhs[1,] %>% bind_rows() 
# class(variog_perm_results_df)

## plot variograms
# setting colors
variog_colors <- c("Sph" = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[4],
                   "Ste" = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[5])

# mean_hg - use Sph for this
plot_mean.hh <- ggplot() +
  geom_ribbon(data = variog_perm_resultshhs_df %>% filter(formula == "mean_hg~1"), 
              # originally temp_permute
              aes(x = dist, ymin = q0.025, ymax = q0.975, fill = fill_flag),
              alpha = 0.3) +
  geom_point(data = variog_expvarhhs_df %>% filter(formula == "mean_hg~1"), 
             aes(x = dist, y = gamma)) + # replaced temp_expvar with variog_expvar_df
  geom_text(data = variog_expvarhhs_df %>% filter(formula == "mean_hg~1"), 
            aes(x = dist, y = 0.06*0.95, label = np, size = 8),
            size = 2) + # empirical variog, does not depend on model
  geom_line(data = variog_linehhs_df %>% filter(formula == "mean_hg~1" & model == "Sph"), 
            aes(x = dist, y = gamma, color = model)) +
  geom_vline(data = variog_statshhs_df %>% filter(formula == "mean_hg~1" & model == "Sph"), 
             aes(xintercept = effrange, color = model),
             lty = "dashed", show.legend = FALSE) +
  scale_color_manual(values = variog_colors,
                     labels = c("Sph" = "Spherical", "Ste" = "Matern")) +
  scale_fill_manual(values = c("1" = "grey"), labels = c("1" = "Monte Carlo\nenvelope")) +
  labs(x = "Distance (km)",
       y = "Semivariance",
       color = "Fitted model",
       fill = "") +
  coord_cartesian(x = c(0, max((variog_expvarhhs_df%>% filter(formula == "mean_hg~1")) %>% 
                                 pull(dist))),
                  y = c(0, 1.5)) +
  # max((variog_expvarhhs_df%>% filter(formula == "mean_hg~1")) %>% pull(gamma))
  theme_classic() +
  theme(title = element_text(size = 10),
        panel.grid.major = element_line(color = "grey95"),
  ) + # legend.position = "none"
  guides(fill  = guide_legend(order = 2),
         color = guide_legend(order = 1)) + 
  labs(title = "A) Mean Hemoglobin")
plot_mean.hh

# median hg - use Sph for this
plot_median.hh <- ggplot() +
  geom_ribbon(data = variog_perm_resultshhs_df %>% filter(formula == "median_hg~1"), 
              aes(x = dist, ymin = q0.025, ymax = q0.975, fill = fill_flag),
              alpha = 0.3) +
  geom_point(data = variog_expvarhhs_df %>% filter(formula == "median_hg~1"), 
             aes(x = dist, y = gamma)) + # replaced temp_expvar with variog_expvar_df
  geom_text(data = variog_expvarhhs_df %>% filter(formula == "median_hg~1"), 
            aes(x = dist, y = 0.06*0.95, label = np),
            size = 2) + # empirical variog, does not depend on model
  geom_line(data = variog_linehhs_df %>% filter(formula == "median_hg~1" & model == "Sph"), 
            aes(x = dist, y = gamma, color = model)) +
  geom_vline(data = variog_statshhs_df %>% filter(formula == "median_hg~1" & model == "Sph"), 
             aes(xintercept = effrange, color = model),
             lty = "dashed", show.legend = FALSE) +
  scale_color_manual(values = variog_colors,
                     labels = c("Sph" = "Spherical", "Ste" = "Matern")) +
  scale_fill_manual(values = c("1" = "grey"), labels = c("1" = "Monte Carlo\nenvelope")) +
  labs(x = "Distance (km)",
       y = "Semivariance",
       color = "Fitted model",
       fill = "") +
  coord_cartesian(x = c(0, max((variog_expvarhhs_df%>%filter(formula == "median_hg~1"))%>% pull(dist))),
                  y = c(0, 1.5)) +
  # max((variog_expvarhhs_df%>% filter(formula == "median_hg~1")) %>% pull(gamma)))
  theme_classic() +
  theme(title = element_text(size = 10),
        panel.grid.major = element_line(color = "grey95"),
  ) + # legend.position = "none"
  guides(fill  = guide_legend(order = 2),
         color = guide_legend(order = 1)) + 
  labs(title = "B) Median Hemoglobin")
plot_median.hh

# rand hg -- use Ste for this
plot_rand.hh <- ggplot() +
  geom_ribbon(data = variog_perm_resultshhs_df %>% filter(formula == "rand_hg~1"), 
              # originally temp_permute
              aes(x = dist, ymin = q0.025, ymax = q0.975, fill = fill_flag),
              alpha = 0.3) +
  geom_point(data = variog_expvarhhs_df %>% filter(formula == "rand_hg~1"), 
             aes(x = dist, y = gamma)) + # replaced temp_expvar with variog_expvar_df
  geom_text(data = variog_expvarhhs_df %>% filter(formula == "rand_hg~1"), 
            aes(x = dist, y = 0.06*0.95, label = np),
            size = 2) + # empirical variog, does not depend on model
  geom_line(data = variog_linehhs_df %>% filter(formula == "rand_hg~1" & model == "Ste"), 
            aes(x = dist, y = gamma, color = model)) +
  geom_vline(data = variog_statshhs_df %>% filter(formula == "rand_hg~1" & model == "Ste"), 
             aes(xintercept = effrange, color = model),
             lty = "dashed", show.legend = FALSE) +
  scale_color_manual(values = variog_colors,
                     labels = c("Sph" = "Spherical", "Ste" = "Matern")) +
  scale_fill_manual(values = c("1" = "grey"), labels = c("1" = "Monte Carlo\nenvelope")) +
  labs(x = "Distance (km)",
       y = "Semivariance",
       color = "Fitted model",
       fill = "") +
  coord_cartesian(x = c(0, max((variog_expvarhhs_df%>% filter(formula == "rand_hg~1")) %>% 
                                 pull(dist))),
                  y = c(0, 1.5)) +
  # max((variog_expvarhhs_df%>% filter(formula == "rand_hg~1")) %>% pull(gamma))
  theme_classic() +
  theme(title = element_text(size = 10),
        panel.grid.major = element_line(color = "grey95"),
  ) + # legend.position = "none"
  guides(fill  = guide_legend(order = 2),
         color = guide_legend(order = 1)) + 
  labs(title = "C) Random Hemoglobin")
plot_rand.hh

# putting them all together
semivar_plots.hh <- ggpubr::ggarrange(plot_mean.hh, plot_median.hh, plot_rand.hh,
                                      nrow = 1, common.legend = T, legend = "right")
semivar_plots.hh

# saving it
# setwd("~/Library/Mobile Documents/com~apple~CloudDocs/R projects/aa-anemia")
ggsave(here("figures","hh_semivariogram.eps"), semivar_plots.hh, 
       device = cairo_ps, # this allows me to save eps with alpha
       height = 6, width = 15)

# getting fitted values for the semi-variogram
# will just do this manually for each type of model and formula

# mean hg
automap::autofitVariogram(formula = mean_hg ~ 1,
                          input_data = data.hhs %>% 
                            st_as_sf(coords = c("lon","lat"),
                                     remove = F,
                                     crs = crs) %>% as("Spatial"),
                          miscFitOptions = list(merge.small.bins = T,
                                                min.np.bin = 150),
                          model = "Sph") %>% plot(.)
automap::autofitVariogram(formula = mean_hg ~ 1,
                          input_data = data.hhs %>% 
                            st_as_sf(coords = c("lon","lat"),
                                     remove = F,
                                     crs = crs) %>% as("Spatial"),
                          miscFitOptions = list(merge.small.bins = T,
                                                min.np.bin = 150),
                          model = "Ste") %>% plot(.)

# median hg
automap::autofitVariogram(formula = median_hg ~ 1,
                          input_data = data.hhs %>% 
                            st_as_sf(coords = c("lon","lat"),
                                     remove = F,
                                     crs = crs) %>% as("Spatial"),
                          miscFitOptions = list(merge.small.bins = T,
                                                min.np.bin = 150)) %>% plot(.)
automap::autofitVariogram(formula = median_hg ~ 1,
                          input_data = data.hhs %>% 
                            st_as_sf(coords = c("lon","lat"),
                                     remove = F,
                                     crs = crs) %>% as("Spatial"),
                          miscFitOptions = list(merge.small.bins = T,
                                                min.np.bin = 150),
                          model = "Sph") %>% plot(.)
automap::autofitVariogram(formula = median_hg ~ 1,
                          input_data = data.hhs %>% 
                            st_as_sf(coords = c("lon","lat"),
                                     remove = F,
                                     crs = crs) %>% as("Spatial"),
                          miscFitOptions = list(merge.small.bins = T,
                                                min.np.bin = 150),
                          model = "Ste") %>% plot(.)

# random hg
automap::autofitVariogram(formula = rand_hg ~ 1,
                          input_data = data.hhs %>% 
                            st_as_sf(coords = c("lon","lat"),
                                     remove = F,
                                     crs = crs) %>% as("Spatial"),
                          miscFitOptions = list(merge.small.bins = T,
                                                min.np.bin = 150),
                          model = "Sph") %>% plot(.)
automap::autofitVariogram(formula = rand_hg ~ 1,
                          input_data = data.hhs %>% 
                            st_as_sf(coords = c("lon","lat"),
                                     remove = F,
                                     crs = crs) %>% as("Spatial"),
                          miscFitOptions = list(merge.small.bins = T,
                                                min.np.bin = 150),
                          model = "Ste") %>% plot(.)


# #### OLD CODE - village level semivariogram ####
# #### semi-variograms ####
# 
# # create a spatial points dataframe to allow for variogram building
# data_spatial <- data %>%
#   st_as_sf(coords = c("lon","lat"), crs = crs) %>%
#   as("Spatial")
# 
# # trying to see what model the autofitVariogram function chooses
# automap::autofitVariogram(formula = any.logit ~ 1,
#                           input_data = data %>% 
#                             st_as_sf(coords = c("lon","lat"),
#                                      remove = F,
#                                      crs = crs) %>% as("Spatial"),
#                           miscFitOptions = list(merge.small.bins = T,
#                                                 min.np.bin = 10)) %>% plot(.)
# # Gau
# 
# automap::autofitVariogram(formula = modsev.logit ~ 1,
#                           input_data = data %>% 
#                             st_as_sf(coords = c("lon","lat"),
#                                      remove = F,
#                                      crs = crs) %>% as("Spatial"),
#                           miscFitOptions = list(merge.small.bins = T,
#                                                 min.np.bin = 10)) %>% plot(.)
# # Sph -- no spatial structure to the data
# 
# automap::autofitVariogram(formula = mean_hg ~ 1,
#                           input_data = data %>% 
#                             st_as_sf(coords = c("lon","lat"),
#                                      remove = F,
#                                      crs = crs) %>% as("Spatial"),
#                           miscFitOptions = list(merge.small.bins = T,
#                                                 min.np.bin = 10),
#                           model = "Sph") %>% plot(.)
# # Sph
# 
# automap::autofitVariogram(formula = median_hg ~ 1,
#                           input_data = data %>% 
#                             st_as_sf(coords = c("lon","lat"),
#                                      remove = F,
#                                      crs = crs) %>% as("Spatial"),
#                           miscFitOptions = list(merge.small.bins = T,
#                                                 min.np.bin = 10),
#                           model = "Sph") %>% plot(.)
# # Sph
# 
# # rule of thumb: limit max distance to hald of max interpoint distance
# # estimate distance from each point to all other points
# clu_distances <- sapply(1:nrow(data_spatial), 
#                         function(x) geosphere::distGeo(p1 = data_spatial, 
#                                                        p2 = data_spatial[x,]))
# # divide by 1000 to turn meters into KM
# max_dist <- max(clu_distances) / 2 / 1000
# 
# # getvariogram function estimates variograms, fitted models, and monte carlo envelopes
# # input_df: input dataframe
# # input_formula: formula for variogram estimation
# # input_survey: survey time point to concentrate on
# # input_model: model to fit to variogram
# # return line, fit (sserr), effective range, and empirical values (exp_var) for variograms
# get_variogram <- function(input_df, input_formula, input_model = "Exp"){
#   
#   # create temporary spatial dataframe based on inputs
#   temp_spatial <- input_df %>%
#     st_as_sf(coords = c("lon", "lat"), remove = FALSE, crs = crs) %>%
#     as("Spatial")
#   
#   ## get variogram model fits
#   # note that automap creates fewer, but larger bins than gstat 
#   # (gstat results were much more unstable)
#   # note: `dist` in results is the average distance of all pairs in that bin 
#   # (see `gstat` documentation)
#   temp_fit <- automap::autofitVariogram(formula = as.formula(input_formula),
#                                         input_data = temp_spatial,
#                                         model = input_model, 
#                                         # ability to fit multiple models is deprecated?
#                                         kappa = c(1.5, 2.5, 3.5, 5, 10, 15, 20), 
#                                         # only used for Matern
#                                         miscFitOptions = list(min.np.bin = 10,
#                                                               merge.small.bins = T)) 
#   # automap has this feature, gstat does not 
#   
#   temp_var_model <- temp_fit$var_model
#   
#   # save values for fitted model
#   temp_line <- gstat::variogramLine(temp_var_model, maxdist = max_dist)
#   variog_line <- temp_line %>% mutate(formula = input_formula, model = input_model)
#   
#   # estimate effective range using nugget + 95% of sill
#   temp_sill <- temp_var_model[which(temp_var_model$model == input_model), "psill"]
#   temp_nugget <- temp_var_model[which(temp_var_model$model == "Nug"), "psill"]
#   temp_effrange_y <- temp_sill*0.95 + temp_nugget
#   temp_effrange <- ifelse(max(temp_line$gamma) < temp_effrange_y,
#                           NA,
#                           temp_line$dist[min(which(temp_line$gamma >= temp_effrange_y))])
#   variog_stats <- data.frame(formula = input_formula, model = input_model,
#                              sserr = temp_fit$sserr, effrange = temp_effrange)
#   
#   # save values for empirical variogram
#   variog_expvar <- temp_fit$exp_var %>% dplyr::select(np, dist, gamma) %>%
#     mutate(formula = input_formula, model = input_model)
#   
#   ret <- list(variog_line = variog_line,
#               variog_stats = variog_stats,
#               variog_expvar = variog_expvar)
#   
#   return(ret)
# }
# 
# ## `geo_perm`: helper function to conduct permutations
# # input_df: input dataframe
# # input_seed: for replication of permutations
# # return permuted df - fixed lat/lon values, permuted outcomes
# geo_perm <- function(input_df, input_seed){
#   fix_locations <- input_df %>% 
#     dplyr::select(lon, lat)
#   
#   set.seed(input_seed)
#   temp_permute <- input_df %>% 
#     dplyr::select(-c(lon, lat)) %>% 
#     sample_n(size = nrow(.), replace = FALSE)
#   
#   return(bind_cols(fix_locations, temp_permute))
# }
# 
# ## `get_variog_perms`: function to conduct Monte Carlo variogram permutations
# # input_df: input dataframe
# # input_formula: formula for variogram estimation
# # input_survey: survey time point to concentrate on
# # n_permute: number of permutations, must be >0
# # return empirical values (exp_var) and 95% interval for permutations
# get_variog_perms <- function(input_df, input_formula, n_permute){
#   
#   # keep track of seeds; return list for replication
#   rand_seed_list <- sample(1:10000, size = n_permute, replace = TRUE)
#   
#   permute_results <- lapply(
#     c(1:n_permute),
#     function(x){
#       rand_seed <- rand_seed_list[x]
#       
#       temp_permute_df <- geo_perm(input_df = input_df,
#                                   input_seed = rand_seed)
#       
#       temp_permute_variog <- get_variogram(input_df = temp_permute_df,
#                                            input_formula = input_formula)
#       
#       temp_permute_variog$variog_expvar %>% rename(!!paste0("gamma_", x) := gamma)
#     }) %>%
#     reduce(left_join, by = c("np", "dist", "formula", "model"))
#   
#   # add observed results
#   variog_results <- get_variogram(input_df = input_df,
#                                   input_formula = input_formula)
#   
#   permute_results <- permute_results %>% 
#     left_join(variog_results$variog_expvar %>% rename(gamma_0 = gamma),
#               by = c("np", "dist", "formula", "model"))
#   
#   # get upper and lower bounds for 95% envelope
#   permute_bounds <- apply(permute_results %>% dplyr::select(starts_with("gamma")), 
#                           1, quantile, c(0.025, 0.975)) %>%
#     t() %>% 
#     as.data.frame() %>% 
#     rename(q0.025 = "2.5%", q0.975 = "97.5%") %>% 
#     mutate(fill_flag = "1") # used for creating ggplot legend
#   
#   ret <- list(permute_results = permute_results %>%
#                 bind_cols(permute_bounds) %>%
#                 dplyr::select(-c(model)), 
#               # note that model is not important for permutations, uses only empirical
#               rand_seed_list = rand_seed_list)
#   
#   return(ret)
# }  
# 
# ## iterate through formulas, surveys, and models
# # 200 permutations takes ~1.1 minutes per formula
# # 1000 permutations takes ~6.1 minutes per formula
# n_permute_variog <- 200
# variog_formula_list <- c("any.logit~1", "modsev.logit~1", "mean_hg~1", "median_hg~1", 
#                          # JN: I think this is the null hypothesis, i.e., prev unrelated to location
#                          "any.logit~lat+lon", "modsev.logit~lat+lon", "mean_hg~lat+lon", 
#                          "median_hg~lat+lon") # and this is the alternative
# 
# variog_results <- mapply(get_variogram,
#                          input_formula = rep(variog_formula_list, each = 8),
#                          # input_survey = rep(rep(survey_list, 2), 6), 
#                          input_model = rep(rep(c("Gau", "Sph"), each = 4), 6),
#                          MoreArgs = list(input_df = data),
#                          SIMPLIFY = TRUE)
# 
# variog_perm_results <- mapply(get_variog_perms,
#                               input_formula = rep(variog_formula_list, each = 4),
#                               # input_survey = rep(survey_list, 6),
#                               n_permute = n_permute_variog,
#                               MoreArgs = list(input_df = data),
#                               SIMPLIFY = TRUE)
# 
# variog_line_df <- variog_results[1,] %>% bind_rows() # class(variog_line_df)
# variog_stats_df <- variog_results[2,] %>% bind_rows() # class(variog_stats_df)
# variog_expvar_df <- variog_results[3,] %>% bind_rows() # class(variog_expvar_df)
# variog_perm_results_df <- variog_perm_results[1,] %>% bind_rows() # class(variog_perm_results_df)
# 
# ## plot variograms
# # setting colors
# variog_colors <- c("Sph" = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[4],
#                    "Gau" = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[5])
# # any anemia
# plot_any <- ggplot() +
#   geom_ribbon(data = variog_perm_results_df %>% filter(formula == "any.logit~1"), 
#               # originally temp_permute
#               aes(x = dist, ymin = q0.025, ymax = q0.975, fill = fill_flag),
#               alpha = 0.3) +
#   geom_point(data = variog_expvar_df %>% filter(formula == "any.logit~1"), 
#              aes(x = dist, y = gamma)) + # replaced temp_expvar with variog_expvar_df
#   geom_text(data = variog_expvar_df %>% filter(formula == "any.logit~1"), 
#             aes(x = dist, y = 0.06*0.95, label = np),
#             size = 2) + # empirical variog, does not depend on model
#   geom_line(data = variog_line_df %>% filter(formula == "any.logit~1" & model == "Gau"), 
#             aes(x = dist, y = gamma, color = model)) +
#   geom_vline(data = variog_stats_df %>% filter(formula == "any.logit~1" & model == "Gau"), 
#              aes(xintercept = effrange, color = model),
#              lty = "dashed", show.legend = FALSE) +
#   scale_color_manual(values = variog_colors,
#                      labels = c("Sph" = "Spherical", "Gau" = "Gaussian")) +
#   scale_fill_manual(values = c("1" = "grey"), labels = c("1" = "Monte Carlo\nenvelope")) +
#   labs(x = "Distance (km)",
#        y = "Semivariance",
#        color = "Fitted model",
#        fill = "") +
#   coord_cartesian(x = c(0, max((variog_expvar_df%>% filter(formula == "any.logit~1")) %>% 
#                                  pull(dist))),
#                   y = c(0, 1)) +
#   # max((variog_expvar_df%>% filter(formula == "any.logit~1")) %>% pull(gamma))
#   theme_classic() +
#   theme(title = element_text(size = 9),
#         panel.grid.major = element_line(color = "grey95"),
#   ) + # legend.position = "none"
#   guides(fill  = guide_legend(order = 2),
#          color = guide_legend(order = 1)) + 
#   labs(title = "A) Prevalence of mild-or-worse anemia")
# plot_any
# 
# # moderate to severe anemia
# plot_modsev <- ggplot() +
#   geom_ribbon(data = variog_perm_results_df %>% filter(formula == "modsev.logit~1"), 
#               # originally temp_permute
#               aes(x = dist, ymin = q0.025, ymax = q0.975, fill = fill_flag),
#               alpha = 0.3) +
#   geom_point(data = variog_expvar_df %>% filter(formula == "modsev.logit~1"), 
#              aes(x = dist, y = gamma)) + # replaced temp_expvar with variog_expvar_df
#   geom_text(data = variog_expvar_df %>% filter(formula == "modsev.logit~1"), 
#             aes(x = dist, y = 0.06*0.95, label = np),
#             size = 2) + # empirical variog, does not depend on model
#   geom_line(data = variog_line_df %>% filter(formula == "modsev.logit~1" &model=="Sph"), 
#             aes(x = dist, y = gamma, color = model)) +
#   geom_vline(data = variog_stats_df %>% filter(formula == "modsev.logit~1" &model=="Sph"), 
#              aes(xintercept = effrange, color = model),
#              lty = "dashed", show.legend = FALSE) +
#   scale_color_manual(values = variog_colors,
#                      labels = c("Sph" = "Spherical", "Gau" = "Gaussian")) +
#   scale_fill_manual(values = c("1" = "grey"), labels = c("1" = "Monte Carlo\nenvelope")) +
#   labs(x = "Distance (km)",
#        y = "Semivariance",
#        color = "Fitted model",
#        fill = "") +
#   coord_cartesian(x = c(0, max((variog_expvar_df%>%filter(formula=="modsev.logit~1")) %>% 
#                                  pull(dist))),
#                   y = c(0, 1.0)) +
#   # max((variog_expvar_df%>%filter(formula=="modsev.logit~1")) %>% pull(gamma))
#   theme_classic() +
#   theme(title = element_text(size = 9),
#         panel.grid.major = element_line(color = "grey95"),
#   ) + # legend.position = "none"
#   guides(fill  = guide_legend(order = 2),
#          color = guide_legend(order = 1)) +
#   labs(title = "B) Prevalence of moderate-or-worse anemia")
# plot_modsev
# 
# # mean Hg
# plot_meanhg <- ggplot() +
#   geom_ribbon(data = variog_perm_results_df %>% filter(formula == "mean_hg~1"), 
#               # originally temp_permute
#               aes(x = dist, ymin = q0.025, ymax = q0.975, fill = fill_flag),
#               alpha = 0.3) +
#   geom_point(data = variog_expvar_df %>% filter(formula == "mean_hg~1"), 
#              aes(x = dist, y = gamma)) + # replaced temp_expvar with variog_expvar_df
#   geom_text(data = variog_expvar_df %>% filter(formula == "mean_hg~1"), 
#             aes(x = dist, y = 0.06*0.95, label = np),
#             size = 2) + # empirical variog, does not depend on model
#   geom_line(data = variog_line_df %>% filter(formula == "mean_hg~1" &model=="Sph"), 
#             aes(x = dist, y = gamma, color = model)) +
#   geom_vline(data = variog_stats_df %>% filter(formula == "mean_hg~1"&model=="Sph"), 
#              aes(xintercept = effrange, color = model),
#              lty = "dashed", show.legend = FALSE) +
#   scale_color_manual(values = variog_colors,
#                      labels = c("Sph" = "Spherical", "Gau" = "Gaussian")) +
#   scale_fill_manual(values = c("1" = "grey"), labels = c("1" = "Monte Carlo\nenvelope")) +
#   labs(x = "Distance (km)",
#        y = "Semivariance",
#        color = "Fitted model",
#        fill = "") +
#   coord_cartesian(x = c(0, max((variog_expvar_df%>% filter(formula == "mean_hg~1")) %>% 
#                                  pull(dist))),
#                   y = c(0, 1)) +
#   # max((variog_expvar_df%>% filter(formula == "mean_hg~1")) %>% pull(gamma))
#   theme_classic() +
#   theme(title = element_text(size = 9),
#         panel.grid.major = element_line(color = "grey95"),
#   ) + # legend.position = "none"
#   guides(fill  = guide_legend(order = 2),
#          color = guide_legend(order = 1)) +
#   labs(title = "C) Mean Hemoglobin")
# plot_meanhg
# 
# # median Hg
# plot_medianhg <- ggplot() +
#   geom_ribbon(data = variog_perm_results_df %>% filter(formula == "median_hg~1"), 
#               # originally temp_permute
#               aes(x = dist, ymin = q0.025, ymax = q0.975, fill = fill_flag),
#               alpha = 0.3) +
#   geom_point(data = variog_expvar_df %>% filter(formula == "median_hg~1"), 
#              aes(x = dist, y = gamma)) + # replaced temp_expvar with variog_expvar_df
#   geom_text(data = variog_expvar_df %>% filter(formula == "median_hg~1"), 
#             aes(x = dist, y = 0.06*0.95, label = np),
#             size = 2) + # empirical variog, does not depend on model
#   geom_line(data = variog_line_df %>% filter(formula == "median_hg~1"&model=="Sph"), 
#             aes(x = dist, y = gamma, color = model)) +
#   geom_vline(data = variog_stats_df %>% filter(formula == "median_hg~1" &model=="Sph"), 
#              aes(xintercept = effrange, color = model),
#              lty = "dashed", show.legend = FALSE) +
#   scale_color_manual(values = variog_colors,
#                      labels = c("Sph" = "Spherical", "Gau" = "Gaussian")) +
#   scale_fill_manual(values = c("1" = "grey"), labels = c("1" = "Monte Carlo\nenvelope")) +
#   labs(x = "Distance (km)",
#        y = "Semivariance",
#        color = "Fitted model",
#        fill = "") +
#   coord_cartesian(x = c(0, max((variog_expvar_df%>% filter(formula == "median_hg~1")) %>% 
#                                  pull(dist))),
#                   y = c(0, 1)) +
#   #max((variog_expvar_df%>% filter(formula == "median_hg~1")) %>% pull(gamma))
#   theme_classic() +
#   theme(title = element_text(size = 9),
#         panel.grid.major = element_line(color = "grey95"),
#   ) + # legend.position = "none"
#   guides(fill  = guide_legend(order = 2),
#          color = guide_legend(order = 1)) + 
#   labs(title = "D) Median Hemoglobin")
# plot_medianhg
# 
# # putting them all together
# semivar_plots <- ggpubr::ggarrange(plot_any, plot_modsev, plot_meanhg, plot_medianhg,
#                                    nrow = 1, common.legend = T, legend = "right")
# semivar_plots
# 
# # saving
# # setwd("~/Library/Mobile Documents/com~apple~CloudDocs/R projects/aa-anemia")
# # ggsave(here("figures","village_semivariogram.eps"), semivar_plots, device = cairo_ps, # this allows me to save eps with alpha
# #        height = 6, width = 15)
# 
# # now getting the nug, range, and sill values for each semivariogram
# # trying to see what model the autofitVariogram function chooses
# automap::autofitVariogram(formula = any.logit ~ 1,
#                           input_data = data %>% 
#                             st_as_sf(coords = c("lon","lat"),
#                                      remove = F,
#                                      crs = crs) %>% as("Spatial"),
#                           miscFitOptions = list(merge.small.bins = T,
#                                                 min.np.bin = 10)) %>% plot(.)
# # Gau
# 
# automap::autofitVariogram(formula = modsev.logit ~ 1,
#                           input_data = data %>% 
#                             st_as_sf(coords = c("lon","lat"),
#                                      remove = F,
#                                      crs = crs) %>% as("Spatial"),
#                           miscFitOptions = list(merge.small.bins = T,
#                                                 min.np.bin = 10)) %>% plot(.)
# # Sph -- no spatial structure to the data
# 
# automap::autofitVariogram(formula = mean_hg ~ 1,
#                           input_data = data %>% 
#                             st_as_sf(coords = c("lon","lat"),
#                                      remove = F,
#                                      crs = crs) %>% as("Spatial"),
#                           miscFitOptions = list(merge.small.bins = T,
#                                                 min.np.bin = 10),
#                           model = "Sph") %>% plot(.)
# # Sph
# 
# automap::autofitVariogram(formula = median_hg ~ 1,
#                           input_data = data %>% 
#                             st_as_sf(coords = c("lon","lat"),
#                                      remove = F,
#                                      crs = crs) %>% as("Spatial"),
#                           miscFitOptions = list(merge.small.bins = T,
#                                                 min.np.bin = 10),
#                           model = "Sph") %>% plot(.)
# # Sph