# ALTO AMAZONAS ANEMIA SPATIAL ANALYSIS
# JOHN NESEMANN

#### packages ####
library(tidyverse)
library(haven)
library(foreign)
library(sf)
library(gridExtra)

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
  select(unique_id:hh_numb, pplinhh:gps_alt, dbs_hg, dbs_hg_na, anemia, anemia_level_sa, anemia_orig) %>%
  dplyr::rename(anemia.f=anemia_orig) %>%
  group_by(hh_numb) %>%
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
data.hhs

#### descriptive code ####
ggplot(data, aes(x=p_any)) + geom_histogram()
ggplot(data, aes(x=p_mild)) + geom_histogram()
ggplot(data, aes(x=p_modsev)) + geom_histogram() # none of these look normally distributed

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
data_spatial

# rule of thumb: limit max distance to hald of max interpoint distance
# estimate distance from each point to all other points
clu_distances <- sapply(1:nrow(data_spatial), function(x) geosphere::distGeo(p1 = data_spatial, p2 = data_spatial[x,]))
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
    st_as_sf(coords = c("lon", "lat"), remove = FALSE, crs = swift_crs) %>%
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

variog_line_df <- variog_results[1,] %>% bind_rows()
variog_stats_df <- variog_results[2,] %>% bind_rows()
variog_expvar_df <- variog_results[3,] %>% bind_rows()
variog_perm_results_df <- variog_perm_results[1,] %>% bind_rows()

## plot variograms
variog_colors <- c("Sph" = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[7],
                   "Exp" = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[4],
                   "Mat" = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[5])

## `get_variog_plot`: plots variograms across 4 survey time points
# variog_results
# permute_results
# returns figure containing 4 variogram plots, one for each survey point
get_variog_plot <- function(permute_results, variog_line, variog_stats, variog_expvar, input_title){
  
  # initialize container to save plots
  variog_plot_list <- list()
  
  #for(s in survey_list){
    
    temp_permute <- permute_results #%>% filter(survey == s)
    temp_line <- variog_line #%>% filter(survey == s)
    temp_stats <- variog_stats #%>% filter(survey == s)
    temp_expvar <- variog_expvar #%>% filter(survey == s)
    
    # create variogram plot of points and each model fit
    temp_plot <- ggplot() +
      geom_ribbon(data = temp_permute,
                  aes(x = dist, ymin = q0.025, ymax = q0.975, fill = fill_flag),
                  alpha = 0.3) +
      geom_point(data = temp_expvar, aes(x = dist, y = gamma)) + # empirical variog, does not depend on model
      geom_text(data = temp_expvar, aes(x = dist, y = 0.06*0.95, label = np),
                size = 2) + # empirical variog, does not depend on model
      geom_line(data = temp_line, aes(x = dist, y = gamma, color = model)) +
      geom_vline(data = temp_stats, aes(xintercept = effrange, color = model),
                 lty = "dashed", show.legend = FALSE) +
      scale_color_manual(values = variog_colors,
                         labels = c("Exp" = "Exponential", "Mat" = "Matern", "Sph" = "Spherical")) +
      scale_fill_manual(values = c("1" = "grey"), labels = c("1" = "Monte Carlo\nenvelope")) +
      labs(x = "Distance (km)",
           y = "Semivariance",
           color = "Fitted model",
           fill = "") +
      coord_cartesian(x = c(0, max(temp_expvar %>% pull(dist))),
                      y = c(0, 0.062)) +
      theme_classic() +
      theme(title = element_text(size = 10),
            panel.grid.major = element_line(color = "grey95"),
            legend.position = "none") +
      guides(fill  = guide_legend(order = 2),
             color = guide_legend(order = 1))
    
    # if(s != 0){
      temp_plot <- temp_plot +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank())
    # }
    
    # if(s == 0){
    #   temp_plot <- temp_plot +
    #     geom_text(aes(x = 8, y = 0.062), label = "Bin sample sizes", size = 3) +
    #     theme(plot.margin = unit(c(0.3,0,0.3,1), 'lines'))
    # }
    
    # if(s == 36){
    #   temp_plot <- temp_plot + 
    #     theme(legend.position = "right",
    #           legend.title = element_text(size = 10),
    #           legend.text = element_text(size = 8))
    # }
    
    variog_plot_list <- c(variog_plot_list, list(temp_plot))
  # }
  
  variog_plot <- arrangeGrob(grobs = variog_plot_list, nrow = 1, widths = c(1.2,1,1,1.6),
                             top = input_title)
  return(variog_plot)
}

variog_plots <- lapply(c("p_any~1", "p_mild~1", "p_modsev~1"),
                       function(x){
                         get_variog_plot(
                           permute_results = variog_perm_results_df %>% filter(formula == x),
                           variog_expvar = variog_expvar_df %>% filter(formula == x),
                           variog_line = variog_line_df %>% filter(formula == x),
                           variog_stats = variog_stats_df %>% filter(formula == x),
                           input_title = paste0("Variograms by study month: ", x)
                         )})

# not working just yet but will come back and try to plot these.


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
