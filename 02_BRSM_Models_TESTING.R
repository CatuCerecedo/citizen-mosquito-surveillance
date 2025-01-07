############################## BG Model ########################################

# if (!require("tidyverse")) install.packages("tidyverse")
# if (!require("rstanarm")) install.packages("rstanarm")
# if (!require("brms")) install.packages("brms")
# if (!require("cmdstanr")) install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# # WARNING: installing cmdstanr in windows 10 requires patience, surely you will
# # search for help
# if (!require("shinystan")) install.packages("shinystan")
# if (!require("loo")) install.packages("loo")

library(tidyverse)
library(rstanarm)
library(brms)
library(shinystan)
library(loo)
library(cmdstanr)
library(sf)

check_cmdstan_toolchain()
cmdstan_version()

# remotes::install_github("stan-dev/cmdstanr")

rm(list = ls())
# Directories ------------------------------------------------------------------

# In local
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")

# In cluster
loc.output <- paste0(getwd(), "/Spain_tiger/OUTPUT/")
loc.data <- paste0(getwd(), "/Spain_tiger/DATA/")

# Loading ----------------------------------------------------------------------
tiger <- readRDS(paste0(loc.output, "bg_tiger_spain_daily.rds"))
# tiger$year <- as.factor(lubridate::year(tiger$end_date))

tiger <- tiger %>% 
  filter(trapping_effort != 0) %>%
  mutate(
    y = as.factor(year(end_date)),
    # urban = industrial_transport + urban_fabric,
    occu = case_when(
      females > 0 ~ 1,
      TRUE ~ 0),
    FH = case_when(mean_relative_humidity < 40~0, mean_relative_humidity >95~0, 
                   (mean_relative_humidity >=40 & mean_relative_humidity <= 95)~((mean_relative_humidity/55)-(40/55))),
    FT = case_when(mean_temperature<=15~0, mean_temperature>30~0, (mean_temperature>15 & mean_temperature <=20)~(.2*mean_temperature)-3, (mean_temperature>20 & mean_temperature<=25)~1, (mean_temperature>25 & mean_temperature <= 30)~(-.2*mean_temperature)+6),
    mwi = FH*FT,
    urban = other_artificial + roads_rails + sports_leisure + discont_urban_fabric
  ) %>%
  filter(mean_relative_humidity < 101) # There is a wrong data --> ERA5 has no sense

# print(names(tiger))

# Filtering by municpilaties with more than 2 samples:
mun_selected <- tiger %>% group_by(id) %>% summarize(n = n()) %>% filter(n > 2)
length(unique(mun_selected$id))

tiger <- tiger %>% filter(id %in% mun_selected$id)

# Modeling process -------------------------------------------------------------
nchains = 4
threads_per_chain = 1

iteret = 400
wup = 200

# NOT RUN
# Only when you are working in local
# Sys.setenv(PATH = "C:/Users/Catu/Documents/.cmdstan/cmdstan-2.31.0/stan/lib/stan_math/lib/tbb")

mtiger1 <- brm(females ~ poly(mean_temperature, 2) +  precipitation +
                 log(trapping_effort) + log(n_traps) +
                 # offset(log(trapping_effort)) + offset(log(n_traps)) +
                 (1 | prov_name)  + (1 | y),
  data = tiger,
  prior = set_prior("cauchy(0,2.5)", class="b"),
  family = negbinomial(link = "log"),
  iter = iteret,
  chains = nchains,
  cores = nchains,
  backend = "cmdstanr",
  threads = threading(threads_per_chain),
  control = list(adapt_delta = 0.99))

mtiger2 <- brm(females ~ poly(mean_temperature, 2) + poly(precipitation, 2) +
                 log(trapping_effort) + log(n_traps) +
                 # offset(log(trapping_effort)) + offset(log(n_traps)) +
                 (1 | prov_name)  + (1 | y),
               data = tiger,
               prior = set_prior("cauchy(0,2.5)", class="b"),
               family = negbinomial(link = "log"),
               iter = iteret,
               chains = nchains,
               cores = nchains,
               backend = "cmdstanr",
               threads = threading(threads_per_chain),
               control = list(adapt_delta = 0.99))

mtiger3 <- brm(females ~ poly(mean_temperature, 2) + poly(precipitation, 2) +
                 wind_speed +
                 log(trapping_effort) + log(n_traps) +
                 # offset(log(trapping_effort)) + offset(log(n_traps)) +
                 (1 | prov_name) + (1 | y),
                 data = tiger,
                 prior = set_prior("cauchy(0,2.5)", class="b"),
                 family = negbinomial(link = "log"),
                 iter = iteret,
                 chains = nchains,
                 cores = nchains,
                 backend = "cmdstanr",
                 threads = threading(threads_per_chain),
                 control = list(adapt_delta = 0.99))

mtiger4 <- brm(females ~ poly(mean_temperature, 2) + poly(precipitation, 2) +
                 discont_urban_fabric + forests_scrub +
                 log(trapping_effort) + log(n_traps) +
                 # offset(log(trapping_effort)) + offset(log(n_traps)) +
                 (1 | prov_name) + (1 | y),
               data = tiger,
               prior = set_prior("cauchy(0,2.5)", class="b"),
               family = negbinomial(link = "log"),
               iter = iteret,
               chains = nchains,
               cores = nchains,
               backend = "cmdstanr",
               threads = threading(threads_per_chain),
               control = list(adapt_delta = 0.99))

mtiger5 <-  brm(females ~ poly(mean_temperature, 2) + poly(precipitation, 2) +
                  green_urban +
                  log(trapping_effort) + log(n_traps) +
                  # offset(log(trapping_effort)) + offset(log(n_traps)) +
                  (1 | prov_name)  + (1 | y),
               data = tiger,
               prior = set_prior("cauchy(0,2.5)", class="b"),
               family = negbinomial(link = "log"),
               iter = iteret,
               chains = nchains,
               cores = nchains,
               backend = "cmdstanr",
               threads = threading(threads_per_chain),
               control = list(adapt_delta = 0.99))

mtiger6 <-  brm(females ~ poly(mean_temperature, 2) + poly(precipitation, 2) +
                  mean_relative_humidity +
                  log(trapping_effort) + log(n_traps) +
                  # offset(log(trapping_effort)) + offset(log(n_traps)) +
                  (1 | prov_name)  + (1 | y),
               data = tiger,
               prior = set_prior("cauchy(0,2.5)", class="b"),
               family = negbinomial(link = "log"),
               iter = iteret,
               chains = nchains,
               cores = nchains,
               backend = "cmdstanr",
               threads = threading(threads_per_chain),
               control = list(adapt_delta = 0.99))

mtiger7 <-  brm(females ~ poly(mean_temperature, 2) + poly(precipitation, 2) +
                  green_urban + agricultural +
                  log(trapping_effort) + log(n_traps) +
                  # offset(log(trapping_effort)) + offset(log(n_traps)) +
                  (1 | prov_name)  + (1 | y),
               data = tiger,
               prior = set_prior("cauchy(0,2.5)", class="b"),
               family = negbinomial(link = "log"),
               iter = iteret,
               chains = nchains,
               cores = nchains,
               backend = "cmdstanr",
               threads = threading(threads_per_chain),
               control = list(adapt_delta = 0.99))


mtiger8 <-  brm(females ~ poly(mean_temperature, 2) +  poly(precipitation, 2) +
                  green_urban + agricultural +
                  log(trapping_effort) + log(n_traps) +
                  # offset(log(trapping_effort)) + offset(log(n_traps)) +
                  (1 | prov_name),
               data = tiger,
               prior = set_prior("cauchy(0,2.5)", class="b"),
               family = negbinomial(link = "log"),
               iter = iteret,
               chains = nchains,
               cores = nchains,
               backend = "cmdstanr",
               threads = threading(threads_per_chain),
               control = list(adapt_delta = 0.99))

mtiger9 <- brm(females ~ poly(mean_temperature, 2) +  poly(precipitation, 2) +
                 discont_urban_fabric +
                 log(trapping_effort) + log(n_traps) +
                 # offset(log(trapping_effort)) + offset(log(n_traps)) +
                 (1 | prov_name),
               data = tiger,
               prior = set_prior("cauchy(0,2.5)", class="b"),
               family = negbinomial(link = "log"),
               iter = iteret,
               chains = nchains,
               cores = nchains,
               backend = "cmdstanr",
               threads = threading(threads_per_chain),
               control = list(adapt_delta = 0.99))

mtiger10 <- brm(females ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                    agricultural + cont_urban_fabric +
                    offset(log(trapping_effort)) + offset(log(n_traps)) +
                  (1 | prov_name) +  (1 | y),
                  data = tiger,
                  prior = set_prior("cauchy(0,2.5)", class="b"),
                  family = negbinomial(link = "log"),
                  iter = iteret,
                  warmup = wup,
                  chains = nchains,
                  cores = nchains,
                  backend = "cmdstanr",
                  threads = threading(threads_per_chain),
                  control = list(adapt_delta = 0.99))

mtiger11 <- brm(females ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                  forests_scrub + cont_urban_fabric + discont_urban_fabric +
                  offset(log(trapping_effort)) + offset(log(n_traps)) +
                  (1 | prov_name) +  (1 | y),
                data = tiger,
                prior = set_prior("cauchy(0,2.5)", class="b"),
                family = negbinomial(link = "log"),
                iter = iteret,
                warmup = wup,
                chains = nchains,
                cores = nchains,
                backend = "cmdstanr",
                threads = threading(threads_per_chain),
                control = list(adapt_delta = 0.99))

mtiger12 <- brm(females ~ poly(min_temperature, 2) + precipitation + 
                  discont_urban_fabric  + cont_urban_fabric + agricultural +
                  forests_scrub + green_urban + inland_water + inland_wetlands +
                  marine_water + marine_wetlands + open + other_artificial +
                  roads_rails + sports_leisure +
                  log(trapping_effort) + log(n_traps) +
                  (1 | id) +  (1 | y),
                data = tiger,
                prior = set_prior("cauchy(0,2.5)", class="b"),
                family = negbinomial(link = "log"),
                iter = iteret,
                warmup = wup,
                chains = nchains,
                cores = nchains,
                backend = "cmdstanr",
                threads = threading(threads_per_chain),
                control = list(adapt_delta = 0.99))

mtiger13 <- brm(females ~ poly(min_temperature, 2) + precipitation + 
                  cont_urban_fabric  +
                  log(trapping_effort) + log(n_traps) +
                  (1 | id) +  (1 | y),
                data = tiger,
                prior = set_prior("cauchy(0,2.5)", class="b"),
                family = negbinomial(link = "log"),
                iter = iteret,
                warmup = wup,
                chains = nchains,
                cores = nchains,
                backend = "cmdstanr",
                threads = threading(threads_per_chain),
                control = list(adapt_delta = 0.99))

mtiger14 <- brm(females ~ poly(mean_temperature, 2) + precipitation + 
                  other_artificial + sports_leisure + roads_rails +
                  offset(log(trapping_effort)) + offset(log(n_traps)) +
                  (1 | id) + (1 | y),
                data = tiger,
                prior = set_prior("cauchy(0,2.5)", class="b"),
                family = negbinomial(link = "log"),
                iter = iteret,
                chains = nchains,
                cores = nchains,
                backend = "cmdstanr",
                threads = threading(threads_per_chain),
                control = list(adapt_delta = 0.99))

mtiger16b <- brm(females ~ poly(mean_temperature, 2) + precipitation + log(pop) + 
                  offset(log(trapping_effort)) +
                  (1 | id) + (1 | y),
                data = tiger,
                prior = set_prior("cauchy(0,2.5)", class="b"),
                family = negbinomial(link = "log"),
                iter = iteret,
                chains = nchains,
                cores = nchains,
                backend = "cmdstanr",
                save_pars = save_pars(all = TRUE),
                threads = threading(threads_per_chain),
                control = list(adapt_delta = 0.999))

mtiger17 <- brm(females ~ poly(min_temperature, 2) + precipitation + 
                  offset(log(trapping_effort)) +
                  (1 | id) + (1 | y),
                data = tiger,
                prior = set_prior("cauchy(0,2.5)", class="b"),
                family = negbinomial(link = "log"),
                iter = iteret,
                chains = nchains,
                cores = nchains,
                backend = "cmdstanr",
                save_pars = save_pars(all = TRUE),
                threads = threading(threads_per_chain),
                control = list(adapt_delta = 0.999))

mtiger1_occu <- brm(occu ~  poly(min_temperature, 2) + precipitation + 
                      cont_urban_fabric + discont_urban_fabric + 
                      log(trapping_effort) + log(n_traps) +
                      (1 | prov_name) + (1 | y),
                  data = tiger,
                  prior = set_prior("cauchy(0,2.5)", class="b"),
                  family = bernoulli(link = "logit"),
                  iter = iteret,
                  warmup = wup,
                  chains = nchains,
                  cores = nchains,
                  backend = "cmdstanr",
                  threads = threading(threads_per_chain),
                  control = list(adapt_delta = 0.99))

mtiger2_occu <- brm(occu ~  poly(mean_temperature, 2) + precipitation + 
                      cont_urban_fabric + discont_urban_fabric + 
                      log(trapping_effort) + log(n_traps) +
                      (1 | prov_name) + (1 | y),
                    data = tiger,
                    prior = set_prior("cauchy(0,2.5)", class="b"),
                    family = bernoulli(link = "logit"),
                    iter = iteret,
                    warmup = wup,
                    chains = nchains,
                    cores = nchains,
                    backend = "cmdstanr",
                    threads = threading(threads_per_chain),
                    control = list(adapt_delta = 0.99))

mtiger3_occu <- brm(occu ~  poly(mean_temperature, 2) + poly(precipitation, 2) + 
                      offset(log(trapping_effort)) + offset(log(n_traps)) +
                      (1 | prov_name)  + (1 | y),
                    data = tiger,
                    prior = set_prior("cauchy(0,2.5)", class="b"),
                    family = bernoulli(link = "logit"),
                    iter = iteret,
                    warmup = wup,
                    chains = nchains,
                    cores = nchains,
                    backend = "cmdstanr",
                    threads = threading(threads_per_chain),
                    control = list(adapt_delta = 0.99))

mtiger4_occu <- brm(occu ~  poly(mean_temperature, 2) + poly(precipitation, 2) + 
                      wind_speed +
                      offset(log(trapping_effort)) + 
                      (1 | prov_name) + (1 | y),
                    data = tiger,
                    prior = set_prior("cauchy(0,2.5)", class="b"),
                    family = bernoulli(link = "logit"),
                    iter = iteret,
                    warmup = wup,
                    chains = nchains,
                    cores = nchains,
                    backend = "cmdstanr",
                    threads = threading(threads_per_chain),
                    control = list(adapt_delta = 0.99))

mtiger5_occu <- brm(occu ~  poly(mean_temperature, 2) + poly(precipitation, 2) + 
                      wind_speed +
                      green_urban +
                      offset(log(trapping_effort)) + 
                      (1 | prov_name) +(1 | y),
                    data = tiger,
                    prior = set_prior("cauchy(0,2.5)", class="b"),
                    family = bernoulli(link = "logit"),
                    iter = iteret,
                    warmup = wup,
                    chains = nchains,
                    cores = nchains,
                    backend = "cmdstanr",
                    threads = threading(threads_per_chain),
                    control = list(adapt_delta = 0.99))

mtiger6_occu <- brm(occu ~  poly(mean_temperature, 2) + poly(precipitation, 2) + 
                      wind_speed +
                      cont_urban_fabric + 
                      offset(log(trapping_effort)) +
                      (1 | prov_name) + (1 | y),
                    data = tiger,
                    prior = set_prior("cauchy(0,2.5)", class="b"),
                    family = bernoulli(link = "logit"),
                    iter = iteret,
                    warmup = wup,
                    chains = nchains,
                    cores = nchains,
                    backend = "cmdstanr",
                    threads = threading(threads_per_chain),
                    control = list(adapt_delta = 0.99))

mtiger7_occu <- brm(occu ~  poly(mean_temperature, 2) + poly(precipitation, 2) + 
                      wind_speed +
                      discont_urban_fabric +  +
                      offset(log(trapping_effort)) + 
                      (1 | prov_name) + (1 | y),
                    data = tiger,
                    prior = set_prior("cauchy(0,2.5)", class="b"),
                    family = bernoulli(link = "logit"),
                    iter = iteret,
                    warmup = wup,
                    chains = nchains,
                    cores = nchains,
                    backend = "cmdstanr",
                    threads = threading(threads_per_chain),
                    control = list(adapt_delta = 0.99))

mtiger8_occu <- brm(occu ~  poly(mean_temperature, 2) + poly(precipitation, 2) + 
                      wind_speed +
                      agricultural + urban_fabric +
                      offset(log(trapping_effort)) + 
                      (1 | prov_name),
                    data = tiger,
                    prior = set_prior("cauchy(0,2.5)", class="b"),
                    family = bernoulli(link = "logit"),
                    iter = iteret,
                    warmup = wup,
                    chains = nchains,
                    cores = nchains,
                    backend = "cmdstanr",
                    threads = threading(threads_per_chain),
                    control = list(adapt_delta = 0.99))

mtiger9_occu <- brm(occu ~ mean_temperature + precipitation + 
                      wind_speed +
                      agricultural + urban_fabric +
                      offset(log(trapping_effort)) + 
                      (1 | prov_name),
                    data = tiger,
                    prior = set_prior("cauchy(0,2.5)", class="b"),
                    family = bernoulli(link = "logit"),
                    iter = iteret,
                    warmup = wup,
                    chains = nchains,
                    cores = nchains,
                    backend = "cmdstanr",
                    threads = threading(threads_per_chain),
                    control = list(adapt_delta = 0.99))

mtiger10_occu <- brm(occu ~ poly(mean_temperature, 2) +
                      wind_speed + urban_fabric +
                      offset(log(trapping_effort)) + 
                      (1 | prov_name),
                    data = tiger,
                    prior = set_prior("cauchy(0,2.5)", class="b"),
                    family = bernoulli(link = "logit"),
                    iter = iteret,
                    warmup = wup,
                    chains = nchains,
                    cores = nchains,
                    backend = "cmdstanr",
                    threads = threading(threads_per_chain),
                    control = list(adapt_delta = 0.99))

# Mosquito Alert models --------------------------------------------------------
ma_df <- readRDS(paste0(loc.output, "ma_tiger_spain_daily.rds")) %>%
  mutate(
    id = as.factor(id),
    y = as.factor(year(date))
  ) %>%
  filter(SE > 0) %>%
  drop_na()

# Filtering by municpilaties with more than 2 samples:
mun_selected <- ma_df %>% group_by(id) %>% summarize(n = n()) %>% filter(n > 3)
length(unique(mun_selected$id))

ma_df <- ma_df %>% filter(id %in% mun_selected$id)

# these_priors = c(prior(student_t(4, 0, 6), class = "b"))
mtiger1_ma <- brm(any_reps ~ poly(mean_temperature, 2) + mean_relative_humidity + dens +
                    (1 | id) + (1 | y) + (1 | app_phase) +
                    offset(log(SE)),
                     data = ma_df,
                     prior = set_prior("cauchy(0,2.5)", class="b"),
                     family = bernoulli(link = "logit"),
                     iter = iteret,
                     warmup = wup,
                     chains = nchains,
                     cores = nchains,
                     backend = "cmdstanr",
                     threads = threading(threads_per_chain),
                     control = list(adapt_delta = 0.99))

mtiger1a_ma <- brm(any_reps ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                    (1 | id) + (1 | y) + (1 | app_phase) +
                    offset(log(SE)),
                  data = ma_df,
                  prior = set_prior("cauchy(0,2.5)", class="b"),
                  family = bernoulli(link = "logit"),
                  iter = iteret,
                  warmup = wup,
                  chains = nchains,
                  cores = nchains,
                  backend = "cmdstanr",
                  threads = threading(threads_per_chain),
                  control = list(adapt_delta = 0.99))

mtiger2_ma <- brm(any_reps ~ poly(mean_temperature, 2) + 
                    (1 | id) + (1 | y) + (1 | app_phase) +
                    offset(log(SE)),
                  data = ma_df,
                  prior = set_prior("cauchy(0,2.5)", class="b"),
                  family = bernoulli(link = "logit"),
                  iter = iteret,
                  warmup = wup,
                  chains = nchains,
                  cores = nchains,
                  backend = "cmdstanr",
                  threads = threading(threads_per_chain),
                  control = list(adapt_delta = 0.99))

mtiger3_ma <- brm(any_reps ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                    precipitation +
                    (1 | id) + (1 | y) + (1 | app_phase) +
                    offset(log(SE)),
                  data = ma_df,
                  prior = set_prior("cauchy(0,2.5)", class="b"),
                  family = bernoulli(link = "logit"),
                  iter = iteret,
                  warmup = wup,
                  chains = nchains,
                  cores = nchains,
                  backend = "cmdstanr",
                  threads = threading(threads_per_chain),
                  control = list(adapt_delta = 0.99))

mtiger4_ma <-  brm(any_reps ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                     wind_speed +
                     (1 | id) + (1 | y) + (1 | app_phase) +
                     offset(log(SE)),
                  data = ma_df,
                  prior = set_prior("cauchy(0,2.5)", class="b"),
                  family = bernoulli(link = "logit"),
                  iter = iteret,
                  warmup = wup,
                  chains = nchains,
                  cores = nchains,
                  backend = "cmdstanr",
                  threads = threading(threads_per_chain),
                  control = list(adapt_delta = 0.99))

mtiger5_ma <- brm(any_reps ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                    agricultural + artificial_green_urban + forests + urban_fabric +
                    (1 | id) + (1 | y) + (1 | app_phase) +
                    offset(log(SE)),
                  data = ma_df,
                  prior = set_prior("cauchy(0,2.5)", class="b"),
                  family = bernoulli(link = "logit"),
                  iter = iteret,
                  warmup = wup,
                  chains = nchains,
                  cores = nchains,
                  backend = "cmdstanr",
                  threads = threading(threads_per_chain),
                  control = list(adapt_delta = 0.99))

mtiger6_ma <- brm(any_reps ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                    urban_fabric + industrial_transport + 
                    (1 | id) + (1 | y) + (1 | app_phase) +
                    offset(log(SE)),
                  data = ma_df,
                  prior = set_prior("cauchy(0,2.5)", class="b"),
                  family = bernoulli(link = "logit"),
                  iter = iteret,
                  warmup = wup,
                  chains = nchains,
                  cores = nchains,
                  backend = "cmdstanr",
                  threads = threading(threads_per_chain),
                  control = list(adapt_delta = 0.99))

mtiger7_ma <- brm(any_reps ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                    urban_fabric +
                    (1 | id) + (1 | y) + (1 | app_phase) +
                    offset(log(SE)),
                  data = ma_df,
                  prior = set_prior("cauchy(0,2.5)", class="b"),
                  family = bernoulli(link = "logit"),
                  iter = iteret,
                  warmup = wup,
                  chains = nchains,
                  cores = nchains,
                  backend = "cmdstanr",
                  threads = threading(threads_per_chain),
                  control = list(adapt_delta = 0.99))

mtiger8_ma <- brm(any_reps ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                    urban +
                    (1 | id) + (1 | y) + (1 | app_phase) +
                    offset(log(SE)),
                  data = ma_df,
                  prior = set_prior("cauchy(0,2.5)", class="b"),
                  family = bernoulli(link = "logit"),
                  iter = iteret,
                  warmup = wup,
                  chains = nchains,
                  cores = nchains,
                  backend = "cmdstanr",
                  threads = threading(threads_per_chain),
                  control = list(adapt_delta = 0.99))

mtiger9_ma <- brm(any_reps ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                    agricultural + cont_urban_fabric + discont_urban_fabric +
                    (1 | id) + (1 | y) + (1 | app_phase) +
                    offset(log(SE)),
                  data = ma_df,
                  prior = set_prior("cauchy(0,2.5)", class="b"),
                  family = bernoulli(link = "logit"),
                  iter = iteret,
                  warmup = wup,
                  chains = nchains,
                  cores = nchains,
                  backend = "cmdstanr",
                  threads = threading(threads_per_chain),
                  control = list(adapt_delta = 0.99))

mtiger11_ma <- brm(any_reps ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                     urban_fabric + inland_water +
                     (1 | prov_name) + 
                     offset(log(SE)),
                   data = ma_df,
                   prior = set_prior("cauchy(0,2.5)", class="b"),
                   family = bernoulli(link = "logit"),
                   iter = iteret,
                   warmup = wup,
                   chains = nchains,
                   cores = nchains,
                   backend = "cmdstanr",
                   threads = threading(threads_per_chain),
                   control = list(adapt_delta = 0.99))

mtiger12_ma <- brm(any_reps ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                     urban_fabric + inland_water + forests +
                     (1 | prov_name) + 
                     offset(log(SE)),
                   data = ma_df,
                   prior = set_prior("cauchy(0,2.5)", class="b"),
                   family = bernoulli(link = "logit"),
                   iter = iteret,
                   warmup = wup,
                   chains = nchains,
                   cores = nchains,
                   backend = "cmdstanr",
                   threads = threading(threads_per_chain),
                   control = list(adapt_delta = 0.99))

# Models using number of report as SE (avoiding SE calculated from historical data)
ma_df <- ma_df %>% 
  filter(n_total_reports > 0) # Si no hay esfuerzo de muestreo el dato no debería de estar.

mtiger1_ma_rep <- brm(any_reps ~ poly(mean_temperature, 2) + 
                        precipitation + 
                        (1 | pixel_id) + (1 | year) + offset(log(n_total_reports)),
                      data = ma_df,
                      prior = set_prior("cauchy(0,2.5)", class="b"),
                      family = bernoulli(link = "logit"),
                      iter = iteret,
                      warmup = wup,
                      chains = nchains,
                      cores = nchains,
                      backend = "cmdstanr",
                      threads = threading(threads_per_chain),
                      control = list(adapt_delta = 0.99))
# saveRDS(mtiger1_ma_rep, file = paste0(loc.output, "mtiger1_ma_rep.rds"))

mtiger2_ma_rep <- brm(any_reps ~ poly(mean_temperature, 2) + 
                        mean_relative_humidity + 
                        (1 | pixel_id) + + (1 | year) offset(log(n_total_reporters)),
                      data = ma_df,
                      prior = set_prior("cauchy(0,2.5)", class="b"),
                      family = bernoulli(link = "logit"),
                      iter = iteret,
                      warmup = wup,
                      chains = nchains,
                      cores = nchains,
                      backend = "cmdstanr",
                      threads = threading(threads_per_chain),
                      control = list(adapt_delta = 0.99))
# saveRDS(mtiger2_ma_rep, file = paste0(loc.output, "mtiger2_ma_rep.rds"))

mtiger3_ma_rep <- brm(any_reps ~ poly(mean_temperature, 2) + 
                        mean_relative_humidity + 
                        (1 | pixel_id) + (1 | year) + offset(log(n_total_reports/n_total_reporters)),
                      data = ma_df,
                      prior = set_prior("cauchy(0,2.5)", class="b"),
                      family = bernoulli(link = "logit"),
                      iter = iteret,
                      warmup = wup,
                      chains = nchains,
                      cores = nchains,
                      backend = "cmdstanr",
                      threads = threading(threads_per_chain),
                      control = list(adapt_delta = 0.99))
# saveRDS(mtiger3_ma_rep, file = paste0(loc.output, "mtiger3_ma_rep.rds"))

mtiger4_ma_rep <- brm(any_reps ~ poly(mean_temperature, 2) + 
                        mean_relative_humidity + 
                        (1 | pixel_id) + (1 | year) + 
                        offset(log(n_total_reports)) + offset(log(n_total_reporters)),
                      data = ma_df,
                      prior = set_prior("cauchy(0,2.5)", class="b"),
                      family = bernoulli(link = "logit"),
                      iter = iteret,
                      warmup = wup,
                      chains = nchains,
                      cores = nchains,
                      backend = "cmdstanr",
                      threads = threading(threads_per_chain),
                      control = list(adapt_delta = 0.99))
# saveRDS(mtiger4_ma_rep, file = paste0(loc.output, "mtiger4_ma_rep.rds"))

# Bites models -----------------------------------------------------------------
bites <- readRDS(file = paste0(loc.output, "bites_spain_daily.rds"))
nchains = 4
threads_per_chain = 1

iteret = 400
wup = 200

# Following the variables selected by dredge and glmmTMB
mbites_1 <- brm(any_reps ~ l21mean_relative_humidity + l21min_temperature + l21precipitation +
                  agricultural + inland_wetlands + other_artificial + marine_water +
                  inland_water + inland_water + cont_urban_fabric + sports_leisure +
                        (1 | id) + (1 | y) + 
                        offset(log(SE)),
                      data = bites,
                      prior = set_prior("cauchy(0,2.5)", class="b"),
                      family = bernoulli(link = "logit"),
                      iter = iteret,
                      warmup = wup,
                      chains = nchains,
                      cores = nchains,
                      backend = "cmdstanr",
                      threads = threading(threads_per_chain),
                      control = list(adapt_delta = 0.99))
summary(mbites_1)
saveRDS(mbites_1, file = paste0(loc.output, "testing_","mbites_1.rds"))

mbites_1 <- brm(any_reps ~ l21mean_relative_humidity + l21min_temperature + l21precipitation +
                  agricultural + inland_wetlands + other_artificial + marine_water +
                  inland_water + inland_water + cont_urban_fabric + sports_leisure +
                  (1 | id) + (1 | y) + 
                  offset(log(SE)),
                data = bites,
                prior = set_prior("cauchy(0,2.5)", class="b"),
                family = bernoulli(link = "logit"),
                iter = iteret,
                warmup = wup,
                chains = nchains,
                cores = nchains,
                backend = "cmdstanr",
                threads = threading(threads_per_chain),
                control = list(adapt_delta = 0.99))
summary(mbites_1)
loo(mbites_1)
saveRDS(mbites_1, file = paste0(loc.output, "testing_","mbites_1.rds"))

mbites_2 <- brm(any_reps ~ poly(l21min_temperature, 2) +
                  agricultural + other_artificial + marine_water +
                  inland_water + inland_water + cont_urban_fabric + sports_leisure +
                  (1 | id) + (1 | y) + 
                  offset(log(SE)),
                data = bites,
                prior = set_prior("cauchy(0,2.5)", class="b"),
                family = bernoulli(link = "logit"),
                iter = iteret,
                warmup = wup,
                chains = nchains,
                cores = nchains,
                backend = "cmdstanr",
                threads = threading(threads_per_chain),
                control = list(adapt_delta = 0.99))
summary(mbites_2)
loo(mbites_2)
saveRDS(mbites_2, file = paste0(loc.output, "testing_","mbites_2.rds"))

mbites_3 <- brm(any_reps ~ poly(l21min_temperature, 2) +
                  agricultural + other_artificial +
                  (1 | id) + (1 | y) + 
                  offset(log(SE)),
                data = bites,
                prior = set_prior("cauchy(0,2.5)", class="b"),
                family = bernoulli(link = "logit"),
                iter = iteret,
                warmup = wup,
                chains = nchains,
                cores = nchains,
                backend = "cmdstanr",
                threads = threading(threads_per_chain),
                control = list(adapt_delta = 0.99))
summary(mbites_3)
loo(mbites_3)
saveRDS(mbites_3, file = paste0(loc.output, "testing_","mbites_3.rds"))

mbites_4 <- brm(any_reps ~ l21mean_relative_humidity + l21min_temperature + l21precipitation +
                  agricultural + other_artificial +
                  (1 | id) + (1 | y) + 
                  offset(log(SE)),
                data = bites,
                prior = set_prior("cauchy(0,2.5)", class="b"),
                family = bernoulli(link = "logit"),
                iter = iteret,
                warmup = wup,
                chains = nchains,
                cores = nchains,
                backend = "cmdstanr",
                threads = threading(threads_per_chain),
                control = list(adapt_delta = 0.99))
summary(mbites_4)
loo(mbites_4)
saveRDS(mbites_4, file = paste0(loc.output, "testing_","mbites_4.rds"))

mbites_5 <- brm(any_reps ~ l21mean_relative_humidity + poly(l21min_temperature, 2) + 
                  agricultural + other_artificial +
                  (1 | id) + (1 | y) + 
                  offset(log(SE)),
                data = bites,
                prior = set_prior("cauchy(0,2.5)", class="b"),
                family = bernoulli(link = "logit"),
                iter = iteret,
                warmup = wup,
                chains = nchains,
                cores = nchains,
                backend = "cmdstanr",
                threads = threading(threads_per_chain),
                control = list(adapt_delta = 0.99))
summary(mbites_5)
loo(mbites_5)
saveRDS(mbites_5, file = paste0(loc.output, "testing_","mbites_5.rds"))

mbites_6 <- brm(any_reps ~ poly(l21min_temperature, 2) + 
                  agricultural + other_artificial + green_urban + discont_urban_fabric +
                  (1 | id) + (1 | y) + 
                  offset(log(SE)),
                data = bites,
                prior = set_prior("cauchy(0,2.5)", class="b"),
                family = bernoulli(link = "logit"),
                iter = iteret,
                warmup = wup,
                chains = nchains,
                cores = nchains,
                backend = "cmdstanr",
                threads = threading(threads_per_chain),
                control = list(adapt_delta = 0.99))
summary(mbites_6)
loo(mbites_6)
saveRDS(mbites_6, file = paste0(loc.output, "testing_","mbites_6.rds"))

mbites_7 <- brm(any_reps ~ poly(l21min_temperature, 2) + 
                  agricultural + other_artificial + discont_urban_fabric +
                  (1 | id) + (1 | y) + 
                  offset(log(SE)),
                data = bites,
                prior = set_prior("cauchy(0,2.5)", class="b"),
                family = bernoulli(link = "logit"),
                iter = iteret,
                warmup = wup,
                chains = nchains,
                cores = nchains,
                backend = "cmdstanr",
                threads = threading(threads_per_chain),
                control = list(adapt_delta = 0.99))
summary(mbites_7)
loo(mbites_7)
saveRDS(mbites_7, file = paste0(loc.output, "testing_","mbites_7.rds"))
