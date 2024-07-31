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
tiger <- readRDS(paste0(loc.output, "bg_tiger_spain.rds"))
# tiger$year <- as.factor(lubridate::year(tiger$end_date))

tiger <- tiger %>% 
  mutate(
    urban = discont_urban_fabric + sports_leisure + green_urban + cont_urban_fabric,
    occu = case_when(
      females > 0 ~ 1,
      TRUE ~ 0),
    FH = case_when(mean_relative_humidity < 40~0, mean_relative_humidity >95~0, 
                   (mean_relative_humidity >=40 & mean_relative_humidity <= 95)~((mean_relative_humidity/55)-(40/55))),
    FT = case_when(mean_temperature<=15~0, mean_temperature>30~0, (mean_temperature>15 & mean_temperature <=20)~(.2*mean_temperature)-3, (mean_temperature>20 & mean_temperature<=25)~1, (mean_temperature>25 & mean_temperature <= 30)~(-.2*mean_temperature)+6),
    mwi = FH*FT,
  ) %>%
  filter(mean_relative_humidity < 101) # There is a wrong data --> ERA5 has no sense

# print(names(tiger))

# Modeling process -------------------------------------------------------------
nchains = 4
threads_per_chain = 1

iteret = 400
wup = 200

# NOT RUN
# Only when you are working in local
# Sys.setenv(PATH = "C:/Users/Catu/Documents/.cmdstan/cmdstan-2.31.0/stan/lib/stan_math/lib/tbb")

mtiger1 <- brm(females ~ poly(mean_temperature, 2) +
                 agricultural + forests_scrub +
                 offset(log(trapping_effort)) + offset(log(n_traps)) +
                 (1 | pixel_id) + (1 | y),
  data = tiger,
  prior = set_prior("cauchy(0,2.5)", class="b"),
  family = negbinomial(link = "log"),
  iter = iteret,
  chains = nchains,
  cores = nchains,
  backend = "cmdstanr",
  threads = threading(threads_per_chain),
  control = list(adapt_delta = 0.99))
# saveRDS(mtiger1, file = paste0(loc.output, "mtiger1.rds"))

mtiger2 <- brm(females ~ poly(mean_temperature, 2) + precipitation +
                 agricultural + forests_scrub + discont_urban_fabric + roads_rails +
                 offset(log(trapping_effort)) + offset(log(n_traps)) +
                 (1 | pixel_id) + (1 | y),
               data = tiger,
               prior = set_prior("cauchy(0,2.5)", class="b"),
               family = negbinomial(link = "log"),
               iter = iteret,
               chains = nchains,
               cores = nchains,
               backend = "cmdstanr",
               threads = threading(threads_per_chain),
               control = list(adapt_delta = 0.99))
# saveRDS(mtiger2, file = paste0(loc.output, "mtiger2.rds"))

mtiger3 <- brm(females ~ poly(mean_temperature, 2) + precipitation + wind_speed +
                 cont_urban_fabric +
                 log(trapping_effort) + log(n_traps) +
                 (1 | pixel_id) + (1 | y),
                 data = tiger,
                 prior = set_prior("cauchy(0,2.5)", class="b"),
                 family = negbinomial(link = "log"),
                 iter = iteret,
                 chains = nchains,
                 cores = nchains,
                 backend = "cmdstanr",
                 threads = threading(threads_per_chain),
                 control = list(adapt_delta = 0.99))
# saveRDS(mtiger3, file = paste0(loc.output, "mtiger3.rds"))

mtiger4 <- brm(females ~ poly(mean_temperature, 2) + precipitation + 
                cont_urban_fabric +
                log(trapping_effort) + log(n_traps) +
                (1 | pixel_id) + (1 | y),
               data = tiger,
               prior = set_prior("cauchy(0,2.5)", class="b"),
               family = negbinomial(link = "log"),
               iter = iteret,
               chains = nchains,
               cores = nchains,
               backend = "cmdstanr",
               threads = threading(threads_per_chain),
               control = list(adapt_delta = 0.99))
# saveRDS(mtiger4, file = paste0(loc.output, "mtiger4.rds"))

mtiger5 <- brm(females ~ poly(mean_temperature, 2) + precipitation + 
                 roads_rails + green_urban +
                log(trapping_effort) + log(n_traps) +
                (1 | pixel_id) + (1 | y),
               data = tiger,
               prior = set_prior("cauchy(0,2.5)", class="b"),
               family = negbinomial(link = "log"),
               iter = iteret,
               chains = nchains,
               cores = nchains,
               backend = "cmdstanr",
               threads = threading(threads_per_chain),
               control = list(adapt_delta = 0.99))
# saveRDS(mtiger5, file = paste0(loc.output, "mtiger5.rds"))

mtiger6 <- brm(females ~ discont_urban_fabric + agricultural + forests_scrub + open +
                 inland_wetlands + inland_water + other_artificial + sports_leisure +
                 marine_water + roads_rails + cont_urban_fabric + marine_wetlands + green_urban +
                 offset(log(trapping_effort)) + offset(log(n_traps)) +
                 (1 | pixel_id) + (1 | y),
               data = tiger,
               prior = set_prior("cauchy(0,2.5)", class="b"),
               family = negbinomial(link = "log"),
               iter = iteret,
               chains = nchains,
               cores = nchains,
               backend = "cmdstanr",
               threads = threading(threads_per_chain),
               control = list(adapt_delta = 0.99))
# saveRDS(mtiger6, file = paste0(loc.output, "mtiger6.rds"))

mtiger7 <- brm(females ~ poly(min_temperature, 2) +
                 precipitation +
                 urban +
                 log(trapping_effort) + log(n_traps) +
                 (1 | pixel_id) + (1 | y),
              data = tiger,
              prior = set_prior("cauchy(0,2.5)", class="b"),
              family = negbinomial(link = "log"),
              iter = iteret,
              chains = nchains,
              cores = nchains,
              backend = "cmdstanr",
              threads = threading(threads_per_chain),
              control = list(adapt_delta = 0.99))
# saveRDS(mtiger7, file = paste0(loc.output, "mtiger7.rds"))

mtiger8 <- brm(females ~ poly(mean_temperature, 2) + precipitation + 
                 cont_urban_fabric + discont_urban_fabric + sports_leisure + green_urban +
                 log(trapping_effort) + log(n_traps) +
                 (1 | pixel_id) + (1 | y),
               data = tiger,
               prior = set_prior("cauchy(0,2.5)", class="b"),
               family = negbinomial(link = "log"),
               iter = iteret,
               warmup = wup,
               chains = nchains,
               cores = nchains,
               backend = "cmdstanr",
               threads = threading(threads_per_chain),
               control = list(adapt_delta = 0.999))
# saveRDS(mtiger8, file = paste0(loc.output, "mtiger8.rds"))

mtiger9 <- brm(females ~ poly(mean_temperature, 2) + poly(precipitation, 2) + 
                cont_urban_fabric +
                log(trapping_effort) + log(n_traps) +
                (1 | pixel_id) + (1 | y),
               data = tiger,
               prior = set_prior("cauchy(0,2.5)", class="b"),
               family = negbinomial(link = "log"),
               iter = iteret,
               chains = nchains,
               cores = nchains,
               backend = "cmdstanr",
               threads = threading(threads_per_chain),
               control = list(adapt_delta = 0.99))
# saveRDS(mtiger9, file = paste0(loc.output, "mtiger9.rds"))

mtiger10 <- brm(females ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                    agricultural + cont_urban_fabric +
                     log(trapping_effort) + log(n_traps) +
                    (1 | id) + (1 | year),
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
# saveRDS(mtiger10, file = paste0(loc.output, "mtiger10.rds"))

mtiger11 <- brm(females ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                  agricultural + cont_urban_fabric +
                  offset(log(trapping_effort)) + offset(log(n_traps)) +
                  (1 | id) + (1 | year),
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
# saveRDS(mtiger11, file = paste0(loc.output, "mtiger11.rds"))

mtiger1_occu <- brm(occu ~ poly(mean_temperature, 2) + poly(precipitation, 2) + 
                      cont_urban_fabric +
                      log(trapping_effort) + log(n_traps) +
                      (1 | pixel_id) + (1 | y),
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
# saveRDS(mtiger1_occu, file = paste0(loc.output, "mtiger1_occu.rds"))

mtiger2_occu <- brm(occu ~ poly(min_temperature, 2) +
                      mean_relative_humidity + 
                      agricultural + forests_scrub + discont_urban_fabric +
                      offset(log(trapping_effort)) + offset(log(n_traps)) +
                      (1 | pixel_id) + (1 | year),
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
# saveRDS(mtiger2_occu, file = paste0(loc.output, "mtiger2_occu.rds"))

mtiger3_occu <- brm(occu ~ poly(min_temperature, 2) +
                      agricultural + forests_scrub + 
                      offset(log(trapping_effort)) + offset(log(n_traps)) +
                      (1 | pixel_id) + (1 | year),
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
# saveRDS(mtiger3_occu, file = paste0(loc.output, "mtiger3_occu.rds"))

mtiger4_occu <- brm(occu ~ poly(mean_temperature, 2) +
                      agricultural + forests_scrub + 
                      offset(log(trapping_effort)) + offset(log(n_traps)) +
                      (1 | pixel_id) + (1 | year),
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
# saveRDS(mtiger4_occu, file = paste0(loc.output, "mtiger4_occu.rds"))

mtiger5_occu <- brm(occu ~ poly(min_temperature, 2) + precipitation +
                      agricultural + forests_scrub + 
                      offset(log(trapping_effort)) + offset(log(n_traps)) +
                      (1 | pixel_id) + (1 | year),
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
# saveRDS(mtiger5_occu, file = paste0(loc.output, "mtiger5_occu.rds"))

mtiger6_occu <- brm(occu ~ poly(min_temperature, 2) + mean_relative_humidity +
                      agricultural + 
                      offset(log(trapping_effort)) + offset(log(n_traps)) +
                      (1 | pixel_id) + (1 | year),
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
# saveRDS(mtiger6_occu, file = paste0(loc.output, "mtiger6_occu.rds"))

mtiger7_occu <- brm(occu ~ poly(min_temperature, 2) + mean_relative_humidity +
                      agricultural + 
                      offset(log(trapping_effort)) + offset(log(n_traps)) +
                      (1 | pixel_id) + (1 | year),
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
# saveRDS(mtiger7_occu, file = paste0(loc.output, "mtiger7_occu.rds"))

mtiger8_occu <- brm(occu ~ poly(mean_temperature, 2) +
                      mean_relative_humidity + 
                      agricultural + forests_scrub + discont_urban_fabric +
                      log(trapping_effort) + log(n_traps) +
                      (1 | pixel_id) + (1 | year),
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
# saveRDS(mtiger8_occu, file = paste0(loc.output, "mtiger8_occu.rds"))

mtiger9_occu <- brm(occu ~ poly(mean_temperature, 2) +
                      mean_relative_humidity + 
                      agricultural + forests_scrub + 
                      log(trapping_effort) + log(n_traps) +
                      (1 | pixel_id) + (1 | year),
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
# saveRDS(mtiger9_occu, file = paste0(loc.output, "mtiger9_occu.rds"))

mtiger10_occu <- brm(occu ~ poly(mean_temperature, 2) +
                       agricultural + forests_scrub + 
                       log(trapping_effort) + log(n_traps) +
                       (1 | pixel_id) + (1 | year),
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
# saveRDS(mtiger10_occu, file = paste0(loc.output, "mtiger10_occu.rds"))

mtiger11_occu <- brm(occu ~ poly(mean_temperature, 2) +
                      agricultural + forests_scrub + discont_urban_fabric +
                      log(trapping_effort) + log(n_traps) +
                      (1 | pixel_id) + (1 | year),
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
# saveRDS(mtiger11_occu, file = paste0(loc.output, "mtiger11_occu.rds"))

loo_compare(loo(mtiger1_occu), loo(mtiger2_occu), loo(mtiger3_occu), loo(mtiger4_occu), 
            loo(mtiger5_occu), loo(mtiger6_occu), loo(mtiger7_occu), loo(mtiger8_occu), 
            loo(mtiger9_occu), loo(mtiger10_occu), loo(mtiger11_occu))
loo(mtiger10_occu)

# Mosquito Alert models --------------------------------------------------------
ma_df <- readRDS(paste0(loc.output, "ma_tiger_spain.rds"))

# these_priors = c(prior(student_t(4, 0, 6), class = "b"))
mtiger1_ma <- brm(any_reps ~ poly(max_temperature, 2) + 
                    mean_relative_humidity + 
                    (1 | pixel_id) + (1 | year) +  offset(log(SE)),
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
# saveRDS(mtiger1_ma, file = paste0(loc.output, "mtiger1_ma.rds"))

mtiger2_ma <- brm(any_reps ~ poly(mean_temperature, 2) + 
                    mean_relative_humidity + 
                    (1 | municipality) + (1 | year) + offset(log(SE)),
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
# saveRDS(mtiger2_ma, file = paste0(loc.output, "mtiger2_ma.rds"))

mtiger3_ma <- brm(any_reps ~ poly(min_temperature, 2) + 
                    mean_relative_humidity + 
                    (1 | pixel_id) + (1 | year) + offset(log(SE)),
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
# saveRDS(mtiger3_ma, file = paste0(loc.output, "mtiger3_ma.rds"))

mtiger4_ma <- brm(any_reps ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                    agricultural + green_urban +
                    (1 | municipality) + (1 | year) + offset(log(SE)),
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
# saveRDS(mtiger4_ma, file = paste0(loc.output, "mtiger4_ma.rds"))

mtiger5_ma <- brm(any_reps ~ poly(mean_temperature, 2) + mean_relative_humidity + 
                    agricultural + cont_urban_fabric +
                    (1 | id) + (1 | year) + offset(log(SE)),
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
# saveRDS(mtiger5_ma, file = paste0(loc.output, "mtiger5_ma.rds"))

loo_compare(loo(mtiger2_ma), loo(mtiger4_ma))

# Models using number of report as SE (avoiding SE calculated from historical data)
ma_df <- ma_df %>% 
  filter(n_total_reports > 0) # Si no hay esfuerzo de muestreo el dato no debería de estar.

mtiger1_ma_rep <- brm(any_reps ~ poly(mean_temperature, 2) + 
                        mean_relative_humidity + 
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
