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
set_cmdstan_path("/home/usuaris/ccerecedo/cmdstan")

check_cmdstan_toolchain()
cmdstan_version()
 
# remotes::install_github("stan-dev/cmdstanr")

rm(list = ls())
# Directories ------------------------------------------------------------------
# In local
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")

# In Cluster
loc.output <- paste0(getwd(), "/Spain_Tiger/OUTPUT/")
loc.data <- paste0(getwd(), "/Spain_Tiger/DATA/")

# Loading ---------------------------------------------------------------------
# Count data
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

# Mosquito Alert data
ma_df <- readRDS(paste0(loc.output, "ma_tiger_spain.rds")) %>%
  mutate(
    id = as.factor(id),
    year = as.factor(year)
  ) %>%
  filter(n_total_reports > 0)

# Modeling process -------------------------------------------------------------
nchains = 4
threads_per_chain = 1

iteret = 5000
wup = 2000

# NOT RUN
# Only when you are working in local
# Sys.setenv(PATH = "C:/Users/Catu/Documents/.cmdstan/cmdstan-2.31.0/stan/lib/stan_math/lib/tbb")

mtiger5 <- brm(females ~ poly(mean_temperature, 2) + precipitation + 
                 agricultural + green_urban +
                 offset(log(trapping_effort)) + offset(log(n_traps)) +
                 (1 | id) + (1 | year),
               data = tiger,
               prior = set_prior("cauchy(0,2.5)", class="b"),
               family = negbinomial(link = "log"),
               iter = iteret,
               chains = nchains,
               cores = nchains,
               backend = "cmdstanr",
               threads = threading(threads_per_chain),
               control = list(adapt_delta = 0.99))
saveRDS(mtiger10, file = paste0(loc.output, "mtiger10.rds"))

mtiger1_occu <- brm(occu ~ poly(mean_temperature, 2) + precipitation + 
                      agricultural + green_urban +
                      offset(log(trapping_effort)) + offset(log(n_traps)) +
                      (1 | id) + (1 | year),
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
saveRDS(mtiger1_occu, file = paste0(loc.output, "mtiger1_occu.rds"))

mtiger5_ma <- brm(any_reps ~ poly(mean_temperature, 2) + precipitation + 
                    agricultural + green_urban +
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
saveRDS(mtiger5_ma, file = paste0(loc.output, "mtiger5_ma.rds"))

# Models using number of report as SE (avoiding SE calculated from historical data)
ma_df <- ma_df %>% 
  filter(n_total_reports > 0) # Si no hay esfuerzo de muestreo el dato no debería de estar.

mtiger2_ma_rep <- brm(any_reps ~ poly(mean_temperature, 2) + precipitation + 
                        agricultural + green_urban +
                        (1 | id) + (1 | year) +  offset(log(n_total_reports + 0.0001)),
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
saveRDS(mtiger2_ma_rep, file = paste0(loc.output, "mtiger2_ma_rep.rds"))
