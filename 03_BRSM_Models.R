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
# set_cmdstan_path("/home/usuaris/ccerecedo/cmdstan")

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
tiger <- readRDS(paste0(loc.output, "bg_tiger_spain_daily.rds")) 
# tiger$year <- as.factor(lubridate::year(tiger$end_date))

tiger <- tiger %>% 
  mutate(
    # urban = artificial_green_urban + industrial_transport + urban_fabric,
    y = as.factor(year(end_date)),
    occu = case_when(
      females > 0 ~ 1,
      TRUE ~ 0),
    FH = case_when(mean_relative_humidity < 40~0, mean_relative_humidity >95~0, 
                   (mean_relative_humidity >=40 & mean_relative_humidity <= 95)~((mean_relative_humidity/55)-(40/55))),
    FT = case_when(mean_temperature<=15~0, mean_temperature>30~0, (mean_temperature>15 & mean_temperature <=20)~(.2*mean_temperature)-3, (mean_temperature>20 & mean_temperature<=25)~1, (mean_temperature>25 & mean_temperature <= 30)~(-.2*mean_temperature)+6),
    mwi = FH*FT,
  ) %>%
  filter(mean_relative_humidity < 101) # There is a wrong data --> ERA5 has no sense

# Filtering by municpilaties with more than 2 samples:
mun_selected <- tiger %>% group_by(id) %>% summarize(n = n()) %>% filter(n > 3)
length(unique(mun_selected$id))

tiger <- tiger %>% filter(id %in% mun_selected$id)

# Mosquito Alert data
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

# Modeling process -------------------------------------------------------------
nchains = 4
threads_per_chain = 1

iteret = 5000
wup = 2000

# NOT RUN
# Only when you are working in local
# Sys.setenv(PATH = "C:/Users/Catu/Documents/.cmdstan/cmdstan-2.31.0/stan/lib/stan_math/lib/tbb")

# mtiger16 <- brm(females ~ poly(l21mean_temperature, 2) + l21precipitation + 
#                   offset(log(trapping_effort)) +
#                   (1 | id) + (1 | y),
#                 data = tiger,
#                 prior = set_prior("cauchy(0,2.5)", class="b"),
#                 family = negbinomial(link = "log"),
#                 iter = iteret,
#                 chains = nchains,
#                 cores = nchains,
#                 backend = "cmdstanr",
#                 save_pars = save_pars(all = TRUE),
#                 threads = threading(threads_per_chain),
#                 control = list(adapt_delta = 0.999))
# loo(mtiger16)
# loo(mtiger16, moment_match = TRUE, recompile = TRUE, reloo = TRUE)
# bayes_R2(mtiger16)
# saveRDS(mtiger16, file = paste0(loc.output, "mtiger16.rds"))
# 
# mtiger3_occu <- brm(occu ~ poly(l21mean_temperature, 2) + 
#                       discont_urban_fabric + agricultural +
#                       offset(log(trapping_effort)) +
#                       (1 | id) + (1 | y),
#                     data = tiger,
#                     prior = set_prior("cauchy(0,2.5)", class="b"),
#                     family = bernoulli(link = "logit"),
#                     iter = iteret,
#                     warmup = wup,
#                     chains = nchains,
#                     cores = nchains,
#                     backend = "cmdstanr",
#                     threads = threading(threads_per_chain),
#                     save_pars = save_pars(all = TRUE),
#                     control = list(adapt_delta = 0.99))
# loo(mtiger3_occu)
# loo(mtiger3_occu, moment_match = TRUE, recompile = TRUE, reloo = TRUE)
# bayes_R2(mtiger3_occu)
# saveRDS(mtiger3_occu, file = paste0(loc.output, "mtiger3_occu.rds"))
# 
# mtiger7_ma <- brm(any_reps ~ poly(mean_temperature, 2) + mean_relative_humidity + 
#                      (1 | id) + (1 | y) +
#                      offset(log(SE)),
#                   data = ma_df,
#                   prior = set_prior("cauchy(0,2.5)", class="b"),
#                   family = bernoulli(link = "logit"),
#                   iter = iteret,
#                   warmup = wup,
#                   chains = nchains,
#                   cores = nchains,
#                   backend = "cmdstanr",
#                   threads = threading(threads_per_chain),
#                   save_pars = save_pars(all = TRUE),
#                   control = list(adapt_delta = 0.99))
# loo(mtiger7_ma)
# loo(mtiger7_ma, moment_match = TRUE, recompile = TRUE, reloo = TRUE)
# bayes_R2(mtiger7_ma)
# saveRDS(mtiger7_ma, file = paste0(loc.output, "mtiger7_ma.rds"))
# 
# # Models using number of report as SE (avoiding SE calculated from historical data)
# ma_df <- ma_df %>% 
#   filter(n_total_reports > 0) # Si no hay esfuerzo de muestreo el dato no debería de estar.
# 
# mtiger2_ma_rep <- brm(any_reps ~ poly(mean_temperature, 2) + precipitation + 
#                         agricultural + green_urban +
#                         (1 | prov_name) + (1 | year) +
#                         offset(log(n_total_reports + 0.0001)),
#                       data = ma_df,
#                       prior = set_prior("cauchy(0,2.5)", class="b"),
#                       family = bernoulli(link = "logit"),
#                       iter = iteret,
#                       warmup = wup,
#                       chains = nchains,
#                       cores = nchains,
#                       backend = "cmdstanr",
#                       threads = threading(threads_per_chain),
#                       save_pars = save_pars(all = TRUE),
#                       control = list(adapt_delta = 0.99))
# saveRDS(mtiger2_ma_rep, file = paste0(loc.output, "mtiger2_ma_rep.rds"))

bites <- readRDS(file = paste0(loc.output, "bites_spain_daily.rds"))

# Reduce he number of FALSE
bites_t <- bites %>% filter(any_reps == TRUE)
bites_f <- bites %>% filter(any_reps == FALSE) %>% sample_n(8000)

bites <- rbind(bites_t, bites_f)

mbites_6_reduce <- brm(any_reps ~ poly(l21min_temperature, 2) + 
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
saveRDS(mbites_6_reduce, file = paste0(loc.output,"mbites_6_reduce.rds"))
summary(mbites_6_reduce)
loo(mbites_6_reduce)
