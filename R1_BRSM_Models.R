############################## BG Model ########################################
#' Mosquito Distribution Modeling:
#' Hierarchical TRAP and CITSCI model using Bayesian approach

#' Cleaned data is available on DATA folder of CatuCerecedo/citizen-mosquito-surveillance 
#' repository

################################################################################

# Dependencies

if (!require("tidyverse")) install.packages("tidyverse")
if (!require("rstanarm")) install.packages("rstanarm")
if (!require("brms")) install.packages("brms")
if (!require("cmdstanr")) install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
if (!require("loo")) install.packages("loo")

library(tidyverse)
library(rstanarm)
library(brms)
library(shinystan)
library(loo)
library(cmdstanr)

# Checking cmdstanr library
check_cmdstan_toolchain()
cmdstan_version()

rm(list = ls())
# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")

# Loading data -----------------------------------------------------------------
tiger <- readRDS(paste0(loc.data, "trap_surveillance_data.rds")) 
ma_df <- readRDS(paste0(loc.data, "citizen_observations.rds")) 

# Modeling process -------------------------------------------------------------
# Parameters of brm models
nchains = 4
threads_per_chain = 1
iteret = 5000
wup = 2000

mtiger16_trap <- brm(females ~ poly(l21mean_temperature, 2) + l21precipitation +
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
saveRDS(mtiger16_trap, file = paste0(loc.output, "mtiger16_trap.rds"))

mtiger7_ma <- brm(any_reps ~ poly(mean_temperature, 2) + mean_relative_humidity +
                     (1 | id) + (1 | y) +
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
                  save_pars = save_pars(all = TRUE),
                  control = list(adapt_delta = 0.99))
saveRDS(mtiger7_ma, file = paste0(loc.output, "mtiger7_ma.rds"))

# Visual check of models
plot(mtiger16_trap) # Converge
plot(mtiger7_ma) # Converge

# Other checks (Pareto and R2_bayes)
loo(mtiger16_trap, moment_match = TRUE, recompile = TRUE, reloo = TRUE)
loo(mtiger7_ma, moment_match = TRUE, recompile = TRUE, reloo = TRUE)

bayes_R2(mtiger16_trap)
bayes_R2(mtiger7_ma)


