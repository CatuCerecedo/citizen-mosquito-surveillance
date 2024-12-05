################### Cross-validation of counts models ##########################
#' We calculate the accuracy of the model following the cross-validation of 
#' machine learning. Using Pareto: 20/80

library(tidyverse)
library(rstanarm)
library(brms)
library(shinystan)
library(loo)
library(cmdstanr)
library(caret)

# Directories ------------------------------------------------------------------

# In local
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")
loc.era5 <- "/home/catuxa/Documents/Mosquito_Models/EU_Culex/ERA5_Download/"

# In cluster
loc.output <- paste0(getwd(), "/Spain_Tiger/OUTPUT/")
loc.data <- paste0(getwd(), "/Spain_Tiger/DATA/")
loc.era5 <- paste0(getwd(), "/EU_Culex/ERA5_Download/")

# Load data ----------------------------------------------------------------
tiger_integrating <- readRDS(file = paste0(loc.output, "tiger_integrating_daily.rds")) %>%
  drop_na(ma) %>%
  mutate(
    log_ma = log(ma)
  )

# Plotting relationship bg ~ ma
ggplot(tiger_integrating, aes(y = females, x = log(ma))) + 
  geom_point() +
  geom_smooth() +
  theme_classic() 

ggplot(tiger_integrating, aes(x = females > 0, y = log(ma))) + 
  geom_boxplot() +
  theme_classic() 

# BG MODELS---------------------------------------------------------------------
# Split the data: train/test ---------------------------------------------------

set.seed(12123)
sample <- sample(c(TRUE, FALSE), nrow(tiger_integrating), replace=TRUE, prob=c(0.8, 0.2))
train  <- tiger_integrating[sample, ]
test   <- tiger_integrating[!sample, ]

nchains = 4
threads_per_chain = 1

iteret = 2500
wup = 1000
# Count model ---------------------------------------------------------------
mtiger16 <- brm(females ~ poly(l21mean_temperature, 2) + l21precipitation + 
                  offset(log(trapping_effort)) +
                  (1 | id) + (1 | y),
                data = train,
                prior = set_prior("cauchy(0,2.5)", class="b"),
                family = negbinomial(link = "log"),
                iter = iteret,
                chains = nchains,
                cores = nchains,
                backend = "cmdstanr",
                save_pars = save_pars(all = TRUE),
                threads = threading(threads_per_chain),
                control = list(adapt_delta = 0.999))
saveRDS(mtiger16, file = paste0(loc.output, "mtiger16_cv.rds"))

mtiger16 <- readRDS(file = paste0(loc.output, "mtiger16_cv.rds"))
# Calculating the predictions: test --------------------------------------------

pp <- apply(posterior_predict(mtiger16, 
                              newdata = test,
                              allow_new_levels = TRUE, 
                              re_formula = NA,
                              ndraws = 1000), 
            2, function(x) mean(x)) %>%
  as.data.frame()

test$predicted <- pp$.

ggplot(test, aes(x = predicted, y = females)) +
  geom_point() +
  # annotate("text", x = 110, y = 500, label = "y = -6.13 + 1.42x") +
  # annotate("text", x = 110, y = 400, label = "S = 0.63") +
  geom_smooth(method = "lm") +
  theme_bw() 
summary(lm(females ~ predicted, data = test))
summary(lm(I(females - 1) ~ predicted - 1, data = test))

library(caret)
cv_metrics <- data.frame(Scor = cor(test$predicted, test$females, method = "spearman"),
                         R2 = R2(test$predicted, test$females),
                         RMSE = RMSE(test$predicted, test$females),
                         MAE = MAE(test$predicted, test$females),
                         model = "Count model")

# Integrated model -------------------------------------------------------------
mtiger_integrated <- brm(bf(females ~ poly(l21mean_temperature, 2) + l21precipitation + 
                              offset(log(trapping_effort))  + 
                              (1 | id) + (1 | y),
                            zi ~ log_ma + (1 | id) + (1 | y)),
                         data = train,
                         prior = set_prior("cauchy(0,2.5)", class="b"),
                         family = zero_inflated_negbinomial(link = "log"),
                         iter = iteret,
                         chains = nchains,
                         cores = nchains,
                         backend = "cmdstanr",
                         save_pars = save_pars(all = TRUE),
                         threads = threading(threads_per_chain),
                         control = list(adapt_delta = 0.999))
loo::loo(mtiger_integrated)
# plot(conditional_effects(mtiger_integrated, effects = "log_ma", dpar = "zi"))
summary(mtiger_integrated)

saveRDS(mtiger_integrated, file = paste0(loc.output, "mtiger_integrated_cv.rds"))

mtiger_integrated <- readRDS(paste0(loc.output, "mtiger_integrated_cv.rds"))

# Calculating the predictions: test --------------------------------------------

pp <- apply(posterior_predict(mtiger_integrated, 
                              newdata = test,
                              allow_new_levels = TRUE, 
                              re_formula = NA,
                              ndraws = 1000), 
            2, function(x) mean(x)) %>%
  as.data.frame()

test$predicted <- pp$.
ggplot(test, aes(x = predicted, y = females)) +
  geom_point() +
  # annotate("text", x = 110, y = 500, label = "y = -6.13 + 1.42x") +
  # annotate("text", x = 110, y = 400, label = "S = 0.63") +
  geom_smooth(method = "lm") +
  theme_bw() 
summary(lm(I(females - 1) ~ predicted - 1, data = test))

cv_metrics <- cv_metrics %>% 
  bind_rows(data.frame(Scor = cor(test$predicted, test$females, method = "spearman"),
                       R2 = R2(test$predicted, test$females),
                       RMSE = RMSE(test$predicted, test$females),
                       MAE = MAE(test$predicted, test$females),
                       model = "Integrated model (zi on)"))

# Integrated model -------------------------------------------------------------
mtiger_integrated_2 <- brm(bf(females ~ log_ma + 
                              (1 | id) + (1 | y),
                            zi ~ poly(l21mean_temperature, 2) + l21precipitation + 
                              offset(log(trapping_effort))  + (1 | id) + (1 | y)),
                         data = train,
                         prior = set_prior("cauchy(0,2.5)", class="b"),
                         family = zero_inflated_negbinomial(link = "log"),
                         iter = iteret,
                         chains = nchains,
                         cores = nchains,
                         backend = "cmdstanr",
                         save_pars = save_pars(all = TRUE),
                         threads = threading(threads_per_chain),
                         control = list(adapt_delta = 0.999))
loo::loo(mtiger_integrated_2)
# plot(conditional_effects(mtiger_integrated, effects = "log_ma", dpar = "zi"))
summary(mtiger_integrated_2)

loo_compare(loo(mtiger_integrated), loo(mtiger_integrated_2), loo(mtiger16))

saveRDS(mtiger_integrated_2, file = paste0(loc.output, "mtiger_integrated_2_cv.rds"))

mtiger_integrated_2 <- readRDS(paste0(loc.output, "mtiger_integrated_2_cv.rds"))

# Calculating the predictions: test --------------------------------------------

pp <- apply(posterior_predict(mtiger_integrated_2, 
                              newdata = test,
                              allow_new_levels = TRUE, 
                              re_formula = NA,
                              ndraws = 1000), 
            2, function(x) mean(x)) %>%
  as.data.frame()

test$predicted <- pp$.
ggplot(test, aes(x = predicted, y = females)) +
  geom_point() +
  # annotate("text", x = 110, y = 500, label = "y = -6.13 + 1.42x") +
  # annotate("text", x = 110, y = 400, label = "S = 0.63") +
  geom_smooth(method = "lm") +
  theme_bw() 
summary(lm(I(females - 1) ~ predicted - 1, data = test))

cv_metrics <- cv_metrics %>% 
  bind_rows(data.frame(Scor = cor(test$predicted, test$females, method = "spearman"),
                       R2 = R2(test$predicted, test$females),
                       RMSE = RMSE(test$predicted, test$females),
                       MAE = MAE(test$predicted, test$females),
                       model = "Integrated model (zi off)"))

# Quicky way to get meaures ----------------------------------------------------

iteret = 250
wup = 100
nchains = 4
threads_per_chain = 1

cv_metrics <- data.frame()

for (i in 1:10){

  # Preparing test/training datasets
  set.seed(i)
  sample <- sample(c(TRUE, FALSE), nrow(tiger_integrating), replace=TRUE, prob=c(0.8, 0.2))
  train  <- tiger_integrating[sample, ]
  test   <- tiger_integrating[!sample, ]
  
  # Count model
  mtiger16 <- brm(bf(females ~ poly(l21mean_temperature, 2) + l21precipitation + 
                       offset(log(trapping_effort))  + 
                       (1 | id) + (1 | y),
                     zi ~ 1),
                  data = train,
                  prior = set_prior("cauchy(0,2.5)", class="b"),
                  family = zero_inflated_negbinomial(link = "log"),
                  iter = iteret,
                  chains = nchains,
                  cores = nchains,
                  backend = "cmdstanr",
                  save_pars = save_pars(all = TRUE),
                  threads = threading(threads_per_chain),
                  control = list(adapt_delta = 0.999))
    pp <- apply(posterior_predict(mtiger16, 
                                newdata = test,
                                allow_new_levels = TRUE, 
                                re_formula = NA,
                                ndraws = 100), 
              2, function(x) mean(x)) %>%
    as.data.frame()
  test$predicted <- pp$.
  cv_metrics <- cv_metrics %>% 
    bind_rows(data.frame(iter = i,
                         Scor = cor(test$predicted, test$females, method = "spearman"),
                         R2 = R2(test$predicted, test$females),
                         RMSE = RMSE(test$predicted, test$females),
                         MAE = MAE(test$predicted, test$females),
                         model = "Count model"))
  
  # Integrated model -------------------------------------------------------------
  mtiger_integrated <- brm(bf(females ~ poly(l21mean_temperature, 2) + l21precipitation + 
                                offset(log(trapping_effort))  + 
                                (1 | id) + (1 | y),
                              zi ~ log_ma + (1 | id) + (1 | y)),
                           data = train,
                           prior = set_prior("cauchy(0,2.5)", class="b"),
                           family = zero_inflated_negbinomial(link = "log"),
                           iter = iteret,
                           chains = nchains,
                           cores = nchains,
                           backend = "cmdstanr",
                           save_pars = save_pars(all = TRUE),
                           threads = threading(threads_per_chain),
                           control = list(adapt_delta = 0.999))
  
  pp <- apply(posterior_predict(mtiger_integrated, 
                                newdata = test,
                                allow_new_levels = TRUE, 
                                re_formula = NA,
                                ndraws = 100), 
              2, function(x) mean(x)) %>%
    as.data.frame()
  test$predicted <- pp$.
  cv_metrics <- cv_metrics %>% 
    bind_rows(data.frame(iter = i,
                         Scor = cor(test$predicted, test$females, method = "spearman"),
                         R2 = R2(test$predicted, test$females),
                         RMSE = RMSE(test$predicted, test$females),
                         MAE = MAE(test$predicted, test$females),
                         model = "Integrated model (zi on)"))
  
  # Integrated model -----------------------------------------------------------
  mtiger_integrated_2 <- brm(bf(females ~ log_ma + 
                                  (1 | id) + (1 | y),
                                zi ~ poly(l21mean_temperature, 2) + l21precipitation + 
                                  offset(log(trapping_effort))  + (1 | id) + (1 | y)),
                             data = train,
                             prior = set_prior("cauchy(0,2.5)", class="b"),
                             family = zero_inflated_negbinomial(link = "log"),
                             iter = iteret,
                             chains = nchains,
                             cores = nchains,
                             backend = "cmdstanr",
                             save_pars = save_pars(all = TRUE),
                             threads = threading(threads_per_chain),
                             control = list(adapt_delta = 0.999))
  pp <- apply(posterior_predict(mtiger_integrated_2, 
                                newdata = test,
                                allow_new_levels = TRUE, 
                                re_formula = NA,
                                ndraws = 100), 
              2, function(x) mean(x)) %>%
    as.data.frame()
  test$predicted <- pp$.
  
  cv_metrics <- cv_metrics %>% 
    bind_rows(data.frame(iter = i,
                         Scor = cor(test$predicted, test$females, method = "spearman"),
                         R2 = R2(test$predicted, test$females),
                         RMSE = RMSE(test$predicted, test$females),
                         MAE = MAE(test$predicted, test$females),
                         model = "Integrated model (zi off)"))
}

# saveRDS(cv_metrics, file = paste0(loc.output, "cv_metrics_10.rds"))
cv_metrics <- readRDS(paste0(loc.output, "cv_metrics_10.rds"))

cv_metrics <- cv_metrics %>% mutate(model = case_when(model == "Count model" ~ "ZI-MA Null",
                                                      model ==  "Integrated model (zi on)" ~ "ZI-MA Presence", 
                                                      model ==  "Integrated model (zi off)" ~ "ZI-MA Counts"))

cv_metrics$model <- factor(cv_metrics$model, levels = c("ZI-MA Null", "ZI-MA Counts", "ZI-MA Presence"))

a <- ggplot(cv_metrics) +
  geom_boxplot(aes(y = model, x = Scor )) +
  labs(
    y = "Model\n",
    x = "\nSpearman Spatial Correlation (S)\n"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",    # La leyenda se posiciona abajo
    legend.title = element_text(size = 22, face = "bold", margin = margin(b=20)),
    legend.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    # axis.text.x = element_text(angle = 35, vjust = 0.5),
    legend.key.spacing.y = unit(0.5, 'cm')
  ) 

b <- ggplot(cv_metrics) +
  geom_boxplot(aes(y = model, x = R2 )) +
  labs(
    y = "",
    x = "\nGoodness-of-fit (R²)\n"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",    # La leyenda se posiciona abajo
    legend.title = element_text(size = 22, face = "bold", margin = margin(b=20)),
    legend.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    # axis.text.x = element_text(angle = 35, vjust = 0.5),
    legend.key.spacing.y = unit(0.5, 'cm')
  ) 

c <- ggplot(cv_metrics) +
  geom_boxplot(aes(y = model, x = RMSE )) +
  labs(
    y = "Model\n",
    x = "\nRoot Mean Squared Error (RMSE)\n"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",    # La leyenda se posiciona abajo
    legend.title = element_text(size = 22, face = "bold", margin = margin(b=20)),
    legend.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    # axis.text.x = element_text(angle = 35, vjust = 0.5),
    legend.key.spacing.y = unit(0.5, 'cm')
  ) 

d <- ggplot(cv_metrics) +
  geom_boxplot(aes(y = model, x = MAE )) +
  labs(
    y = "",
    x = "\nMean Absolute Error (MAE)\n"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",    # La leyenda se posiciona abajo
    legend.title = element_text(size = 22, face = "bold", margin = margin(b=20)),
    legend.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    # axis.text.x = element_text(angle = 35, vjust = 0.5),
    legend.key.spacing.y = unit(0.5, 'cm')
  ) 

ggpubr::ggarrange(a,b,c,d)
ggsave(file = paste0(loc.fig, "CV_zero_inflated.png"), 
       units = "cm", height = 25, width = 42, bg = "white")

loo_compare(loo(mtiger16), loo(mtiger_integrated), loo(mtiger_integrated_2))

# Ploting the elpd -------------------------------------------------------------
elpd <- data.frame(
  model = factor(c("ZI-MA NULL", "ZI-MA Counts", "ZI-MA Presence"), levels = c("ZI-MA NULL", "ZI-MA Counts", "ZI-MA Presence")),
  elpd = c(-6542.6, -6616.2, -6492.1 ),
  elpd_se = c(74.2, 91.6, 75.7)
)

ggplot(elpd, aes(x = elpd, y = model)) + 
  geom_errorbar(aes(xmin =  elpd - elpd_se, xmax = elpd + elpd_se)) +
  geom_point(shape = 21, fill = "#3288bd", size = 3) +
  labs(
    x = "elpd",
    y = ""
  ) +
  theme_classic() +
  theme(
    legend.position = "right",    # La leyenda se posiciona abajo
    legend.title = element_text(size = 22, face = "bold", margin = margin(b=20)),
    legend.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    legend.key.spacing.y = unit(0.5, 'cm')
  ) 
ggsave(file = paste0(loc.fig, "CV_elpd_zero_inflated.png"), 
       units = "cm", height = 15, width = 20, bg = "white")
