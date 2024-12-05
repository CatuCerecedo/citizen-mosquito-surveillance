########################## Creating MA predictor ###############################
library(tidyverse)
library(parallel)
library(data.table)
library(ggplot2)
library(loo)
library(rstanarm)
library(brms)

rm(list = ls())
# Directories ------------------------------------------------------------------

# In local
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.era5 <- "/home/catuxa/Documents/Mosquito_Models/EU_Culex/ERA5_Download/"

# In cluster
loc.output <- paste0(getwd(), "/Spain_Tiger/OUTPUT/")
loc.data <- paste0(getwd(), "/Spain_Tiger/DATA/")
loc.era5 <- paste0(getwd(), "/EU_Culex/ERA5_Download/")

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

# Creating MA predictor --------------------------------------------------------
tiger_integrating <- bind_rows(mclapply(1:nrow(tiger), function(i){
  
  cat(paste0("Number row:", i, "\n"))
  
  data_row <- tiger[i, ]
  
  start <- as.Date(data_row$start_date)
  end <- as.Date(data_row$end_date)
  id_row <- data_row$id
  
  ma_vector <- vector()
  
  if (year(end) == "2023"){
    data_row$ma <- NA
  } else {for (d in seq.Date(start, end, 1)){
    d <- as.Date(d)
    # print(d)
    # load the predictions by month and year
    m <- sprintf("%02d", (month(d)))
    y <- as.character(year(d))
    pred <- readRDS(paste0(loc.output, "PREDICTIONS/MA/mtiger7_ma/tiger_", m, "_", y, "_ma.rds"))
    
    # Filter the data that we need
    pred <- pred %>% filter(id == id_row)
    ma <- pred[,as.character(d)]
    
    ma_vector <- append(ma_vector, ma)
  }
    data_row$ma <- prod(ma_vector, na.rm = TRUE)
    }
  return(data_row)
}
, mc.cores = 8))

saveRDS(tiger_integrating, file = paste0(loc.output, "tiger_integrating_daily.rds"))

# plotting ----------------------------------------------------------------------

ggplot(tiger_integrating, aes(y = females, x = log(ma))) + 
  geom_point() +
  geom_smooth() +
  theme_classic() 

ggplot(tiger_integrating, aes(x = females > 0, y = log(ma))) + 
  geom_boxplot() +
  theme_classic() 

# Zero-inflated model ----------------------------------------------------------
library(brms)
library(shinystan)
library(loo)
library(cmdstanr)

nchains = 4
threads_per_chain = 1

iteret = 5000
wup = 2000

tiger_integrating <- tiger_integrating %>% mutate(
  log_ma = log(ma)
)
mtiger_integrated <- brm(bf(females ~ poly(l21mean_temperature, 2) + l21precipitation + 
                                offset(log(trapping_effort))  + 
                                (1 | id) + (1 | y),
                    zi ~ log_ma + (1 | id) + (1 | y)),
                 data = tiger_integrating,
                 prior = set_prior("cauchy(0,2.5)", class="b"),
                 family = zero_inflated_negbinomial(link = "log"),
                 iter = iteret,
                 chains = nchains,
                 cores = nchains,
                 backend = "cmdstanr",
                 save_pars = save_pars(all = TRUE),
                 threads = threading(threads_per_chain),
                 control = list(adapt_delta = 0.999))

loo(mtiger_integrated)
# loo(mtiger_integrated, moment_match = TRUE, recompile = TRUE, reloo = TRUE)
bayes_R2(mtiger_integrated)
saveRDS(mtiger_integrated, file = paste0(loc.output, "mtiger_integrated.rds"))

# Plotting the marginal and conditional effect
a <- conditional_effects(mtiger_integrated, effects = "log_ma", dpar = "zi")
conditional_data <- a[["log_ma"]] %>%
  mutate(
    ma = exp(effect1__),
    females = estimate__)

ggplot(conditional_data) +
  geom_point(aes(x = ma, y = females)) +  
  scale_x_continuous(trans = "log") +
  labs(
    x = "Mosquito Alert Probability",
    y = "P(y = 0)"
  ) +
  theme_classic()

b <- marginal_effects(mtiger_integrated)
marginal_data <- b[["log_ma"]] %>%
  mutate(
    ma = exp(effect1__),
    females = estimate__)

scientific_format_custom <- function(x) {
  parse(text = ifelse(x == 1, 1, gsub("e", " %*% 10^", scales::scientific_format(digits = 2)(x))))
}

ggplot(marginal_data, aes(x = ma, y = females)) + 
  geom_line(size = 1, color = "#2c7fb8") +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#2c7fb8", alpha = 0.2) + 
  scale_x_continuous(trans = "log", labels = scientific_format_custom) +  
  labs(
    x = "Mosquito Alert Probability",
    y = "Estimated Counts of Females"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),  # Texto de ejes en negro
    axis.title = element_text(face = "bold", size = 14),  # Títulos de los ejes en negrita
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)  # Centrar el título y aumentar tamaño
  )

ggplot(marginal_data, aes(x = females, y = ma)) + 
  geom_line(size = 1, color = "#2c7fb8") +  
  # geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#2c7fb8", alpha = 0.2) + 
  scale_y_continuous(trans = "log", labels = scientific_format_custom) +  
  labs(
    x = "Mosquito Alert Probability",
    y = "Estimated Counts of Females"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),  # Texto de ejes en negro
    axis.title = element_text(face = "bold", size = 14),  # Títulos de los ejes en negrita
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)  # Centrar el título y aumentar tamaño
  )
# Comparing with the non-zero one
mtiger16 <- readRDS(file = paste0(loc.output, "mtiger16.rds"))
loo_compare(loo(mtiger16), loo(mtiger_integrated)) 


