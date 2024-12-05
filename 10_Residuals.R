############################## Residuals analysis ##############################
library(ggplot2)
library(tidyverse)
library(sf)
library(parallel)
# library(mgcv)

rm(list = ls())

# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.heavy <- paste0(getwd(), "/Heavy_files/")
loc.clc <- paste0(getwd(), "/OUTPUT/land_cover_raster/")

# Loading data -----------------------------------------------------------------
all <- readRDS(file = paste0(loc.output, "all_monthly_pred_mun.rds"))

# Fit the gam curve ------------------------------------------------------------
# We fit the same curve that ggplot2 returns
# fit_gam <- mgcv::gam(bg ~ s(ma, bs = "cs"), data = all)
fit_pol <- lm(bg ~ poly(ma, 2), data = all)
summary(fit_pol)

# residuals = observed - fitted values
all$fitted_values <- predict(fit_pol)
all$residuals <- all$bg - all$fitted_values

# The same plot as ggplot2
ggplot(all %>% filter(m == "08")) +
  geom_point(aes(x = ma, y = bg), alpha = 0.6, size = 0.5) +
  geom_smooth(aes(x = ma, y = bg), linetype = "dashed", color = "purple") +
  # geom_smooth(aes(x = ma, y = bg), method = "lm", formula = y ~ poly(x, 2), linetype = "dashed", color = "purple") +
  labs(
    x = "Mosquito Alert predicted probabilities",
    y = "Count predicted values"
  ) +
  theme_classic()

# Load variables data ----------------------------------------------------------
clc <- readRDS(paste0(loc.output, "clc_surface_mun.rds"))

all <- merge(all %>% filter(year == 2020), clc, by = "id", all.x = TRUE)
clc_names <- names(all)[10:23]

clc_residuals_plots <- list()
for (i in 1:length(clc_names)){
  print(paste(clc_names[i]))
  clc_residuals_plots[[i]] <- ggplot(all, aes_string(x = "residuals", y = clc_names[i])) +
    geom_point(alpha = 0.5, color = "#3288bd") +
    geom_smooth(method = "lm", color = "#d53e4f", se = FALSE) +
    labs(
      x = "Residuals",
      y = clc_names[i]
    ) +
    theme_classic()
}

clc_residuals_plots[[1]]
clc_residuals_plots[[2]]
clc_residuals_plots[[3]]
clc_residuals_plots[[4]]
clc_residuals_plots[[5]]
clc_residuals_plots[[6]]
clc_residuals_plots[[7]]
clc_residuals_plots[[8]]
clc_residuals_plots[[9]]
clc_residuals_plots[[10]]
clc_residuals_plots[[11]]
clc_residuals_plots[[12]]

ggplot(all, aes_string(x = "residuals", y = "m")) +
  geom_point(alpha = 0.5, color = "#3288bd") +
  geom_smooth(method = "lm", color = "#d53e4f", se = FALSE) +
  labs(
    x = "Residuals",
    y = "Months"
  ) +
  theme_classic()

# Joining temperature and precipitation ----------------------------------------
all_2020 <- all %>% filter(year == "2020")

all_data <- mclapply(1:nrow(all), function(i){
  
  message(paste0("Number row:", i, "\n"))
  
  data_row <- all[i, ]
  
  wth <- readRDS(paste0(loc.output, "monthly_weather_data/prep_predictions_", 
                        data_row$m, "_",data_row$year, ".rds"))
  
  data_row <- merge(data_row, wth, by = "id", all.x = TRUE)
}
, mc.cores = 1)

all_data <- do.call(rbind, all_data) 
