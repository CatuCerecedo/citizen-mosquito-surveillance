############################## Residuals analysis ##############################
library(ggplot2)
library(tidyverse)
library(sf)
library(parallel)
library(mapSpain) # Load library mapSpain to have shapefile with NATCODE municipalities Spain
# library(mgcv)

rm(list = ls())

# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.heavy <- paste0(getwd(), "/Heavy_files/")
loc.clc <- paste0(getwd(), "/OUTPUT/land_cover_raster/")

# Loading data -----------------------------------------------------------------
all <- readRDS(file = paste0(loc.output, "all_monthly_pred_mun.rds"))
all <- all %>% filter(m == c("07", "08", "09"))

# Fit the gam curve ------------------------------------------------------------
# We fit the same curve that ggplot2 returns
# fit_gam <- mgcv::gam(bg ~ s(ma, bs = "cs"), data = all)
fit_pol <- lm(bg ~ poly(ma, 2), data = all)
summary(fit_pol)

# residuals = observed - fitted values
all$fitted_values <- predict(fit_pol)
all$residuals <- all$bg - all$fitted_values

# The same plot as ggplot2
ggplot(all) +
  geom_point(aes(x = ma, y = bg), alpha = 0.6, size = 0.5) +
  geom_smooth(aes(x = ma, y = bg), linetype = "dashed", color = "purple") +
  # geom_smooth(aes(x = ma, y = bg), method = "lm", formula = y ~ poly(x, 2), linetype = "dashed", color = "purple") +
  labs(
    x = "Mosquito Alert predicted probabilities",
    y = "Count predicted values"
  ) +
  theme_classic()

# Filtering by the critics months
ggplot(all) +
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

all_2020 <- merge(all %>% filter(year == 2020), clc, by = "id", all.x = TRUE)
clc_names <- names(all_2020)[10:23]

clc_residuals_plots <- list()
for (i in 1:length(clc_names)){
  print(paste(clc_names[i]))
  clc_residuals_plots[[i]] <- ggplot(all_2020, aes_string(x = "residuals", y = clc_names[i])) +
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

# Joining human density --------------------------------------------------------
# Read File with the number of people per year
pop_full <- data.frame()
for(n in c("20", "21", "22", "23")){
  pop <- read.csv(paste0(loc.data, "pobmun", n, ".csv"))
  
  pop$cmun <- ifelse(pop$CMUN<10, paste0("00",pop$CMUN),
                     ifelse(pop$CMUN<100, paste0("0",pop$CMUN),as.character(pop$CMUN)))
  pop$pop <- as.numeric(pop$POB) # This variable name, number of inhabitants, might change depending on the file. Check the name before should be POB or POB22 !
  pop$cpro <- ifelse(pop$CPRO < 10, paste0("0",pop$CPRO), pop$CPRO)
  
  pop <- pop %>% dplyr::select(cpro, cmun, pop)
  
  esp_can <- esp_get_munic_siane() %>% janitor::clean_names() %>%
    mutate(area = as.double(st_area(st_transform(., 3035)) / 1000000))
  esp_can$id <- as.factor(paste0(esp_can$codauto,
                                 esp_can$cpro,
                                 esp_can$cmun,
                                 esp_can$lau_code)) 
  esp_can <- esp_can %>% dplyr::select(id, cpro, cmun, area)
  esp_can <- merge(esp_can, pop,  by = c("cpro","cmun"), all.x = TRUE) # Join population density with esp_can to obtain NATCODE
  esp_can <- esp_can %>% mutate(dens = pop/area)
  
  esp_can <- esp_can %>% dplyr::select(id, dens, pop) %>% st_drop_geometry() %>%
    mutate(y = as.factor(paste0("20", n)))
  
  pop_full <- rbind(pop_full, esp_can)
}

all <- merge(all %>% mutate(y = year), pop_full, by = c("id", "y"), all.x = TRUE)

ggplot(all, aes_string(x = "residuals", y = "dens")) +
  geom_point(alpha = 0.5, color = "#3288bd") +
  geom_smooth(method = "lm", color = "#d53e4f", se = FALSE) +
  labs(
    x = "Residuals",
    y = "Human density"
  ) +
  theme_classic()

ggplot(all, aes_string(x = "residuals", y = "pop")) +
  geom_point(alpha = 0.5, color = "#3288bd") +
  geom_smooth(method = "lm", color = "#d53e4f", se = FALSE) +
  labs(
    x = "Residuals",
    y = "Human density"
  ) +
  theme_classic()

# Joining temperature and precipitation ----------------------------------------
all_2020 <- all %>% filter(year == "2020")

all_data <- mclapply(1:nrow(all_2020), function(i){
  
  cat(paste0("Number row:", i, "\n"))
  
  data_row <- all_2020[i, ]
  
  wth <- readRDS(paste0(loc.output, "monthly_weather_data/prep_", 
                        data_row$m, "-",data_row$year, ".rds"))[, c(2,9:12)]
  
  data_row <- merge(data_row, wth, by = "id", all.x = TRUE)
}
, mc.cores = 8)

all_data <- do.call(rbind, all_data) 

ggplot(all_data, aes_string(x = "residuals", y = "mean_temperature")) +
  geom_point(alpha = 0.5, color = "#3288bd") +
  geom_smooth(method = "lm", color = "#d53e4f", se = FALSE) +
  labs(
    x = "Residuals",
    y = "Mean temperature"
  ) +
  theme_classic()

ggplot(all_data, aes_string(x = "residuals", y = "min_temperature")) +
  geom_point(alpha = 0.5, color = "#3288bd") +
  geom_smooth(method = "lm", color = "#d53e4f", se = FALSE) +
  labs(
    x = "Residuals",
    y = "Min temperature"
  ) +
  theme_classic()

ggplot(all_data, aes_string(x = "residuals", y = "mean_relative_humidity")) +
  geom_point(alpha = 0.5, color = "#3288bd") +
  geom_smooth(method = "lm", color = "#d53e4f", se = FALSE) +
  labs(
    x = "Residuals",
    y = "Mean Relative humidity"
  ) +
  theme_classic()

ggplot(all_data, aes_string(x = "residuals", y = "precipitation")) +
  geom_point(alpha = 0.5, color = "#3288bd") +
  geom_smooth(method = "lm", color = "#d53e4f", se = FALSE) +
  labs(
    x = "Residuals",
    y = "Mean temperature"
  ) +
  theme_classic()

ggplot(all_data, aes(x = abs(residuals), y = mean_relative_humidity)) +
  geom_point(alpha = 0.5, color = "#3288bd") +
  geom_smooth(method = "lm", color = "#d53e4f", se = FALSE) +
  labs(
    x = "Residuals",
    y = "Mean Relative humidity"
  ) +
  theme_classic()

ggplot(all_data) +
  geom_point(aes(x = ma, y = bg, color = mean_relative_humidity), alpha = 0.6, size = 0.5) +
  # geom_smooth(aes(x = ma, y = bg), method = "lm", formula = y ~ poly(x, 2), linetype = "dashed", color = "purple") +
  labs(
    x = "Mosquito Alert predicted probabilities",
    y = "Count predicted values"
  ) +
  theme_classic()
