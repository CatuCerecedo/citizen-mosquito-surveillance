###################### Spatial analysis:logistic interactions ##################
library(ggplot2)
library(patchwork)
library(tidyverse)
library(sf)
library(parallel)

rm(list = ls())
# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")
loc.heavy <- paste0(getwd(), "/Heavy_files/")
loc.clc <- paste0(getwd(), "/OUTPUT/land_cover_raster/")

sf::sf_use_s2(FALSE)
# Load data --------------------------------------------------------------------
# Country
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  )

# Mosquitoes
mdl <- "mtiger3_occu"
mdl_name <- "/mtiger3_occu"
fldr <- "Suitability"
sub <- "_occu" # with _

months = c("01", "02", "03", "07", "08", "09", "04", "05", "06", "10", "11", "12")

y = "2020"

calc_pred <- function(){for(i in 1:length(months)){
  print(i)
  if (i == 1){
    pred <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                           "/tiger_", months[i], "_", y, sub, ".rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    
    pred$count <- rowMeans(pred[,6:ncol(pred)], na.rm = TRUE)
    
    pred <- pred[c("municipality","id", "prov_name", "lon", "lat", "count")]
  } else {
    p <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                        "/tiger_", months[i], "_", y, sub, ".rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    
    p$count <- rowMeans(p[,6:ncol(p)], na.rm = TRUE)
    
    p <- p[c("municipality","id", "prov_name", "lon", "lat", "count")]
    
    pred <- cbind(pred, p["count"])
  }
}
  return(pred)
}
pred_suit <- calc_pred()
pred_suit$pred_count <- rowMeans(pred_suit[,6:ncol(pred_suit)], na.rm = TRUE)
pred_suit <- pred_suit[c("municipality","id", "prov_name", "lon", "lat", "pred_count")]

mdl <- "mtiger16"
mdl_name <- "/mtiger16"
fldr <- "Counts"
sub <- "" # with _

pred_count <- calc_pred()
pred_count$pred_count <- rowMeans(pred_count[,6:ncol(pred_count)], na.rm = TRUE)
pred_count <- pred_count[c("municipality","id", "prov_name", "lon", "lat", "pred_count")]

mdl <- "mbites_7_reduce"
mdl_name <- "/mbites_7_reduce"
fldr <- "Bites"
sub <- "bts" # with _

pred_bite <- calc_pred()
pred_bite$pred_count <- rowMeans(pred_bite[,6:ncol(pred_bite)], na.rm = TRUE)
pred_bite <- pred_bite[c("municipality","id", "prov_name", "lon", "lat", "pred_count")]

mdl <- "mtiger7_ma"
mdl_name <- "/mtiger7_ma"
fldr <- "MA"
sub <- "_ma" # with _

pred_ma <- calc_pred()
pred_ma$pred_count <- rowMeans(pred_ma[,6:ncol(pred_ma)], na.rm = TRUE)
pred_ma <- pred_ma[c("municipality","id", "prov_name", "lon", "lat", "pred_count")]

rank_order_conversion <- function(vctr) {
  rank = rank(vctr, na.last = "keep", ties.method = "first")       
  rank_normalized = rank / max(rank, na.rm = TRUE) 
  return(rank_normalized)
}

pred_suit$pred_count_rnk <- rank_order_conversion(pred_suit$pred_count)
pred_count$pred_count_rnk <- rank_order_conversion(pred_count$pred_count)
pred_bite$pred_count_rnk <- rank_order_conversion(pred_bite$pred_count)
pred_ma$pred_count_rnk <- rank_order_conversion(pred_ma$pred_count)

# Building the data set
pred_suit <- pred_suit %>% mutate(
  Model_comp = "SM",
  diff = pred_ma$pred_count_rnk - pred_count_rnk
)

pred_count <- pred_count %>% mutate(
  Model_comp = "CM",
  diff = pred_ma$pred_count_rnk - pred_count_rnk
)

pred_bite <- pred_bite %>% mutate(
  Model_comp = "BM",
  diff = pred_ma$pred_count_rnk - pred_count_rnk
)

pred <- rbind(pred_suit, pred_count, pred_bite) %>%
  mutate(
    CSM_great = case_when(diff > 0 ~ 1, diff <= 0 ~ 0)
  )

spain_wth <- mclapply(1:nrow(spain), function(i){
  
  cat(paste0("Number row:", i, "\n"))
  
  data_row <- spain[i, ]
  
  data_point <- st_centroid(data_row)
  
  for (m in 1:length(months)){
    wth <- readRDS(paste0(loc.output, "monthly_weather_data/prep_", 
                          months[m], "-", y, ".rds"))[, c(1,9:12, 32)]
    data_row <- merge(data_row, wth, by = "id", all.x = TRUE)
  }
  
  data_row <- data_row %>%
    mutate( 
      min_temperature = rowMeans(across(starts_with("min_temperature")), na.rm = TRUE),
      mean_temperature = rowMeans(across(starts_with("mean_temperature")), na.rm = TRUE),
      max_temperature = rowMeans(across(starts_with("max_temperature")), na.rm = TRUE),
      precipitation = rowMeans(across(starts_with("precipitation")), na.rm = TRUE),
      mean_relative_humidity = rowMeans(across(starts_with("mean_relative_humidity")), na.rm = TRUE) 
    ) %>% 
    dplyr::select(id, municipality, prov_name, min_temperature, max_temperature,
                  mean_temperature, precipitation, mean_relative_humidity)
}
, mc.cores = 8)
spain_wth <- do.call(rbind, spain_wth)
saveRDS(spain_wth, paste0(loc.output, "wth_mun_2020.rds"))

cluster_wth <- merge(pred, 
                     spain_wth %>% st_drop_geometry(),
                     by = c("id", "municipality", "prov_name"))
cluster_wth <- cluster_wth %>% mutate(
  ext_t = case_when(min_temperature < 13 & max_temperature > 31 ~ TRUE,
                     .default = FALSE)
)

custom_colors <- c("SM" = "#abdda4", "CM" = "#d53e4f", "BM" = "#3288bd")

ggplot(cluster_wth, aes(Model_comp, diff, color = Model_comp)) + 
  geom_boxplot() 

logistic_model <- glm(CSM_great ~ Model_comp*min_temperature, 
                      data = cluster_wth, 
                      family = binomial(link = "logit"))
summary(logistic_model)
# e <- predict_response(logistic_model, terms = c("min_temperature", "Model_comp"))
# a <- plot(e, show_ci = TRUE) +
#   scale_color_manual("Model\ncomparison", values = custom_colors) +
#   scale_fill_manual("Model\ncomparison", values = custom_colors) +
#   labs(title = "",
#        y = "Pr(MA > other model prediction)",
#        x = "Min Temperature (ºC)") +
#   theme_classic()
e <- effects::effect("Model_comp*min_temperature", logistic_model) %>% as.data.frame()
a <- ggplot(e, aes(min_temperature, fit, color=Model_comp, group = Model_comp)) + 
  scale_color_manual("Model\ncomparison", values = custom_colors) +
  geom_point(size = 1.5) + 
  geom_line(linewidth = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.01) +
  labs(x = "Min Temperature (ºC)", y = "Pr(MA > other model prediction)") + 
  theme_classic()

logistic_model <- glm(CSM_great ~ Model_comp*max_temperature, 
                      data = cluster_wth, 
                      family = binomial(link = "logit"))
summary(logistic_model)
# e <- predict_response(logistic_model, terms = c("max_temperature", "Model_comp"))
# b <- plot(e, show_ci = TRUE) +
#   scale_color_manual("Model\ncomparison", values = custom_colors) +
#   scale_fill_manual("Model\ncomparison", values = custom_colors) +
#   labs(title = "",
#        y = "Pr(MA > other model prediction)",
#        x = "Max Temperature (ºC)") +
#   theme_classic()
e <- effects::effect("Model_comp*max_temperature", logistic_model) %>% as.data.frame()
b <- ggplot(e, aes(max_temperature, fit, color=Model_comp, group = Model_comp)) + 
  scale_color_manual("Model\ncomparison", values = custom_colors) +
  geom_point(size = 1.5) + 
  geom_line(linewidth = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.01) +
  labs(x = "Max Temperature (ºC)", y = "Pr(MA > other model prediction)") + 
  theme_classic()

logistic_model <- glm(CSM_great ~ Model_comp*mean_temperature, 
                      data = cluster_wth, 
                      family = binomial(link = "logit"))
summary(logistic_model)
# e <- predict_response(logistic_model, terms = c("mean_temperature", "Model_comp"))
# c <- plot(e, show_ci = TRUE) +
#   scale_color_manual("Model\ncomparison", values = custom_colors) +
#   scale_fill_manual("Model\ncomparison", values = custom_colors) +
#   labs(title = "",
#        y = "Pr(MA > other model prediction)",
#        x = "Mean Temperature (ºC)") +
#   theme_classic()
e <- effects::effect("Model_comp*mean_temperature", logistic_model) %>% as.data.frame()
c <- ggplot(e, aes(mean_temperature, fit, color=Model_comp, group = Model_comp)) + 
  scale_color_manual("Model\ncomparison", values = custom_colors) +
  geom_point(size = 1.5) + 
  geom_line(linewidth = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.01) +
  labs(x = "Mean Temperature (ºC)", y = "Pr(MA > other model prediction)") + 
  theme_classic()

logistic_model <- glm(CSM_great ~ Model_comp*mean_relative_humidity, 
                      data = cluster_wth, 
                      family = binomial(link = "logit"))
summary(logistic_model)

# e <- predict_response(logistic_model, terms = c("mean_relative_humidity", "Model_comp"))
# d <- plot(e, show_ci = TRUE) +
#   scale_color_manual("Model\ncomparison", values = custom_colors) +
#   scale_fill_manual("Model\ncomparison", values = custom_colors) +
#   labs(title = "",
#        y = "Pr(MA > other model prediction)",
#        x = "Mean Relative Humidity (%)") +
#   theme_classic()
e <- effects::effect("Model_comp*mean_relative_humidity", logistic_model) %>% as.data.frame()
d <- ggplot(e, aes(mean_relative_humidity, fit, color=Model_comp, group = Model_comp)) + 
  scale_color_manual("Model\ncomparison", values = custom_colors) +
  geom_point(size = 1.5) + 
  geom_line(linewidth = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.01) +
  labs(x = "Mean Relative Humidity (%)", y = "Pr(MA > other model prediction)") + 
  theme_classic() 

logistic_model <- glm(CSM_great ~ Model_comp*precipitation, 
                      data = cluster_wth, 
                      family = binomial(link = "logit"))
summary(logistic_model)
# e <- predict_response(logistic_model, terms = c("precipitation", "Model_comp"))
# f <- plot(e, show_ci = TRUE) +
#   scale_color_manual("Model\ncomparison", values = custom_colors) +
#   scale_fill_manual("Model\ncomparison", values = custom_colors) +
#   labs(title = "",
#        y = "Pr(MA > other model prediction)",
#        x = "Precipitation (m)") +
#   theme_classic()
e <- effects::effect("Model_comp*precipitation", logistic_model) %>% as.data.frame()
f <- ggplot(e, aes(precipitation, fit, color = Model_comp, group = Model_comp)) + 
  scale_color_manual("Model\ncomparison", values = custom_colors) +
  geom_point(size = 1.5) + 
  geom_line(linewidth = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.01) +
  labs(x = "Precipitation (m)", y = "Pr(MA > other model prediction)") + 
  theme_classic()

a + b + c + d + f +
  plot_annotation(tag_levels = c("A"), tag_suffix = ")") 

ggsave(paste0(loc.fig, "Spatial_clusters/glm/interactions/mar_effect_pred2.png"),
       width = 35,  height = 22, units = "cm")
