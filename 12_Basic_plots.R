############################### Basic plots ####################################
library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(tidyverse)
library(parallel)
library(patchwork)

rm(list = ls())

# Directories ------------------------------------------------------------------
# In local
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")
loc.era5 <- paste0("/home/catuxa/Documents/Mosquito_Models/EU_Culex/ERA5_Download/")

sf::sf_use_s2(FALSE)
# Loading spain municipality map -----------------------------------------------
# Country
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  )

# Plot of samples location (traps and reports) ---------------------------------
tiger <- readRDS(paste0(loc.output, "bg_tiger_spain_daily.rds")) 
  # filter(females > 0) %>%
  dplyr::select(id, municipality, prov_name) %>%
  distinct()
tiger <- merge(tiger, spain, by = c("id", "municipality", "prov_name"))
st_geometry(tiger) <- "geometry"

ma_df <- readRDS(paste0(loc.output, "ma_tiger_spain_daily.rds")) %>%
  # filter(any_reps == TRUE) %>%
  dplyr::select(id, municipality) %>%
  distinct()
ma_df <- merge(ma_df, spain, by = c("id", "municipality"))
st_geometry(ma_df) <- "geometry"

can_prov <- mapSpain::esp_get_prov() %>% 
  filter(nuts2.name != "Canarias") %>%
  st_transform(4236)

a <- ggplot() +
  geom_sf(data = ma_df, fill = "#d53e4f", color = "white", alpha = 0.4) +
  geom_sf(data = tiger, fill = "black", color = "white") +
  geom_sf(data = can_prov, fill = "transparent", color = "black") +
  labs(
    color = "Province",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_classic() 

a + b + 
  plot_annotation(tag_levels = c("A"), tag_suffix = ")") 

ggsave(paste0(loc.fig, "study_area_samples.png"),
       width = 30,  height = 30, units = "cm")

# Temporal resolution of surveillance ------------------------------------------
tiger <- readRDS(paste0(loc.output, "bg_tiger_spain_daily.rds")) %>%
  filter(y != 2023) %>%
  filter(females > 0) %>%
  mutate(m = month(end_date), method = "TRAPS") %>%
  dplyr::select(y, m, method) %>%
  distinct()

ma_df <- readRDS(paste0(loc.output, "ma_tiger_spain_daily.rds")) %>%
  filter(any_reps == TRUE) %>%
  mutate(m = month(date), method = "CITSI") %>%
  dplyr::select(y, m, method) %>%
  distinct()

merged_df <- full_join(tiger, ma_df, by = c("y", "m")) %>%
  mutate(
    method = case_when(
      !is.na(method.x) & !is.na(method.y) ~ "BOTH", 
      !is.na(method.x) ~ "TRAPS",  
      !is.na(method.y) ~ "CITSI"  
    )
  ) %>%
  dplyr::select(y, m, method) %>%
  distinct() 

ggplot(merged_df, aes(x = factor(m), y = factor(y), fill = method)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("CITSI" = "#d53e4f", "TRAPS" = "#abdda4", "BOTH" = "#3288bd")) +
  labs(x = "Months",
       y = "Year",
       fill = "Surveillance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(loc.fig, "temporal_coverage_samples.png"),
       width = 20,  height = 10, units = "cm")


# Plotting weather variables with raw data -------------------------------------
# Count data
tiger <- readRDS(paste0(loc.output, "bg_tiger_spain_daily.rds"))

# Mosquito Alert data
ma_df <- readRDS(paste0(loc.output, "ma_tiger_spain_daily.rds")) %>%
  mutate(
    id = as.factor(id),
    y = as.factor(year(date))
  ) %>%
  filter(SE > 0) %>%
  drop_na()

tiger <- tiger %>% dplyr::select(females, id, y, municipality, min_temperature, 
                                 max_temperature, mean_temperature, precipitation,
                                 mean_relative_humidity) %>%
  mutate(
    presence = ifelse(females > 0, 1, 0),
    method = "COUNT"
  ) %>% 
  dplyr::select(-females)
ma_df <- ma_df %>% dplyr::select(any_reps, id, y, municipality, min_temperature, 
                                 max_temperature, mean_temperature, precipitation,
                                 mean_relative_humidity) %>%
  mutate(
    presence = ifelse(any_reps == TRUE, 1, 0),
    method = "CITSI"
  ) %>% 
  dplyr::select(-any_reps)
combine <- bind_rows(tiger, ma_df) %>%
  filter(presence == 1)

climate_vars <- c("min_temperature", "max_temperature", "mean_temperature", 
                  "precipitation", "mean_relative_humidity")

combine_long <- combine %>%
  pivot_longer(cols = c("min_temperature", "max_temperature", "mean_temperature", 
                        "precipitation", "mean_relative_humidity"), 
               names_to = "climate_var", values_to = "value") 
custom_labels <- c(
  min_temperature = "Minimum Temperature (ºC)",
  max_temperature = "Maximum Temperature (ºC)",
  mean_temperature = "Mean Temperature (ºC)",
  precipitation = "Precipitation (m)",
  mean_relative_humidity = "Mean Relative Humidity (%)"
)

ggplot(combine_long, aes(x = value, y = method, fill = method)) +
  scale_fill_manual(values = c("COUNT" = "#abdda4", "CITSI" = "#d53e4f")) +
  ggridges::geom_density_ridges(alpha = 0.5, scale = 1.2, rel_min_height = 0.01) +  
  geom_vline(data = combine_long %>% 
               group_by(climate_var) %>% 
               summarise(mean_val = mean(value, na.rm = TRUE)), 
             aes(xintercept = mean_val), linetype = "dashed", size = 1) +
  facet_wrap(~climate_var, scales = "free_x",
             labeller = as_labeller(custom_labels)) +  
  theme_classic() +
  labs(
    x = "\nClimatic Variables",
    y = "Method\n"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold")
  )
ggsave(file = paste0(loc.fig, "Spatial_clusters/glm/weather_dist_raw_data.png"), 
       units = "cm", height = 20, width = 35, bg = "white")

# Plotting spatial correlation along variables ---------------------------------
# Yearly correlations between BG and MA
cor_predictions <- readRDS(file = paste0(loc.output, "temporal_corr.rds"))

# # Load max temperature values
# years <- c("2020", "2021", "2022")
# months = c("01", "02", "03", "07", "08", "09", "04", "05", "06", "10", "11", "12")
# 
# y = "2022"
# ncname <- paste0(loc.era5, "ERA5_EU_hourly_", y)
# ncfname <- paste(ncname, ".nc", sep="")
# 
# tmp_raster <- terra::rast(ncfname)["t2m"]
# tmp_raster <- tapp(tmp_raster, index = "months", fun = "max")
# 
# vls_max_temp <- terra::extract(tmp_raster, st_centroid(spain))[,-1] - 273.15
# vls_max_temp$id <- spain$id
# 
# for (i in 1:length(months)){
#   wth <- readRDS(paste0(loc.output, "monthly_weather_data/prep_",
#                         months[i], "-", y, ".rds"))
#   vls_max <- vls_max_temp[, c(i, 13)]
#   colnames(vls_max) <- c("max_temperature", "id")
#   wth <- merge(wth, vls_max, by = "id", all.x = TRUE)
#   saveRDS(wth, paste0(loc.output, "monthly_weather_data/prep_",
#                  months[i], "-", y, ".rds"))
# }
months = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
years = c("2020", "2021", "2022")

cor_table_wth <- data.frame()
for(y in years){
  print(y)
  for (i in 1:length(months)){
    print(months[i])
    wth <- readRDS(paste0(loc.output, "monthly_weather_data/prep_", 
                          months[i], "-", y, ".rds"))[, c(1, 9:12, 32)]
    wth$year <- y
    wth$month <- as.numeric(months[i])
    
    wth <- wth %>% group_by(year, month) %>% summarise(
      min_temperature = mean(min_temperature, na.rm = TRUE),
      mean_temperature = mean(mean_temperature, na.rm = TRUE),
      precipitation = sum(precipitation, na.rm = TRUE),
      mean_relative_humidity = mean(mean_relative_humidity, na.rm = TRUE),
      max_temperature = mean(max_temperature, na.rm = TRUE)
    ) %>% as.data.frame()
    
    cor_sample <- cor_predictions %>%
      filter(year == y & month == as.numeric(months[i]))
    cor_sample <- merge(cor_sample, wth, by = c("year", "month"),
                               all.x = TRUE)
    cor_table_wth <- rbind(cor_table_wth, cor_sample)
  }
}
rm(y, i, cor_sample)

cor_table_wth <- cor_table_wth %>% 
  filter(comparison == c("suit_ma", "count_ma"))
custom_colors <- c("suit_ma" = "#abdda4", "count_ma" = "#d53e4f") #, "bites_ma" = "#3288bd")

vars = c("min_temperature", "max_temperature", 
         "mean_temperature", "precipitation", "mean_relative_humidity")
nms = c("Min Temperature (ºC)", "Max Temperature (ºC)", "Mean Temperature (ºC)",
        "Precipitation (m", "Mean Relative Humidity (%)")

plot_wth <- list()
for(i in 1:length(vars)){
  p <- ggplot(cor_table_wth, aes_string(x = "month", y = vars[i])) +
    geom_col(size = 7, alpha = 0.6) +
    # geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = custom_colors) +
    labs(
      y = nms[i],
      x = "Month",
      color = "Comparison",      # Aquí se especifica "Comparison" para la leyenda unificada
      shape = "Comparison",      # Se especifica también para shape
      linetype = "Comparison",   # Y para linetype
      # title = "Spanish Tiger Mosquito Models: Spatial correlation over time"
    ) +
    scale_x_continuous(breaks = seq(1, 12, 1)) +
    geom_point(aes_string(y = "value*100", color = "comparison"), size = 7, alpha = 0.6) +
    geom_line(aes_string(y = "value*100", color = "comparison"), size = 2, alpha = 0.6) +
    scale_y_continuous(
      sec.axis = sec_axis(~./100, name = "Spearman Correlation (S)")
    ) +
    theme_classic() +
    facet_wrap(~year) +
    theme(
      legend.position = "bottom",    # La leyenda se posiciona abajo
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 20),
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
    ) +
    guides(                       # Unificación de la leyenda
      color = guide_legend(title = "Comparison"),
      shape = guide_legend(title = "Comparison"),
      linetype = guide_legend(title = "Comparison")
    )
  plot_wth[[i]] <- p
}

ggpubr::ggarrange(plot_wth[[1]] ,
plot_wth[[2]],
plot_wth[[3]],
plot_wth[[4]],
plot_wth[[5]], nrow = 5)

ggsave(paste0(loc.fig, "Spatial_clusters/glm/corr_wth.png"),
       width = 40,  height = 60, units = "cm")

# Plotting weather on decorrelation period -------------------------------------
decor_time <- list()
decor_time[[1]] <-  c("2020", "2021", "2022")
decor_time[[2]] <- list(list("07", "08"), list("08", "09"), list("07", "08"))

wth_table <-data.frame()
for (i in 1:3){
  print(decor_time[[1]][i])
  
  for (j in decor_time[[2]][[i]]){
    print(j)
    
    wth <- readRDS(paste0(loc.output, "monthly_weather_data/prep_", 
                         j, "-", decor_time[[1]][i], ".rds"))[, c(1, 9:12, 32)]
    wth$year <- decor_time[[1]][i]
    wth$month <- j

    wth_table <- rbind(wth_table, wth)
  }
}
  
wth_table <- wth_table %>% 
  group_by(id) %>%
  summarise(
    min_temperature = mean(min_temperature, na.rm = TRUE),
    max_temperature = mean(max_temperature, na.rm = TRUE),
    mean_temperature = mean(mean_temperature, na.rm = TRUE),
    precipitation = sum(precipitation, na.rm = TRUE),
    mean_relative_humidity = mean(mean_relative_humidity, na.rm = TRUE)
  )
wth_table <- merge(spain, wth_table, by = "id")

variables <- c("min_temperature", "max_temperature", "mean_temperature", 
               "precipitation", "mean_relative_humidity")
nms = c("Min Temperature (ºC)", "Max Temperature (ºC)", "Mean Temperature (ºC)",
        "Precipitation (m)", "Mean Relative Humidity (%)")
wth_plot <- list()
for (i  in 1:length(variables)){
  wth_plot[[i]] <- ggplot() +
    geom_sf(data = wth_table, aes_string(fill = variables[i]), color = "transparent") +
    scale_fill_distiller(nms[i], palette = "Spectral") +
    theme_classic()
}

wth_plot[[1]] + wth_plot[[2]] + wth_plot[[3]] + wth_plot[[4]] + wth_plot[[5]] +
  plot_annotation(tag_levels = c("A"), tag_suffix = ")") 
ggsave(paste0(loc.fig, "Spatial_clusters/glm/weather_suring_secorrelation.png"), units = "cm", bg = "white",
       height = 45, width = 35)

# Plot MRH vs temperature ------------------------------------------------------
# La MRH aumenta con la cantidad de vapor de agua con la temperatura es constante.
# La MRH aumenta con la temperatura cuando hay más cantidad de aire
# Qué pasa en la cornisa cantábrica (temperaturas bajas y mucha vapor de agua), las
# temperaturas bajas disminuyen los valores de MRH, pero como el vapor de agua es 
# mayor en esta zona, los valores de MRH se compensan.
# Las zonas más humedas en España son las costas, espeialmente, en el norte de la península.
months = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
years = c("2020", "2021", "2022")

wth_plot <- list()
for (y in 1:3){
  wth_year <- wth["id"]
  
  for (m in 1:12){
    wth <- readRDS(paste0(loc.output, "monthly_weather_data/prep_", 
                          months[m], "-", years[y], ".rds"))
    colnames(wth)[colnames(wth) == "mean_temperature"] <- paste0("mean_temperature_", m)
    colnames(wth)[colnames(wth) == "mean_relative_humidity"] <- paste0("mean_relative_humidity_", m)
    wth_year <- bind_cols(wth_year, wth[c(paste0("mean_temperature_", m), 
                                          paste0("mean_relative_humidity_", m))])
  }
  wth_year <- wth_year %>%
    mutate(
      mean_temperature = rowMeans(across(starts_with("mean_temperature"))),
      mean_relative_humidity = rowMeans(across(starts_with("mean_relative_humidity")))
    ) %>% 
    dplyr::select(id, mean_relative_humidity, mean_temperature)
  
  wth_plot[[y]] <- ggplot(wth_year, aes(y = mean_relative_humidity, x = mean_temperature)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_classic() +
    labs(
      x = "Mean Temperature (ºC)",
      y = "Mean Relative Humidity (%)"
    )
}
wth_plot[[1]] +  wth_plot[[2]] + wth_plot[[3]] +
  plot_annotation(tag_levels = "A", tag_suffix = ")")

ggsave(paste0(loc.fig, "Spatial_clusters/glm/annual_temp_vs_humidity.png"),
       width = 25,  height = 10, units = "cm")
