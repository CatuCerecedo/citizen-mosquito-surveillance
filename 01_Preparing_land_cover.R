########################## Preparing land Cover ################################
#' Percentages by grid

library(tidyverse)
library(sf)
library(janitor)
library(parallel)
library(dplyr)
library(lubridate)
library(terra)
library(data.table)
library(exactextractr)

rm(list = ls())

# Directories ------------------------------------------------------------------
# In local
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.era5 <- "/home/catuxa/Documents/Mosquito_Models/EU_Culex/ERA5_Download/"
loc.clc <- "/home/catuxa/Documents/Mosquito_Models/EU_Culex/DATA/u2018_clc2018_v2020_20u1_raster100m/"

# In cluster
loc.output <- paste0(getwd(), "/Spain_Tiger/OUTPUT/")
loc.data <- paste0(getwd(), "/Spain_Tiger/DATA/")
loc.era5 <- paste0(getwd(), "/EU_Culex/ERA5_Download/")
loc.clc <- paste0(getwd(), "/EU_Culex/DATA/u2018_clc2018_v2020_20u1_raster100m/")


sf::sf_use_s2(FALSE)
# Calculate clc by municipalities ----------------------------------------------
# Municipality map
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  )

# clc raster
path <- paste0(loc.clc, "DATA/U2018_CLC2018_V2020_20u1.tif")
landcover <- rast(path)
landcover <- project(landcover, crs(spain))
# plot(landcover) # Plot the raster

# loop to calculate clc percentages
results_list <- list()
for (i in 1:nrow(spain)) {
  
  mun <- spain[i, ]
  
  if (i %% 100 == 0){
   cat("row number: ", i, "/n") 
  }
  
  # Crop the raster by each municipality
  landcover_crop <- crop(landcover, mun)
  landcover_mask <- mask(landcover_crop, mun)
  
  # Calculate the surface of each category
  freq_table <- as.data.frame(landcover_mask) %>%
    group_by(LABEL3) %>%
    summarise(n = n())
  
  # Percentage
  freq_table$per_clc <- (freq_table$n / sum(freq_table$n))
  
  # adding municipality
  freq_table$municipality <- mun$municipality
  freq_table$id <- mun$id
  
  results_list[[i]] <- freq_table
}
# saveRDS(results_list, file = paste0(loc.output, "mun_clc_original_list.rds"))

results_df <- bind_rows(results_list)

landcover_crop <- crop(landcover, spain)
values_id <- levels(landcover_crop)[[1]]

results_df <- merge(results_df, values_id, by = "LABEL3") 
results_df <- results_df %>%
  rename(
    landcover = "Value",
    value = "LABEL3"
  ) %>% mutate(
    value = as.numeric(value)
  )

clc_surface <- results_df %>%
  mutate(
    value = case_when(value == 1 ~ "cont_urban_fabric",
                      value ==2 ~ "discont_urban_fabric",
                      value == 4 ~ "roads_rails",
                      value == 10 ~ "green_urban",
                      value == 11 ~ "sports_leisure",
                      value %in% c(3, 5:9) ~ "other_artificial",
                      value %in% 12:22 ~ "agricultural",
                      value %in% 23:29 ~ "forests_scrub",
                      value %in% 30:34 ~ "open",
                      value %in% 35:36 ~ "inland_wetlands",
                      value %in% 37:39 ~ "marine_wetlands",
                      value %in% 40:41 ~ "inland_water",
                      value %in% 42:44 ~ "marine_water",
                      value >= 45 ~ "no_data")
  ) %>%
  units::drop_units() %>%
  as.data.table()

clc_surface <- clc_surface %>%
  group_by(value, municipality, id) %>%
  summarise(
    per_clc = sum(per_clc, na.rm = TRUE)
  ) %>% 
  pivot_wider(
    names_from = value,
    values_from = per_clc
      )

clc_surface[is.na(clc_surface)] <- 0

saveRDS(clc_surface, file = paste0(loc.output, "clc_surface_mun.rds"))
