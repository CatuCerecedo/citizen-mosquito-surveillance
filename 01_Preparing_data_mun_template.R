###################### Preparing Mosquito Alert Data ###########################

library(tidyverse)
library(sf)
library(terra)
library(parallel)
library(data.table)

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

sf::sf_use_s2(FALSE)
# Loading mosquito data --------------------------------------------------------
tiger <- readRDS(paste0(loc.output, "bg_traps_spain_2020_2022.rds"))
tiger <- st_as_sf(tiger, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
# names(tiger)

# Loading Spain map ------------------------------------------------------------
# spain <- mapSpain::esp_get_munic_siane() %>%
#   st_transform(4326) %>%
#   filter(!(ine.ccaa.name %in% c("Canarias", "Ceuta", "Melilla"))) %>%
#   mutate(
#     id = paste0(codauto, cpro, cmun, LAU_CODE)
#   ) %>%
#   ungroup() %>%
#   dplyr::select(name, id, ine.prov.name) %>%
#   rename(municipality = "name", prov_name = "ine.prov.name")
# saveRDS(spain, file = paste0(loc.output, "spain_mun.rds"))
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  )

# Checking all municipalities
ggplot() +
  geom_sf(data = spain, fill = "black", color = "transparent",
          size = 0.01, alpha = 0.8, na.rm = TRUE) 

# We have to add a new random effect: the municipality!
tiger <- st_join(tiger, spain, join = st_intersects) %>% st_drop_geometry()

# Grouping counts by municipality and month
tiger <- tiger %>%
  mutate(
    m = as.integer(month(end_date)), 
    y = year(end_date)
  ) %>% 
  group_by(y, m, municipality, id) %>% 
  summarize(
    females = sum(females, na.rm = TRUE),
    trapping_effort = sum(trapping_effort, na.rm = TRUE),
    n_traps = n(),
    .groups ="drop") %>%
  ungroup() # 256 rows
sum(duplicated(tiger))

plot(table(tiger$females)) # Negative binomial

tiger <- merge(tiger, spain, by = c("municipality", "id")) # 256
st_geometry(tiger) <- "geometry"

# Extracting climatic variables from grib raster -------------------------------
# centroids <- st_centroid(tiger) %>% 
#   mutate(longitude = st_coordinates(.)[,1], latitude = st_coordinates(.)[,2]) %>% 
#   st_drop_geometry() %>%
#   st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
# # plot(st_geometry(centoids))

# Full function ----------------------------------------------------------------
# Loading the functions
# In local
# source("read_and_extract_grib_data.R")
# In cluster
# source("/home/usuaris/ccerecedo/EU_Culex/read_and_extract_grib_data.R")

# Loading the functions
# In local
source("read_and_extract_nc_data.R")
# In cluster
# source("/home/usuaris/ccerecedo/EU_Culex/read_and_extract_nc_data.R")
cat("------------------------- Number of rows:", nrow(tiger), "\n")

weather_data <- mclapply(1:nrow(tiger), function(i){
  
  cat(paste0("Number row:", i, "\n"))
  
  data_row <- tiger[i, ]
  
  ex_wt <- extract_weather(data_row)
  
  wth <- calculating_weather(ex_wt)
}
, mc.cores = 4)

bg_wth <- do.call(rbind, weather_data) 

saveRDS(bg_wth, file = paste0(loc.output, "bg_tiger_spain_climatic_variables_mun_monthly.rds"))

# Adding Corine Land Covers ---------------------------------------------------
clc_surface <- readRDS(paste0(loc.output, "clc_surface_mun_level_0.rds"))

bg <- merge(bg_wth, clc_surface, by = c("municipality", "id"), all.x = TRUE)
bg$no_data <- NULL

saveRDS(bg, file = paste0(loc.output, "bg_tiger_spain.rds"))
