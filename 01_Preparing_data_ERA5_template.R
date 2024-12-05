############### Preparing monthly mosquito data: ERA5 template #################
library(tidyverse)
library(sf)
library(mapSpain)
library(terra)
library(mosquitoR)
library(parallel)
library(data.table)

rm(list = ls())
sf::sf_use_s2(FALSE)
# Directories ------------------------------------------------------------------

# In local
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.era5 <- "/home/catuxa/Documents/Mosquito_Models/EU_Culex/ERA5_Download/"

# In cluster
loc.output <- paste0(getwd(), "/Spain_Tiger/OUTPUT/")
loc.data <- paste0(getwd(), "/Spain_Tiger/DATA/")
loc.era5 <- paste0(getwd(), "/EU_Culex/ERA5_Download/")

# Loading mosquito data --------------------------------------------------------
tiger <- readRDS(paste0(loc.output, "bg_traps_spain_2020_2022.rds"))
tiger <- st_as_sf(tiger, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
# names(tiger)

# Reading ERA5 -----------------------------------------------------------------
pname <- "ERA5_EU_monthly"
# NOT RUN:
# Only to know the variables names
# varnames(rast(paste0(loc.era5, pname, "2021.nc")) )

# Getting a template
tmp_temp <- rast(paste0(loc.era5, pname, ".grib")) 
template <- tmp_temp[[1]] # Raster to include the mosquito data

# Check the mosquito raster using a template: summing all counts during years
tiger_rast <- terra::rasterize(tiger, template, "females", fun = sum, na.rm=TRUE)
plot(tiger_rast) # It works

# We have to add a new random effect: the pixel_id!
tiger$pixel_id <- raster::cellFromXY(template, st_coordinates(tiger))
length(unique(tiger$pixel_id)) # 29

# Grouping by pixeles ----------------------------------------------------------
xy <- st_as_sf(as.points(template)) # We need the coordinates of pixels
xy$pixel_id <- as.factor(raster::cellFromXY(template, st_coordinates(xy))) # Including the ID_pixel
xy$`SFC (Ground or water surface); 10 metre u wind component [m/s]` <- NULL

# Selecting pixels only with mosquito data
xy <- xy %>% filter(pixel_id %in% unique(tiger$pixel_id))
xy$longitude <- st_coordinates(xy)[,1]
xy$latitude <- st_coordinates(xy)[,2]
plot(st_geometry(xy), add = TRUE) # It works

# In this step, we remove all information related to location, apart from country
tiger <- st_drop_geometry(tiger)
tiger$latitude <- NULL 
tiger$longitude <- NULL 

tiger <- merge(tiger, xy, by = "pixel_id") # 2188

tiger <- tiger %>% 
  st_drop_geometry() %>%
  mutate(
    m = month(end_date),
    y = year(end_date)
  ) %>%
  group_by(pixel_id, start_date, end_date, 
           longitude, latitude, m, y) %>%
  summarise(
    females = sum(females, na.rm = TRUE),
    trapping_effort = sum(trapping_effort, na.rm = TRUE),
    n_traps = n()
  )
sum(duplicated(tiger))

tiger <- st_as_sf(tiger, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) # N = 1077
# Only RUN the first time
# saveRDS(tiger, file = paste0(loc.output, "tiger_spain_pixel_monthly.rds"))

# Full function ----------------------------------------------------------------
# Loading the functions
# In local
source("read_and_extract_grib_data.R")
# In cluster
# source("/home/usuaris/ccerecedo/EU_Culex/read_and_extract_grib_data.R")
cat("------------------------- Number of rows:", nrow(tiger), "\n")

weather_data <- mclapply(1:nrow(tiger), function(i){
  
  cat(paste0("Number row:", i))
  
  data_row <- tiger[i, ]
  
  ex_wt <- extract_weather(data_row)
  
  wth <- calculating_weather(ex_wt)
}
, mc.cores = 10)

bg_wth <- do.call(rbind, weather_data) 

saveRDS(bg_wth, file = paste0(loc.output, "tiger_spain_climatic_variables_pixel_monthly.rds"))

# Adding Corine Land Covers ---------------------------------------------------
clc_surface <- readRDS(paste0(loc.output, "clc_all_pixels.rds"))

tiger <- merge(bg_wth, clc_surface, by = "pixel_id", all.x = TRUE)
tiger$no_data <- NULL

saveRDS(tiger, file = paste0(loc.output, "tiger_spain.rds"))
