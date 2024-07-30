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
pname <- "ERA5_EU_hourly_"
# NOT RUN:
# Only to know the variables names
# varnames(rast(paste0(loc.era5, pname, "2021.nc")) )

# Getting a template
tmp_temp <- rast(paste0(loc.era5, pname, "2021.nc"), "t2m") 
template <- tmp_temp[[1]] # Raster to include the mosquito data

# Check the mosquito raster using a template: summing all counts during years
tiger_rast <- terra::rasterize(tiger, template, "females", fun = sum, na.rm=TRUE)
plot(tiger_rast) # It works

# We have to add a new random effect: the pixel_id!
tiger$pixel_id <- raster::cellFromXY(template, st_coordinates(tiger))
length(unique(tiger$pixel_id)) # 17

# Grouping by pixeles ----------------------------------------------------------
xy <- st_as_sf(as.points(template)) # We need the coordinates of pixels
xy$pixel_id <- c(1:nrow(xy)) # Including the ID_pixel
xy$t2m_1 <- NULL

# Selecting pixels only with mosquito data
xy <- xy %>% filter(pixel_id %in% unique(tiger$pixel_id))
xy$longitude <- st_coordinates(xy)[,1]
xy$latitude <- st_coordinates(xy)[,2]
plot(st_geometry(xy), add = TRUE) # It works

# In this step, we remove all information related to location, apart from country
tiger <- st_drop_geometry(tiger)
tiger$latitude <- NULL 
tiger$longitude <- NULL 

tiger <- merge(tiger, xy, by = "pixel_id")

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

tiger <- st_as_sf(tiger, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
# Only RUN the first time
# saveRDS(tiger, file = paste0(loc.output, "tiger_spain_pixel_monthly.rds"))

# Extracting climatic variables from nc raster ---------------------------------
# Preparing ERA5 folder
my_nc_20 = paste0(loc.era5, "ERA5_EU_hourly_2020.nc")
my_nc_21 = paste0(loc.era5, "ERA5_EU_hourly_2021.nc")
my_nc_22 = paste0(loc.era5, "ERA5_EU_hourly_2022.nc")
my_nc_23 = paste0(loc.era5, "ERA5_EU_hourly_2023.nc")

extract_weather <- function(data_row){
  
  # Specify desired single point (within the bounds of your .nc file)
  x <- data_row$longitude # Selecting longitude coordinate
  y <- data_row$latitude # Selecting latitude coordinate
  
  if (data_row$y == 2020){
    my_nc = my_nc_20
  } 
  
  if (data_row$y == 2021){
    my_nc = my_nc_21
  } 
  
  if (data_row$y == 2022){
    my_nc = my_nc_22
  } 
  
  if (data_row$y == 2023){
    my_nc = my_nc_23
  } 
  
  # Extraction data for exposure week
  
  # Temperature
  print("Temperature")
  
  tmp_raster <- raster::brick(my_nc, varname=c("t2m"))
  
  values_t2m <- raster::extract(tmp_raster, data_row) %>%
    as.data.frame() %>%
    pivot_longer(
      starts_with("X"),
      names_to = "date",
      values_to = "temperature"
    ) %>%
    mutate(
      temperature = temperature -  273.15
    ) %>%
    st_drop_geometry() %>%
    setDT()
  
  values_t2m <- values_t2m[, date := as.Date(date, format = "X%Y.%m.%d")]
  values_t2m <- values_t2m[, .(max_temperature = max(temperature, na.rm = TRUE),
                               min_temperature = min(temperature, na.rm = TRUE),
                               mean_temperature = (max(temperature, na.rm = TRUE) + min(temperature, na.rm = TRUE))/2), 
                           by = date]
  
  values_t2m <- values_t2m[lubridate::month(date) == data_row$m, ]
  
  values_t2m <- values_t2m[, .(max_temperature = mean(max_temperature, na.rm = TRUE),
                               min_temperature = min(min_temperature, na.rm = TRUE),
                               mean_temperature = mean(mean_temperature, na.rm = TRUE)), 
                           by = lubridate::month(date)] 
  data_row <- data_row %>%
    mutate(
      max_temperature = values_t2m$max_temperature,
      min_temperature = values_t2m$min_temperature,
      mean_temperature = values_t2m$mean_temperature
    )
  
  # Precipitation
  print("Precipitation")
  tmp_raster <- raster::brick(my_nc, varname="tp")
  
  values_tp <- raster::extract(tmp_raster, data_row) %>%
    as.data.frame() %>%
    pivot_longer(
      starts_with("X"),
      names_to = "date",
      values_to = "precipitation"
    ) %>%
    mutate(
      precipitation = precipitation
    ) %>%
    st_drop_geometry() %>%
    setDT()
  
  values_tp <- values_tp[, date := as.Date(date, format = "X%Y.%m.%d")]
  
  values_tp <- values_tp[lubridate::month(date) == data_row$m, ]
  
  values_tp <- values_tp[, .(precipitation = sum(precipitation, na.rm = TRUE)),
                         by = lubridate::month(date)] # Daily grouping values
  
  data_row <- data_row %>%
    mutate(
      precipitation = values_tp$precipitation
    )
  
  # Dewpoint 
  print("Dewpoint")
  tmp_raster <- raster::brick(my_nc, varname="d2m")
  
  values_dp <- raster::extract(tmp_raster, data_row) %>%
    as.data.frame() %>%
    pivot_longer(
      starts_with("X"),
      names_to = "date",
      values_to = "dew_point"
    ) %>%
    mutate(
      dew_point = dew_point - 273.15
    ) %>%
    st_drop_geometry() %>%
    setDT()
  
  values_dp <- values_dp[, date := as.Date(date, format = "X%Y.%m.%d")]
  
  values_dp <- values_dp[lubridate::month(date) == data_row$m, ]
  
  values_dp <- values_dp[, .(dew_point = mean(dew_point, na.rm = TRUE)),
                         by = lubridate::month(date)]
  
  data_row <- data_row %>%
    mutate(
      dew_point = values_dp$dew_point
    )
  
  # Wind speed
  print("v10")
  tmp_raster <- raster::brick(my_nc, varname="v10")
  
  values_v <- raster::extract(tmp_raster, data_row) %>%
    as.data.frame() %>%
    pivot_longer(
      starts_with("X"),
      names_to = "date",
      values_to = "v10"
    ) %>%
    st_drop_geometry() %>%
    setDT()
  
  values_v <- values_v[, date := as.Date(date, format = "X%Y.%m.%d")]
  
  values_v <- values_v[lubridate::month(date) == data_row$m, ]
  
  values_v <- values_v[, .(v10 = mean(v10, na.rm = TRUE)),
                       by = lubridate::month(date)]
  data_row <- data_row %>%
    mutate(
      v10 = values_v$v10
    )
  
  print("u10")
  tmp_raster <- raster::brick(my_nc, varname="u10")
  
  values_u <- raster::extract(tmp_raster, data_row) %>%
    as.data.frame() %>%
    pivot_longer(
      starts_with("X"),
      names_to = "date",
      values_to = "u10"
    ) %>%
    st_drop_geometry() %>%
    setDT()
  
  values_u <- values_u[, date := as.Date(date, format = "X%Y.%m.%d")]
  
  values_u <- values_u[lubridate::month(date) == data_row$m, ]
  
  values_u <- values_u[, .(u10 = mean(u10, na.rm = TRUE)),
                       by = lubridate::month(date)]
  data_row <- data_row %>%
    mutate(
      u10 = values_u$u10
    )
  point_out <- data_row
  
  return(point_out)
}

calculating_weather <- function(point_out, data_row){
  # Dealing with humidity
  # Using Magnus approximation for relative humidity. 
  # From https://web.archive.org/web/20200212215746im_/
  # https://www.vaisala.com/en/system/files?file=documents/Humidity_Conversion_Formulas_B210973EN.pdf
  # valid for temp between -20 and 50 (which is our range here)
  print("calculating variables...")
  
  result <- point_out %>%
    mutate(
      mean_relative_humidity = 100*10^(7.59138*((dew_point/(dew_point + 240.726))  - (mean_temperature/(mean_temperature + 240.726)) )),
      wind_speed = sqrt(u10^2 + v10^2)
    ) 
  result <- st_drop_geometry(result)
  return(result)
}

# Full function ----------------------------------------------------------------
weather_data <- mclapply(1:nrow(tiger), function(i){
  
  cat(paste0("Number row:", i))
  
  data_row <- tiger[i, ]
  
  ex_wt <- extract_weather(data_row)
  
  wth <- calculating_weather(ex_wt, data_row)
}
, mc.cores = 10)

bg_wth <- do.call(rbind, weather_data) 

saveRDS(bg_wth, file = paste0(loc.output, "tiger_spain_climatic_variables_pixel_monthly.rds"))

# Adding Corine Land Covers ---------------------------------------------------
clc_surface <- readRDS(paste0(loc.output, "clc_all_pixels.rds"))

tiger <- merge(bg_wth, clc_surface, by = "pixel_id", all.x = TRUE)
tiger$no_data <- NULL

saveRDS(tiger, file = paste0(loc.output, "tiger_spain.rds"))
