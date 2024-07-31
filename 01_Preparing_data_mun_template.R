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
#   dplyr::select(name, id) %>%
#   rename(municipality = "name")
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
    month = as.integer(month(end_date)), 
    year = year(end_date)
  ) %>% 
  group_by(year, month, municipality, id) %>% 
  summarize(
    females = sum(females, na.rm = TRUE),
    trapping_effort = sum(trapping_effort, na.rm = TRUE),
    n_traps = n(),
    .groups ="drop") %>%
  ungroup() # 256 rows
sum(duplicated(tiger))

plot(table(tiger$females)) # Negative binomial

tiger <- merge(tiger, spain, by = c("municipality", "id")) 
st_geometry(tiger) <- "geometry"

# Extracting climatic variables from nc raster ---------------------------------
centroids <- st_centroid(tiger) %>% 
  mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>% 
  st_drop_geometry() %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
# plot(st_geometry(centoids))

# Preparing ERA5 folder
my_nc_20 = paste0(loc.era5, "ERA5_EU_hourly_2020.nc")
my_nc_21 = paste0(loc.era5, "ERA5_EU_hourly_2021.nc")
my_nc_22 = paste0(loc.era5, "ERA5_EU_hourly_2022.nc")

extract_weather <- function(data_row){
  
  # Specify desired single point (within the bounds of your .nc file)
  x <- data_row$lon # Selecting longitude coordinate
  y <- data_row$lat# Selecting latitude coordinate
  
  if (data_row$year == "2020"){
    my_nc = my_nc_20
  } 
  
  if (data_row$year == "2021"){
    my_nc = my_nc_21
  } 
  
  if (data_row$year == "2022"){
    my_nc = my_nc_22
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
  
  values_t2m <- values_t2m[lubridate::month(date) == data_row$month, ]
  
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
  
  values_tp <- values_tp[lubridate::month(date) == data_row$month, ]
  
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
  
  values_dp <- values_dp[lubridate::month(date) == data_row$month, ]
  
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
  
  values_v <- values_v[lubridate::month(date) == data_row$month, ]
  
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
  
  values_u <- values_u[lubridate::month(date) == data_row$month, ]
  
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
weather_data <- mclapply(1:nrow(centroids), function(i){
  
  cat(paste0("Number row:", i, "\n"))
  
  data_row <- centroids[i, ]
  
  ex_wt <- extract_weather(data_row)
  
  wth <- calculating_weather(ex_wt, data_row)
}
, mc.cores = 8)

bg_wth <- do.call(rbind, weather_data) 

saveRDS(bg_wth, file = paste0(loc.output, "bg_tiger_spain_climatic_variables_mun_monthly.rds"))

# Adding Corine Land Covers ---------------------------------------------------
clc_surface <- readRDS(paste0(loc.output, "clc_surface_mun.rds"))

bg <- merge(bg_wth, clc_surface, by = c("municipality", "id"), all.x = TRUE)
bg$no_data <- NULL

saveRDS(bg, file = paste0(loc.output, "bg_tiger_spain.rds"))
