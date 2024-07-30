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

# Loading Spain map ------------------------------------------------------------
spain <- mapSpain::esp_get_country(moveCAN = TRUE) %>% 
  st_transform(4326)

# Download Sampling effort data -----------------------------------------------
# Functions from mosquitoR package:
decimal_places <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
round_down = function(x, n) round( floor( (x*1000)/ (n*1000))*n, decimal_places(n))

trs_daily <- data.table::fread("https://github.com/Mosquito-Alert/sampling_effort_data/raw/main/sampling_effort_daily_cellres_025.csv.gz") %>% 
  as_tibble() %>% 
  rename(SEev = SE_expected) %>% 
  mutate(date = as_date(date)) %>%
  mutate(
    date = as_date(date), 
    masked_lon_6min = round_down(masked_lon, .1), 
    masked_lat_6min = round_down(masked_lat, .1)
  ) %>%
  st_as_sf(coords = c("masked_lon", "masked_lat"), remove = FALSE, crs = 4326)

index_intersects <- st_intersects(trs_daily, spain) 
index_intersects <- lengths(index_intersects) > 0

trs_daily <- trs_daily[index_intersects, ] # Only selecting Europe data

# # Checking the selected points
# plot(st_geometry(spain))
# plot(st_geometry(trs_daily), add = TRUE)

trs_monthly = trs_daily %>% 
  st_drop_geometry() %>%
  mutate(
    year = year(date), 
    month = month(date),
    n_total_reports = n_reports_albopictus + n_reports_bite + n_reports_culex + n_reports_japonicus + n_reports_aegypti + n_reports_koreicus,
    n_total_reporters = n_reporters_albopictus + n_reporters_bite + n_reporters_culex + n_reporters_japonicus + n_reporters_aegypti + n_reporters_koreicus
  ) %>% 
  group_by(year, month, TigacellID, masked_lon_6min, masked_lat_6min) %>%
  reframe(
    SEev = sum(SEev), 
    SE = 1-prod(1-SE), # This comes directly fom John's script
    n_total_reports = sum(n_total_reports),
    n_total_reporters = sum(n_total_reporters)
  ) %>%
  filter(year %in% c(2020, 2021, 2022, 2023))

# Reading ERA5 -----------------------------------------------------------------
pname <- "ERA5_EU_hourly_"

# Getting a template
tmp_temp <- rast(paste0(loc.era5, pname, "2021.nc"), "t2m") 
template <- tmp_temp[[1]] # Raster to include the mosquito data

# We have to add a new random effect: the pixel_id!
trs_monthly$pixel_id <- as.factor(raster::cellFromXY(template, st_coordinates(trs_monthly %>% 
                                                                                st_as_sf(coords = c("masked_lon_6min", "masked_lat_6min"), 
                                                                                         remove = FALSE, crs = 4326))))
length(unique(trs_monthly$pixel_id)) # 900
# Now we have incorporated the new ID cells

# Grouping sampling effort by pixel ID
trs_monthly = trs_monthly %>% 
  st_drop_geometry() %>%
  group_by(year, month, pixel_id) %>%
  summarise(
    SEev = sum(SEev), 
    SE = 1-prod(1-SE), # This comes directly fom John's script
    n_total_reports = sum(n_total_reports),
    n_total_reporters = sum(n_total_reporters),
    .groups = "drop") 

# Calculating new coordinates:
xy <- st_as_sf(as.points(template)) # We need the coordinates of pixels
xy$pixel_id <- c(1:nrow(xy)) # Including the ID_pixel
xy$t2m_1 <- NULL
xy <- xy %>% filter(pixel_id %in% unique(trs_monthly$pixel_id))
xy$lon <- st_coordinates(xy)[,1]
xy$lat <- st_coordinates(xy)[,2]

trs_monthly <- merge(trs_monthly, xy %>% st_drop_geometry(), by = "pixel_id", all.x = TRUE)
trs_monthly <- trs_monthly %>% drop_na(lon, lat) # Out of the European limits

# # Checking SE selection
# ggplot() +
#   geom_sf(data = spain, color = "black", fill = NA, size = 0.1, alpha = 0.5, na.rm = TRUE) +
#   geom_sf(data = trs_monthly %>%
#             filter(month == 7 & year == 2022) %>%
#             st_as_sf(coords = c("lon", "lat"), crs = 4326),
#           aes(color = n_total_reports), size = 0.5) +
#   scale_color_distiller("", palette = "Spectral") +
#   theme_classic()

# Loading reports data --------------------------------------------------------
vrs_full = readRDS(paste0(loc.data, 'mosquito_alert_cleaned_reports_cs_responses.Rds')) %>% 
  filter(!is.na(lon) & !is.na(lat)) %>% 
  st_as_sf(coords = c("lon", "lat"), crs=4326, remove=FALSE) %>% 
  mutate(package_name_version = paste0(package_name, "_", package_version)) %>%
  filter(movelab_certainty_category_euro_class_label == "aedes-albopictus", movelab_certainty_category_euro_class_value >= 1) 

index_intersects <- st_intersects(vrs_full, spain) 
index_intersects <- lengths(index_intersects) > 0

vrs_full <- vrs_full[index_intersects, ] # Only selecting Europe data

# # Checking the selected points
# plot(st_geometry(spain))
# plot(st_geometry(vrs_full), add = TRUE)

# We have to add a new random effect: the pixel_id!
vrs_full$pixel_id <- as.factor(raster::cellFromXY(template, st_coordinates(vrs_full)))
length(unique(vrs_full$pixel_id)) # 254
# Now we have incorporated the new ID cells

# Grouping reports by pixel ID
vrs_full <- vrs_full %>%
  st_drop_geometry() %>%
  mutate(
    month = as.integer(month(date)), 
    year = year(date)
  ) %>% 
  group_by(year, month, pixel_id) %>% 
  summarize(n_target_reps = n(), .groups ="drop") 

# Joining reports and Sampling Effort
ma_df = trs_monthly %>% 
  left_join(vrs_full, by = c("pixel_id", "year", "month"), multiple = "all") %>%
  replace_na(list(n_target_reps=0)) %>%
  filter(SE>0) %>% 
  mutate(
    reps_per_SEev = n_target_reps/(SEev),
    rep_per_reporter = n_total_reports/n_total_reporters
  ) %>% 
  mutate(SE_log_odds = log( (SE+.0001)/(1-SE+.0001)))

# adding app phase 
app_phases = tibble(app_phase = paste0("p", 1:7), 
                    date = c(as_date("2014-06-01"), 
                             as_date("2014-10-01"), 
                             as_date("2015-07-01"), 
                             as_date("2016-10-01"), 
                             as_date("2017-04-01"), 
                             as_date("2019-02-01"), 
                             as_date("2020-10-01"))) %>% 
  complete(date = seq.Date(as_date("2014-06-01"), today(), by="month")) %>% 
  fill(app_phase) %>% 
  mutate(month=month(date), year=year(date))

ma_df = ma_df %>% 
  left_join(app_phases %>% dplyr::select(app_phase, month, year))

ma_df = ma_df %>%
  mutate(
    any_reps = n_target_reps > 0
  )

# Adding other information -----------------------------------------------------

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
ma_df <- ma_df %>% st_as_sf(coords = c("lon", "lat"), remove = FALSE, crs = 4326)

weather_data <- mclapply(1:nrow(ma_df), function(i){
  
  cat(paste0("Number row:", i, "\n"))
  
  data_row <- ma_df[i, ]
  
  ex_wt <- extract_weather(data_row)
  
  wth <- calculating_weather(ex_wt, data_row)
}
, mc.cores = 10)

ma_wth <- do.call(rbind, weather_data) 

saveRDS(ma_wth, file = paste0(loc.output, "ma_tiger_spain_climatic_variables_pixel_monthly.rds"))

# Adding Corine Land Covers ---------------------------------------------------
clc_surface <- readRDS(paste0(loc.output, "clc_all_pixels.rds"))

ma <- merge(ma_df, clc_surface, by = "pixel_id", all.x = TRUE)
ma$no_data <- NULL

saveRDS(ma, file = paste0(loc.output, "ma_tiger_spain.rds"))

# # Martas' landcover from raster
# clc_surface <- readRDS(paste0(loc.data, "df_perc_landcover_EU.Rds"))
# clc_surface$pixel_id <- as.factor(raster::cellFromXY(template, st_coordinates(clc_surface %>% 
#                                                                                 st_as_sf(coords = c("lon", "lat"), 
#                                                                                          remove = FALSE, crs = 4326))))
# ma <- merge(ma_wth, clc_surface, by = "pixel_id", all.x = TRUE)
# ma$no_data <- NULL
# 
# saveRDS(ma, file = paste0(loc.output, "ma_df_europe.rds"))

