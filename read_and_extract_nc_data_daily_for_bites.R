# Extracting climatic variables from nc raster ---------------------------------
# Preparing ERA5 folder
my_nc_20 = paste0(loc.era5, "ERA5_EU_hourly_2020.nc")
my_nc_21 = paste0(loc.era5, "ERA5_EU_hourly_2021.nc")
my_nc_22 = paste0(loc.era5, "ERA5_EU_hourly_2022.nc")
my_nc_23 = paste0(loc.era5, "ERA5_EU_hourly_2023.nc")

extract_weather <- function(data_row){
  
  point_out <- data_row
  # Specify desired single point (within the bounds of your .nc file)
  x <- data_row$longitude # Selecting longitude coordinate
  y <- data_row$latitude # Selecting latitude coordinates
  d <- data_row$end_date
    
  if (year(data_row$end_date) == 2020){
    my_nc = my_nc_20
  } 
  
  if (year(data_row$end_date) == 2021){
    my_nc = my_nc_21
  } 
  
  if (year(data_row$end_date) == 2022){
    my_nc = my_nc_22
  } 
  
  if (year(data_row$end_date) == 2023){
    my_nc = my_nc_23
  } 
  
  # Extraction data for exposure week
  
  # Temperature
  if (month(d) == 1){
    my_nc_2 <- paste0(loc.era5, "ERA5_EU_hourly_", year(d) -1, ".nc")
    
    tmp_raster <- rast(my_nc) 
    tmp_raster_2 <- rast(my_nc_2) 
    
    tmp_raster <- c(tmp_raster, tmp_raster_2)
  } else {
    tmp_raster <- rast(my_nc) 
  }
  
  # Dates from raster
  date_index <- which(format(time(tmp_raster), "%Y-%m-%d") >= (d - days(21)) &
                        format(time(tmp_raster), "%Y-%m-%d") <= d)
  tmp_raster <- tmp_raster[[date_index]]
  
  # Extracting hourly rainfull data
  print("Extracting temperature")
  tmp_raster_t2m <- tmp_raster["t2m"]
  
  fraction <- terra::extract(tmp_raster_t2m, data_row, exact = TRUE) %>% t()
  tm <- format(time(tmp_raster_t2m), "%Y-%m-%d")
  
  # Transform to dataframe
  values <- data.frame(
    fraction = fraction[-1, 1] -  273.15, # Removing the ID column!
    tm = tm
    ) %>%
    group_by(tm) %>%
    summarise(
      max_temperature = max(fraction, na.rm = TRUE),
      min_temperature = min(fraction, na.rm = TRUE),
      mean_temperature = (max(fraction, na.rm = TRUE) + min(fraction, na.rm = TRUE))/2
    )
  
  point_out <- point_out %>%
    mutate(
      max_temperature = values$max_temperature[values$tm == d], 
      min_temperature = values$min_temperature[values$tm == d],
      mean_temperature = values$mean_temperature[values$tm == d]
    ) %>%
    mutate(
      l21max_temperature = mean(values$max_temperature[values$tm >= d - days(21) & values$tm < d], na.rm = TRUE),
      l21min_temperature = mean(values$min_temperature[values$tm >= d - days(21) & values$tm < d], na.rm = TRUE),
      l21mean_temperature = mean(values$mean_temperature[values$tm >= d - days(21) & values$tm < d], na.rm = TRUE),
    ) 
    
  # Precipitation
  print("Extracting precipitation")
  tmp_raster_tp <- tmp_raster["tp"]
  
  fraction <- terra::extract(tmp_raster_tp, data_row, exact = TRUE) %>% t()
  tm <- format(time(tmp_raster_tp), "%Y-%m-%d")
  
  # Transform to dataframe
  values <- data.frame(
    fraction = fraction[-1, 1], # Removing the ID column!
    tm = tm) %>%
    group_by(tm) %>%
    summarise(
      precipitation = sum(fraction, na.rm = TRUE)
    )
  
  point_out <- point_out %>%
    mutate(
      precipitation = values$precipitation[values$tm == d],
      l21precipitation = sum(values$precipitation[values$tm >= d - days(21) & values$tm < d], na.rm = TRUE)
    )
  
  # Dewpoint 
  print("Extracting dewpoint")
  tmp_raster_d2m <- tmp_raster["d2m"]
  
  fraction <- terra::extract(tmp_raster_d2m, data_row, exact = TRUE) %>% t()
  tm <- format(time(tmp_raster_d2m), "%Y-%m-%d")
  
  # Transform to dataframe
  values <- data.frame(
    fraction = fraction[-1, 1] - 273.15, # Removing the ID column!
    tm = tm
    ) %>%
    group_by(tm) %>%
    summarise(
      dew_point = (max(fraction, na.rm = TRUE) + min(fraction, na.rm = TRUE))/2
    )
  
  point_out <- point_out %>%
    mutate(
      dew_point = values$dew_point[values$tm == d],
      l21dew_point = mean(values$dew_point[values$tm >= d - days(21) & values$tm < d], na.rm = TRUE)
    )
  
  # Wind speed
  print("Extracting v10")
  tmp_raster_v10 <- tmp_raster["v10"]
  
  fraction <- terra::extract(tmp_raster_v10, data_row, exact = TRUE) %>% t()
  tm <- format(time(tmp_raster_v10), "%Y-%m-%d")
  
  # Transform to dataframe
  values <- data.frame(
    fraction = fraction[-1, 1], # Removing the ID column!
    tm = tm
  ) %>%
    group_by(tm) %>%
    summarise(
      v10 = mean(fraction, na.rm = TRUE)
    )
  
  point_out <- point_out %>%
    mutate(
      v10 = values$v10[values$tm == d],
      l21v10 = mean(values$v10[values$tm >= d - days(21) & values$tm < d], na.rm = TRUE)
    )
  
  print("Extracting u10")
  tmp_raster_u10 <- tmp_raster["u10"]
  
  fraction <- terra::extract(tmp_raster_u10, data_row, exact = TRUE) %>% t()
  tm <- format(time(tmp_raster_u10), "%Y-%m-%d")
  
  # Transform to dataframe
  values <- data.frame(
    fraction = fraction[-1, 1], # Removing the ID column!
    tm = tm
  ) %>%
    group_by(tm) %>%
    summarise(
      u10 = mean(fraction, na.rm = TRUE)
    )
  
  point_out <- point_out %>%
    mutate(
      u10 = values$u10[values$tm == d],
      l21u10 = mean(values$u10[values$tm >= d - days(21) & values$tm < d], na.rm = TRUE)
    )
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
      mean_relative_humidity = 100*10^(7.59138*((dew_point/(dew_point + 240.726))  - (mean_temperature/(mean_temperature + 240.726)))),
      l21mean_relative_humidity = 100*10^(7.59138*((l21dew_point/(l21dew_point + 240.726))  - (l21mean_temperature/(l21mean_temperature + 240.726)))),
      wind_speed = sqrt(u10^2 + v10^2),
      l21wind_speed = sqrt(l21u10^2 + l21v10^2)
    ) 
  result <- st_drop_geometry(result)
  return(result)
}