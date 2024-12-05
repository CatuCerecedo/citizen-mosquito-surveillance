# Extracting climatic variables from nc raster ---------------------------------
# Preparing ERA5 folder
my_nc_20 = paste0(loc.era5, "ERA5_EU_hourly_2020.nc")
my_nc_21 = paste0(loc.era5, "ERA5_EU_hourly_2021.nc")
my_nc_22 = paste0(loc.era5, "ERA5_EU_hourly_2022.nc")
my_nc_23 = paste0(loc.era5, "ERA5_EU_hourly_2023.nc")

extract_weather <- function(data_row){
  
  # Specify desired single point (within the bounds of your .nc file)
  x <- data_row$longitude # Selecting longitude coordinate
  y <- data_row$latitude # Selecting latitude coordinates
    
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
  values_t2m <- values_t2m[date >= data_row$start_date & date <= data_row$end_date]
  values_t2m <- values_t2m[, .(max_temperature = max(temperature, na.rm = TRUE),
                               min_temperature = min(temperature, na.rm = TRUE),
                               mean_temperature = (max(temperature, na.rm = TRUE) + min(temperature, na.rm = TRUE))/2), 
                           by = date]
  
  values_t2m <- values_t2m[, .(max_temperature = mean(max_temperature, na.rm = TRUE),
                               min_temperature = mean(min_temperature, na.rm = TRUE),
                               mean_temperature = mean(mean_temperature, na.rm = TRUE))] 
  
  data_row <- data_row %>%
    mutate(
      max_temperature = values_t2m$max_temperature,
      min_temperature = values_t2m$min_temperature,
      mean_temperature = values_t2m$mean_temperature
    )
  
  # Temperature 21 days lags
  print("Temperature 21")
  
  tmp_raster <- raster::brick(my_nc, varname=c("t2m"))
  
  values_t2m_21 <- raster::extract(tmp_raster, data_row) %>%
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
  
  values_t2m_21 <- values_t2m_21[, date := as.Date(date, format = "X%Y.%m.%d")]
  values_t2m_21 <- values_t2m_21[date >= (data_row$start_date - days(21))& date < data_row$start_date]
  values_t2m_21 <- values_t2m_21[, .(max_temperature = max(temperature, na.rm = TRUE),
                               min_temperature = min(temperature, na.rm = TRUE),
                               mean_temperature = (max(temperature, na.rm = TRUE) + min(temperature, na.rm = TRUE))/2), 
                           by = date]
  
  values_t2m_21 <- values_t2m_21[, .(max_temperature = mean(max_temperature, na.rm = TRUE),
                               min_temperature = mean(min_temperature, na.rm = TRUE),
                               mean_temperature = mean(mean_temperature, na.rm = TRUE))] 
  
  data_row <- data_row %>%
    mutate(
      l21max_temperature = values_t2m_21$max_temperature,
      l21min_temperature = values_t2m_21$min_temperature,
      l21mean_temperature = values_t2m_21$mean_temperature
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
  
  values_tp <- values_tp[date >= data_row$start_date & date <= data_row$end_date]
  
  values_tp <- values_tp[, .(precipitation = sum(precipitation, na.rm = TRUE))] # Daily grouping values
  
  data_row <- data_row %>%
    mutate(
      precipitation = values_tp$precipitation
    )
  
  # Precipitation
  print("Precipitation 21")
  tmp_raster <- raster::brick(my_nc, varname="tp")
  
  values_tp_21 <- raster::extract(tmp_raster, data_row) %>%
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
  
  values_tp_21 <- values_tp_21[, date := as.Date(date, format = "X%Y.%m.%d")]
  
  values_tp_21 <- values_tp_21[date >= (data_row$start_date - days(21))& date < data_row$start_date]
  
  values_tp_21 <- values_tp_21[, .(precipitation = sum(precipitation, na.rm = TRUE))] # Daily grouping values
  
  data_row <- data_row %>%
    mutate(
      l21precipitation = values_tp_21$precipitation
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
  
  values_dp <- values_dp[date >= data_row$start_date & date <= data_row$end_date]
  
  values_dp <- values_dp[, .(dew_point = mean(dew_point, na.rm = TRUE))]
  
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
  
  values_v <- values_v[date >= data_row$start_date & date <= data_row$end_date]
  
  values_v <- values_v[, .(v10 = mean(v10, na.rm = TRUE))]
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
  
  values_u <- values_u[date >= data_row$start_date & date <= data_row$end_date]
  
  values_u <- values_u[, .(u10 = mean(u10, na.rm = TRUE))]
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