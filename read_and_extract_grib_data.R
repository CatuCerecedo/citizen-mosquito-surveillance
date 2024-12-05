# Adding other information -----------------------------------------------------

# Extracting climatic variables from grib raster ---------------------------------
# Preparing ERA5 folder
# In local
my_grib = "/home/catuxa/Documents/Mosquito_Models/EU_Culex/ERA5_Download/ERA5_EU_monthly.grib"

# In cluster
# my_grib = paste0(getwd(), "/EU_Culex/ERA5_Download/ERA5_EU_monthly.grib")


extract_weather <- function(data_row){
  
  # Specify desired single point (within the bounds of your .grib file)
  mnth <- sprintf("%02d", data_row$m)
  yr <- as.character(data_row$y)
  
  ## Selecting the specific date of the sample ---------------------------------
  tmp_raster <- rast(my_grib)
  target_date <- paste0(yr, "-", mnth, "-01")
  cat("----- date", target_date, "\n")
  
  date_index <- which(format(time(tmp_raster), "%Y-%m-%d") == target_date)
  
  tmp_raster <- tmp_raster[[date_index]]
  
  # Temperature
  print("Temperature")
  
  temp_layers <- grep("2.metre.temperature", names(tmp_raster), value = TRUE)
  tmp_raster_date <- tmp_raster[[temp_layers]] - 273.15
  # plot(tmp_raster_date)
  
  fraction <- terra::extract(tmp_raster_date, data_row, exact = TRUE)
  data_row$mean_temperature <- mean(fraction[,2], na.rm = TRUE)
  
  # Precipitation
  print("Precipitation")
  temp_layers <- grep("precipitation", names(tmp_raster), value = TRUE)
  tmp_raster_date <- tmp_raster[[temp_layers]] 
  
  fraction <- terra::extract(tmp_raster_date, data_row, exact = TRUE)
  data_row$precipitation <- mean(fraction[,2], na.rm = TRUE)
  
  # Dewpoint 
  print("Dewpoint")
  temp_layers <- grep("dewpoint", names(tmp_raster), value = TRUE)
  tmp_raster_date <- tmp_raster[[temp_layers]] - 273.15
  
  fraction <- terra::extract(tmp_raster_date, data_row, exact = TRUE)
  data_row$dew_point <- mean(fraction[,2], na.rm = TRUE)
  
  # Wind speed
  print("v10")
  temp_layers <- grep("v wind", names(tmp_raster), value = TRUE)
  tmp_raster_date <- tmp_raster[[temp_layers]]
  
  fraction <- terra::extract(tmp_raster_date, data_row, exact = TRUE)
  data_row$v10 <- mean(fraction[,2], na.rm = TRUE)
  
  print("u10")
  temp_layers <- grep("u wind", names(tmp_raster), value = TRUE)
  tmp_raster_date <- tmp_raster[[temp_layers]]
  
  fraction <- terra::extract(tmp_raster_date, data_row, exact = TRUE)
  data_row$u10 <- mean(fraction[,2], na.rm = TRUE)
  
  data_row <- data_row %>% st_drop_geometry() 
}

calculating_weather <- function(point_out){
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
