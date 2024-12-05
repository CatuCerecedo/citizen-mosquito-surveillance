######################### Prediction: era5, clc, model #########################

library(tidyverse)
library(sf)
library(parallel)
library(dplyr)
library(lubridate)
library(data.table)
library(ncdf4) 
library(brms)
library(terra)

rm(list = ls())

# Directories ------------------------------------------------------------------
# In local
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.era5 <- paste0("/home/catuxa/Documents/Mosquito_Models/EU_Culex/ERA5_Download/")

# In Cluster
loc.output <- paste0(getwd(), "/Spain_Tiger/OUTPUT/")
loc.data <- paste0(getwd(), "/Spain_Tiger/DATA/")
loc.era5 <- paste0(getwd(), "/EU_Culex/ERA5_Download/")

sf::sf_use_s2(FALSE)

mdl <- "mtiger3_occu.rds"
mdl_name <- "/mtiger3_occu"
fldr <- "Suitability"
sub <- "_occu" # with _

# Obtaining centroids from municipalities --------------------------------------
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  )

xy <- st_centroid(spain) %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
    ) %>%
  st_drop_geometry() %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# plot(st_geometry(spain))
# plot(st_geometry(xy["id"]), add = TRUE, fill = "red", size = 2)

# Loading CLC data ------------------------------------------------------------- 
clc_surface <- readRDS(paste0(loc.output, "clc_surface_mun_level_0.rds"))

# Function to calculate daily average values------------------------------------
extract_acc_values <- function(tmp_raster_day_delay, new_points, var = "t", func = c("min", "max", "mean"), d = 21){
  
  if (var == "t") {
    # Depends on the max/min average we make the following calculations:
    if ("min" %in% func){
      # Calculate the min value per day
      tmp_raster_day_delay_value <-  raster::stackApply(tmp_raster_day_delay, indices = rep(seq(1:d), each = 23), fun = min)
      
      # Averge the min values of all days
      l21min <-  raster::stackApply(tmp_raster_day_delay_value, indices = rep(seq(1:1), 
                                                                              each = d), 
                                    fun = mean) %>%
        raster::extract(., new_points) %>%
        as.data.frame() %>%
        mutate(
          l21min = . - 273.15 # Our should be in degrees Celsius
        ) %>%
        dplyr::select(l21min)
    } else {l21min <- NULL}
    
    if ("max" %in% func){
      tmp_raster_day_delay_value <-  raster::stackApply(tmp_raster_day_delay, indices = rep(seq(1:d), each = 23), fun = max)
      l21max <-  raster::stackApply(tmp_raster_day_delay_value, indices = rep(seq(1:1), 
                                                                              each = d), 
                                    fun = mean) %>%
        raster::extract(., new_points) %>%
        as.data.frame() %>%
        mutate(
          l21max = . - 273.15
        ) %>%
        dplyr::select(l21max)
    } else {l21max <- NULL}
    
    if ("mean" %in% func){
      l21mean <-  raster::stackApply(tmp_raster_day_delay, indices = rep(seq(1:1), 
                                                                         each = d), 
                                     fun = mean) %>%
        raster::extract(., new_points) %>%
        as.data.frame() %>%
        mutate(
          l21mean = . - 273.15
        ) %>%
        dplyr::select(l21mean)
    } else {l21mean <- NULL}
    
    values <- bind_cols(new_points, l21min, l21max, l21mean) %>%
      mutate(
        date = sel_date
      ) %>%
      st_drop_geometry() %>%
      setDT()
  }
  
  if(var == "p"){
    
    # Precipitation has no min or max values on our study
    l21mean <-  raster::stackApply(tmp_raster_day_delay, indices = rep(seq(1:1), 
                                                                       each = d), 
                                   fun = sum) %>%
      raster::extract(., new_points) %>%
      as.data.frame() %>%
      rename(
        l21mean = "."
      ) %>%
      dplyr::select(l21mean)
    
    values <- bind_cols(new_points, l21mean) %>%
      mutate(
        date = sel_date
      ) %>%
      st_drop_geometry() %>%
      setDT()
  }
  
  if(var == "dp"){
    
    # Precipitation has no min or max values on our study
    l21mean <-  raster::stackApply(tmp_raster_day_delay, indices = rep(seq(1:1), 
                                                                       each = d), 
                                   fun = mean) %>%
      raster::extract(., new_points) %>%
      as.data.frame() %>%
      mutate(
        l21mean = . - 273.15
      ) %>%
      dplyr::select(l21mean)
    
    values <- bind_cols(new_points, l21mean) %>%
      mutate(
        date = sel_date
      ) %>%
      st_drop_geometry() %>%
      setDT()
  }
  
  return(values)
  
}

# NC to raster and extraction --------------------------------------------------
ncores = 1

new_points <- xy

# Loading model
model <- readRDS(paste0(loc.output, mdl))

# Preparing ERA5 folder
recalculate_variables = FALSE

for (year in c("2020", "2021", "2022")){
  # set path and file name
  ncname <- paste0(loc.era5, "ERA5_EU_hourly_", year)
  ncfname <- paste(ncname, ".nc", sep="")
  
  months <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
  
  for (m in months){
    
    # Empty box
    pred_points <- xy %>% st_drop_geometry()
    pred_points_sd <- xy %>% st_drop_geometry()
    
    print(paste("month: ", m, "year: ", year))
    
    # Calculate the first date of the month
    first_date <- lubridate::as_date(paste(year, m, "01", sep = "-"))
    
    # Calculate the last date of the month
    last_day <- ceiling_date(lubridate::as_date(paste(year, m, "01", sep = "-")), "month") - 1
    
    # Then, we select the date which we want filter
    sel_date <- as.Date(last_day)
    
    # Function to extract daily data -----------------------------------------------
    if (recalculate_variables == TRUE){
      
      data_prep_day <- bind_rows(mclapply(last_day, function(last_day){
        
        # Temperature ------------------------------------------------------------------
        # The tmp raster is a rasterBrick with a lot of layers
        tmp_raster <- raster::brick(ncfname, varname=c("t2m"))
        
        # We will join raster info and new points info
        new_points <- new_points %>% st_transform(st_crs(tmp_raster))
        
        # To extract the info from the rasterBrick by day.
        # First, er have to transform names attribute to dates
        dates <- as.Date(names(tmp_raster), format = "X%Y.%m.%d.%H.%M.%S") 
        
        print(month(sel_date))
        
        # We search the index number
        # index <- which(dates == sel_date)
        index_delay <- which(dates %in% seq(sel_date - days(sel_date - first_date), sel_date, by = "day"))
        
        # We filter on rasterBrick
        tmp_raster_day_delay <- tmp_raster[[index_delay]]
        
        # Extraction data for one day
        values_t2m <- extract_acc_values(tmp_raster_day_delay, new_points, var = "t", func = c("min", "mean"), d = 2)
        values_t2m <- setnames(values_t2m, old = c("l21min", "l21mean"), new = c("min_temperature", "mean_temperature")) %>% unique()
        values_t2m <- values_t2m %>% 
          mutate(
            m = lubridate::month(date),
            y = lubridate::year(date)
          ) 
        
        # Precipitation ------------------------------------------------------------------
        # The tmp raster is a rasterBrick with a lot of layers
        tmp_raster <- raster::brick(ncfname, varname=c("tp"))
        
        # To extract the info from the rasterBrick by day.
        # First, er have to transform names attribute to dates
        dates <- as.Date(names(tmp_raster), format = "X%Y.%m.%d.%H.%M.%S") 
        
        print(month(sel_date))
        
        # We search the index number
        # index <- which(dates == sel_date)
        index_delay <- which(dates %in% seq(sel_date - days(sel_date - first_date), sel_date, by = "day"))
        
        # We filter on rasterBrick
        tmp_raster_day_delay <- tmp_raster[[index_delay]]
        
        # Extraction data for one day
        values_tp <- extract_acc_values(tmp_raster_day_delay, new_points, var = "p", d = 1)
        values_tp <- setnames(values_tp, old = c("l21mean"), new = c("precipitation")) %>% unique()
        values_tp <- values_tp %>% 
          mutate(
            m = lubridate::month(date),
            y = lubridate::year(date)
          ) 
        
        # Dewpoint ---------------------------------------------------------------------
        tmp_raster <- raster::brick(ncfname, varname="d2m")
        
        # To extract the info from the rasterBrick by day.
        # First, er have to transform names attribute to dates
        dates <- as.Date(names(tmp_raster), format = "X%Y.%m.%d.%H.%M.%S") 
        
        # We search the index number
        index_delay <- which(dates %in% seq(sel_date - days(sel_date - first_date), sel_date, by = "day"))
        
        # We filter on rasterBrick
        tmp_raster_day_delay <- tmp_raster[[index_delay]]
        
        # Extraction data for 1 accumulated days
        values_dp <- extract_acc_values(tmp_raster_day_delay, new_points, var = "t", func = c("mean"), d = 2)
        values_dp <- setnames(values_dp, old = c('l21mean'), new = c('dew_point')) %>% unique()
        values_dp <- values_dp %>% 
          mutate(
            m = lubridate::month(date),
            y = lubridate::year(date)
          ) 
        
        # Gathering information ----------------------------------------------------
        prep_data_day <- merge(values_t2m, values_dp, 
                               by = c("municipality", "id", "prov_name", "lon", "lat", "date", "m", "y"))
        prep_data_day <- merge(prep_data_day, values_tp, 
                               by = c("municipality", "id", "prov_name", "lon", "lat", "date", "m", "y"))
        
        
        prep_data_day <- prep_data_day %>%
          mutate(
            mean_relative_humidity = 100*10^(7.59138*( (dew_point/(dew_point + 240.726))  - (mean_temperature/(mean_temperature+240.726)) )),
            # l21mean_relative_humidity = 100*10^(7.59138*( (l21dew_point/(l21dew_point + 240.726))  - (l21mean_temperature/(l21mean_temperature+240.726)) )),
            # l21FH = case_when(l21mean_relative_humidity < 40~0, l21mean_relative_humidity >95~0, (l21mean_relative_humidity >=40 & l21mean_relative_humidity <= 95)~((l21mean_relative_humidity/55)-(40/55))),
            # l21FT = case_when(l21mean_temperature<=15~0, l21mean_temperature>30~0, (l21mean_temperature>15 & l21mean_temperature <=20)~(.2*l21mean_temperature)-3, (l21mean_temperature>20 & l21mean_temperature<=25)~1, (l21mean_temperature>25 & l21mean_temperature <= 30)~(-.2*l21mean_temperature)+6),
            # l21mwi = l21FH*l21FT,
            trapping_effort = 30,
            attractor = "CO2",
            n_traps = 1,
            month = m,
            year = year,
          ) %>%
          dplyr::select(-dew_point) %>% 
          setDT()
        
        # Adding clc data
        prep_data_day <- merge(prep_data_day, 
                               clc_surface, 
                               by = c("municipality", "id"), all.x = TRUE)
        
        # # Adding MA probability
        # ma_raster <- rast(paste0(loc.output, "MA_data/ma_raster_", year, "_", m, ".tif"))
        # 
        # prep_data_day_sf <- prep_data_day %>%
        #   st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = TRUE)
        # 
        # prep_data_day <- terra::extract(ma_raster, prep_data_day_sf) %>%
        #   cbind(prep_data_day, .) %>%
        #   dplyr::select(-ID) %>%
        #   rename(ma_prob = "est")
        
        saveRDS(prep_data_day, file = paste0(loc.output, "monthly_weather_data/prep_",
                                             m, "-", year, ".rds"))
        
        return(prep_data_day)
        
      }, mc.cores = ncores))
    } else {
      data_prep_day <- readRDS(paste0(loc.output, "monthly_weather_data/prep_",
                                      m, "-", year, ".rds"))
    }
    
    print("Calculating variables done")
    
    # # Only cat
    # 
    # data_prep_day <- prep_data_day 
    
    # Function to predict daily data -----------------------------------------------
    
    # It requires a lot of space
    # We split the memory
    data_prep_day <- data_prep_day %>% 
      mutate(
        trapping_effort = 30,
        attractor = "CO2",
        n_traps = 1,
        SE = 1,
        n_total_reporters = 1,
        n_total_reports = 1
        # min_temperature = scale(min_temperature),
        # agricultural = scale(agricultural),
        # forests_scrub = scale(forests_scrub)
      )
    
    nrow_these_pred_points = nrow(data_prep_day)
    max_chunksize = 300000
    chunksize = min(as.integer((nrow_these_pred_points/ncores)), max_chunksize)
    
    pred <- bind_rows(mclapply(seq(1, nrow(data_prep_day), chunksize), function(i){
      print(i)
      
      data_chunk = data_prep_day[i:min(nrow(data_prep_day), (i+(chunksize-1))), ]
      flush.console()
      
      pp <- apply(posterior_predict(model,
                                    newdata = data_chunk,
                                    allow_new_levels = TRUE,
                                    re_formula =  NA,
                                    ndraws = 1000),
                  2, function(x) mean(x)) %>% # Or mean
        as.data.frame()
      
      colnames(pp) <- as.Date(sel_date)
      
      data_chunk <- data_chunk %>% dplyr::select(municipality, id, prov_name, lon, lat)
      
      pp <- bind_cols(data_chunk, pp)
      
      return(pp)
      
    }, mc.cores = ncores))
    
    pred_points <- merge(pred_points, pred, by = c("municipality", "id", "prov_name", "lon", "lat"))
    
    print(paste("Saving pred: ", m))
    saveRDS(pred_points %>% st_drop_geometry(), file = paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                                                              "/tiger_", m, "_", year, sub, ".rds"))
    pred_sd <- bind_rows(mclapply(seq(1, nrow(data_prep_day), chunksize), function(i){
      print(i)
      
      data_chunk = data_prep_day[i:min(nrow(data_prep_day), (i+(chunksize-1))), ]
      flush.console()
      
      pp <- apply(posterior_predict(model, 
                                    newdata = data_chunk,
                                    allow_new_levels = TRUE, 
                                    re_formula = NA,
                                    ndraws = 1000), 
                  2, function(x) sd(x)) %>% 
        as.data.frame()
      
      colnames(pp) <- as.Date(sel_date)
      
      data_chunk <- data_chunk %>% dplyr::select(municipality, id, prov_name, lon, lat)
      
      pp <- bind_cols(data_chunk, pp)
      
      return(pp)
      
    }, mc.cores = ncores))
    
    pred_points_sd <- merge(pred_points_sd, pred_sd, by =  c("municipality", "id", "prov_name", "lon", "lat"))
    print(paste("Saving pred sd: ", m))
    saveRDS(pred_points_sd %>% st_drop_geometry(), file = paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                                                                 "/tiger_", m, "_", year, sub, "_sd.rds"))
  }
}
  
  
  
  


