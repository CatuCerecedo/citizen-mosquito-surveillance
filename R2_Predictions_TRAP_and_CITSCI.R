############################## TRAP AND CITSCI PREDICTIONS #####################
#' Daily predictions of TRAP and CITSCI models
#' 
#' The script contains extraction nc files download directly from ERA5 dataset
#' and calculation functions of covariates, and Bayesian prediction (using brms library)

################################################################################
library(tidyverse)
library(sf)
library(parallel)
library(dplyr)
library(lubridate)
library(data.table)
library(brms)
library(terra)

rm(list = ls())

sf::sf_use_s2(FALSE)
# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.era5 <- paste0(getwd(), "/ERA5_Download/") # Directory where ERA5 nc files are stores (one file per month)

# Obtaining centroids from municipalities --------------------------------------
spain <- readRDS(paste0(loc.data, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  ) %>%
  st_transform(4326) # This map was download from mapSpain R librery:  https://ropenspain.github.io/mapSpain/articles/x02_mapasesp.html

xy <- st_centroid(spain) %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry() %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Loading CLC data ------------------------------------------------------------- 
clc_surface <- readRDS(paste0(loc.output, "clc_surface_mun_level_0.rds"))

# Model loading ----------------------------------------------------------------

# Select one of the models to calculate the predictions

# mdl <- "mtiger16_trap.rds"
# mdl_name <- "/mtiger16_trap"
# fldr <- "COUNTS"
# sub <- "_trap" # with _

mdl <- "mtiger7_ma.rds"
mdl_name <- "/mtiger7"
fldr <- "MA"
sub <- "_ma"

# NC to raster and extraction --------------------------------------------------
ncores = 8

new_points <- xy

# Loading model
model <- readRDS(paste0(loc.output, mdl))

# Preparing ERA5 folder
recalculate_variables = FALSE

for (year in c("2020", "2021", "2022")){
  # set path and file name
  ncname <- paste0(loc.era5, "ERA5_ES_hourly_", year) 
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
    
    # Function to extract daily data -----------------------------------------------
    for (sel_date in seq.Date(first_date, last_day, 1)){
      sel_date <- as_date(sel_date)
      print(sel_date)
      
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
          
          # We search the index number
          index <- which(dates == sel_date)
          
          # We filter on rasterBrick
          tmp_raster_day_delay <- tmp_raster[[index]]
          
          # Extraction data for one day
          values_t2m <-  raster::calc(tmp_raster_day_delay, fun = function(x) {
            (max(x, na.rm = TRUE) + min(x, na.rm = TRUE)) / 2
          }) %>%
            raster::extract(., new_points) %>%
            as.data.frame() %>%
            mutate(
              mean_temperature = . - 273.15 # Our should be in degrees Celsius
            ) %>%
            dplyr::select(mean_temperature)
          
          # Temperature 21
          # We search the index number
          if (year(sel_date) == year(sel_date - days(21))){
            index <- which(dates == (sel_date - days(21)))
          } else { 
            index <- which(dates >= first_date & dates <= sel_date)
          }
          
          # We filter on rasterBrick
          tmp_raster_day_delay <- tmp_raster[[index]]
          
          # Extraction data for one day
          values_t2m_21 <-  raster::calc(tmp_raster_day_delay, fun = function(x) {
            (max(x, na.rm = TRUE) + min(x, na.rm = TRUE)) / 2
          }) %>%
            raster::extract(., new_points) %>%
            as.data.frame() %>%
            mutate(
              l21mean_temperature = . - 273.15 # Our should be in degrees Celsius
            ) %>%
            dplyr::select(l21mean_temperature)
          
          values_t2m_21_min <-  raster::calc(tmp_raster_day_delay, fun = function(x) {
            min(x, na.rm = TRUE)
          }) %>%
            raster::extract(., new_points) %>%
            as.data.frame() %>%
            mutate(
              l21min_temperature = . - 273.15 # Our should be in degrees Celsius
            ) %>%
            dplyr::select(l21min_temperature)
          
          # Precipitation ------------------------------------------------------------------
          # The tmp raster is a rasterBrick with a lot of layers
          tmp_raster <- raster::brick(ncfname, varname=c("tp"))
          
          # We search the index number
          index <- which(dates == sel_date)
          
          # We filter on rasterBrick
          tmp_raster_day_delay <- tmp_raster[[index]]
          
          # Extraction data for one day
          values_tp <-  raster::calc(tmp_raster_day_delay, fun = sum) %>%
            raster::extract(., new_points) %>%
            as.data.frame() %>%
            mutate(
              precipitation = . + 0
            ) %>%
            dplyr::select(precipitation)
          
          # Precipitation 21
          # We search the index number
          if (year(sel_date) == year(sel_date - days(21))){
            index <- which(dates == (sel_date - days(21)))
          } else { 
            index <- which(dates >= first_date & dates <= sel_date)
          }
          
          # We filter on rasterBrick
          tmp_raster_day_delay <- tmp_raster[[index]]
          
          # Extraction data for one day
          values_tp_21 <-  raster::calc(tmp_raster_day_delay, fun = sum) %>%
            raster::extract(., new_points) %>%
            as.data.frame() %>%
            mutate(
              l21precipitation = . + 0
            ) %>%
            dplyr::select(l21precipitation)
          
          # Dewpoint ---------------------------------------------------------------------
          tmp_raster <- raster::brick(ncfname, varname="d2m")
          
          # We filter on rasterBrick
          tmp_raster_day_delay <- tmp_raster[[index]]
          
          # Extraction data for one day
          values_dp <-  raster::calc(tmp_raster_day_delay, fun = mean) %>%
            raster::extract(., new_points) %>%
            as.data.frame() %>%
            mutate(
              dew_point = . - 273.15
            ) %>%
            dplyr::select(dew_point)
          
          # Gathering information ----------------------------------------------------
          prep_data_day <- new_points %>% mutate(
            mean_temperature = values_t2m$mean_temperature,
            l21mean_temperature = values_t2m_21$l21mean_temperature,
            l21min_temperature = values_t2m_21_min$l21min_temperature,
            precipitation = values_tp$precipitation,
            l21precipitation = values_tp_21$l21precipitation,
            dew_point = values_dp$dew_point,
          )
          
          prep_data_day <- prep_data_day %>%
            mutate(
              mean_relative_humidity = 100*10^(7.59138*( (dew_point/(dew_point + 240.726))  - (mean_temperature/(mean_temperature+240.726)) )),
              trapping_effort = 1,
              n_traps = 1,
              date = sel_date
            ) %>%
            dplyr::select(-dew_point) %>% 
            setDT()
          
          # Adding clc data
          prep_data_day <- merge(prep_data_day, 
                                 clc_surface, 
                                 by = c("municipality", "id"), all.x = TRUE)
          
          saveRDS(prep_data_day, file = paste0(loc.output, "daily_weather_data/prep_",
                                               sel_date, ".rds"))
          
          return(prep_data_day)
          
        }, 
        mc.cores = ncores))
      } else {
        data_prep_day <- readRDS(paste0(loc.output, "daily_weather_data/prep_",
                                        sel_date, ".rds"))
      }
      print("Calculating variables done")
      
      # Function to predict daily data -----------------------------------------------
      
      data_prep_day <- data_prep_day %>% 
        mutate(
          attractor = "CO2",
          n_traps = 1,
          SE = 1,
          n_total_reporters = 1,
          n_total_reports = 1,
          trapping_effort = 7
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
                    2, function(x) mean(x)) %>%
          as.data.frame()
        
        colnames(pp) <- as.Date(sel_date)
        
        data_chunk <- data_chunk %>% dplyr::select(municipality, id, prov_name, lon, lat)
        
        pp <- bind_cols(data_chunk, pp)
        
        return(pp)
        
      }, mc.cores = ncores))
      
      pred_points <- merge(pred_points, pred,  by = c("municipality", "id", "prov_name", "lon", "lat"))
      
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
      
      pred_points_sd <- merge(pred_points_sd, pred_sd, by = c("municipality", "id", "prov_name", "lon", "lat"))
      
      
      } # Days
    
    print(paste("Saving pred: ", m))
    saveRDS(pred_points %>% st_drop_geometry(), file = paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                                                       "/tiger_", m, "_", year,  sub, ".rds"))
    
    print(paste("Saving pred sd: ", m))
    saveRDS(pred_points_sd %>% st_drop_geometry(), file = paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                                                          "/tiger_", m, "_", year, sub, "_sd.rds"))
  } # months
} # years
