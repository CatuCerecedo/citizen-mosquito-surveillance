###################### Preparing Mosquito Alert Data ###########################
library(tidyverse)
library(sf)
library(terra)
library(parallel)
library(data.table)
library(readxl)

rm(list = ls())
# Directories ------------------------------------------------------------------
# In cluster
loc.output <- paste0(getwd(), "/Spain_Tiger/OUTPUT/")
loc.data <- paste0(getwd(), "/Spain_Tiger/DATA/")
loc.era5 <- paste0(getwd(), "/EU_Culex/ERA5_Download/")

# # In local
# loc.output <- paste0(getwd(), "/OUTPUT/")
# loc.data <- paste0(getwd(), "/DATA/")
# loc.era5 <- "/home/catuxa/Documents/Mosquito_Models/EU_Culex/ERA5_Download/"

sf::sf_use_s2(FALSE)
# Loading Spain map ------------------------------------------------------------
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  ) %>% st_transform(4326)

# Checking all municipalities
ggplot() +
  geom_sf(data = spain, fill = "black", color = "transparent",
          size = 0.01, alpha = 0.8, na.rm = TRUE) 

# Download Bites data ----------------------------------------------------------
bites <- read_excel(paste0(loc.data, "bites.xlsx"), col_names = TRUE) %>%
  janitor::clean_names() %>%
  mutate(
    date = as.Date(date), # Extract date
    time = format(as.POSIXct(date, tz = "UTC"), "%H:%M:%S") # Extrat hour
  ) %>%
  dplyr::select(date, time, longitude, latitude) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) 
  
index_intersects <- st_intersects(bites, spain) 
index_intersects <- lengths(index_intersects) > 0

bites <- bites[index_intersects, ] # Only selecting Spain data

# We have to add a new random effect: the municipality!
bites <- st_join(bites, spain, join = st_intersects) 
bites <- merge(bites, spain %>% st_drop_geometry() , by = c("municipality", "id", "prov_name")) # N = 5870

bites$n_target_reps <- 1
# # Checking the selected points
# plot(st_geometry(spain))
# st_geometry(bites) <- "geometry"
# plot(st_geometry(bites), add = TRUE)
bites$geometry <- NULL

# loading sampling_effort data -------------------------------------------------
trs_daily <- data.table::fread("https://github.com/Mosquito-Alert/sampling_effort_data/raw/main/sampling_effort_daily_cellres_025.csv.gz") %>%
  as_tibble() %>%
  mutate(
    date = as_date(date),
    y = year(as_date(date)),
    m = month(as_date(date)),
    n_total_reports = n_reports_albopictus + n_reports_bite + n_reports_culex + n_reports_japonicus + n_reports_aegypti + n_reports_koreicus,
    n_total_reporters = n_reporters_albopictus + n_reporters_bite + n_reporters_culex + n_reporters_japonicus + n_reporters_aegypti + n_reporters_koreicus
  ) %>%
  filter(y %in% c("2020", "2021", "2022")) %>%
  st_as_sf(coords = c("masked_lon", "masked_lat"), remove = FALSE, crs = 4326)

index_intersects <- st_intersects(trs_daily, spain) 
index_intersects <- lengths(index_intersects) > 0

trs_daily <- trs_daily[index_intersects, ] # Only selecting Spain data

# We have to add a new random effect: the municipality!
trs_daily <- st_join(trs_daily, spain, join = st_intersects) %>% st_drop_geometry() 
trs_daily <- merge(trs_daily, spain, by = c("municipality", "id", "prov_name")) # N = 55030

trs_daily = trs_daily %>%
  rename(SEev = SE_expected) %>%
  st_drop_geometry() %>%
  group_by(date, municipality, id, prov_name, geometry) %>%
  summarize(
    SEev = sum(SEev),
    SE = 1-prod(1-SE), # This comes directly fom John's script: one user at least send one report
    n_total_reports = sum(n_total_reports),
    n_total_reporters = sum(n_total_reporters)
  ) %>%
  ungroup() # n = 42306
# st_geometry(trs_daily) <- "geometry"

bites_df = trs_daily %>% 
  left_join(bites, by = c("municipality", "id", "date", "prov_name")) %>%
  replace_na(list(n_target_reps=0)) %>%
  filter(SE>0) %>% 
  mutate(
    any_reps = n_target_reps > 0,
    reps_per_SEev = n_target_reps/(SEev),
    rep_per_reporter = n_total_reports/n_total_reporters
  ) %>% 
  mutate(SE_log_odds = log( (SE+.0001)/(1-SE+.0001))) # 43443

bites_df <-  bites_df %>% mutate(
  app_phase = case_when(
    (date >= as_date("2014-06-01") & date < as_date("2014-10-01")) ~ "p1",
    (date >= as_date("2014-10-01") & date < as_date("2015-07-01")) ~ "p2",
    (date >= as_date("2015-07-01") & date < as_date("2016-10-01")) ~ "p3",
    (date >= as_date("2016-10-01") & date < as_date("2017-04-01")) ~ "p4",
    (date >= as_date("2017-04-01") & date < as_date("2019-02-01")) ~ "p5",
    (date >= as_date("2019-02-01") & date < as_date("2020-10-01")) ~ "p6",
    (date >= as_date("2020-10-01")) ~ "p7")
  )

summary(as.factor(bites_df$n_target_reps))

# bites_df$geometry <- NULL
st_geometry(bites_df) <- "geometry"
# Adding other information -----------------------------------------------------
# Extracting climatic variables from nc raster ---------------------------------
bites_df <- st_centroid(bites_df) %>%
  mutate(longitude = st_coordinates(.)[,1], latitude = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
# plot(st_geometry(centoids))

# Full function ----------------------------------------------------------------
# Loading the functions
# In local
# source("read_and_extract_grib_data.R")
# In cluster
# source("/home/usuaris/ccerecedo/EU_Culex/read_and_extract_grib_data.R")

# In local
# source("read_and_extract_nc_data_daily_for_bites.R")
source("Spain_Tiger/read_and_extract_nc_data_daily_for_bites.R")

bites_df <- bites_df %>%
  mutate(
    start_date = date,
    end_date = date
  )

# In cluster
# source("/home/usuaris/ccerecedo/EU_Culex/read_and_extract_grib_data.R")
cat("------------------------- Number of rows:", nrow(bites_df), "\n")
bites_df_1 <- bites_df[1:1371, ]
bites_df_2 <- bites_df[1372:nrow(bites_df), ]
weather_data <- mclapply(1:nrow(bites_df_2), function(i){
  
  cat(paste0("Number row:", i, "\n"))
  
  data_row <- bites_df_2[i, ]
  
  ex_wt <- extract_weather(data_row)
  
  wth <- calculating_weather(ex_wt)
}
, mc.cores = 10)

bites_wth <- do.call(rbind, weather_data) 

saveRDS(ma_wth, file = paste0(loc.output, "bites_spain_climatic_variables_pixel_daily.rds"))

# Adding Corine Land Covers ---------------------------------------------------
clc_surface <- readRDS(paste0(loc.output, "clc_surface_mun_level_0.rds"))

ma <- merge(ma_wth, clc_surface, by = c("municipality", "id"), all.x = TRUE) # N = 12566 // 42025
ma$no_data <- NULL

saveRDS(ma, file = paste0(loc.output, "bites_spain_daily.rds"))

# Adding population density ----------------------------------------------------
ma <- ma %>% mutate(y = as.factor(year(end_date)))

library(mapSpain) # Load library mapSpain to have shapefile with NATCODE municipalities Spain
# Read File with the number of people per year
pop_full <- data.frame()
for(n in c("20", "21", "22", "23")){
  pop <- read.csv(paste0(loc.data, "pobmun", n, ".csv"))
  
  pop$cmun <- ifelse(pop$CMUN<10, paste0("00",pop$CMUN),
                     ifelse(pop$CMUN<100, paste0("0",pop$CMUN),as.character(pop$CMUN)))
  pop$pop <- as.numeric(pop$POB) # This variable name, number of inhabitants, might change depending on the file. Check the name before should be POB or POB22 !
  pop$cpro <- ifelse(pop$CPRO < 10, paste0("0",pop$CPRO), pop$CPRO)
  
  pop <- pop %>% dplyr::select(cpro, cmun, pop)
  
  esp_can <- esp_get_munic_siane() %>% janitor::clean_names() %>%
    mutate(area = as.double(st_area(st_transform(., 3035)) / 1000000))
  esp_can$id <- as.factor(paste0(esp_can$codauto,
                                 esp_can$cpro,
                                 esp_can$cmun,
                                 esp_can$lau_code)) 
  esp_can <- esp_can %>% dplyr::select(id, cpro, cmun, area)
  esp_can <- merge(esp_can, pop,  by = c("cpro","cmun"), all.x = TRUE) # Join population density with esp_can to obtain NATCODE
  esp_can <- esp_can %>% mutate(dens = pop/area)
  
  esp_can <- esp_can %>% dplyr::select(id, dens, pop) %>% st_drop_geometry() %>%
    mutate(y = as.factor(paste0("20", n)))
  
  pop_full <- rbind(pop_full, esp_can)
}

ma <- merge(ma, pop_full, by = c("id", "y"), all.x = TRUE)

saveRDS(ma, file = paste0(loc.output, "bites_spain_daily.rds"))
