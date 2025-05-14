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
# spain <- mapSpain::esp_get_munic_siane() %>%
#   st_transform(4326) %>%
#   filter(!(ine.ccaa.name %in% c("Canarias", "Ceuta", "Melilla"))) %>%
#   mutate(
#     id = paste0(codauto, cpro, cmun, LAU_CODE)
#   ) %>%
#   ungroup() %>%
#   dplyr::select(name, id, ine.prov.name) %>%
#   rename(municipality = "name", prov_name = "ine.prov.name")
# saveRDS(spain, file = paste0(loc.output, "spain_mun.rds"))
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  )

# Checking all municipalities
ggplot() +
  geom_sf(data = spain, fill = "black", color = "transparent",
          size = 0.01, alpha = 0.8, na.rm = TRUE) 

# Download Sampling effort data -----------------------------------------------
trs_daily <- read_csv("https://github.com/Mosquito-Alert/sampling_effort_data/raw/main/sampling_effort_daily_cellres_025.csv.gz") %>%
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

# # Checking the selected points
# plot(st_geometry(spain))
# plot(st_geometry(trs_daily), add = TRUE)

# BORRAR: EN CASO DE QUE SE CALCULE DIARIAMENTE
# trs_monthly = trs_daily %>%
#   rename(SEev = SE_expected) %>%
#   st_drop_geometry() %>%
#   group_by(m, y, municipality, id, geometry) %>%
#   summarize(
#     SEev = sum(SEev),
#     SE = 1-prod(1-SE), # This comes directly fom John's script: one usser at least send one report
#     n_total_reports = sum(n_total_reports),
#     n_total_reporters = sum(n_total_reporters)
#   ) %>%
#   ungroup() # n = 12624
# st_geometry(trs_monthly) <- "geometry"

trs_daily = trs_daily %>%
  rename(SEev = SE_expected) %>%
  st_drop_geometry() %>%
  group_by(date, municipality, id, geometry) %>%
  summarize(
    SEev = sum(SEev),
    SE = 1-prod(1-SE), # This comes directly fom John's script: one usser at least send one report
    n_total_reports = sum(n_total_reports),
    n_total_reporters = sum(n_total_reporters)
  ) %>%
  ungroup() # n = 42306
st_geometry(trs_daily) <- "geometry"

trs_random <- trs_daily %>% 
  mutate(
    date =  sample(seq(as.Date("2020-01-01"), as.Date("2022-12-31"), by = "day"), 
                   nrow(trs_daily), replace = TRUE),
    id = sample(spain$id, nrow(trs_daily), replace = TRUE),
    municipality = sample(spain$municipality, nrow(trs_daily), replace = TRUE)
    ) %>% # random dates and SE
  st_drop_geometry()
trs_random <- merge(trs_random, spain %>% dplyr::select(id), by = "id")

# Checking SE selection
# ggplot() +
#   geom_sf(data = spain, color = "grey", fill = NA, size = 0.005, alpha = 0.05, na.rm = TRUE) +
#   geom_sf(data = trs_monthly %>%
#             filter(month == 3 & year == 2020), 
#           aes(fill = SE), color = "transparent", size = 0.1, alpha = 0.5, na.rm = TRUE) +
#   scale_fill_distiller("", palette = "Spectral") +
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

# We have to add a new random effect: the municipality!
vrs_full <- st_join(vrs_full, spain, join = st_intersects) %>% st_drop_geometry()

# Grouping reports by pixel ID
vrs_full <- vrs_full %>%
  st_drop_geometry() %>%
  # mutate(
  #   m = as.integer(month(date)), 
  #   y = year(date)
  # ) %>% 
  # group_by(y, m, municipality, id) %>% 
  group_by(date, municipality, id) %>% 
  summarize(n_target_reps = n(), .groups ="drop") %>%
  ungroup() # monthly 5841 rows // daily 12071

# Joining reports and Sampling Effort
# BORRAR: PARA MES
# ma_df = trs_monthly %>% 
#   left_join(vrs_full, by = c("municipality", "id", "y", "m")) %>%
#   replace_na(list(n_target_reps=0)) %>%
#   filter(SE>0) %>% 
#   mutate(
#     any_reps = n_target_reps > 0,
#     reps_per_SEev = n_target_reps/(SEev),
#     rep_per_reporter = n_total_reports/n_total_reporters
#   ) %>% 
#   mutate(SE_log_odds = log( (SE+.0001)/(1-SE+.0001))) # 12566

ma_df = trs_daily %>% 
  left_join(vrs_full, by = c("municipality", "id", "date")) %>%
  replace_na(list(n_target_reps=0)) %>%
  filter(SE>0) %>% 
  mutate(
    any_reps = n_target_reps > 0,
    reps_per_SEev = n_target_reps/(SEev),
    rep_per_reporter = n_total_reports/n_total_reporters
  ) %>% 
  mutate(SE_log_odds = log( (SE+.0001)/(1-SE+.0001))) # 42025

ma_df = trs_random %>%
  left_join(vrs_full, by = c("municipality", "id", "date")) %>%
  replace_na(list(n_target_reps=0)) %>%
  filter(SE>0) %>%
  mutate(
    any_reps = n_target_reps > 0,
    reps_per_SEev = n_target_reps/(SEev),
    rep_per_reporter = n_total_reports/n_total_reporters
  ) %>%
  mutate(SE_log_odds = log( (SE+.0001)/(1-SE+.0001))) # 42025

# adding app phase 
# app_phases = tibble(app_phase = paste0("p", 1:7), 
#                     date = c(as_date("2014-06-01"), 
#                              as_date("2014-10-01"), 
#                              as_date("2015-07-01"), 
#                              as_date("2016-10-01"), 
#                              as_date("2017-04-01"), 
#                              as_date("2019-02-01"), 
#                              as_date("2020-10-01"))) %>% 
#   complete(date = seq.Date(as_date("2014-06-01"), today(), by="month")) %>% 
#   fill(app_phase) %>% 
#   mutate(m=month(date), y=year(date))
# 
# ma_df = ma_df %>% 
#   left_join(app_phases %>% dplyr::select(app_phase, m, y))

ma_df <-  ma_df %>% mutate(
  app_phase = case_when(
    (date >= as_date("2014-06-01") & date < as_date("2014-10-01")) ~ "p1",
    (date >= as_date("2014-10-01") & date < as_date("2015-07-01")) ~ "p2",
    (date >= as_date("2015-07-01") & date < as_date("2016-10-01")) ~ "p3",
    (date >= as_date("2016-10-01") & date < as_date("2017-04-01")) ~ "p4",
    (date >= as_date("2017-04-01") & date < as_date("2019-02-01")) ~ "p5",
    (date >= as_date("2019-02-01") & date < as_date("2020-10-01")) ~ "p6",
    (date >= as_date("2020-10-01")) ~ "p7")
  )

# Adding other information -----------------------------------------------------
# Extracting climatic variables from nc raster ---------------------------------
# centroids <- st_centroid(ma_df) %>%
#   mutate(longitude = st_coordinates(.)[,1], latitude = st_coordinates(.)[,2]) %>%
#   st_drop_geometry() %>%
#   st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
# plot(st_geometry(centoids))

# Full function ----------------------------------------------------------------
# Loading the functions
# In local
# source("read_and_extract_grib_data.R")


# In local
source("read_and_extract_nc_data_daily.R")
# In cluster
source("/home/usuaris/ccerecedo/Spain_Tiger/read_and_extract_nc_data_daily.R")

ma_df <- ma_df %>%
  mutate(
    start_date = date,
    end_date = date
  )
st_geometry(ma_df) <- "geometry"

# In cluster
# source("/home/usuaris/ccerecedo/EU_Culex/read_and_extract_grib_data.R")
cat("------------------------- Number of rows:", nrow(ma_df), "\n")

weather_data <- mclapply(1:nrow(ma_df), function(i){
  
  cat(paste0("Number row:", i, "\n"))
  
  data_row <- ma_df[i, ]
  
  ex_wt <- extract_weather(data_row)
  
  wth <- calculating_weather(ex_wt)
}
, mc.cores = 8)

ma_wth <- do.call(rbind, weather_data) 

# saveRDS(ma_wth, file = paste0(loc.output, "ma_tiger_spain_climatic_variables_pixel_daily.rds"))
saveRDS(ma_wth, file = paste0(loc.output, "ma_tiger_spain_random.rds"))

# Adding Corine Land Covers ---------------------------------------------------
clc_surface <- readRDS(paste0(loc.output, "clc_surface_mun_level_0.rds"))

ma <- merge(ma_wth, clc_surface, by = c("municipality", "id"), all.x = TRUE) # N = 12566 // 42025
ma$no_data <- NULL

saveRDS(ma, file = paste0(loc.output, "ma_tiger_spain_daily_random.rds"))

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

saveRDS(ma, file = paste0(loc.output, "ma_tiger_spain_daily.rds"))
