################# Preparing Daily Mosquito Alert Data ##########################

library(tidyverse)
library(sf)
library(terra)
library(parallel)
library(data.table)

rm(list = ls())
# Directories ------------------------------------------------------------------

# In cluster
loc.output <- paste0(getwd(), "/Spain_Tiger/OUTPUT/")
loc.data <- paste0(getwd(), "/Spain_Tiger/DATA/")
loc.era5 <- paste0(getwd(), "/EU_Culex/ERA5_Download/")

# In local
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.era5 <- "/home/catuxa/Documents/Mosquito_Models/EU_Culex/ERA5_Download/"

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

# We have to add a new random effect: the municipality!
tiger <- st_join(tiger, spain, join = st_intersects) %>% st_drop_geometry()

# BORRAR:
# Grouping counts by municipality and month
# tiger <- tiger %>%
#   mutate(
#     m = as.integer(month(end_date)), 
#     y = year(end_date)
#   ) %>% 
#   group_by(y, m, municipality, id) %>% 
#   summarize(
#     females = sum(females, na.rm = TRUE),
#     trapping_effort = sum(trapping_effort, na.rm = TRUE),
#     n_traps = n(),
#     .groups ="drop") %>%
#   ungroup() # 256 rows
sum(duplicated(tiger))

# removing duplicates
tiger <- tiger[!duplicated(tiger), ]

plot(table(tiger$females)) # Negative binomial

tiger <- merge(tiger, spain, by = c("municipality", "id")) # 256 // 2816
st_geometry(tiger) <- "geometry"

# WARNING there are trapping effor = 0
tiger <- tiger %>% filter(trapping_effort != 0)

# Full function ----------------------------------------------------------------
# Loading the functions
# In local
# source("read_and_extract_grib_data.R")
# In cluster
# source("/home/usuaris/ccerecedo/EU_Culex/read_and_extract_grib_data.R")

# Loading the functions
# In local
source("read_and_extract_nc_data_daily.R")
# In cluster
# source("/home/usuaris/ccerecedo/EU_Culex/read_and_extract_nc_data.R")
cat("------------------------- Number of rows:", nrow(tiger), "\n")

weather_data <- mclapply(1:nrow(tiger), function(i){
  
  message(paste0("Number row:", i, "\n"))
  
  data_row <- tiger[i, ]
  
  ex_wt <- extract_weather(data_row)
  
  wth <- calculating_weather(ex_wt)
}
, mc.cores = 8)

bg_wth <- do.call(rbind, weather_data) 

saveRDS(bg_wth, file = paste0(loc.output, "bg_tiger_spain_climatic_variables_mun_daily.rds"))

# Adding Corine Land Covers ---------------------------------------------------
clc_surface <- readRDS(paste0(loc.output, "clc_surface_mun_level_0.rds"))

bg <- merge(bg_wth, clc_surface, by = c("municipality", "id"), all.x = TRUE)
bg$no_data <- NULL

saveRDS(bg, file = paste0(loc.output, "bg_tiger_spain_daily.rds"))

# Adding population density ----------------------------------------------------
bg <- bg %>% mutate(y = as.factor(year(end_date)))

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

bg <- merge(bg, pop_full, by = c("id", "y"), all.x = TRUE)

saveRDS(bg, file = paste0(loc.output, "bg_tiger_spain_daily.rds"))
