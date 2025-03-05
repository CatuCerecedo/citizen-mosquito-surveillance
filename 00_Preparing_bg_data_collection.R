######################## Historical weather collection #########################

#' Here we will clean the raw data, following Simone and Laura recommendations
#' WARNING: the orinal csv must not be changed

# if (!require("tidyverse")) install.packages("tidyverse")
# if (!require("dplyr")) install.packages("dplyr")
# if (!require("lubridate")) install.packages("lubridate")
# if (!require("janitor")) install.packages("janitor")
# if (!require("readxl")) install.packages("readxl")
# if (!require("units")) install.packages("units")

library(tidyverse)
library(dplyr)
library(lubridate)
library(janitor)
library(readxl)
library(units)
library(data.table)
library(sf)
library(ggplot2)
library(mapSpain)

rm(list = ls())
# Directories ------------------------------------------------------------------

loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")

# Girona -----------------------------------------------------------------------

gi_202 <- read_excel(paste0(loc.data, 'Surveillance2020_2022_CEAB.xlsx'),
                       skip = 3, col_names = TRUE , col_types = c("guess", "guess","guess",
                                                                  "guess", "guess", "guess", "guess","guess", "guess",
                                                                  "guess", "guess", "guess", "guess", "guess", "guess",
                                                                  "date", "guess", "guess", "date","guess", "guess",
                                                                  "guess", "guess", "guess", "guess", "guess","numeric",
                                                                  "guess", "guess", "guess", "guess")) %>%
  clean_names() %>%
  filter(!(comment_sampling %in% c("ants, no bodies",
                                   "Ants, no mosquitoes",
                                   "Ants; insect  bodies",
                                   "Bag not collected (no mosquitoes)",
                                   "Diptera, no culcidae",
                                   "Many mosquito wings and legs",
                                   "No data",
                                   "tapa bloqueada"))) %>%
  filter(!(trap_status %in% c("Bag not collected", 
                              "Bag punctured", 
                              "data temporarily unavailable", 
                              "trap not working",
                              "trap fallen",
                              "Funnel door closed",
                              "Funnel upside down",
                              "Funnel out",
                              "Working 12 h/day only", 
                              "Unplugged",
                              "Trap upside down", 
                              "Trap unplugged several hours",
                              "Trap unplugged",
                              "Trap not working",
                              "trap not working",
                              "Trap off (cable cut)",
                              "Trap not connected",
                              "Trap lying down",
                              "Trap laying down",
                              "Trap fallen",
                              "Trap door faulty",
                              "Rubbish in the funnel",
                              "No trap door until june 6th",
                              "No trap",
                              "No bag",
                              "lost, won't be replaced",
                              "Last week; trap disconnected"
                              ))) %>% 
  filter(!(comment_trap %in% c("Heavy rain; Trap unistalled",
                               "no trap funnel",
                               "Switched off",
                               "trap's door closed",
                               "trap's door half-open (heavy rain)",
                               "Trap unplugged"
                               ))) %>% # Laura and Simone told me to do this filtration
  filter(!(city %in% c("Badalona"))) %>% # Simone: eliminate the Badalona traps
  filter(trapping_effort > 0) %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, land_use,
                start_date, end_date, trapping_effort, species, females) %>%
  transform(latitude = as.numeric(latitude),
            longitude = as.numeric(longitude)) %>%
  pivot_wider(id_cols = c(city, province, trap_name, latitude, longitude, trap_type, land_use,
                          start_date, end_date, trapping_effort),
              names_from = species, values_from = females) %>%
  clean_names() %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, aedes_albopictus) %>%
  rename(females = aedes_albopictus) %>% 
  filter(females != "Not available") %>%
  mutate(
    females = as.numeric(ifelse(is.na(females), 0, females))
  ) %>% 
  group_by(trap_name) %>% 
  arrange(end_date) %>% 
  ungroup() %>%
  mutate(
    trapping_effort = as.numeric(end_date - start_date)
  ) 

gi_2023 <- read_excel(paste0(loc.data, 'VEO_Surveilance2023_Data report sheet_CEAB.xlsx'),
                     skip = 3, sheet=1, col_names = TRUE , col_types = c("guess", "guess","guess",
                                                                "guess", "guess", "guess", "guess","guess", "guess",
                                                                "guess", "guess", "guess", "guess", "guess", "guess",
                                                                "date", "guess", "guess", "date","guess", "guess",
                                                                "guess", "guess", "guess", "guess", "guess","numeric",
                                                                "guess", "guess", "guess")) %>%
  clean_names() %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date) %>%
  mutate(females = NA) %>%
  group_by(trap_name) %>% 
  arrange(end_date) %>% 
  ungroup %>%
  mutate(
    trapping_effort = as.numeric(end_date - start_date)
  ) 

# Barcelona --------------------------------------------------------------------
bcn_2020 <- read_excel(paste0(loc.data, '2020_BG_BarcelonaCity.xlsx'),
                       skip = 3,col_names = TRUE , col_types = c("guess", "guess","guess",
                                                                 "guess", "guess", "guess", "guess","guess", "guess",
                                                                 "guess", "guess", "guess", "guess", "guess", "guess",
                                                                 "date", "guess", "guess", "date","guess", "guess",
                                                                 "guess", "guess", "guess", "guess", "guess","numeric",
                                                                 "guess", "guess", "guess", "guess", "guess")) %>%
  clean_names() %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, species, females, comments_classification) %>%
  pivot_wider(id_cols = c(city, province, trap_name, latitude, longitude, trap_type, 
                          start_date, end_date, comments_classification),
              names_from = species, values_from = females, values_fill = 0) %>% # Adding ceros when NAs
  clean_names() %>% 
  filter(comments_classification > 1) %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
           start_date, end_date, aedes_albopictus) %>%
  rename(females = aedes_albopictus) %>% 
  group_by(trap_name) %>% 
  arrange(end_date) %>% 
  ungroup() %>%
  mutate(
    trapping_effort = as.numeric(end_date - start_date)
  ) 

bcn_2021 <- read_excel(paste0(loc.data, '2021_BG_BarcelonaCity_Volunteers.xlsx'),
                       skip = 3,col_names = TRUE , col_types = c("guess", "guess","guess",
                                                                 "guess", "guess", "guess", "guess","guess", "guess",
                                                                 "guess", "guess", "guess", "guess", "guess", "guess",
                                                                 "date", "guess", "guess", "date","guess", "guess",
                                                                 "guess", "guess", "guess", "guess", "guess","numeric",
                                                                 "guess", "guess", "guess", "guess", "guess")) %>%
  clean_names() %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, species, females, comments_classification) %>%
  pivot_wider(id_cols = c(city, province, trap_name, latitude, longitude, trap_type, 
                          start_date, end_date, comments_classification),
              names_from = species, values_from = females, values_fill = 0) %>%
  clean_names() %>% 
  filter(comments_classification > 1) %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, aedes_albopictus) %>%
  rename(females = aedes_albopictus) %>% 
  group_by(trap_name) %>% 
  arrange(end_date) %>% 
  ungroup() %>%
  mutate(
    trapping_effort = as.numeric(end_date - start_date)
  ) 

bcn_2021b <- read_excel(paste0(loc.data, '2021_BG_BarcelonaCity.xlsx'),
                        skip = 3,col_names = TRUE , col_types = c("guess", "guess","guess",
                                                                  "guess", "guess", "guess", "guess","guess", "guess",
                                                                  "guess", "guess", "guess", "guess", "guess", "guess",
                                                                  "date", "guess", "guess", "date","guess", "guess",
                                                                  "guess", "guess", "guess", "guess", "guess","numeric",
                                                                  "guess", "guess", "guess", "guess", "guess")) %>%
  clean_names() %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type,
                start_date, end_date, species, females, comments_classification) %>%
  pivot_wider(id_cols = c(city, province, trap_name, latitude, longitude, trap_type,
                          start_date, end_date, comments_classification),
              names_from = species, values_from = females, values_fill = 0) %>%
  clean_names() %>%
  filter(comments_classification > 1) %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type,
                start_date, end_date, aedes_albopictus) %>%
  rename(females = aedes_albopictus) %>%
  group_by(trap_name) %>%
  arrange(end_date) %>%
  ungroup() %>%
  mutate(
    trapping_effort = as.numeric(end_date - start_date)
  )

bcn_2022 <- read_excel(paste0(loc.data, '2022_BG_Mtego_BarcelonaCity.xlsx'),
                       skip = 3,col_names = TRUE , col_types = c("guess", "guess","guess",
                                                                 "guess", "guess", "guess", "guess","guess", "guess",
                                                                 "guess", "guess", "guess", "guess", "guess", "guess",
                                                                 "date", "guess", "guess", "date","guess", "guess",
                                                                 "guess", "guess", "guess", "guess", "guess","numeric",
                                                                 "guess", "guess", "guess", "guess", "guess")) %>%
  clean_names() %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, species, females, comments_classification) %>%
  pivot_wider(id_cols = c(city, province, trap_name, latitude, longitude, trap_type, 
                          start_date, end_date, comments_classification),
              names_from = species, values_from = females, values_fill = 0) %>%
  clean_names() %>% 
  filter(comments_classification > 1, trap_type != "Mtego") %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, aedes_albopictus) %>%
  rename(females = aedes_albopictus) %>% 
  group_by(trap_name) %>% 
  arrange(end_date) %>% 
  ungroup() %>%
  mutate(
    trapping_effort = as.numeric(end_date - start_date)
  ) 

bcn_2023 <- read_excel(paste0(loc.data, 'BCN_BGs_AIMCOST_Surveilance2023_Data report sheet_sampling_effort_sent.xlsx'),
                       sheet = 2, skip = 3,col_names = TRUE , col_types = c("guess", "guess","guess",
                                                                            "guess", "guess", "guess", "guess","guess", "guess",
                                                                            "guess", "guess", "guess", "guess", "guess", "guess",
                                                                            "date", "guess", "guess", "date","guess", "guess",
                                                                            "guess", "guess", "guess", "guess", "guess","numeric",
                                                                            "guess", "guess", "guess", "guess", "guess")) %>%
  clean_names() %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, species, females, comments_classification) %>%
  mutate(females = NA) %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, females) %>%
  group_by(trap_name) %>% 
  arrange(end_date) %>% 
  ungroup() %>%
  mutate(
    trapping_effort = as.numeric(end_date - start_date)
  ) 
# Barlearic islands (CB) -------------------------------------------------------
  
cb <- read_excel(paste0(loc.data, 'Plantilla_datos_BG_CB.xlsx'),
                       skip = 3, col_names = TRUE , col_types = c("guess", "guess","guess",
                                                                 "guess", "guess", "guess", "guess","guess", "guess",
                                                                 "guess", "guess", "guess", "guess", "guess", "guess",
                                                                 "date", "guess", "guess", "date","guess", "guess",
                                                                 "guess", "guess", "guess", "guess", "guess","numeric",
                                                                 "guess", "guess", "guess", "guess", "guess", "guess")) %>%
  clean_names() %>%
  mutate(
    females = females + bf # WARNING: two columns with females values. Isis me lo tiene que confirmar
  ) %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, species, females) %>%
  pivot_wider(id_cols = c(city, province, trap_name, latitude, longitude, trap_type, 
                          start_date, end_date), # Warning: data validity has different values on the same day for different species. No data valididty < 1
              names_from = species, values_from = females, values_fill = 0) %>%
  clean_names() %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, aedes_albopictus) %>%
  rename(females = aedes_albopictus) %>% 
  group_by(trap_name) %>% 
  arrange(end_date) %>% 
  ungroup() %>%
  mutate(
    trapping_effort = as.numeric(end_date - start_date)
  ) 

mb <- read_excel(paste0(loc.data, 'Plantilla_datos_BG_MB.xlsx'),
                 skip = 3, col_names = TRUE , col_types = c("guess", "guess","guess",
                                                            "guess", "guess", "guess", "guess","guess", "guess",
                                                            "guess", "guess", "guess", "guess", "guess", "guess",
                                                            "date", "guess", "guess", "date","guess", "guess",
                                                            "guess", "guess", "guess", "guess", "guess","numeric",
                                                            "guess", "guess", "guess", "guess", "guess")) %>%
  clean_names() %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, species, females) %>% # Warning: data validity has different values on the same day for different species. No data valididty < 1
  pivot_wider(id_cols = c(city, province, trap_name, latitude, longitude, trap_type, 
                          start_date, end_date),
              names_from = species, values_from = females, values_fill = 0) %>%
  clean_names()  %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, aedes_albopictus) %>%
  rename(females = aedes_albopictus) %>% 
  group_by(trap_name) %>% 
  arrange(end_date) %>% 
  ungroup() %>%
  mutate(
    trapping_effort = as.numeric(end_date - start_date)
  ) 

mg <- read_excel(paste0(loc.data, 'Plantilla_datos_BG_MG.xlsx'),
                 skip = 3, col_names = TRUE , col_types = c("guess", "guess","guess",
                                                            "guess", "guess", "guess", "guess","numeric", "numeric",
                                                            "guess", "guess", "guess", "guess", "guess", "guess",
                                                            "date", "guess", "guess", "date","guess", "guess",
                                                            "guess", "guess", "guess", "guess", "guess","numeric",
                                                            "guess", "guess", "guess", "guess", "guess")) %>%
  clean_names() %>%
  slice(1:68) %>% # There are a lot of rows at the end of the df with NAs
  mutate(
    trap_type = ifelse(trap_type != "BG sentinel + BG lure7", "BG sentinel + BG lure7", trap_type) # WARNING: there are different trap_types, seems a mistake
  ) %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, species, females, data_validity) %>% # Warning: data validity has different values on the same day for different species. No data valididty < 1
  pivot_wider(id_cols = c(city, province, trap_name, latitude, longitude, trap_type, 
                          start_date, end_date, data_validity),
              names_from = species, values_from = females, values_fill = 0) %>%
  clean_names()  %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, aedes_albopictus) %>%
  rename(females = aedes_albopictus) %>% 
  group_by(trap_name) %>% 
  arrange(end_date) %>% 
  ungroup() %>%
  mutate(
    trapping_effort = as.numeric(as_date(end_date) - as_date(start_date)),
    latitude = ifelse(is.na(latitude), 39.5604, latitude),
    longitude = ifelse(is.na(longitude), 4.08055, longitude)
  ) 

# Extremadura (db) -------------------------------------------------------------

db <- read_excel(paste0(loc.data, 'Plantilla_datos_BG_DB.xlsx'),
                 skip = 3, col_names = TRUE , col_types = c("guess", "guess","guess",
                                                            "guess", "guess", "guess", "guess","guess", "guess",
                                                            "guess", "guess", "guess", "guess", "guess", "guess",
                                                            "date", "guess", "guess", "date","guess", "guess",
                                                            "guess", "guess", "guess", "guess", "guess","numeric",
                                                            "guess", "guess", "guess", "guess", "guess")) %>%
  clean_names() %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, species, females) %>% # Warning: data validity has different values on the same day for different species. No data valididty < 1
  pivot_wider(id_cols = c(city, province, trap_name, latitude, longitude, trap_type, 
                          start_date, end_date),
              names_from = species, values_from = females, values_fill = 0) %>%
  clean_names() %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date) %>%
  mutate(females = 0) %>% # There are not A. albopictus
  group_by(trap_name) %>% 
  arrange(end_date) %>% 
  ungroup %>%
  mutate(
    trapping_effort = as.numeric(end_date - start_date)
  ) 

# Logroño ----------------------------------------------------------------------

ir <- read_excel(paste0(loc.data, 'Plantilla_datos_BG_IR.xlsx'),
                 skip = 3, col_names = TRUE , col_types = c("guess", "guess","guess",
                                                            "guess", "guess", "guess", "guess","guess", "guess",
                                                            "guess", "guess", "guess", "guess", "guess", "guess",
                                                            "date", "guess", "guess", "date","guess", "guess",
                                                            "guess", "guess", "guess", "guess", "guess","numeric",
                                                            "guess", "guess", "guess", "guess", "guess")) %>%
  clean_names() %>%
  mutate(
    trap_type = ifelse(trap_type == "BG sentinel + BG lure7", "BG sentinel + BG lure 7", trap_type)
  ) %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, species, females) %>% # Warning: data validity has different values on the same day for different species. No data valididty < 1
  pivot_wider(id_cols = c(city, province, trap_name, latitude, longitude, trap_type, 
                          start_date, end_date),
              names_from = species, values_from = females, values_fill = 0) %>%
  clean_names() %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date) %>%
  mutate(females = 0) %>% # There are not A. albopictus
  group_by(trap_name) %>% 
  arrange(end_date) %>% 
  ungroup %>%
  mutate(
    trapping_effort = as.numeric(end_date - start_date)
  ) 

# Valencia ---------------------------------------------------------------------

sd <- read_excel(paste0(loc.data, 'Plantilla_datos_BG_SD .xlsx'),
                 skip = 3, col_names = TRUE , col_types = c("guess", "guess","guess",
                                                            "guess", "guess", "guess", "guess","guess", "guess",
                                                            "guess", "guess", "guess", "guess", "guess", "guess",
                                                            "date", "guess", "guess", "date","guess", "guess",
                                                            "guess", "guess", "guess", "guess", "guess","numeric",
                                                            "guess", "guess", "guess", "guess", "guess")) %>%
  clean_names() %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
              start_date, end_date, species, females, data_validity) %>% # Warning: data validity has different values on the same day for different species. No data valididty < 1
  pivot_wider(id_cols = c(city, province, trap_name, latitude, longitude, trap_type, 
                          start_date, end_date, data_validity),
              names_from = species, values_from = females, values_fill = 0) %>%
  clean_names()  %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, aedes_albopictus) %>%
  rename(females = aedes_albopictus) %>% 
  group_by(trap_name) %>% 
  arrange(end_date) %>% 
  ungroup %>%
  mutate(
    trapping_effort = as.numeric(as_date(end_date) - as_date(start_date))
  ) 

pa <- read_excel(paste0(loc.data, 'Plantilla_datos_BG_PA.xlsx'),
                 skip = 3, col_names = TRUE , col_types = c("guess", "guess","guess",
                                                            "guess", "guess", "guess", "guess","guess", "guess",
                                                            "guess", "guess", "guess", "guess", "guess", "guess",
                                                            "date", "guess", "guess", "date","guess", "guess",
                                                            "guess", "guess", "guess", "guess", "guess","numeric",
                                                            "guess", "guess", "guess", "guess", "guess")) %>%
  clean_names() %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, species, females, data_validity) %>% # Warning: data validity has different values on the same day for different species. No data valididty < 1
  pivot_wider(id_cols = c(city, province, trap_name, latitude, longitude, trap_type, 
                          start_date, end_date, data_validity),
              names_from = species, values_from = females, values_fill = 0) %>%
  clean_names()  %>%
  dplyr::select(city, province, trap_name, latitude, longitude, trap_type, 
                start_date, end_date, aedes_albopictus) %>%
  rename(females = aedes_albopictus) %>% 
  group_by(trap_name) %>% 
  arrange(end_date) %>% 
  ungroup %>%
  mutate(
    trapping_effort = as.numeric(as_date(end_date) - as_date(start_date))
  ) 


bg_traps <- rbind(gi_202, gi_2023, bcn_2020, bcn_2021, bcn_2021b, bcn_2022, bcn_2023, cb, mb, mg, db, ir, sd, pa)
rm(gi_202, bcn_2020, bcn_2021, bcn_2021b, bcn_2022, cb, mb, mg, db, ir, sd, pa, bcn_2023, gi_2023)

# bg_traps = rbind(bg_traps, bcn)
# There are some NAs on start and enda dates

bg_traps <- bg_traps %>%
  drop_na(start_date, end_date)


