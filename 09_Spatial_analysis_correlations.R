################ Spatial analysis - spatial ecology ############################

library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(tidyverse)
library(ggpubr)

rm(list = ls())

# Directories ------------------------------------------------------------------
# In local
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")
loc.era5 <- "/home/catuxa/Documents/Mosquito_Models/EU_tiger/ERA5_Download/"

# Loading spatial correlations -------------------------------------------------
# In the case of Count vs MA predictions
spatial_corr <- readRDS(paste0(loc.output, "spatial_corr_bg_ma.rds"))

# In the case of Suit vs MA predictions
# spatial_corr <- readRDS(paste0(loc.output, "spatial_corr_suit_ma.rds"))

# Loading spatial variables ----------------------------------------------------
# Country
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  )

## CLC data --------------------------------------------------------------------
clc_surface <- readRDS(paste0(loc.output, "clc_surface_mun_level_0.rds"))

# Relationships with each clc category
df <- merge(spatial_corr[["2022"]], clc_surface, by = c("municipality", "id"))

# Selecting columns to plot
column_names <- names(clc_surface)[3:16]

df_long <- df %>%
  select(rho, all_of(column_names)) %>%
  pivot_longer(cols = all_of(column_names), names_to = "CLC", values_to = "Value") %>%
  filter(CLC != "no_data")

ggplot(df_long, aes(x = Value, y = rho)) +
  geom_point(size = 2) +
  geom_smooth(aes(x = Value, y = rho), method = "loess", se = TRUE, color = "darkgrey", size = 0.8) +
  facet_wrap(~ CLC, scales = "free_x") +
  labs(
    x = "Percentage of CORINE Land Covers",
    y = "Spatial Spearman Rank Correlation (S)",
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

##### Maybe with a PCA???

## Population data -------------------------------------------------------------
# Loading data

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

# Plotting relationships 
# Creating data with values from all years
df_full <- data.frame() 
for (year in c("2020", "2021", "2022")){
  # Filter spatial correlation by year
  df <- merge(spatial_corr[[year]], spain, by = c("municipality", "id", "prov_name"))
  df$y <- year
  st_geometry(df) <- "geometry"
  
  # Filter population data by year
  pop <- pop_full %>% filter(y == year)
  
  # Marge by year
  df <- merge(df, pop, by = c("id", "y"), all.x = TRUE)
  
  # Joining tables
  df_full <- rbind(df_full, df) 
}

custom_colors <- c("2020" = "#abdda4", "2021" = "#d53e4f", "2022" = "#3288bd")
  
ggplot(df_full, aes(x = dens, y = rho,  color = y, fill = y)) +
  geom_point(size = 4, alpha = 0.6) +
  geom_smooth(method = lm) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "\nPopulation density (Nº of people/km²)",
    y = "Spearman Spatial Correlation (S)\n",
    color = "Year",
    fill = "Year"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",    
    legend.title = element_text(size = 22, face = "bold", margin = margin(b=20)),
    legend.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(hjust = 2, size = 20),
    legend.key.spacing.y = unit(0.5, 'cm')
  ) 




