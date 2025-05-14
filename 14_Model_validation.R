###################### Model Validation: MA on BG model ########################
library(tidyverse)
library(sf)
library(parallel)
library(dplyr)
library(lubridate)

rm(list = ls())

sf::sf_use_s2(FALSE)
# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")

# Load information -------------------------------------------------------------
val_data <- read_csv(paste0(loc.data, "ES_albo_info.csv")) %>%
  janitor::clean_names() %>%
  mutate(
    presence = TRUE
  ) %>%
  # filter(quien_fue_antes %in% c(-2023, -2024, 2023, 2024) & meth == "FSUR") %>%
  dplyr::select(municipio_nombre, natcode, presence, meth, quien_fue_antes)
spain <- mapSpain::esp_get_munic() %>%
  janitor::clean_names() %>%
  mutate(
    # natcode = 34 + codauto + cpro  + lau_code
    natcode = paste0("34", codauto, cpro, lau_code), 
    id = paste0(codauto, cpro, cmun, lau_code) 
  ) %>%
  st_drop_geometry()
val_data <- merge(val_data, spain, by = "natcode", all.x = TRUE) %>%
  dplyr::select(id, meth, presence, quien_fue_antes)
rm(spain)

val_data <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  ) %>% dplyr::select(id) %>%
  left_join(val_data, by = "id", keep = FALSE) %>%
  mutate(
    presence = ifelse(is.na(presence), FALSE, presence)
  )

# Load count predictions from COUNT model
mdl_name <- "/mtiger16"
fldr <- "Counts"
sub <- "" # with _

years <- c("2020", "2021", "2022")
months <- c("01", "02", "03", "04", "05", "06" ,"07", "08", "09", "10", "11", "12")

month_values <- list()
iter <- 0
for(y in years){
  for (m in months){
    print(paste0("Plotting: ", m, "-", y))
    iter <- iter + 1
    
    pred <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr, mdl_name, "/tiger_", m ,"_", y, sub, ".rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    
    pred$pred_count <- rowMeans(pred[,6:ncol(pred)], na.rm = TRUE)
    
    pred <- pred[c("municipality","id", "prov_name", "lon", "lat", "pred_count")]
    pred <- pred %>% mutate(
      month = m,
      year = y
    )
    month_values[[iter]] <- pred
  }
}
pred_counts <- do.call(rbind, month_values) 
rm(month_values, pred, iter, m, y, years, months)

# Load count predictions from INTE model
fldr = "Monthly_Integration/"
mdl_name = "tiger_inter_monthly_only_ma"
sub = ""

years <- c("2020", "2021", "2022")
months <- c("01", "02", "03", "04", "05", "06" ,"07", "08", "09", "10", "11", "12")

month_values <- list()
iter <- 0

for(y in years){
  for (m in months){
    print(paste0("Procesing: ", m, "-", y))
    iter <- iter + 1
    
    pred <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr, mdl_name, "/tiger_", m ,"_", y, sub, ".rds")) %>%
      janitor::clean_names() %>%
      mutate(
        month = m,
        year = y
      ) %>%
      st_drop_geometry()
    
    colnames(pred) <- c("municipality","id", "prov_name", "lon", "lat", "pred_count", "month", "year") 
    
    month_values[[iter]] <- pred 
  }
}
pred_inte <- do.call(rbind, month_values) 
rm(month_values, pred, iter, m, y, years, months)

# Validation with Rogers' information ------------------------------------------
# Global validation
pred_counts_global <- pred_counts %>%
  group_by(municipality, id, prov_name, lon, lat) %>%
  summarise(
    counts = mean(pred_count),
    model = "Counts"
  ) %>% 
  ungroup()
pred_inte_global <- pred_inte %>%
  group_by(municipality, id, prov_name, lon, lat) %>%
  summarise(
    counts = mean(pred_count),
    model = "Integrated"
  ) %>% ungroup()

global <- rbind(pred_counts_global, pred_inte_global)
global <- merge(global, val_data, by = "id", all.x = TRUE)

# global <- global %>% filter(quien_fue_antes %in% c(-2024, 2024))

custom_colors <- c("Integrated" = "#abdda4", "Counts" = "#3288bd")
ggplot(global) +
  geom_boxplot(aes(y = counts, x = presence, fill = model), outlier.shape = NA) +
  scale_fill_manual(values = custom_colors, labels = c("COUNT", "INTEGRATED")) +
  labs(
    x = "Presence of Ae. albopictus",
    y = "Predicted Abundance Estimates",
    fill = "Model") +
  theme_classic() 

