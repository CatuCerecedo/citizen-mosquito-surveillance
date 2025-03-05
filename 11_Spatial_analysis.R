############################## Spatial analysis ################################
library(ggplot2)
library(tidyverse)
library(sf)
library(parallel)
library(spdep)
library(sfdep)

rm(list = ls())

# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")
loc.heavy <- paste0(getwd(), "/Heavy_files/")
loc.clc <- paste0(getwd(), "/OUTPUT/land_cover_raster/")

# Model names ------------------------------------------------------------------
mdl <- "mtiger16"
mdl_name <- "/mtiger16"
fldr <- "Counts"
sub <- "" # with _

# mdl <- "mtiger3_occu"
# mdl_name <- "/mtiger3_occu"
# fldr <- "Suitability"
# sub <- "_occu" # with _

# mdl <- "mbites_7_reduce"
# mdl_name <- "/mbites_7_reduce"
# fldr <- "Bites"
# sub <- "bts" # with _

mdl_ma <- "mtiger7_ma"
mdl_name_ma <- "/mtiger7_ma"
fldr_ma <- "MA"
sub_ma <- "_ma" # with _

m = "02"
y = "2020"
# Load files -------------------------------------------------------------------
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  )

# Load the predictionin which you are interested:
pred <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                       "/tiger_", m ,"_", y, sub, ".rds")) %>%
  janitor::clean_names() %>%
  st_drop_geometry()

pred$pred_count <- rowMeans(pred[,6:ncol(pred)], na.rm = TRUE)

pred <- pred[c("municipality","id", "prov_name", "lon", "lat", "pred_count")]
# print(summary(pred$pred_count))

pred <- merge(pred, spain, by = c("municipality", "id", "prov_name"))
st_geometry(pred) <- "geometry"

a <- ggplot() +
  geom_sf(data = pred, aes(fill = pred_count), color = "transparent",
          size = 0.01, alpha = 0.8, na.rm = TRUE) +
  scale_fill_distiller("", palette = "Spectral", na.value = "transparent") + 
  ggtitle(paste(fldr, m, "-", y)) +
  theme_classic()

pred_ma <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr_ma, mdl_name_ma,
                       "/tiger_", m ,"_", y, sub_ma, ".rds")) %>%
  janitor::clean_names() %>%
  st_drop_geometry()

pred_ma$pred_count <- rowMeans(pred_ma[,6:ncol(pred_ma)], na.rm = TRUE)

pred_ma <- pred_ma[c("municipality","id", "prov_name", "lon", "lat", "pred_count")]
# print(summary(pred$pred_count))

pred_ma <- merge(pred_ma, spain, by = c("municipality", "id", "prov_name"))
st_geometry(pred_ma) <- "geometry"

b <- ggplot() +
  geom_sf(data = pred_ma, aes(fill = pred_count), color = "transparent",
          size = 0.01, alpha = 0.8, na.rm = TRUE) +
  scale_fill_distiller("", palette = "Spectral", na.value = "transparent") + 
  ggtitle(paste(fldr_ma, m, "-", y)) +
  theme_classic()

ggpubr::ggarrange(a, b)

# Spatial autocorrelation analysis ---------------------------------------------
sf::sf_use_s2(FALSE)

# Create a list of neighbors for each polygon
list_nb <- poly2nb(pred, queen = TRUE) 
list_nb_ma <- poly2nb(pred_ma, queen = TRUE) 
# he relationship will be based on a queen contiguity criterion, which considers 
# polygons to be neighbors if they share a boundary or a vertex.

# Check for empty neighbor sets
# card() calculates number of neighbors for each polygon in the list
# which() finds polygons with 0 neighbors
empty_nb <- which(card(list_nb) == 0)
empty_nb_ma <- which(card(list_nb_ma) == 0)

# Remove polygons with empty neighbor sets from the data
pred_subset <- pred[-empty_nb, ]# Remove polygons with empty neighbor sets from the data
pred_subset_ma <- pred_ma[-empty_nb_ma, ]

# # Subset 'tes_data' to extract polygons with empty neighbor sets
# empty_polygons <- pred[empty_nb, ]
# empty_polygons$municipality  # print municipality names

## Global G test ---------------------------------------------------------------
# Now that we removed empty neighbor sets (tes_subset)
# Identify neighbors with queen contiguity (edge/vertex touching)
tes_nb <- poly2nb(pred_subset, queen = TRUE)
tes_nb_ma <- poly2nb(pred_subset_ma, queen = TRUE)

# Binary weighting assigns a weight of 1 to all neighboring features 
# and a weight of 0 to all other features
tes_w_binary <- nb2listw(tes_nb, style="B", zero.policy=T )
tes_w_binary_ma <- nb2listw(tes_nb_ma, style="B", zero.policy=T )

# Calculate spatial lag of predicted values
tes_lag <- lag.listw(tes_w_binary, pred_subset$pred_count)
tes_lag_ma <- lag.listw(tes_w_binary_ma, pred_subset_ma$pred_count)

# Test for global G statistic of predicted values
globalG.test(pred_subset$pred_count, tes_w_binary)
globalG.test(pred_subset_ma$pred_count, tes_w_binary_ma)

## Local G test ----------------------------------------------------------------
# Identify neighbors, create weights, calculate spatial lag
tes_nbs <- pred_subset %>% 
  mutate(
    nb = st_contiguity(geometry),        # neighbors share border/vertex
    wt = st_weights(nb),                 # row-standardized weights
    tes_lag = st_lag(pred_count, nb, wt)    # calculate spatial lag of TreEqty
  ) 
tes_nbs_ma <- pred_subset_ma %>% 
  mutate(
    nb = st_contiguity(geometry),        
    wt = st_weights(nb),                 
    tes_lag = st_lag(pred_count, nb, wt)    
  ) 

# Calculate the Gi using local_g_perm
tes_hot_spots <- tes_nbs %>%
  mutate(
    Gi = local_g_perm(pred_count, nb, wt, nsim = 999)
  ) %>%
  unnest(Gi) 
tes_hot_spots_ma <- tes_nbs_ma %>%
  mutate(
    Gi = local_g_perm(pred_count, nb, wt, nsim = 999)
  ) %>%
  unnest(Gi) 

tes_hot_spots %>% 
  ggplot((aes(fill = gi))) +
  geom_sf(color = "black", lwd = 0.15) +
  scale_fill_gradient2()
tes_hot_spots_ma %>% 
  ggplot((aes(fill = gi))) +
  geom_sf(color = "black", lwd = 0.15) +
  scale_fill_gradient2()

# Create a new data frame called 'tes_hot_spots"
cluster <- tes_hot_spots %>%
  select(gi, p_folded_sim, id) %>% 
  mutate(
    classification = case_when(
      # Classify based on the following p-value criteria:
      # gi > 0 & p_folded_sim <= 0.01 ~ "Very hot",
      gi > 0 & p_folded_sim <= 0.05 ~ "Hot",
      # # gi > 0 & p_folded_sim <= 0.1 ~ "Somewhat hot",
      # # gi < 0 & p_folded_sim <= 0.01 ~ "Very cold",
      gi < 0 & p_folded_sim <= 0.05 ~ "Cold",
      # gi < 0 & p_folded_sim <= 0.01 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    classification = factor(
      classification,
      levels = c(
                 # "Very hot", 
                 "Hot", 
                 # "Somewhat hot",
                 "Insignificant",
                 # "Somewhat cold", 
                 "Cold" 
                 # "Very cold"
                 )
    )
  ) 

cluster_ma <- tes_hot_spots_ma %>%
  select(gi, p_folded_sim, id) %>% 
  mutate(
    classification = case_when(
      # Classify based on the following p-value criteria:
      # gi > 0 & p_folded_sim <= 0.01 ~ "Very hot",
      gi > 0 & p_folded_sim <= 0.05 ~ "Hot",
      # # gi > 0 & p_folded_sim <= 0.1 ~ "Somewhat hot",
      # # gi < 0 & p_folded_sim <= 0.01 ~ "Very cold",
      gi < 0 & p_folded_sim <= 0.05 ~ "Cold",
      # gi < 0 & p_folded_sim <= 0.01 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    classification = factor(
      classification,
      levels = c(
        # "Very hot", 
        "Hot", 
        # "Somewhat hot",
        "Insignificant",
        # "Somewhat cold", 
        "Cold"
        # "Very cold"
        )
    )
  )

ggpubr::ggarrange(
  ggplot(cluster, aes(fill = classification)) +
    geom_sf(color = "black", lwd = 0.1) +
    scale_fill_brewer(type = "div", palette = 5) +
    theme_void() +
    labs(
      fill = "Hot Spot Classification",
      title = paste("Predicted Hot Spots of", fldr)
    ),
  ggplot(cluster_ma, aes(fill = classification)) +
    geom_sf(color = "black", lwd = 0.1) +
    scale_fill_brewer(type = "div", palette = 5) +
    theme_void() +
    labs(
      fill = "Hot Spot Classification",
      title = paste("Predicted Hot Spots of", fldr_ma)
    )
)
  
diff_cluster <- cluster %>%
  mutate(
    classification_ma = st_drop_geometry(cluster_ma$classification),
    diff = case_when(
      classification == "Insignificant" & classification_ma =="Insignificant" ~ "Overlap",
      classification == "Hot" & classification_ma =="Hot" ~ "Overlap",
      classification == "Cold" & classification_ma =="Cold" ~ "Overlap",
      classification == "Hot" & classification_ma == "Cold" ~ paste(fldr, ">", fldr_ma),
      classification == "Hot" & classification_ma == "Insignificant" ~ paste(fldr, ">", fldr_ma),
      classification == "Cold" & classification_ma == "Hot" ~ paste(fldr, "<", fldr_ma),
      classification == "Cold" & classification_ma == "Insignificant" ~ paste(fldr, "<", fldr_ma),
      classification == "Insignificant" & classification_ma == "Hot" ~ paste(fldr, "<", fldr_ma),
      classification == "Insignificant" & classification_ma == "Cold" ~ paste(fldr, ">", fldr_ma),
    ),
    diff = factor(
      diff,
      levels = c(paste(fldr, ">", fldr_ma), "Overlap",  paste(fldr, "<", fldr_ma)),
    )
  ) 

custom_colors <- setNames(
  c("#abdda4", "#f6f6f6", "#d53e4f"),  
  c(paste(fldr, ">", fldr_ma),        
    "Overlap",                         
    paste(fldr, "<", fldr_ma))         
)
ggplot(diff_cluster, aes(fill = diff)) +
  geom_sf(color = "black", lwd = 0.1) +
  scale_fill_manual(values = custom_colors) +
  theme_void() +
  labs(
    fill = " Classification",
    title = paste("Difference on Spatial Cluster of", fldr, "and", fldr_ma)
  )

# Analysing the climatic variables----------------------------------------------
spain_wth <- mclapply(1:nrow(spain), function(i){
  
  cat(paste0("Number row:", i, "\n"))
  
  data_row <- spain[i, ]
  
  data_point <- st_centroid(data_row)
  
  wth <- readRDS(paste0(loc.output, "monthly_weather_data/prep_", 
                        m, "-", y, ".rds"))[, c(2,9:12)]
  
  data_row <- merge(data_row, wth, by = "id", all.x = TRUE)
}
, mc.cores = 8)
spain_wth <- do.call(rbind, spain_wth)

diff_cluster_wth <- merge(diff_cluster %>% st_drop_geometry(), 
           spain_wth %>% st_drop_geometry(),
           by = "id")

ggpubr::ggarrange(
  ggplot(diff_cluster_wth) +
    geom_boxplot(aes(x = diff, y = min_temperature)) +
    theme_classic(),
  ggplot(diff_cluster_wth) +
    geom_boxplot(aes(x = diff, y = mean_temperature)) +
    theme_classic(),
  ggplot(diff_cluster_wth) +
    geom_boxplot(aes(x = diff, y = precipitation)) +
    theme_classic(),
  ggplot(diff_cluster_wth) +
    geom_boxplot(aes(x = diff, y = mean_relative_humidity)) +
    theme_classic()
)


