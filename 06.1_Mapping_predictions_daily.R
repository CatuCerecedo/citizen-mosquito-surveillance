########################### Mapping the predictions ############################

library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(tidyverse)

rm(list = ls())

# Directories ------------------------------------------------------------------
# In local
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")
loc.era5 <- "/home/catuxa/Documents/Mosquito_Models/EU_tiger/ERA5_Download/"

# # In cluster
# loc.output <- paste0(getwd(), "/Spain_Tiger/OUTPUT/")
# loc.data <- paste0(getwd(), "/Spain_Tiger/DATA/")
# loc.figures <- paste0(getwd(), "/Spain_Tiger/FIGURES/")
# loc.era5 <- paste0(getwd(), "/EU_tiger/ERA5_Download/")

sf::sf_use_s2(FALSE)

mdl <- "mtiger16"
mdl_name <- "/mtiger16"
fldr <- "Counts"
sub <- "" # with _

# mdl <- "mtiger3_occu"
# mdl_name <- "/mtiger3_occu"
# fldr <- "Suitability"
# sub <- "_occu" # with _

# mdl <- "mtiger7_ma"
# mdl_name <- "/mtiger7_ma"
# fldr <- "MA"
# sub <- "_ma" # with _

# mdl <- "mbites_7_reduce"
# mdl_name <- "/mbites_7_reduce"
# fldr <- "Bites"
# sub <- "bts" # with _

fldr <- "Monthly_Integration"
mdl_name <- "/tiger_inte"
mdl <- "tiger_inte"
sub <- "" # with _

# Some checks ------------------------------------------------------------------
# Mapping climatic variables
wthr <- readRDS(paste0(loc.output, "daily_weather_data/prep_2020-07-20.rds"))
# st_geometry(wthr) <- "geometry"
wthr <- wthr %>%
  janitor::clean_names() %>%
  st_as_sf(coords = c( "lon", "lat"), crs = 4326, remove = FALSE)

ggplot() + 
  geom_sf(data = wthr, aes(color = l21mean_temperature), size = 0.2) +
  scale_color_distiller("", palette = "Spectral") +
  theme_classic()

# Loading spain municipality map -----------------------------------------------
# Country
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  )

# Plot preidctions: raster/tif -------------------------------------------------

years <- c("2020", "2021", "2022")
months <- c("01", "02", "03", "04", "05", "06" ,"07", "08", "09", "10", "11", "12")

tiger_maps <- list()
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
    # print(summary(pred$pred_count))
    
    pred <- merge(pred, spain, by = c("municipality", "id", "prov_name"))
    st_geometry(pred) <- "geometry"
    
    plt <- ggplot() +
      geom_sf(data = pred, aes(fill = pred_count), color = "transparent",
              size = 0.01, alpha = 0.8, na.rm = TRUE) +
      scale_fill_distiller("", palette = "Spectral", na.value = "transparent") + 
      # xlim(-13, 5) +
      # ylim(34, 44) +
      ggtitle(paste0(m, "-", y)) +
      theme_classic()
    
    tiger_maps[[iter]] <- plt
    month_values[[iter]] <- pred$pred_count %>% st_drop_geometry()
  }
}

a1 <- tiger_maps[[1]]
a2 <- tiger_maps[[2]]
a3 <- tiger_maps[[3]]
a4 <- tiger_maps[[4]]
a5 <- tiger_maps[[5]]
a6 <- tiger_maps[[6]]
a7 <- tiger_maps[[7]]
a8 <- tiger_maps[[8]]
a9 <- tiger_maps[[9]]
a10 <- tiger_maps[[10]]
a11 <- tiger_maps[[11]]
a12 <- tiger_maps[[12]]

ggpubr::ggarrange(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12,
                  common.legend = TRUE, nrow = 3, ncol = 4, legend = "right")
ggsave(paste0(loc.fig, "monthly_tiger_2020_", mdl, ".png"), units = "cm", bg = "white",
       height = 25, width = 25)

a1 <- tiger_maps[[13]]
a2 <- tiger_maps[[14]]
a3 <- tiger_maps[[15]]
a4 <- tiger_maps[[16]]
a5 <- tiger_maps[[17]]
a6 <- tiger_maps[[18]]
a7 <- tiger_maps[[19]]
a8 <- tiger_maps[[20]]
a9 <- tiger_maps[[21]]
a10 <- tiger_maps[[22]]
a11 <- tiger_maps[[23]]
a12 <- tiger_maps[[24]]

ggpubr::ggarrange(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12,
                  common.legend = TRUE, nrow = 3, ncol = 4, legend = "right")
ggsave(paste0(loc.fig, "monthly_tiger_2021_", mdl, ".png"), units = "cm", bg = "white",
       height = 25, width = 25)

a1 <- tiger_maps[[25]]
a2 <- tiger_maps[[26]]
a3 <- tiger_maps[[27]]
a4 <- tiger_maps[[28]]
a5 <- tiger_maps[[29]]
a6 <- tiger_maps[[30]]
a7 <- tiger_maps[[31]]
a8 <- tiger_maps[[32]]
a9 <- tiger_maps[[33]]
a10 <- tiger_maps[[34]]
a11 <- tiger_maps[[35]]
a12 <- tiger_maps[[36]]

ggpubr::ggarrange(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12,
                  common.legend = TRUE, nrow = 3, ncol = 4, legend = "right")
ggsave(paste0(loc.fig, "monthly_tiger_2022_", mdl, ".png"), units = "cm", bg = "white",
       height = 25, width = 25)

rm(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12)

# Predicted plot ---------------------------------------------------------------
bg_pred <- ggpubr::ggarrange(a4, a11, common.legend = TRUE, nrow = 1, ncol = 2, legend = "right")
suit_pred <- ggpubr::ggarrange(a4, a11, common.legend = TRUE, nrow = 1, ncol = 2, legend = "right")
ma_pred <- ggpubr::ggarrange(a4, a11, common.legend = TRUE, nrow = 1, ncol = 2, legend = "right")
bite_pred <- ggpubr::ggarrange(a4, a11, common.legend = TRUE, nrow = 1, ncol = 2, legend = "right")

ggpubr::ggarrange(suit_pred, bg_pred, ma_pred, bite_pred,
                  nrow = 4, ncol = 1)
ggsave(file = paste0(loc.fig, "predicted_plot_Spain.png"), 
       units = "cm", height = 25, width = 25, bg = "white")

# Plotting spatial correlations ------------------------------------------------
years <- c("2020", "2021", "2022")
months <- c("01", "02", "03", "04", "05", "06" ,"07", "08", "09", "10", "11", "12")

spatial_corr <- vector(mode = "list")
for(y in years){
  pred_bg_all <- readRDS(paste0(loc.output, "PREDICTIONS/Counts/mtiger16/tiger_", "07" ,"_", y, ".rds")) %>%
    janitor::clean_names() %>%
    dplyr::select(municipality, id, lon, lat, prov_name)
  # Suitability measure
  # pred_bg_all <- readRDS(paste0(loc.output, "PREDICTIONS/Suitability/mtiger3_occu/tiger_", "07" ,"_", y, "_occu.rds")) %>%
  #   janitor::clean_names() %>%
  #   dplyr::select(municipality, id, lon, lat, prov_name)
  # Bites
  # pred_bg_all <- readRDS(paste0(loc.output, "PREDICTIONS/Bites/mbites_7_reduce/tiger_", "07" ,"_", y, "bts.rds")) %>%
  #   janitor::clean_names() %>%
  #   dplyr::select(municipality, id, lon, lat, prov_name)
  
  pred_ma_all <- readRDS(paste0(loc.output, "PREDICTIONS/MA/mtiger7_ma/tiger_", "07" ,"_", y, "_ma.rds")) %>%
    janitor::clean_names() %>%
    dplyr::select(municipality, id, lon, lat, prov_name)
  
  for (m in months){
    print(paste0("Plotting: ", m, "-", y))
    
    pred_bg <- readRDS(paste0(loc.output, "PREDICTIONS/Counts/mtiger16/tiger_", m ,"_", y, ".rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    # Si quiero plotear la correlación con suitability
    # pred_bg <- readRDS(paste0(loc.output, "PREDICTIONS/Suitability/mtiger3_occu/tiger_", m ,"_", y, "_occu.rds")) %>%
    #   janitor::clean_names() %>%
    #   st_drop_geometry()
    # Bites
    # pred_bg <- readRDS(paste0(loc.output, "PREDICTIONS/Bites/mbites_7_reduce/tiger_", m ,"_", y, "bts.rds")) %>%
    #   janitor::clean_names() %>%
    #   st_drop_geometry()
    pred_bg$pred_count <- rowMeans(pred_bg[,6:ncol(pred_bg)], na.rm = TRUE)
    pred_bg <- pred_bg %>% dplyr::select(municipality, id, lon, lat, prov_name, pred_count)
    pred_bg_all <- merge(pred_bg_all, pred_bg, by = c("municipality", "id", "lon", "lat", "prov_name"))
    
    pred_ma <- readRDS(paste0(loc.output, "PREDICTIONS/MA/mtiger7_ma/tiger_", m ,"_", y, "_ma.rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    pred_ma$pred_count <- rowMeans(pred_ma[,6:ncol(pred_ma)], na.rm = TRUE)
    pred_ma <- pred_ma %>% dplyr::select(municipality, id, lon, lat, prov_name, pred_count)
    
    pred_ma_all <- merge(pred_ma_all, pred_ma, by = c("municipality", "id", "lon", "lat", "prov_name"))
    
  }
  colnames(pred_bg_all) <- c("municipality", "id", "lon", "lat", "prov_name", "x01", "x02", "x03", "x04", "x05", "x06" ,"x07", "x08", "x09", "x10", "x11", "x12")
  colnames(pred_ma_all) <- c("municipality", "id", "lon", "lat", "prov_name", "x01", "x02", "x03", "x04", "x05", "x06" ,"x07", "x08", "x09", "x10", "x11", "x12")
  
  cor_df <- data.frame()
  for (i in 1:nrow(pred_bg_all)){
    if(i %% 1000 == 0){
      print(i)
    }
    cor_row <- data.frame(
      municipality = pred_bg_all[i, 1],
      id = pred_bg_all[i, 2],
      lon = pred_bg_all[i, 3],
      lat = pred_bg_all[i, 4],
      prov_name = pred_bg_all[i, 5],
      rho = cor(as.numeric(pred_bg_all[i, 6:17]), as.numeric(pred_ma_all[i, 6:17]), method = "spearman"),
      year = y
    )
    cor_df <- rbind(cor_df, cor_row)
  }
  
  spatial_corr[[y]] <- cor_df
}

# saveRDS(spatial_corr, file = paste0(loc.output, "spatial_corr_suit_ma.rds"))
# saveRDS(spatial_corr, file = paste0(loc.output, "spatial_corr_bg_ma.rds"))

df <- merge(spatial_corr[["2020"]], spain, by = c("municipality", "id"))
st_geometry(df) <- "geometry"
a <- ggplot() +
  geom_sf(data = df, aes(fill = rho), color = "transparent", size = 1) +
  scale_fill_distiller("", palette = "Spectral", limits = c(0.7, 1)) + 
  ggtitle("2020") +
  theme_classic() +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
print(mean(df$rho, na.rm = TRUE))
print(sd(df$rho, na.rm = TRUE))

df <- merge(spatial_corr[["2021"]], spain, by = c("municipality", "id"))
st_geometry(df) <- "geometry"
b <-  ggplot() +
  geom_sf(data = df, aes(fill = rho), color = "transparent", size = 1) +
  scale_fill_distiller("", palette = "Spectral", limits = c(0.7, 1)) + 
  ggtitle("2021") +
  theme_classic() +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
print(mean(df$rho, na.rm = TRUE))
print(sd(df$rho, na.rm = TRUE))

df <- merge(spatial_corr[["2022"]], spain, by = c("municipality", "id"))
st_geometry(df) <- "geometry"
c <- ggplot() +
  geom_sf(data = df, aes(fill = rho), color = "transparent", size = 1) +
  scale_fill_distiller("", palette = "Spectral", limits = c(0.7, 1)) + 
  ggtitle("2022") +
  theme_classic() +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
print(mean(df$rho, na.rm = TRUE))
print(sd(df$rho, na.rm = TRUE))

p <- ggpubr::ggarrange(a, b, c,
                  common.legend = TRUE, nrow = 1, ncol = 3, legend = "right")
bg_ma <- ggdraw(p) +
  theme(
    plot.margin = margin(t = 3, r = 10, b = 3, l = 10), 
    plot.background = element_rect(color = "black", linewidth = 2, fill = "white") 
  )

ggsave(file = paste0(loc.fig, "spatial_correlation_bg_ma_Spain.png"), 
       units = "cm", height = 10, width = 20, bg = "white")

# Plotting temporal correlations -----------------------------------------------
years <- c("2020", "2021", "2022")
months <- c("01", "02", "03", "04", "05", "06" ,"07", "08", "09", "10", "11", "12")

cor_predictions <- data.frame()
for(y in years){
  for (m in months){
    print(paste0("Plotting: ", m, "-", y))
    
    pred_bg <- readRDS(paste0(loc.output, "PREDICTIONS/Counts/mtiger16/tiger_", m ,"_", y, ".rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    pred_bg$bg <- rowMeans(pred_bg[,6:ncol(pred_bg)], na.rm = TRUE)
    pred_bg <- pred_bg %>% dplyr::select(municipality, id, lon, lat, prov_name, bg)
    # pred_bg <- pred_bg %>% filter((pixel_id %in% tiger$pixel_id))
    pred_suit <- readRDS(paste0(loc.output, "PREDICTIONS/Suitability/mtiger3_occu/tiger_", m ,"_", y, "_occu.rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    pred_suit$suit <- rowMeans(pred_suit[,6:ncol(pred_suit)], na.rm = TRUE)
    pred_suit <- pred_suit %>% dplyr::select(municipality, id, lon, lat, prov_name, suit)
    # pred_suit <- pred_suit %>% filter((pixel_id %in% tiger$pixel_id))
    pred_ma <- readRDS(paste0(loc.output, "PREDICTIONS/MA/mtiger7_ma/tiger_", m ,"_", y, "_ma.rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    pred_ma$ma <- rowMeans(pred_ma[,6:ncol(pred_ma)], na.rm = TRUE)
    pred_ma <- pred_ma %>% dplyr::select(municipality, id, lon, lat, prov_name, ma)
    # pred_ma <- pred_ma %>% filter((pixel_id %in% tiger$pixel_id))
    pred_bites <- readRDS(paste0(loc.output, "PREDICTIONS/Bites/mbites_7_reduce/tiger_", m ,"_", y, "bts.rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    pred_bites$bites <- rowMeans(pred_bites[,6:ncol(pred_bites)], na.rm = TRUE)
    pred_bites <- pred_bites %>% dplyr::select(municipality, id, lon, lat, prov_name, bites)
    # pred_bites <- pred_bites %>% filter((pixel_id %in% tiger$pixel_id))
    
    # Joining tables
    pred <- merge(pred_bg, pred_suit, by = c("municipality", "id", "lon", "lat", "prov_name"))
    pred <- merge(pred, pred_ma, by = c("municipality", "id", "lon", "lat", "prov_name"))
    pred <- merge(pred, pred_bites, by = c("municipality", "id", "lon", "lat", "prov_name"))
    
    rm(pred_bg, pred_suit, pred_ma, pred_bites)
    
    cor_row <- data.frame(
      year = y,
      month = as.numeric(m), 
      count_suit = as.numeric(cor.test(pred$bg, pred$suit, method="spearman")$estimate),
      count_suit_p = cor.test(pred$bg, pred$suit, method="spearman")$p.value,
      count_ma = as.numeric(cor.test(pred$bg, pred$ma, method="spearman")$estimate),
      count_ma_p = cor.test(pred$bg, pred$ma, method="spearman")$p.value,
      count_bites = as.numeric(cor.test(pred$bg, pred$bites, method="spearman")$estimate),
      count_bites_p = cor.test(pred$bg, pred$bites, method="spearman")$p.value,
      suit_ma = as.numeric(cor.test(pred$suit, pred$ma, method="spearman")$estimate),
      suit_ma_p = cor.test(pred$suit, pred$ma, method="spearman")$p.value,
      suit_bites = as.numeric(cor.test(pred$suit, pred$bites, method="spearman")$estimate),
      suit_bites_p = cor.test(pred$suit, pred$bites, method="spearman")$p.value,
      bites_ma = as.numeric(cor.test(pred$ma, pred$bites, method="spearman")$estimate),
      bites_ma_p = cor.test(pred$ma, pred$bites, method="spearman")$p.value
    )
    cor_predictions <- rbind(cor_predictions, cor_row)
  }
}

cor_predictions <- tidyr::pivot_longer(
  cor_predictions,
  cols = c(count_suit, count_ma, count_bites, suit_ma, suit_bites, bites_ma),
  names_to = "comparison",
  values_to = "value"
)

saveRDS(cor_predictions, file = paste0(loc.output, "temporal_corr.rds"))

# Define una paleta de colores personalizada
custom_colors <- c("suit_ma" = "#abdda4", "count_ma" = "#d53e4f", "bites_ma" = "#3288bd")

# # the mean annual values
# mean_values <- cor_predictions %>%
#   group_by(year) %>%
#   summarise(mean_value = mean(value, na.rm = TRUE))

ggplot(cor_predictions, aes(x = month, y = value, color = comparison, shape = comparison, linetype = comparison)) +
  geom_point(size = 7, alpha = 0.6) +
  geom_line(size = 2, alpha = 0.6) +
  # geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c("suit_ma" = 15, "count_ma" = 16, "bites_ma" = 17)) +
  scale_linetype_manual(values = c("suit_ma" = "dashed", "count_ma" = "solid", "bites_ma" = "dotdash")) +
  labs(
    y = "Spearman Correlation (S)",
    x = "Month",
    color = "Comparison",      # Aquí se especifica "Comparison" para la leyenda unificada
    shape = "Comparison",      # Se especifica también para shape
    linetype = "Comparison",   # Y para linetype
    # title = "Spanish Tiger Mosquito Models: Spatial correlation over time"
  ) +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  scale_y_continuous(breaks = seq(-0.5, 0.9, 0.2)) +
  theme_classic() +
  facet_wrap(~year) +
  theme(
    legend.position = "bottom",    # La leyenda se posiciona abajo
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
  ) +
  guides(                       # Unificación de la leyenda
    color = guide_legend(title = "Comparison"),
    shape = guide_legend(title = "Comparison"),
    linetype = guide_legend(title = "Comparison")
  )

ggsave(file = paste0(loc.fig, "temporal_correlation_Spain.png"), 
       units = "cm", height = 20, width = 35, bg = "white")

# Plotting Spearman correlation against abundances
tiger <- readRDS(file = paste0(loc.output, "tiger_integrating_daily.rds")) %>%
  drop_na(ma) %>%
  mutate(year = as.factor(year(end_date)),
         month = month(end_date)) %>%
  group_by(year, month) %>%
  summarize(females = mean(females))

df <- merge(cor_predictions, tiger, by = c("year", "month"))
ggplot(df, aes(x = females, y = value, group = comparison, color = comparison)) +
  geom_point(size = 7, alpha = 0.6) +
  geom_smooth(method = lm, aes(fill = comparison)) +
  scale_color_manual(values = custom_colors) +
  labs(
    x = "\nAverage Counts of Tiger Mosquitoes from Traps",
    y = "Spearman Spatial Correlation (S)\n",
    color = "Comparison",
    fill = "Comparison"
  ) +
  # facet_wrap(~year) +
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

ggsave(file = paste0(loc.fig, "spearman_vs_obs_counts_Spain.png"), 
       units = "cm", height = 20, width = 35, bg = "white")

# Plotting the same but with reports
reports <- readRDS("/home/catuxa/Documents/Mosquito_Models/EU_Culex/OUTPUT/mosquito_alert_cleaned_reports_cs_responses.Rds") %>% 
  filter(!is.na(lon) & !is.na(lat)) %>% 
  st_as_sf(coords = c("lon", "lat"), crs=4326, remove=FALSE) %>% 
  mutate(package_name_version = paste0(package_name, "_", package_version)) %>%
  filter(movelab_certainty_category_euro_class_label == "aedes-albopictus", movelab_certainty_category_euro_class_value >= 1) 

reports <- reports %>%
  mutate(year = as.factor(year(date)),
         month = month(date)) %>%
  group_by(year, month) %>%
  summarize(n_reports = n())

df <- merge(cor_predictions, reports, by = c("year", "month"))
ggplot(df, aes(x = n_reports, y = value, group = comparison, color = comparison)) +
  geom_point() +
  geom_smooth(method = lm, aes(fill = comparison), alpha = 0.1) +
  scale_color_manual(values = custom_colors) +
  labs(
    x = "Number of validated Reports of Mosquito Alert",
    y = "Spearman Spatial Correlation",
    color = "Comparison",
    fill = "Comparison"
  ) +
  # theme_classic()
theme_classic() +
facet_wrap(~year)

# Plotting uncertainty against abundances 
years <- c("2020", "2021", "2022")
months <- c("01", "02", "03", "04", "05", "06" ,"07", "08", "09", "10", "11", "12")

uncertainty <- data.frame()
for(y in years){
  for (m in months){
    print(paste0("Plotting: ", m, "-", y))
    
    pred_bg <- readRDS(paste0(loc.output, "PREDICTIONS/Counts/mtiger16/tiger_", m ,"_", y, "_sd.rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    pred_bg$bg <- rowMeans(pred_bg[,6:ncol(pred_bg)], na.rm = TRUE)
    pred_bg <- pred_bg %>% dplyr::select(municipality, id, lon, lat, prov_name, bg)
    # pred_bg <- pred_bg %>% filter((pixel_id %in% tiger$pixel_id))
    pred_suit <- readRDS(paste0(loc.output, "PREDICTIONS/Suitability/mtiger3_occu/tiger_", m ,"_", y, "_occu_sd.rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    pred_suit$suit <- rowMeans(pred_suit[,6:ncol(pred_suit)], na.rm = TRUE)
    pred_suit <- pred_suit %>% dplyr::select(municipality, id, lon, lat, prov_name, suit)
    # pred_suit <- pred_suit %>% filter((pixel_id %in% tiger$pixel_id))
    pred_ma <- readRDS(paste0(loc.output, "PREDICTIONS/MA/mtiger7_ma/tiger_", m ,"_", y, "_ma_sd.rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    pred_ma$ma <- rowMeans(pred_ma[,6:ncol(pred_ma)], na.rm = TRUE)
    pred_ma <- pred_ma %>% dplyr::select(municipality, id, lon, lat, prov_name, ma)
    # pred_ma <- pred_ma %>% filter((pixel_id %in% tiger$pixel_id))
    
    # Joining tables
    pred <- merge(pred_bg, pred_suit, by = c("municipality", "id", "lon", "lat", "prov_name"))
    pred <- merge(pred, pred_ma, by = c("municipality", "id", "lon", "lat", "prov_name"))
    
    rm(pred_bg, pred_suit, pred_ma)
    
    uncer <- data.frame(
      year = y,
      month = as.numeric(m), 
      suit_sd = mean(pred$suit, na.rm = TRUE),
      count_sd = mean(pred$bg, na.rm = TRUE),
      ma_sd = mean(pred$ma, na.rm = TRUE)
    )
    uncertainty <- rbind(uncertainty, uncer)
  }
}

uncertainty <- tidyr::pivot_longer(
  uncertainty,
  cols = c(suit_sd, count_sd, ma_sd),
  names_to = "Uncertainty",
  values_to = "value"
)

# Loading tiger database
tiger <- readRDS(file = paste0(loc.output, "tiger_integrating_daily.rds")) %>%
  drop_na(ma) %>%
  mutate(year = as.factor(year(end_date)),
         month = month(end_date)) %>%
  group_by(year, month) %>%
  summarize(females = mean(females))

df <- merge(uncertainty, tiger, by = c("year", "month"))
custom_colors <- c("count_sd" = "#abdda4", "ma_sd" = "#d53e4f", "suit_sd" = "#3288bd")
ggplot(df, aes(x = females, y = value, color = Uncertainty)) +
  geom_smooth(size = 2, alpha = 0.6) +
  geom_point(size = 7, alpha = 0.6) +
  scale_color_manual(values = custom_colors) +
  labs(
    y = "Uncertainty (sd)\n",
    x = "\nAverage Counts of Tiger Mosquitoes from Traps",
    color = "Model",      # Aquí se especifica "Comparison" para la leyenda unificada
    shape = "Model",      # Se especifica también para shape
    linetype = "Model",   # Y para linetype
  ) +
  theme_classic() +
  facet_wrap(~Uncertainty, scales = "free", labeller = as_labeller(c(count_sd = "Count Model Uncertainty",
                                                                     suit_sd = "Suitability Model Uncertainty",
                                                                     ma_sd = "MA Model Uncertainty"))) +
  theme(
    legend.position = "Right",
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
  ) +
  guides(                       # Unificación de la leyenda
    color = guide_legend(title = "Model Uncertainty")
  )

ggsave(file = paste0(loc.fig, "uncertainty_vs_obs_counts_Spain.png"), 
       units = "cm", height = 20, width = 35, bg = "white")

# Plotting uncertainty against abundances 
reports <- readRDS("/home/catuxa/Documents/Mosquito_Models/EU_Culex/OUTPUT/mosquito_alert_cleaned_reports_cs_responses.Rds") %>% 
  filter(!is.na(lon) & !is.na(lat)) %>% 
  st_as_sf(coords = c("lon", "lat"), crs=4326, remove=FALSE) %>% 
  mutate(package_name_version = paste0(package_name, "_", package_version)) %>%
  filter(movelab_certainty_category_euro_class_label == "aedes-albopictus", movelab_certainty_category_euro_class_value >= 1) 

reports <- reports %>%
  mutate(year = as.factor(year(date)),
         month = month(date)) %>%
  group_by(year, month) %>%
  summarize(n_reports = n())

df <- merge(uncertainty, reports, by = c("year", "month"))
custom_colors <- c("count_sd" = "#abdda4", "ma_sd" = "#d53e4f", "suit_sd" = "#3288bd")
ggplot(df, aes(x = n_reports, y = value, color = Uncertainty)) +
  geom_smooth(size = 2, alpha = 0.6) +
  geom_point(size = 7, alpha = 0.6) +
  scale_color_manual(values = custom_colors) +
  labs(
    y = "Uncertainty (sd)\n",
    x = "\nNumber of validated Reports of Mosquito Alert"
  ) +
  theme_classic() +
  facet_wrap(~Uncertainty, scales = "free", labeller = as_labeller(c(count_sd = "Count Model Uncertainty",
                                                                     suit_sd = "Suitability Model Uncertainty",
                                                                     ma_sd = "MA Model Uncertainty"))) +
  theme(
    legend.position = "Right",
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
  ) +
  guides(                       # Unificación de la leyenda
    color = guide_legend(title = "Model Uncertainty")
  )

ggsave(file = paste0(loc.fig, "uncertainty_vs_n_reports_Spain.png"), 
       units = "cm", height = 20, width = 35, bg = "white")
# Plotting the reports of tiger ------------------------------------------------
reports <- readRDS(paste0(loc.output, "mosquito_alert_cleaned_reports_cs_responses.Rds")) %>% 
  filter(!is.na(lon) & !is.na(lat)) %>% 
  st_as_sf(coords = c("lon", "lat"), crs=4326, remove=FALSE) %>% 
  mutate(package_name_version = paste0(package_name, "_", package_version)) %>%
  filter(movelab_certainty_category_euro_class_label == "tiger-sp", movelab_certainty_category_euro_class_value >= 1) 
tiger <- readRDS(paste0(loc.output, "tiger_clc_clima_ma_monthly.rds")) %>%
  dplyr::select(pixel_id, longitude, latitude) %>%
  unique()

library(leaflet)
leaflet() %>%
  addTiles() %>%
  addCircles(lng = st_coordinates(reports)[,1], lat = st_coordinates(reports)[,2])

# Sampling size for months -----------------------------------------------------
reports <- readRDS(paste0(loc.output, "mosquito_alert_cleaned_reports_cs_responses.Rds")) %>% 
  filter(!is.na(lon) & !is.na(lat)) %>% 
  mutate(package_name_version = paste0(package_name, "_", package_version)) %>%
  filter(movelab_certainty_category_euro_class_label == "aedes-albopictus", movelab_certainty_category_euro_class_value >= 1) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  

index_intersects <- st_intersects(reports, spain) 
index_intersects <- lengths(index_intersects) > 0

reports <- reports[index_intersects, ] # Only selecting Spain data
reports <- reports %>% 
  mutate(
    year = as.factor(lubridate::year(date)),
    month = as.numeric(lubridate::month(date)),
    type = "reports"
  ) %>% 
  filter(year %in% c("2020", "2021", "2022")) %>%
  group_by(year, month, type) %>%
  summarise(n = dplyr::n()) %>%
  mutate(n_prop = n/434)

tiger <- readRDS(paste0(loc.output, "bg_tiger_spain_daily.rds")) %>%
  mutate(
    year = as.factor(lubridate::year(end_date)),
    month = as.numeric(lubridate::month(end_date)),
    type = "traps"
  ) %>% 
  group_by(year, month, type) %>%
  summarise(n = sum(females)) 

sampling_size <- rbind(tiger, reports)

custom_colors <- c("traps" =  "#004949", "reports" = "#924900")
ggplot(sampling_size %>% 
         filter(year != "2023"), 
       aes(x = factor(month), y = n, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(breaks=seq(1, 12, 1)) +
  facet_wrap(~ year, scales = "free_x") +
  labs(
    x = "Month",
    y = "Absolute Number of samples",
    fill = "Type",
    title = "Monthly absolute Sample Size "
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )
ggsave(file = paste0(loc.fig, "Absolute_sample_size_Spain.png"), 
       units = "cm", height = 10, width = 20, bg = "white")

# Comparing with RM  -----------------------------------------------------------
# WARNING: Marta only has RM estimates for 2020 year (no temporal matching)

# Comparing with Bites  --------------------------------------------------------
bites <- readxl::read_excel(paste0(loc.output, "bites/observations.xlsx")) %>%
  dplyr::select(date, longitude, latitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(crs(template))

month_values <- do.call(cbind, month_values)
pred <- cbind(pred[,1:4], month_values)
pred <- pred %>% st_transform(crs(template))

n_columns <- c("pixel_id", "longitude", "latitude", "country")
for (y in years){
  for (m in months){
    n_columns <- append(n_columns, paste0("X", m, "_", y))
  }
}
colnames(pred) <- c(n_columns, "geometry")

a14 <- rasterize(vect(pred), template, field = "X07_2022", fun = mean)
a15 <- rasterize(vect(pred), template, field = "X08_2022", fun = mean)
a16 <- rasterize(vect(pred), template, field = "X09_2022", fun = mean)
a17 <- rasterize(vect(pred), template, field = "X10_2022", fun = mean)

bites$pred_bg_17 <- extract(a17, bites)[,2]
hist(bites$pred_bg_15)
