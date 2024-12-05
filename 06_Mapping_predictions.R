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
# Some checks ------------------------------------------------------------------
# Mapping climatic variables
wthr <- readRDS(paste0(loc.output, "monthly_weather_data/prep_08-2020.rds"))
# st_geometry(wthr) <- "geometry"
wthr <- wthr %>%
  janitor::clean_names() %>%
  st_as_sf(coords = c( "lon", "lat"), crs = 4326, remove = FALSE)

ggplot() + 
  geom_sf(data = wthr, aes(color = precipitation), size = 0.2) +
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
    
    colnames(pred) <- c("municipality","id", "prov_name", "lon", "lat", "pred_count") 
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
bg_pred <- ggpubr::ggarrange(a7, a8, common.legend = TRUE, nrow = 1, ncol = 2, legend = "right")
suit_pred <- ggpubr::ggarrange(a7, a8, common.legend = TRUE, nrow = 1, ncol = 2, legend = "right")
ma_pred <- ggpubr::ggarrange(a7, a8, common.legend = TRUE, nrow = 1, ncol = 2, legend = "right")

ggpubr::ggarrange(suit_pred, bg_pred, ma_pred,
                  nrow = 3, ncol = 1)
ggsave(file = paste0(loc.fig, "predicted_plot.png"), 
       units = "cm", height = 25, width = 25, bg = "white")

# Plotting spatial correlations ------------------------------------------------
years <- c("2020", "2021", "2022")
months <- c("01", "02", "03", "04", "05", "06" ,"07", "08", "09", "10", "11", "12")

spatial_corr <- vector(mode = "list")
for(y in years){
  pred_bg_all <- readRDS(paste0(loc.output, "PREDICTIONS/Counts/mtiger10/tiger_", "07" ,"_", y, ".rds")) %>%
    janitor::clean_names() %>%
    dplyr::select(municipality, id, lon, lat, prov_name)
  pred_ma_all <- readRDS(paste0(loc.output, "PREDICTIONS/MA/mtiger7_ma/tiger_", "07" ,"_", y, "_ma.rds")) %>%
    janitor::clean_names() %>%
    dplyr::select(municipality, id, lon, lat, prov_name)
  
  for (m in months){
    print(paste0("Plotting: ", m, "-", y))
    
    # Si quiero ploter la correlación con suitability
    pred_bg <- readRDS(paste0(loc.output, "PREDICTIONS/Counts/mtiger10/tiger_", m ,"_", y, ".rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    # pred_bg <- readRDS(paste0(loc.output, "PREDICTIONS/BG_abundance/Monthly/tiger_", m ,"_", y, ".rds")) %>%
    #   janitor::clean_names() %>%
    #   st_drop_geometry()
    pred_bg_all <- merge(pred_bg_all, pred_bg, by = c("municipality", "id", "lon", "lat", "prov_name"))
    
    pred_ma <- readRDS(paste0(loc.output, "PREDICTIONS/MA/mtiger7_ma/tiger_", m ,"_", y, "_ma.rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
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

df <- merge(spatial_corr[["2020"]], spain, by = c("municipality", "id"))
st_geometry(df) <- "geometry"
a <- ggplot() +
  geom_sf(data = df, aes(fill = rho), color = "transparent", size = 1) +
  scale_fill_distiller("", palette = "Spectral") + 
  ggtitle("2020") +
  theme_classic()

df <- merge(spatial_corr[["2021"]], spain, by = c("municipality", "id"))
st_geometry(df) <- "geometry"
b <-  ggplot() +
  geom_sf(data = df, aes(fill = rho), color = "transparent", size = 1) +
  scale_fill_distiller("", palette = "Spectral") + 
  ggtitle("2021") +
  theme_classic()

df <- merge(spatial_corr[["2022"]], spain, by = c("municipality", "id"))
st_geometry(df) <- "geometry"
c <- ggplot() +
  geom_sf(data = df, aes(fill = rho), color = "transparent", size = 1) +
  scale_fill_distiller("", palette = "Spectral") + 
  ggtitle("2022") +
  theme_classic()

ggpubr::ggarrange(a, b, c,
                  common.legend = TRUE, nrow = 1, ncol = 3, legend = "right")
ggsave(file = paste0(loc.fig, "spatial_correlation_counts_ma_Spain.png"), 
       units = "cm", height = 10, width = 20, bg = "white")

# Plotting temporal correlations -----------------------------------------------
years <- c("2020", "2021", "2022")
months <- c("01", "02", "03", "04", "05", "06" ,"07", "08", "09", "10", "11", "12")

cor_predictions <- data.frame()
for(y in years){
  for (m in months){
    print(paste0("Plotting: ", m, "-", y))
    
    pred_bg <- readRDS(paste0(loc.output, "PREDICTIONS/Counts/mtiger10/tiger_", m ,"_", y, ".rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    colnames(pred_bg) <- c("municipality", "id", "lon", "lat", "prov_name", "bg")
    # pred_bg <- pred_bg %>% filter((pixel_id %in% tiger$pixel_id))
    pred_suit <- readRDS(paste0(loc.output, "PREDICTIONS/Suitability/mtiger3_occu/tiger_", m ,"_", y, "_occu.rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    colnames(pred_suit) <- c("municipality", "id", "lon", "lat", "prov_name", "suit")
    # pred_suit <- pred_suit %>% filter((pixel_id %in% tiger$pixel_id))
    pred_ma <- readRDS(paste0(loc.output, "PREDICTIONS/MA/mtiger7_ma/tiger_", m ,"_", y, "_ma.rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    colnames(pred_ma) <- c("municipality", "id", "lon", "lat", "prov_name", "ma")
    # pred_ma <- pred_ma %>% filter((pixel_id %in% tiger$pixel_id))
    
    # Joining tables
    pred <- merge(pred_bg, pred_suit, by = c("municipality", "id", "lon", "lat", "prov_name"))
    pred <- merge(pred, pred_ma, by = c("municipality", "id", "lon", "lat", "prov_name"))
    
    rm(pred_bg, pred_suit, pred_ma)
    
    cor_row <- data.frame(
      year = y,
      month = as.numeric(m), 
      bg_suit = as.numeric(cor.test(pred$bg, pred$suit, method="spearman")$estimate),
      bg_suit_p = cor.test(pred$bg, pred$suit, method="spearman")$p.value,
      bg_ma = as.numeric(cor.test(pred$bg, pred$ma, method="spearman")$estimate),
      bg_ma_p = cor.test(pred$bg, pred$ma, method="spearman")$p.value,
      suit_ma = as.numeric(cor.test(pred$suit, pred$ma, method="spearman")$estimate),
      suit_ma_p = cor.test(pred$suit, pred$ma, method="spearman")$p.value
    )
    cor_predictions <- rbind(cor_predictions, cor_row)
  }
}

cor_predictions <- tidyr::pivot_longer(
  cor_predictions,
  cols = c(bg_suit, bg_ma, suit_ma),
  names_to = "comparison",
  values_to = "value"
)

# Define una paleta de colores personalizada
custom_colors <- c("2020" = "#490092", "2021" = "#004949", "2022" = "#924900")

# Crear el gráfico con ggplot2
ggplot(cor_predictions, aes(x = month, y = value, color = factor(year), shape = comparison, linetype = comparison)) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c("bg_suit" = 1, "bg_ma" = 2, "suit_ma" = 5)) +
  scale_linetype_manual(values = c("bg_suit" = 5, "bg_ma" = 1, "suit_ma" = 10)) +
  labs(
    y = "Spearman Correlation (S)",
    x = "Month",
    color = "Year",
    shape = "Comparison",
    linetype = "Comparison",
    title = "Monthly Correlation Coefficients by Year"
  ) +
  scale_x_continuous(breaks=seq(3, 11, 1)) +
  scale_y_continuous(breaks=seq(-0.5, 0.9, 0.2)) +
  theme_classic() +
  facet_wrap(~year) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

ggsave(file = paste0(loc.fig, "temporal_correlation_Spain.png"), 
       units = "cm", height = 10, width = 20, bg = "white")

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
  filter(movelab_certainty_category_euro_class_label == "tiger-sp", movelab_certainty_category_euro_class_value >= 1) %>%
  mutate(
    year = as.factor(lubridate::year(date)),
    month = as.numeric(lubridate::month(date)),
    type = "reports"
  ) %>% 
  group_by(year, month, type) %>%
  summarise(n = dplyr::n()) %>%
  mutate(n_prop = n/3320)

tiger <- readRDS(paste0(loc.output, "tiger_clc_clima_ma_monthly.rds")) %>%
  mutate(
    year = as.factor(lubridate::year(end_date)),
    month = as.numeric(lubridate::month(end_date)),
    type = "traps"
  ) %>% 
  group_by(year, month, type) %>%
  summarise(n = sum(females)) %>%
  mutate(n_prop = n/91345)

sampling_size <- rbind(tiger, reports)

custom_colors <- c("traps" =  "#004949", "reports" = "#924900")
ggplot(sampling_size %>% 
         filter(year != "2020" & !(month %in% c(1,2,12))), 
       aes(x = factor(month), y = n_prop, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ year, scales = "free_x") +
  labs(
    x = "Month",
    y = "Relative Sample Size [Ni/max(Ni)]",
    fill = "Type",
    title = "Monthly Realtive Sample Size "
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
ggsave(file = paste0(loc.fig, "Relative_somple_size.png"), 
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



