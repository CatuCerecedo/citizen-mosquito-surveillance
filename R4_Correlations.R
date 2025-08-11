################# Temporal and spatial correlations ############################
#' Correlations between TRAP and CITSCI model predictions:
#' (1) Temporal variation of spatial patterns
#' (2) Spatial variation of temporal patterns

################################################################################
library(terra)
library(sf)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(cowplot)

rm(list = ls())

# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")
loc.era5 <- paste0(getwd(), "/ERA5_Download/") # Directory where ERA5 nc files are stores (one file per month)

sf::sf_use_s2(FALSE)

# Temporal variation of spatial patterns ---------------------------------------
years <- c("2020", "2021", "2022")
months <- c("01", "02", "03", "04", "05", "06" ,"07", "08", "09", "10", "11", "12")

spatial_corr <- vector(mode = "list")
for(y in years){ # Baseline dataframes
  pred_bg_all <- readRDS(paste0(loc.output, "PREDICTIONS/COUNTS/mtiger16_trap/tiger_", "07" ,"_", y, "_trap.rds")) %>%
    janitor::clean_names() %>%
    dplyr::select(municipality, id, lon, lat, prov_name)
  
  pred_ma_all <- readRDS(paste0(loc.output, "PREDICTIONS/CITSCI/mtiger7_citsci/tiger_", "07" ,"_", y, "_citsci.rds")) %>%
    janitor::clean_names() %>%
    dplyr::select(municipality, id, lon, lat, prov_name)
  
  for (m in months){
    print(paste0("Plotting: ", m, "-", y))
    
    pred_bg <- readRDS(paste0(loc.output, "PREDICTIONS/COUNTS/mtiger16_trap/tiger_", m ,"_", y, "_trap.rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    pred_bg$pred_count <- rowMeans(pred_bg[,6:ncol(pred_bg)], na.rm = TRUE)
    pred_bg <- pred_bg %>% dplyr::select(municipality, id, lon, lat, prov_name, pred_count)
    pred_bg_all <- merge(pred_bg_all, pred_bg, by = c("municipality", "id", "lon", "lat", "prov_name"))
    
    pred_ma <- readRDS(paste0(loc.output, "PREDICTIONS/CITSCI/mtiger7_citsci/tiger_", m ,"_", y, "_citsci.rds")) %>%
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

# Plotting the temporal correlation over space
df <- merge(spatial_corr[["2020"]], spain, by = c("municipality", "id"))
st_geometry(df) <- "geometry"
a <- ggplot() +
  geom_sf(data = df, aes(fill = rho), color = "transparent", size = 1) +
  scale_fill_distiller("", palette = "Spectral", limits = c(0.7, 1)) + 
  theme_classic(base_size = 8, base_family = "Helvetica") +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
print(mean(df$rho, na.rm = TRUE))
print(sd(df$rho, na.rm = TRUE))

df <- merge(spatial_corr[["2021"]], spain, by = c("municipality", "id"))
st_geometry(df) <- "geometry"
b <-  ggplot() +
  geom_sf(data = df, aes(fill = rho), color = "transparent", size = 1) +
  scale_fill_distiller("", palette = "Spectral", limits = c(0.7, 1)) + 
  theme_classic(base_size = 8, base_family = "Helvetica") +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
print(mean(df$rho, na.rm = TRUE))
print(sd(df$rho, na.rm = TRUE))

df <- merge(spatial_corr[["2022"]], spain, by = c("municipality", "id"))
st_geometry(df) <- "geometry"
c <- ggplot() +
  geom_sf(data = df, aes(fill = rho), color = "transparent", size = 1) +
  scale_fill_distiller("", palette = "Spectral", limits = c(0.7, 1)) + 
  theme_classic(base_size = 8, base_family = "Helvetica") +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
print(mean(df$rho, na.rm = TRUE))
print(sd(df$rho, na.rm = TRUE))

ggpubr::ggarrange(a, b, c,
                  common.legend = TRUE, nrow = 1, ncol = 3, legend = "right")

# Global average
df <- rbind(spatial_corr[["2020"]], spatial_corr[["2021"]], spatial_corr[["2022"]])

df <- df %>%
  group_by(municipality, id) %>%
  summarise(rho = mean(rho, na.rm = TRUE))
df <- merge(df, spain, by = c("municipality", "id"))
st_geometry(df) <- "geometry"


a <- ggplot() +
  geom_sf(data = df, aes(fill = rho), color = "transparent", size = 1) +
  scale_fill_distiller("Spearman's rank\ncorrelation (S)\n", palette = "Spectral", limits = c(0.7, 1)) + 
  theme_classic(base_size = 8, base_family = "Helvetica") +
  coord_sf(expand = FALSE) +
  theme_void(base_size = 12, base_family = "Helvetica") 
print(mean(df$rho, na.rm = TRUE))
print(sd(df$rho, na.rm = TRUE))

df <- rbind(spatial_corr[["2020"]], spatial_corr[["2021"]], spatial_corr[["2022"]])

df <- df %>%
  group_by(municipality, id) %>%
  summarise(rho_mean = mean(rho, na.rm = TRUE),
            rho_sd =  sd(rho, na.rm = TRUE),
            rho_se = rho_sd/sqrt(3)) %>%
  arrange(desc(rho_mean)) %>%
  ungroup()
df$id = 1:8123

b <- ggplot(df, aes(y = rho_mean, x = id)) +
  geom_point(alpha = 0.4, size = 0.1, color = "#653496") +
  geom_ribbon(aes(ymin = (rho_mean - rho_se), ymax = (rho_mean + rho_se)),alpha = 0.2) +
  labs (x = "# Municipality", y = "Spearman's rank correlation (S)") +
  theme_classic(base_size = 12, base_family = "Helvetica") 

par(mar = c(0, 0, 0, 0))
a + b + plot_annotation(tag_levels = "a", tag_suffix = "") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(file = paste0(loc.fig, "Fig_2.pdf"), 
       width = 25, height = 18, dpi = 600, units = "cm", device = cairo_pdf)


# Spatial variation of temporal patterns ---------------------------------------
years <- c("2020", "2021", "2022")
months <- c("01", "02", "03", "04", "05", "06" ,"07", "08", "09", "10", "11", "12")

cor_predictions <- data.frame()
for(y in years){
  for (m in months){
    print(paste0("Plotting: ", m, "-", y))
    
    pred_bg <- readRDS(paste0(loc.output, "PREDICTIONS/COUNTS/mtiger16_trap/tiger_", m ,"_", y, "_trap.rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    pred_bg$bg <- rowMeans(pred_bg[,6:ncol(pred_bg)], na.rm = TRUE)
    pred_bg <- pred_bg %>% dplyr::select(municipality, id, lon, lat, prov_name, bg)
    
    pred_ma <- readRDS(paste0(loc.output, "PREDICTIONS/CITSCI/mtiger7_citsci/tiger_", m ,"_", y, "_citsci.rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    pred_ma$ma <- rowMeans(pred_ma[,6:ncol(pred_ma)], na.rm = TRUE)
    pred_ma <- pred_ma %>% dplyr::select(municipality, id, lon, lat, prov_name, ma)
    
    # Joining tables
    pred <- merge(pred_bg, pred_ma, by = c("municipality", "id", "lon", "lat", "prov_name"))
    
    cor_row <- data.frame(
      year = y,
      month = as.numeric(m), 
      correl = as.numeric(cor.test(pred$bg, pred$ma, method="spearman")$estimate),
      correl_p = cor.test(pred$bg, pred$ma, method="spearman")$p.value
    )
    cor_predictions <- rbind(cor_predictions, cor_row)
  }
}

cor_predictions <- tidyr::pivot_longer(
  cor_predictions,
  cols = c(correl, correl_p),
  names_to = "comparison",
  values_to = "value"
)

ggplot(cor_predictions, aes(x = month, y = value)) +
  geom_point(size = 4, alpha = 0.6, color = "#653496") +
  geom_line(size = 1, alpha = 0.6, color = "#653496") +
  geom_hline(yintercept = mean(cor_predictions$value), linetype = "dashed") +
  labs(
    y = "Spearman's rank Correlation (S)",
    x = "Month"
  ) +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  scale_y_continuous(breaks = seq(-0.5, 0.9, 0.2)) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  facet_wrap(~year) 

ggsave(file = paste0(loc.fig, "Fig_3.pdf"), 
       width = 25, height = 18, dpi = 600, units = "cm", device = cairo_pdf) # Without panel b




