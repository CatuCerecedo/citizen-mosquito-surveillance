############################## Prediction plots ################################

#' Predictions plots and maps of TRAP and CITSCI model

# Dependencies -----------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(sf)
library(terra)
library(patchwork)

rm(list = ls())
sf::sf_use_s2(FALSE)

# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")
loc.era5 <- paste0(getwd(), "/ERA5_Download/") # Directory where ERA5 nc files are stores (one file per month)

# Monthly average predictions --------------------------------------------------
monthy_avrg_plot <- function(fldr, mdl_name, which = "TRAP"){
  folder = paste0(loc.output, "PREDICTIONS/", fldr, mdl_name, "/")
  file_list = list.files(folder)
  file_list = file_list[!grepl("_sd", file_list)]
  
  df <- data.frame()
  for (i in file_list){
    print(i)
    
    pred_points <- readRDS(paste0(folder, i)) 
    pred_points$pred_count <- rowMeans(pred_points[,6:ncol(pred_points)], na.rm = TRUE)
    
    pred_points <- pred_points[c("municipality","id", "prov_name", "lon", "lat", "pred_count")] 
    
    parts <- unlist(strsplit(i, "_|\\.rds"))
    
    df_row = data.frame(bg_avrg = mean(pred_points$pred_count, na.rm = TRUE), 
                        bg_se = sd(pred_points$pred_count, na.rm = TRUE)/sqrt(length(pred_points)),
                        m = parts[2], 
                        year = parts[3])
    
    df = rbind(df, df_row)
  }
  
  if (which == "TRAP"){
    avrg_monthly <- ggplot(df, aes(x = as.numeric(m), y = bg_avrg, color = year)) +
      geom_point() +
      geom_ribbon(aes(ymin = (bg_avrg - bg_se), ymax = (bg_avrg + bg_se), fill = year), alpha = 0.2, color = NA) +
      scale_fill_manual(name = "Year", values = c("#653496", "#74B652", "#DAA521")) +
      geom_line(aes(x = as.numeric(m), y = bg_avrg, color = year)) +
      scale_color_manual(name = "Year", values = c("#653496", "#74B652", "#DAA521")) +
      labs(
        y = "Monthly avg. of predicted \ncounts (TRAP model)",
        x = "Month"
      ) +
      scale_x_continuous(breaks=seq(1, 12, 1)) +
      theme_classic(base_size = 10, base_family = "Helvetica") 
  } else if (which == "CITSCI"){
    avrg_monthly <- ggplot(df, aes(x = as.numeric(m), y = bg_avrg, color = year)) +
      geom_point() +
      geom_ribbon(aes(ymin = (bg_avrg - bg_se), ymax = (bg_avrg + bg_se), fill = year), alpha = 0.2, color = NA) +
      scale_fill_manual(name = "Year", values = c("#653496", "#74B652", "#DAA521")) +
      geom_line(aes(x = as.numeric(m), y = bg_avrg, color = year)) +
      scale_color_manual(name = "Year", values = c("#653496", "#74B652", "#DAA521")) +
      labs(
        y = "Monthly avg. of predicted\nMA probability (CITSCI model)\n",
        x = "Month"
      ) +
      scale_x_continuous(breaks=seq(1, 12, 1)) +
      ylim(c(0, 1)) +
      theme_classic(base_size = 10, base_family = "Helvetica") 
  }
  return(avrg_monthly)
}

trap_avrg_monthly <- monthy_avrg_plot(fldr = "COUNTS", mdl_name = "mtiger16_trap", which = "TRAP")
citsci_avrg_monthly <- monthy_avrg_plot(fldr = "COUNTS", mdl_name = "mtiger7_ma", which = "CITSCI")

## Spatial average predictions -------------------------------------------------
years <- c("2020", "2021", "2022")
months <- c("01", "02", "03", "04", "05", "06" ,"07", "08", "09", "10", "11", "12")

tiger_maps <- list()
month_values <- list()
iter <- 0

mdl <- "mtiger16_trap"
mdl_name <- "/mtiger16_trap"
fldr <- "Counts"
sub <- "_trap" 

for(y in years){
  pred_month <- list()
  
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
      scale_fill_distiller("Predicted\nCounts", palette = "Spectral",
                           na.value = "transparent",
      ) +
      theme_void(base_size = 12, base_family = "Helvetica") 
    
    print(max(pred$pred_count, na.rm = TRUE))
    
    pred_month[[iter]] <- pred$pred_count
    
    tiger_maps[[iter]] <- plt
    month_values[[iter]] <- pred$pred_count %>% st_drop_geometry()
  }
  pred_month <- as.data.frame(pred_month)
  colnames(pred_month) <- months
  
  pred_month$total <- rowMeans(pred_month)
  pred_month <- cbind(pred %>% dplyr::select(-pred_count),
                      pred_month["total"])
}

a5_trap <- tiger_maps[[29]] # for Fig. 1
a8_trap <- tiger_maps[[32]] # for Fig. 1

tiger_maps <- list()
month_values <- list()
iter <- 0

mdl <- "mtiger7_citsci"
mdl_name <- "/mtiger7_citsci"
fldr <- "CITSCI"
sub <- "_citsci" 

for(y in years){
  pred_month <- list()
  
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
      scale_fill_distiller("Predicted\nMA_probability", palette = "Spectral",
                           na.value = "transparent",
      ) +
      theme_void(base_size = 12, base_family = "Helvetica") 
    
    print(max(pred$pred_count, na.rm = TRUE))
    
    pred_month[[iter]] <- pred$pred_count
    
    tiger_maps[[iter]] <- plt
    month_values[[iter]] <- pred$pred_count %>% st_drop_geometry()
  }
  pred_month <- as.data.frame(pred_month)
  colnames(pred_month) <- months
  
  pred_month$total <- rowMeans(pred_month)
  pred_month <- cbind(pred %>% dplyr::select(-pred_count),
                      pred_month["total"])
}

a5_citsci <- tiger_maps[[29]] # for Fig. 1
a8_citsci <- tiger_maps[[32]] # for Fig. 1

(trap_avrg_monthly + a5_bg + a8_bg) /
  (citsci_avrg_monthly + a5_citsci + a8_citsci) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag =  element_text(size = 12, face = "bold"))

ggsave(file = paste0(loc.fig, "Fig_1.pdf"), 
       width = 25, height = 18, dpi = 600, units = "cm", device = cairo_pdf)

