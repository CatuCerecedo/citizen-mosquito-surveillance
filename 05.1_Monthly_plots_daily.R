######################## Daily predictions #####################################

#' We will plot the daily predictions of three models: BG, MA and MIX

library(ggplot2)
library(tidyverse)
library(sf)

rm(list = ls())

# Directories ------------------------------------------------------------------
# In cluster
# loc.output <- paste0(getwd(), "/Spain_Culex/OUTPUT/")
# loc.data <- paste0(getwd(), "/Spain_Culex/DATA/")
# loc.era5 <- paste0(getwd(), "/Spain_Culex/ERA5_Download/")

# In local
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")
loc.era5 <- "/home/catuxa/Documents/Mosquito_Models/EU_Culex/ERA5_Download/"

sf::sf_use_s2(FALSE)

mdl_name <- "/tiger_inte"
fldr <- "Monthly_Integration"
sub <- "" # with _

# Checking the daily predictions -----------------------------------------------
folder = paste0(loc.output, "PREDICTIONS/", fldr, mdl_name, "/")
file_list = list.files(folder)
file_list = file_list[!grepl("_sd", file_list)]
# file_list = file_list[grepl("_ma.rds", file_list)]

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

# Saving all values
all <- data.frame(
  m = df$m,
  year = df$year
)

bg <- ggplot(df, aes(x = as.numeric(m), y = bg_avrg, color = year)) +
  geom_point() +
  geom_ribbon(aes(ymin = (bg_avrg - bg_se), ymax = (bg_avrg + bg_se), fill = year), alpha = 0.2, color = NA) +
  geom_line(aes(x = as.numeric(m), y = bg_avrg, color = year)) +
  labs(
    y = "Monthly avg. of predicted \ncounts (Count model)",
    x = "Month"
  ) +
  scale_x_continuous(breaks=seq(1, 12, 1)) +
  theme_classic() 
all$bg <- df$bg_avrg

suitability <- ggplot(df, aes(x = as.numeric(m), y = bg_avrg, color = year)) +
  geom_point() +
  geom_ribbon(aes(ymin = (bg_avrg - bg_se), ymax = (bg_avrg + bg_se), fill = year), alpha = 0.2, color = NA) +
  geom_line() +
  labs(
    y = "Monthly avg. of predicted \nsuitability (Suitability model)",
    x = "Month"
  ) +
  scale_x_continuous(breaks=seq(1, 12, 1)) +
  theme_classic() 
all$suit <- df$bg_avrg

ma <- ggplot(df, aes(x = as.numeric(m), y = bg_avrg, color = year)) +
  geom_point() +
  geom_ribbon(aes(ymin = (bg_avrg - bg_se), ymax = (bg_avrg + bg_se), fill = year), alpha = 0.2, color = NA) +
  geom_line() +
  labs(
    y = "Monthly avg. of predicted \nMA probability (MA model)",
    x = "Month"
  ) +
  scale_x_continuous(breaks=seq(1, 12, 1)) +
  theme_classic() 
all$ma <- df$bg_avrg

bites <- ggplot(df, aes(x = as.numeric(m), y = bg_avrg, color = year)) +
  geom_point() +
  geom_ribbon(aes(ymin = (bg_avrg - bg_se), ymax = (bg_avrg + bg_se), fill = year), alpha = 0.2, color = NA) +
  geom_line() +
  labs(
    y = "Monthly avg. of predicted \nMA probability (Bite model)",
    x = "Month"
  ) +
  scale_x_continuous(breaks=seq(1, 12, 1)) +
  theme_classic() 
all$bites <- df$bg_avrg

ggpubr::ggarrange(suitability, bg, ma, nrow = 3, ncol = 1, common.legend = TRUE)
ggsave(file = paste0(loc.fig, "predicted_values_suit_bg_ma_temporal_plot_Spain.png"), 
       units = "cm", height = 20, width = 15, bg = "white")

# ggpubr::ggarrange(suitability_1, suitability, nrow = 2, ncol = 1)
# ggsave(file = paste0(loc.fig, "suitability_SE.png"), 
#        units = "cm", height = 25, width = 25)

a <- ggplot(all, aes(x = ma, y = bg, color = year)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_regline_equation(aes(label = ..eq.label..), formula = y ~ x, 
                                label.x.npc = "left", label.y.npc = 0.95, show.legend = FALSE) +
  labs(
    title = "a)",
    x = "Monthly avg. of predicted \nMA probability",
    y = "Monthly avg. of predicted \ncounts"
  ) +
  theme_classic() 

b <- ggplot(all, aes(x = ma, y = suit, color = year)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_regline_equation(aes(label = ..eq.label..), formula = y ~ x, 
                                label.x.npc = "left", label.y.npc = 0.95, show.legend = FALSE) +
  labs(
    title = "b)",
    x = "Monthly avg. of predicted \nMA probability",
    y = "Monthly avg. of predicted \nsuitability"
  ) +
  theme_classic() 

ggpubr::ggarrange(a, b, nrow = 2, ncol = 1, common.legend = TRUE)
ggsave(file = paste0(loc.fig, "prediction_linear_regression.png"), 
       units = "cm", height = 20, width = 15, bg = "white")

saveRDS(all, file = paste0(loc.output, "all_monthly_pred.rds"))

# 3D prediction plots ----------------------------------------------------------

# Considering the average of Spain (mean of municipalities)
all <- readRDS(file = paste0(loc.output, "all_monthly_pred.rds"))

ggplot(all) +
  # geom_smooth(aes(x = suit, y = bg), method="lm", formula= y ~ poly(x, 2), linetype = "dashed", color = "black") +
  geom_point(aes(x = suit, y = bg, color = ma), alpha = 0.9, size = 3) +
  scale_color_distiller("Monthly Mosquito Alert \npredicted probabilities", palette = "Spectral", 
                        na.value = "transparent") +
  labs(
    x = "Monthly suitability predicted probabilities",
    y = "Monthly count predicted values"
  ) +
  theme_classic()

ggsave(file = paste0(loc.fig, "prediction_polynomical_regression_3D.png"), 
       units = "cm", height = 10, width = 15, bg = "white")

# Adding one y-axis: bites
ggplot(all) +
  geom_point(aes(x = suit, y = bg, color = ma), alpha = 0.9, size = 3) +
  scale_color_distiller("Monthly Mosquito Alert \npredicted probabilities", 
                        palette = "Spectral", na.value = "transparent") +
  labs(
    x = "Monthly suitability predicted probabilities",
    y = "Monthly count predicted values"
  ) +
  scale_y_continuous(
    name = "Monthly count predicted values",  
    sec.axis = sec_axis(~ . * 0.1, name = "Bites probability") 
  ) +
  theme_classic() +
  theme(
    axis.title.y.right = element_text(color = "blue")  # Diferenciar eje secundario
  )

# Considering each municipalities
mdl_name <- "/mtiger3_occu"
fldr <- "Suitability"
sub <- "_occu" # with _

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
  
  df_row = data.frame(id = pred_points$id,
                      bg_avrg = pred_points$pred_count,
                      m = parts[2], 
                      year = parts[3])
  
  df = rbind(df, df_row)
}

# Saving all values
all <- data.frame(
  id = df$id,
  m = df$m,
  year = df$year
)
all$suit <- df$bg_avrg
all$bg <- df$bg_avrg
all$ma <- df$bg_avrg

# saveRDS(all, file = paste0(loc.output, "all_monthly_pred_mun.rds"))
all <- readRDS(file = paste0(loc.output, "all_monthly_pred_mun.rds"))

ggplot(all) +
  geom_point(aes(x = suit, y = bg, color = ma), alpha = 0.6, size = 0.5) +
  scale_color_distiller("Monthly Mosquito Alert \npredicted probabilities", palette = "Spectral", 
                        na.value = "transparent") +
  # geom_smooth(aes(x = suit, y = bg), linetype = "dashed", color = "black
  geom_smooth(aes(x = suit, y = bg), method="lm", formula= y ~ poly(x, 2), linetype = "dashed", color = "black") +
  labs(
    x = "Monthly Suitability predicted probabilities",
    y = "Monthly Count predicted values"
  ) +
  theme_classic()

ggsave(file = paste0(loc.fig, "prediction_polynomical_regression_3D_municipalities.png"), 
       units = "cm", height = 10, width = 15, bg = "white")
