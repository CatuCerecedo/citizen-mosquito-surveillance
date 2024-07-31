######################## Daily predictions #####################################

#' We will plot the daily predictions of three models: BG, MA and MIX

library(ggplot2)
library(tidyverse)
library(sf)

rm(list = ls())

# Directories ------------------------------------------------------------------
# In local
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")
loc.era5 <- "/home/catuxa/Documents/Mosquito_Models/EU_Culex/ERA5_Download/"

# In cluster
loc.output <- paste0(getwd(), "/Spain_Culex/OUTPUT/")
loc.data <- paste0(getwd(), "/Spain_Culex/DATA/")
loc.era5 <- paste0(getwd(), "/Spain_Culex/ERA5_Download/")

sf::sf_use_s2(FALSE)

# Checking the daily predictions -----------------------------------------------
folder = paste0(loc.output, "PREDICTIONS/MA_rep/")
file_list = list.files(folder)
file_list = file_list[!grepl("_sd", file_list)]
# file_list = file_list[grepl("_ma.rds", file_list)]

df <- data.frame()
for (i in file_list){
  print(i)
  
  pred_points <- readRDS(paste0(folder, i)) 
  colnames(pred_points) <- c("municipality", "id", "lon", "lat", "pred_count") 

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
    x = "Date"
  ) +
  scale_x_continuous(breaks=seq(3, 11, 1)) +
  theme_classic() 
all$bg <- df$bg_avrg

suitability <- ggplot(df, aes(x = as.numeric(m), y = bg_avrg, color = year)) +
  geom_point() +
  geom_ribbon(aes(ymin = (bg_avrg - bg_se), ymax = (bg_avrg + bg_se), fill = year), alpha = 0.2, color = NA) +
  geom_line() +
  labs(
    y = "Monthly avg. of predicted \nsuitability (Suitability model)",
    x = "Date"
  ) +
  scale_x_continuous(breaks=seq(3, 11, 1)) +
  theme_classic() 
all$suit <- df$bg_avrg

ma <- ggplot(df, aes(x = as.numeric(m), y = bg_avrg, color = year)) +
  geom_point() +
  geom_ribbon(aes(ymin = (bg_avrg - bg_se), ymax = (bg_avrg + bg_se), fill = year), alpha = 0.2, color = NA) +
  geom_line() +
  labs(
    y = "Monthly avg. of predicted \nMA probability (MA model)",
    x = "Date"
  ) +
  scale_x_continuous(breaks=seq(1, 12, 1)) +
  theme_classic() 
all$ma <- df$bg_avrg

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
