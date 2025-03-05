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
# loc.output <- paste0(getwd(), "/Spain_Culex/OUTPUT/")
# loc.data <- paste0(getwd(), "/Spain_Culex/DATA/")
# loc.era5 <- paste0(getwd(), "/Spain_Culex/ERA5_Download/")

sf::sf_use_s2(FALSE)

# mdl <- "mtiger16"
# mdl_name <- "/mtiger16"
# fldr <- "Counts"
# sub <- "" # with _

mdl <- "mtiger3_occu"
mdl_name <- "/mtiger3_occu"
fldr <- "Suitability"
sub <- "_occu" # with _

# mdl <- "mtiger7_ma"
# mdl_name <- "/mtiger7_ma"
# fldr <- "MA"
# sub <- "_ma" # with _

# mdl <- "mbites_7_reduce"
# mdl_name <- "/mbites_7_reduce"
# fldr <- "Bites"
# sub <- "bts" # with _


# Checking the daily predictions -----------------------------------------------
folder = paste0(loc.output, "PREDICTIONS/", fldr, mdl_name, "/")
file_list = list.files(folder)
file_list = file_list[!grepl("_sd", file_list)]
# file_list = file_list[grepl("_ma.rds", file_list)]

df <- data.frame()
for (i in file_list){
  print(i)
  
  pred_points <- readRDS(paste0(folder, i)) 
  colnames(pred_points) <- c("municipality","id", "prov_name", "lon", "lat", "pred_count") 

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
    y = "Monthly avg. of predicted \ncounts (CM)",
    x = "Date"
  ) +
  scale_x_continuous(breaks=seq(1, 12, 1)) +
  theme_classic() 
all$bg <- df$bg_avrg

suitability <- ggplot(df, aes(x = as.numeric(m), y = bg_avrg, color = year)) +
  geom_point() +
  geom_ribbon(aes(ymin = (bg_avrg - bg_se), ymax = (bg_avrg + bg_se), fill = year), alpha = 0.2, color = NA) +
  geom_line() +
  labs(
    y = "Monthly avg. of predicted \nsuitability (SM)",
    x = "Date"
  ) +
  scale_x_continuous(breaks=seq(3, 12, 1)) +
  theme_classic() 
all$suit <- df$bg_avrg

ma <- ggplot(df, aes(x = as.numeric(m), y = bg_avrg, color = year)) +
  geom_point() +
  geom_ribbon(aes(ymin = (bg_avrg - bg_se), ymax = (bg_avrg + bg_se), fill = year), alpha = 0.2, color = NA) +
  geom_line() +
  labs(
    y = "Monthly avg. of predicted \nMA probability (CSM)",
    x = "Date"
  ) +
  scale_x_continuous(breaks=seq(1, 12, 1)) +
  theme_classic() 
all$ma <- df$bg_avrg

bites <- ggplot(df, aes(x = as.numeric(m), y = bg_avrg, color = year)) +
  geom_point() +
  geom_ribbon(aes(ymin = (bg_avrg - bg_se), ymax = (bg_avrg + bg_se), fill = year), alpha = 0.2, color = NA) +
  geom_line() +
  labs(
    y = "Monthly avg. of predicted \nbite probability (BM)",
    x = "Date"
  ) +
  scale_x_continuous(breaks=seq(1, 12, 1)) +
  theme_classic() 
all$bites <- df$bg_avrg

ggpubr::ggarrange(suitability, bg, ma, bites, nrow = 4, ncol = 1, common.legend = TRUE)
ggsave(file = paste0(loc.fig, "predicted_values_suit_bg_ma_bites_temporal_plot_Spain.png"), 
       units = "cm", height = 25, width = 15, bg = "white")

saveRDS(all, file = paste0(loc.output, "all_daily_data_4_models.rds"))

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
