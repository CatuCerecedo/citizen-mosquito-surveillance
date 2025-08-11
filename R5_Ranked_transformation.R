################ Ranked Transformation: dissimilarity analysis #################

#' Rank transformation of TRAP and CITSCI model predictions
#' Evaluating the drivers of spatial dissimilarities during summer months (based
#' on GLMMs)

################################################################################
library(glmmTMB)
library(MuMIn)
library(DHARMa)
library(ggplot2)
library(tidyverse)
library(sf)
library(parallel)
library(patchwork)

rm(list = ls())

# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")

# Model names ------------------------------------------------------------------
mdl <- "mtiger16_trap"
mdl_name <- "/mtiger16_trap"
fldr <- "COUNTS"
sub <- "_trap" 

mdl_ma <- "mtiger7_citsci"
mdl_name_ma <- "/mtiger7_citsci"
fldr_ma <- "CITSCI"
sub_ma <- "_citsci" 

months = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
years = c("2020", "2021", "2022")

# Function for rank transformation ---------------------------------------------
rank_order_conversion <- function(vctr) {
  rank = rank(vctr, na.last = "keep", ties.method = "first")       
  rank_normalized = rank / max(rank, na.rm = TRUE) 
  return(rank_normalized)
}

# Rank transformation ----------------------------------------------------------
pred_all <- data.frame()
for(y in years){
  
  print(y)
  for(i in 1:length(months)){
    print(i)
    if (i == 1){
      pred <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                             "/tiger_", months[i], "_", y, sub, ".rds")) %>%
        janitor::clean_names() %>%
        st_drop_geometry()
      
      pred$count <- rowMeans(pred[,6:ncol(pred)], na.rm = TRUE)
      
      pred <- pred[c("municipality","id", "prov_name", "lon", "lat", "count")]
    } else {
      p <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                          "/tiger_", months[i], "_", y, sub, ".rds")) %>%
        janitor::clean_names() %>%
        st_drop_geometry()
      
      p$count <- rowMeans(p[,6:ncol(p)], na.rm = TRUE)
      
      p <- p[c("municipality","id", "prov_name", "lon", "lat", "count")]
      
      pred <- cbind(pred, p["count"])
    }
  }
  colnames(pred) <- c("municipality","id", "prov_name", "lon", "lat", months)
  pred <- pred %>%
    pivot_longer(
      cols = "01":"12",
      names_to = "month",
      values_to = "counts"
    )
  pred$year <- y
  pred_all <- rbind(pred_all, pred)
}
pred_all$rank <- rank_order_conversion(pred_all$counts)
pred_all_c <- pred_all %>% dplyr::select(-rank)
pred_all$counts <-  NULL


pred_all_ma <- data.frame()
for(y in years){
  print(y)
  
  for(i in 1:length(months)){
    print(i)
    if (i == 1){
      pred <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr_ma, mdl_name_ma,
                             "/tiger_", months[i], "_", y, sub_ma, ".rds")) %>%
        janitor::clean_names() %>%
        st_drop_geometry()
      
      pred$count <- rowMeans(pred[,6:ncol(pred)], na.rm = TRUE)
      
      pred <- pred[c("municipality","id", "prov_name", "lon", "lat", "count")]
    } else {
      p <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr_ma, mdl_name_ma,
                          "/tiger_", months[i], "_", y, sub_ma, ".rds")) %>%
        janitor::clean_names() %>%
        st_drop_geometry()
      
      p$count <- rowMeans(p[,6:ncol(p)], na.rm = TRUE)
      
      p <- p[c("municipality","id", "prov_name", "lon", "lat", "count")]
      
      pred <- cbind(pred, p["count"])
    }
  }
  colnames(pred) <- c("municipality","id", "prov_name", "lon", "lat", months)
  pred <- pred %>%
    pivot_longer(
      cols = "01":"12",
      names_to = "month",
      values_to = "counts"
    )
  pred$year <- y
  pred_all_ma <- rbind(pred_all_ma, pred)
}
pred_all_ma$rank <- rank_order_conversion(pred_all_ma$counts)
pred_all_ma_c <- pred_all_ma %>% dplyr::select(-rank)
pred_all_ma$counts <-  NULL

diff_table <- merge(pred_all_ma, pred_all, by = names(pred_all)[1:7])
colnames(diff_table) <- c(names(pred_all)[1:7], "ma_rank", "model_rank")
diff_table$diff <- diff_table$ma_rank - diff_table$model_rank

c_table <- merge(pred_all_ma_c, pred_all_c, by = names(pred_all_c)[c(1:6,8)])
colnames(c_table) <- c(names(pred_all_c)[c(1:6,8)], "ma", "model")

ggplot(c_table) +
  geom_point(aes(y = ma, x = model)) +
  facet_wrap(~ year) +
  theme_classic() +
  labs(
    x = paste("Predictions", fldr),
    y = paste("Predictions Citizen Science")
  )

# Ploting decorrelations maps --------------------------------------------------
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  )

diff_rank_plot <- ggplot() +
  geom_sf(data = diff_table, aes(fill = diff), color = "transparent",
          size = 0.01, alpha = 0.8, na.rm = TRUE) +
  scale_fill_distiller("Diff. ranked\npredictions\n", palette = "Spectral", na.value = "transparent", limits = c(-1, 1)) + 
  ggtitle("CITSCI - TRAP") +
  theme_void(base_size = 8, base_family = "Helvetica") 

# Fig. 3
cor_predictions <- readRDS(file = paste0(loc.output, "temporal_corr.rds"))
cor_predictions <- cor_predictions %>% filter(comparison == "count_ma")

cor_plot <- ggplot(cor_predictions, aes(x = month, y = value)) +
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

maps1 <- (c / d) + plot_layout(guides = 'collect') & theme(legend.position = "right") +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA),
    legend.background = element_rect(fill = NA, color = NA)
  ) 

maps2 <- (maps1 | diff_rank_plot)  + 
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA),
    legend.background = element_rect(fill = NA, color = NA)
  )

ggdraw() +
  draw_plot(cor_plot) +
  draw_plot(maps2, x = 0.26, y = 0, width = 0.55, height = 0.55) +
  draw_grob(
    grid::roundrectGrob(
      x = unit(0.27 + 0.55 / 2, "npc"),  
      y = unit(0 + 0.55 / 2, "npc"),    
      width = unit(0.53, "npc"),         
      height = unit(0.38, "npc"),        
      r = unit(0.05, "snpc"),             
      gp = gpar(fill = NA, col = "black", lwd = 2)
    )
  ) + draw_plot_label(
    label = c("a", "b"),
    x = c(0, 0.29), 
    y = c(1, 0.46), 
    size = 12
  )

ggsave(file = paste0(loc.fig, "Fig_3.pdf"), 
       width = 25, height = 18, dpi = 600, units = "cm", device = cairo_pdf)

# calculate the weather variables-----------------------------------------------
diff_table_wth <- data.frame()
for(y in years){
  print(y)
  for (i in 1:length(months)){
    print(months[i])
    wth <- readRDS(paste0(loc.output, "monthly_weather_data/prep_", 
                          months[i], "-", y, ".rds"))[, c(1, 9:12, 32)] # Monthly average variables calculated from ERA5
    wth$year <- y
    wth$month <- months[i]
    
    diff_sample_table <- diff_table %>%
      filter(year == y & month == months[i])
    diff_sample_table <- merge(diff_sample_table, wth, by = c("id", "year", "month"),
                               all.x = TRUE)
    diff_table_wth <- rbind(diff_table_wth, diff_sample_table)
  }
}
rm(y, i, diff_sample_table)

summary(diff_table_wth)

# decorrelationrelation months 
diff_table_wth <- diff_table_wth %>%
  mutate(
    decorrelation = as.factor(case_when(
      year == "2020" & month %in% c("07", "08", "09") ~ TRUE,
      year == "2021" & month %in% c("08", "09") ~ TRUE,
      year == "2022" & month %in% c("06", "07", "08", "09") ~ TRUE,
      .default = FALSE
    ))
  )

# Modelling the linear regression of differences -------------------------------
diff_table_wth <- diff_table_wth %>%
  filter(decorrelation == TRUE)
hist(diff_table_wth$diff)

diff_model <- glmmTMB(diff ~ scale(min_temperature) + scale(mean_temperature) + scale(max_temperature) +
                     scale(precipitation) + scale(mean_relative_humidity) +
                     (1|id) + (1|month), 
                   data = diff_table_wth)
summary(dmodel1)

options(na.action = "na.fail")
model_dredge <- dredge(diff_model, rank = "AICc", trace = 2,  m.lim = c(2, 5))
options(na.action = "na.omit")
model.avg(model_dredge)

diff_model_d <- glmmTMB(diff ~ scale(mean_temperature) + scale(mean_relative_humidity) +
                     (1|id) + (1|y), 
                   data = diff_table_wth)
summary(diff_model_d)
r.squaredGLMM(dmodel) 

simulationOutput <- simulateResiduals(fittedModel = diff_model_d, plot =T)
testUniformity(simulationOutput)
testDispersion(simulationOutput)
testOutliers(simulationOutput)
testOverdispersion(simulationOutput)


variables <- c("mean_temperature", "mean_relative_humidity")
zero_values <- sapply(variables, function(var) find_zero(var, diff_model_d, diff_table_wth))

me <- as.data.frame(predict_response(diff_model_d, terms = c("mean_temperature")))
a <- ggplot(me, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +                                               
  geom_vline(xintercept = zero_values["mean_temperature"], color = "grey50", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey50", alpha = 0.6) +
  geom_line(color = "#653496", size = 1.2) + 
  labs(
    title = "", 
    y = TeX("$\\Delta$ Predictions [CITSI - TRAP]"),
    x = "Mean temperature (ºC)",
    caption = paste("TN =", round(zero_values["mean_temperature"], 3))
  ) +
  theme_classic(base_size = 12, base_family = "Helvetica")

me <- as.data.frame(predict_response(diff_model_d, terms = c("mean_relative_humidity")))
b <-  ggplot(me, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +   
  geom_vline(xintercept = zero_values["mean_relative_humidity"], color = "grey50", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey50", alpha = 0.6) +
  geom_line(color = "#653496", size = 1.2) + 
  labs(title = "", 
       y = TeX("$\\Delta$ Predictions [CITSI - TRAP]"),
       x = "Mean Relative Humidity (%)",
       caption = paste("MRH =", round(zero_values["mean_relative_humidity"], 3))) +
  theme_classic(base_size = 12, base_family = "Helvetica")

a + b +
  plot_annotation(tag_levels = c("a"), tag_suffix = "") 



