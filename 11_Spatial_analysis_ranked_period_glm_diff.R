########################## Spatial analysis: diff###############################
library(ggplot2)
library(tidyverse)
library(sf)
library(parallel)
library(rstatix)
library(patchwork)
library(lme4)
library(ggeffects)
library(latex2exp)
library(glmmTMB)

rm(list = ls())

# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")

# Model names ------------------------------------------------------------------
# mdl <- "mtiger3_occu"
# mdl_name <- "/mtiger3_occu"
# fldr <- "Suitability"
# sub <- "_occu" # with _

mdl <- "mtiger16"
mdl_name <- "/mtiger16"
fldr <- "Counts"
sub <- "" # with _

# mdl <- "mbites_7_reduce"
# mdl_name <- "/mbites_7_reduce"
# fldr <- "Bites"
# sub <- "bts" # with _

mdl_ma <- "mtiger7_ma"
mdl_name_ma <- "/mtiger7_ma"
fldr_ma <- "MA"
sub_ma <- "_ma" # with _

months = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
years = c("2020", "2021", "2022")

# Load files -------------------------------------------------------------------
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  )

# Load the prediction in which you are interested:
rank_order_conversion <- function(vctr) {
  rank = rank(vctr, na.last = "keep", ties.method = "first")       
  rank_normalized = rank / max(rank, na.rm = TRUE) 
  return(rank_normalized)
}

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
diff_table$diff <- diff_table$ma_rank - diff_table$model_rank # 292428

c_table <- merge(pred_all_ma_c, pred_all_c, by = names(pred_all_c)[c(1:6,8)])
colnames(c_table) <- c(names(pred_all_c)[c(1:6,8)], "ma", "model")

ggplot(c_table) +
  geom_point(aes(y = ma, x = model)) +
  # geom_smooth(aes(y = ma, x = model)) +
  facet_wrap(~ year) +
  theme_classic() +
  labs(
    x = paste("Predictions", fldr),
    y = paste("Predictions Citizen Science")
  )

# calculate the weather variables-----------------------------------------------
diff_table_wth <- data.frame()
for(y in years){
  print(y)
  for (i in 1:length(months)){
    print(months[i])
  wth <- readRDS(paste0(loc.output, "monthly_weather_data/prep_", 
                        months[i], "-", y, ".rds"))[, c(1, 9:12, 32)]
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

# decorrelationrelation var 
diff_table_wth <- diff_table_wth %>%
  mutate(
    decorrelation = as.factor(case_when(
      year == "2020" & month %in% c("07", "08") ~ TRUE,
      year == "2021" & month %in% c("08", "09") ~ TRUE,
      year == "2022" & month %in% c("07", "08") ~ TRUE,
      .default = FALSE
    ))
  )

# Modelling the linear regression of differences -------------------------------
hist(diff_table_wth$diff)

# Distribution between -1 to 1
# diff_table_wth$diff_beta <- (diff_table_wth$diff + 1) / 2
# hist(diff_table_wth$diff_beta)
# 
# dmodel_beta <- glmmTMB(diff_beta ~ scale(min_temperature) + scale(mean_temperature) + scale(max_temperature) +
#                     scale(precipitation) + scale(mean_relative_humidity) +
#                     (1|id) + (1|month),
#                   beta_family(link = "logit"),
#                   data = diff_table_wth)
# summary(dmodel_beta)
# MuMIn::r.squaredGLMM(dmodel_beta) 
#            R2m       R2c
# [1,] 0.3211654 0.9216631

dmodel_unit <- glmmTMB(diff ~ min_temperature + mean_temperature + max_temperature +
                    precipitation + mean_relative_humidity +
                    (1|id) + (1|month), 
                  data = diff_table_wth)
summary(dmodel_unit)
MuMIn::r.squaredGLMM(dmodel_unit) 

dmodel <- glmmTMB(diff ~ scale(min_temperature) + scale(mean_temperature) + scale(max_temperature) +
                 scale(precipitation) + scale(mean_relative_humidity) +
                (1|id) + (1|month), 
              data = diff_table_wth)
summary(dmodel)
MuMIn::r.squaredGLMM(dmodel) 
#            R2m       R2c
# [1,] 0.2475445 0.7549823
# Marginal: the variance explained by the fixed effects
# Conditional: variance explained by the entire model

AIC(dmodel, dmodel_unit)

## Calculate zeros -----
find_zero <- function(var_name, model, data) {
  print(var_name)
  f <- function(x) {
    temp_data <- data
    temp_data[[var_name]] <- x 
    pred <- predict(model, newdata = temp_data, type = "response", re.form = NA) 
    mean(pred) 
  }
  root <- uniroot(f, interval = c(min(data[[var_name]]), max(data[[var_name]])))
  return(root$root)
}
variables <- c("min_temperature", "max_temperature", "mean_temperature", 
               "precipitation", "mean_relative_humidity")
zero_values <- sapply(variables, function(var) find_zero(var, dmodel_unit, diff_table_wth))

# values_cos <- c(
#   "#abdda4", "#abdda4", "#abdda4", "#abdda4", "#abdda4", "#abdda4", "#d53e4f", "#d53e4f",
#   "#d53e4f", "#abdda4", "#abdda4", "#abdda4")

me <- predict_response(dmodel_unit, terms = c("min_temperature")) #, "month"), type = "random")
a <- plot(me, show_ci = TRUE) +
  # scale_color_manual(values = values_cos) +
  geom_vline(xintercept = zero_values["min_temperature"], color = "grey33", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey33", alpha = 0.6) +
  labs(title = "", 
       y = TeX("$\\Delta$ Predictions [CSM - SM]"),
       x = "Min temperature (ºC)",
       caption = paste("TN =", round(zero_values["min_temperature"], 3))) +
  theme_classic()

me <- predict_response(dmodel_unit, terms = c("max_temperature")) #, "month"), type = "random")
b <- plot(me, show_ci = TRUE) +
  # scale_color_manual(values = values_cos) +
  # geom_vline(xintercept = zero_values["max_temperature"], color = "grey33", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey33", alpha = 0.6) +
  labs(title = "", 
       y = TeX("$\\Delta$ Predictions [CITSI - COUNT]"),
       x = "Max temperature (ºC)",
       caption = paste("TM =", round(zero_values["max_temperature"], 3))) +
  theme_classic()

me <- predict_response(dmodel_unit, terms = c("mean_temperature")) #, "month"), type = "random")
c <- plot(me, show_ci = TRUE) +
  # scale_color_manual(values = values_cos) +
  geom_vline(xintercept = zero_values["mean_temperature"], color = "grey33", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey33", alpha = 0.6) +
  labs(title = "", 
       y = TeX("$\\Delta$ Predictions [CITSI - COUNT]"),
       x = "Mean temperature (ºC)",
       caption = paste("TMn =", round(zero_values["mean_temperature"], 3))) +
  theme_classic()

me <- predict_response(dmodel_unit, terms = c("precipitation")) #, "month"), type = "random")
d <- plot(me, show_ci = TRUE) +
  # scale_color_manual(values = values_cos) +
  geom_vline(xintercept = zero_values["precipitation"], color = "grey33", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey33", alpha = 0.6) +
  labs(title = "", 
       y = TeX("$\\Delta$ Predictions [CITSI - COUNT]"),
       x = "Precipitation (m)",
       caption = paste("PPT =", round(zero_values["precipitation"], 3))) +
  theme_classic()

me <- predict_response(dmodel_unit, terms = c("mean_relative_humidity")) #, "month"), type = "random")
e <- plot(me, show_ci = TRUE) +
  # scale_color_manual(values = values_cos) +
  geom_vline(xintercept = zero_values["mean_relative_humidity"], color = "grey33", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey33", alpha = 0.6) +
  labs(title = "", 
       y = TeX("$\\Delta$ Predictions [CITSI - COUNT]"),
       x = "Mean Relative Humidity (%)",
       caption = paste("MRH =", round(zero_values["mean_relative_humidity"], 3))) +
  theme_classic()

a + b + c + d + e +
  plot_annotation(tag_levels = c("A"), tag_suffix = ")") 

ggsave(paste0(loc.fig, "Spatial_clusters/glm/diff_model/month_mar_effect_", fldr, "_", fldr_ma, ".png"),
       width = 35,  height = 22, units = "cm")

# Mapping the ideas
diff_table_wth <- diff_table_wth %>%
  mutate(
    classification = case_when(
      diff > 0 ~ "CSM > CM",
      diff < 0 ~ "CSM < CM",
      .default = NA
    )
  )

diff_table_wth <- merge(diff_table_wth, c_table, by = names(c_table)[1:7], all.x = TRUE)

ggplot(diff_table_wth) +
  geom_point(aes(y = ma, x = model, color = classification)) + 
  facet_wrap(~ year)

# Borrar
diff_table_wth <- diff_table_wth %>%
  mutate(
    class_wth = case_when(
      mean_temperature > 14.357 & mean_relative_humidity > 63.87 ~ "CSM > CM",
      mean_temperature < 14.357 & mean_relative_humidity < 63.87 ~ "CSM < CM",
      .default = NA
    )
  )

ggplot(diff_table_wth) +
  geom_point(aes(x = mean_temperature, y = mean_relative_humidity, color = class_wth)) +
  theme_classic()

# Modelling the linear regression of differences:decorrelation period ----------
diff_table_wth <- diff_table_wth %>%
  filter(decorrelation == TRUE)
hist(diff_table_wth$diff)

dmodel_unit <- glmmTMB(diff ~ min_temperature + mean_temperature + max_temperature +
                         precipitation + mean_relative_humidity +
                         (1|id) + (1|month), 
                       data = diff_table_wth)
summary(dmodel_unit)
MuMIn::r.squaredGLMM(dmodel_unit) 
#            R2m       R2c
# [1,] 0.2876306 0.6450103
dmodel <- glmmTMB(diff ~ scale(min_temperature) + scale(mean_temperature) + scale(max_temperature) +
                    scale(precipitation) + scale(mean_relative_humidity) +
                    (1|id) + (1|month), 
                  data = diff_table_wth)
summary(dmodel)
MuMIn::r.squaredGLMM(dmodel) 
#            R2m       R2c
# [1,] 0.2876236 0.6450191

AIC(dmodel, dmodel_unit)

variables <- c("min_temperature",  "mean_temperature", "mean_relative_humidity" )
zero_values <- sapply(variables, function(var) find_zero(var, dmodel_unit, diff_table_wth))
# Tres variables que la predicción nunca llega al cero

me <- predict_response(dmodel_unit, terms = c("min_temperature")) #, "month"), type = "random")
a <- plot(me, show_ci = TRUE) +
  # scale_color_manual(values = values_cos) +
  geom_vline(xintercept = zero_values["min_temperature"], color = "grey33", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey33", alpha = 0.6) +
  labs(title = "", 
       y = TeX("$\\Delta$ Predictions [CSM - CM]"),
       x = "Min temperature (ºC)",
       # caption = paste("TN =", round(zero_values["min_temperature"], 3))
       ) +
  theme_classic()

me <- predict_response(dmodel_unit, terms = c("max_temperature")) #, "month"), type = "random")
b <- plot(me, show_ci = TRUE) +
  # scale_color_manual(values = values_cos) +
  # geom_vline(xintercept = zero_values["max_temperature"], color = "grey33", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey33", alpha = 0.6) +
  labs(title = "", 
       y = TeX("$\\Delta$ Predictions [CSM - CM]"),
       x = "Max temperature (ºC)",
       caption = paste("TM =", round(zero_values["max_temperature"], 3))) +
  theme_classic()

me <- predict_response(dmodel_unit, terms = c("mean_temperature")) #, "month"), type = "random")
c <- plot(me, show_ci = TRUE) +
  # scale_color_manual(values = values_cos) +
  geom_vline(xintercept = zero_values["mean_temperature"], color = "grey33", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey33", alpha = 0.6) +
  labs(title = "", 
       y = TeX("$\\Delta$ Predictions [CSM - CM]"),
       x = "Mean temperature (ºC)",
       caption = paste("TMn =", round(zero_values["mean_temperature"], 3))) +
  theme_classic()

me <- predict_response(dmodel_unit, terms = c("precipitation")) #, "month"), type = "random")
d <- plot(me, show_ci = TRUE) +
  # scale_color_manual(values = values_cos) +
  # geom_vline(xintercept = zero_values["precipitation"], color = "grey33", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey33", alpha = 0.6) +
  labs(title = "", 
       y = TeX("$\\Delta$ Predictions [CSM - CM]"),
       x = "Precipitation (m)",
       caption = paste("PPT =", round(zero_values["precipitation"], 3))) +
  theme_classic()

me <- predict_response(dmodel_unit, terms = c("mean_relative_humidity")) #, "month"), type = "random")
e <- plot(me, show_ci = TRUE) +
  # scale_color_manual(values = values_cos) +
  geom_vline(xintercept = zero_values["mean_relative_humidity"], color = "grey33", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey33", alpha = 0.6) +
  labs(title = "", 
       y = TeX("$\\Delta$ Predictions [CSM - CM]"),
       x = "Mean Relative Humidity (%)",
       caption = paste("MRH =", round(zero_values["mean_relative_humidity"], 3))) +
  theme_classic()

a + b + c + d + e +
  plot_annotation(tag_levels = c("A"), tag_suffix = ")") 

ggsave(paste0(loc.fig, "Spatial_clusters/glm/diff_model/decor_month_mar_effect_", fldr, "_", fldr_ma, ".png"),
       width = 35,  height = 22, units = "cm")

# Improving the model ----------------------------------------------------------
dmodel1 <- glmmTMB(diff ~ scale(min_temperature) + scale(mean_temperature) + scale(max_temperature) +
                     scale(precipitation) + scale(mean_relative_humidity) +
                     (1|id) + (1|month), 
                   data = diff_table_wth)
summary(dmodel1)

dmodel2 <- glmmTMB(diff ~ scale(min_temperature) + 
                     scale(precipitation) + scale(mean_relative_humidity) +
                     (1|id) + (1|month), 
                   data = diff_table_wth)
summary(dmodel2)

dmodel3 <- glmmTMB(diff ~ scale(mean_temperature) + 
                     scale(precipitation) + scale(mean_relative_humidity) +
                     (1|id) + (1|month), 
                   data = diff_table_wth)
summary(dmodel3)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = dmodel3, plot =T)
residuals(simulationOutput)

dmodel4 <- glmmTMB(diff ~ scale(max_temperature) + 
                     scale(precipitation) + scale(mean_relative_humidity) +
                     (1|id) + (1|month), 
                   data = diff_table_wth)
summary(dmodel4)

dmodel5 <- glmmTMB(diff ~ scale(mean_temperature) + 
                     scale(mean_relative_humidity) +
                     (1|id) + (1|month), 
                   data = diff_table_wth)
summary(dmodel5)
MuMIn::r.squaredGLMM(dmodel5) 

dmodel6 <- glmmTMB(diff ~ scale(max_temperature) + 
                     scale(mean_relative_humidity) +
                     (1|id) + (1|month), 
                   data = diff_table_wth)
summary(dmodel6)

simulationOutput <- DHARMa::simulateResiduals(fittedModel = dmodel5, plot =T)

dmodel5_unit <- glmmTMB(diff ~ mean_temperature + 
                          mean_relative_humidity +
                          (1|id) + (1|month), 
                        data = diff_table_wth)
summary(dmodel5_unit)

AIC(dmodel1, dmodel2, dmodel3, dmodel4)

options(na.action = "na.fail")
dd <- MuMIn::dredge(dmodel1)

variables <- c("mean_temperature", "mean_relative_humidity")
zero_values <- sapply(variables, function(var) find_zero(var, dmodel5_unit, diff_table_wth))

me <- predict_response(dmodel5_unit, terms = c("mean_temperature")) 
a <- plot(me, show_ci = TRUE) +
  # scale_color_manual(values = values_cos) +
  geom_vline(xintercept = zero_values["mean_temperature"], color = "grey33", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey33", alpha = 0.6) +
  labs(title = "", 
       y = TeX("$\\Delta$ Predictions [CITSI - COUNT]"),
       x = "Mean temperature (ºC)",
       caption = paste("TN =", round(zero_values["mean_temperature"], 3))) +
  theme_classic()

me <- predict_response(dmodel5_unit, terms = c("mean_relative_humidity")) #, "month"), type = "random")
b <- plot(me, show_ci = TRUE) +
  # scale_color_manual(values = values_cos) +
  geom_vline(xintercept = zero_values["mean_relative_humidity"], color = "grey33", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey33", alpha = 0.6) +
  labs(title = "", 
       y = TeX("$\\Delta$ Predictions [CITSI - COUNT]"),
       x = "Mean Relative Humidity (%)",
       caption = paste("MRH =", round(zero_values["mean_relative_humidity"], 3))) +
  theme_classic()
a + b +
  plot_annotation(tag_levels = c("A"), tag_suffix = ")") 

ggsave(paste0(loc.fig, "Spatial_clusters/glm/diff_model/imp_mar_effect_", fldr, "_", fldr_ma, ".png"),
       width = 38,  height = 18, units = "cm")

# Differences by landcover -----------------------------------------------------
clc_surface <- readRDS(file = paste0(loc.output, "clc_surface_mun_level_0.rds")) %>% ungroup()

diff_table_clc <- merge(diff_table, clc_surface, by = c("id", "municipality"), all.x = TRUE)

dmodel1 <- glmmTMB(diff ~ scale(cont_urban_fabric) + 
                     scale(discont_urban_fabric) + scale(forests_scrub) + 
                     scale(green_urban) + scale(inland_water) + 
                     scale(inland_wetlands) + scale(marine_water) + 
                     scale(marine_wetlands) + scale(open) + scale(sports_leisure) + 
                     scale(other_artificial) + scale(roads_rails) + 
                     (1|id) + (1|month), 
                   data = diff_table_clc)
summary(dmodel1)
performance::check_collinearity(dmodel1)

options(na.action = "na.fail")
dd <- MuMIn::dredge(dmodel1)
