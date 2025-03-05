############################## Spatial analysis ################################
library(ggplot2)
library(tidyverse)
library(sf)
library(parallel)
library(spdep)
library(sfdep)
library(rstatix)
library(patchwork)
library(overlapping)

rm(list = ls())

# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")
loc.heavy <- paste0(getwd(), "/Heavy_files/")
loc.clc <- paste0(getwd(), "/OUTPUT/land_cover_raster/")

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

decor_time <- list()
decor_time[[1]] <-  c("2020", "2021", "2022")
decor_time[[2]] <- list(list("07", "08"), list("08", "09"), list("07", "08"))

# Load files -------------------------------------------------------------------
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  )

# Load the prediction in which you are interested:
pred <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                           "/tiger_", "07", "_", "2020", ".rds")) %>%
  janitor::clean_names() %>%
  st_drop_geometry()
pred <- pred[,1:5]

for (i in 1:3){
  print(decor_time[[1]][i])
  
  for (j in decor_time[[2]][[i]]){
    print(j)
    
    p <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                           "/tiger_", j, "_", decor_time[[1]][i], ".rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    p$count <- rowMeans(p[,6:ncol(p)], na.rm = TRUE)
    p <- p[c("municipality","id", "prov_name", "lon", "lat", "count")]
    
    pred <- cbind(pred, p["count"])
  }
}

pred$pred_count <- rowMeans(pred[,6:ncol(pred)], na.rm = TRUE)
pred <- pred[c("municipality","id", "prov_name", "lon", "lat", "pred_count")]
pred <- merge(pred, spain, by = c("municipality", "id", "prov_name"))
st_geometry(pred) <- "geometry"

a <- ggplot() +
  geom_sf(data = pred, aes(fill = pred_count), color = "transparent",
          size = 0.01, alpha = 0.8, na.rm = TRUE) +
  scale_fill_distiller("", palette = "Spectral", na.value = "transparent") + 
  ggtitle(paste(fldr, "- Decorrelation")) +
  theme_classic()


pred_ma <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr_ma, mdl_name_ma,
                       "/tiger_", "07", "_", "2020", "_ma.rds")) %>%
  janitor::clean_names() %>%
  st_drop_geometry()
pred_ma <- pred_ma[,1:5]

for (i in 1:3){
  print(decor_time[[1]][i])
  
  for (j in decor_time[[2]][[i]]){
    print(j)
    
    p <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr_ma, mdl_name_ma,
                        "/tiger_", j, "_", decor_time[[1]][i], "_ma.rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    p$count <- rowMeans(p[,6:ncol(p)], na.rm = TRUE)
    p <- p[c("municipality","id", "prov_name", "lon", "lat", "count")]
    
    pred_ma <- cbind(pred_ma, p["count"])
  }
}

pred_ma$pred_count <- rowMeans(pred_ma[,6:ncol(pred_ma)], na.rm = TRUE)
pred_ma <- pred_ma[c("municipality","id", "prov_name", "lon", "lat", "pred_count")]
pred_ma <- merge(pred_ma, spain, by = c("municipality", "id", "prov_name"))
st_geometry(pred_ma) <- "geometry"

b <- ggplot() +
  geom_sf(data = pred_ma, aes(fill = pred_count), color = "transparent",
          size = 0.01, alpha = 0.8, na.rm = TRUE) +
  scale_fill_distiller("", palette = "Spectral", na.value = "transparent") + 
  ggtitle(paste(fldr_ma, "- Decorrelation")) +
  theme_classic()

# Rank ordering convertion -----------------------------------------------------
rank_order_conversion <- function(vctr) {
  rank = rank(vctr, na.last = "keep", ties.method = "first")       
  rank_normalized = rank / max(rank, na.rm = TRUE) 
  return(rank_normalized)
}

pred$pred_count_rnk <- rank_order_conversion(pred$pred_count)
pred_ma$pred_count_rnk <- rank_order_conversion(pred_ma$pred_count)

c <- ggplot() +
  geom_sf(data = pred, aes(fill = pred_count_rnk), color = "transparent",
          size = 0.01, alpha = 0.8, na.rm = TRUE) +
  scale_fill_distiller("", palette = "Spectral", na.value = "transparent") + 
  ggtitle(paste(fldr, "- Decorrelation")) +
  theme_classic()

d <- ggplot() +
  geom_sf(data = pred_ma, aes(fill = pred_count_rnk), color = "transparent",
          size = 0.01, alpha = 0.8, na.rm = TRUE) +
  scale_fill_distiller("", palette = "Spectral", na.value = "transparent") + 
  ggtitle(paste(fldr_ma, "- Decorrelation")) +
  theme_classic()

ggpubr::ggarrange(a, b, c, d)

ggsave(paste0(loc.fig, "Spatial_clusters/period_ranked_", fldr, "_", fldr_ma, "_", period, ".png"),
       width = 25,  height = 25, units = "cm", bg = "white")

pred_ma <- pred_ma %>% 
  rename(pred_ma_rnk = pred_count_rnk) %>%
  select(-pred_count)
pred <- pred %>% 
  select(-pred_count)

diff_pred <- merge(pred, pred_ma %>% st_drop_geometry(),
                   by = c("municipality", "id", "prov_name", "lon", "lat")
                   )

diff_pred$pred_rnk_diff <- diff_pred$pred_ma_rnk - diff_pred$pred_count_rnk

diff_rank_plot <- ggplot() +
  geom_sf(data = diff_pred, aes(fill = pred_rnk_diff), color = "transparent",
          size = 0.01, alpha = 0.8, na.rm = TRUE) +
  scale_fill_distiller("", palette = "Spectral", na.value = "transparent", limits = c(-1, 1)) + 
  ggtitle(paste(fldr_ma, "-", fldr)) +
  theme_classic()

ggpubr::ggarrange(c, d, diff_rank_plot, nrow = 1)
ggsave(paste0(loc.fig, "Spatial_clusters/period_diff_ranked_", fldr, "_", fldr_ma, "_", period, ".png"),
       width = 25,  height = 25, units = "cm", bg = "white")

rm(pred, pred_ma)

# Spatial autocorrelation analysis ---------------------------------------------
sf::sf_use_s2(FALSE)

# Create a list of neighbors for each polygon
list_nb <- poly2nb(diff_pred, queen = TRUE) 

# he relationship will be based on a queen contiguity criterion, which considers 
# polygons to be neighbors if they share a boundary or a vertex.

# Check for empty neighbor sets
# card() calculates number of neighbors for each polygon in the list
# which() finds polygons with 0 neighbors
empty_nb <- which(card(list_nb) == 0)

# Remove polygons with empty neighbor sets from the data
pred_subset <- diff_pred[-empty_nb, ]# Remove polygons with empty neighbor sets from the data

# # Subset 'tes_data' to extract polygons with empty neighbor sets
# empty_polygons <- pred[empty_nb, ]
# empty_polygons$municipality  # print municipality names

## Global G test ---------------------------------------------------------------
# Now that we removed empty neighbor sets (tes_subset)
# Identify neighbors with queen contiguity (edge/vertex touching)
tes_nb <- poly2nb(pred_subset, queen = TRUE)

# Binary weighting assigns a weight of 1 to all neighboring features 
# and a weight of 0 to all other features
pred_subset$area <- st_area(pred_subset)/1000000

wlist <- lapply(seq_along(tes_nb), function(i) {
  nbhs <- tes_nb[[i]]
  as.numeric(pred_subset$area[nbhs])
})
tes_w_binary <- nb2listw(tes_nb, glist = wlist, style="W", zero.policy=T)

# Calculate spatial lag of predicted values
tes_lag <- lag.listw(tes_w_binary, (pred_subset$pred_rnk_diff + 1))

# Test for global G statistic of predicted values
globalG.test((pred_subset$pred_rnk_diff + 1), tes_w_binary)

## Local G test ----------------------------------------------------------------
# Identify neighbors, create weights, calculate spatial lag
tes_nbs <- pred_subset %>% 
  mutate(
    nb = st_contiguity(geometry),        # neighbors share border/vertex
    wt = st_weights(nb),                 # row-standardized weights
    tes_lag = st_lag((pred_rnk_diff + 1), nb, wt)    # calculate spatial lag of TreEqty
  ) 

# Calculate the Gi using local_g_perm
tes_hot_spots <- tes_nbs %>%
  mutate(
    Gi = local_g_perm((pred_rnk_diff + 1), nb, wt, nsim = 999)
  ) %>%
  unnest(Gi) 

tes_hot_spots %>% 
  ggplot((aes(fill = gi))) +
  geom_sf(color = "black", lwd = 0.15) +
  scale_fill_gradient2() +
  theme_classic()

# Create a new data frame called 'tes_hot_spots"
cluster <- tes_hot_spots %>%
  select(gi, p_folded_sim, id) %>% 
  mutate(
    classification = case_when(
      # Classify based on the following p-value criteria:
      # gi > 0 & p_folded_sim <= 0.01 ~ "Very hot",
      gi > 0 & p_folded_sim <= 0.05 ~ paste0(fldr, ">", fldr_ma),
      # gi > 0 & p_folded_sim <= 0.1 ~ "Somewhat hot",
      # gi < 0 & p_folded_sim <= 0.01 ~ "Very cold",
      gi < 0 & p_folded_sim <= 0.05 ~ paste0(fldr, "<", fldr_ma),
      # gi < 0 & p_folded_sim <= 0.01 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    classification = factor(
      classification,
      levels = c(
                 # "Very hot",
                 paste0(fldr, ">", fldr_ma), 
                 # "Somewhat hot",
                 "Insignificant",
                 # "Somewhat cold",
                 paste0(fldr, "<", fldr_ma)
                 # "Very cold"
                 )
    )
  ) 


cluster_plot <-ggplot(cluster, aes(fill = classification)) +
  geom_sf(color = "black", lwd = 0.1) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Cluster Classification",
    title = paste0("Predicted Clusters (", fldr_ma, "-", fldr, ")")
  )

ggpubr::ggarrange(diff_rank_plot, cluster_plot)

ggsave(paste0(loc.fig, "Spatial_clusters/period_clusters_", fldr, "_", fldr_ma, "_", period, ".png"),
       width = 25,  height = 25, units = "cm", bg = "white")

# Analyzing the climatic variables----------------------------------------------
wth_table <- data.frame()
for (i in 1:3){
  print(decor_time[[1]][i])
  
  for (j in decor_time[[2]][[i]]){
    print(j)
    
    wth <- readRDS(paste0(loc.output, "monthly_weather_data/prep_", 
                          j, "-", decor_time[[1]][i], ".rds"))[, c(1, 9:12, 32)]
    wth$year <- decor_time[[1]][i]
    wth$month <- j
    
    wth_table <- rbind(wth_table, wth)
  }
}

wth_table <- wth_table %>% 
  group_by(id) %>%
  summarise(
    min_temperature = mean(min_temperature, na.rm = TRUE),
    max_temperature = mean(max_temperature, na.rm = TRUE),
    mean_temperature = mean(mean_temperature, na.rm = TRUE),
    precipitation = sum(precipitation, na.rm = TRUE),
    mean_relative_humidity = mean(mean_relative_humidity, na.rm = TRUE)
  )

cluster_wth <- merge(cluster %>% st_drop_geometry(), 
                     wth_table,
           by = "id")

saveRDS(cluster_wth, paste0(loc.output, "cluster_", fldr, "_", fldr_ma, "_", 
                            period, ".rds"))

# There are a lot of outliers
remove_outliers <- function(x, na.rm = TRUE) {
  qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  outliers <- which(x < (qnt[1] - H) | x > (qnt[2] + H))
  return(outliers)
}
out <- remove_outliers(cluster_wth$precipitation)

cluster_wth <- cluster_wth  %>%
  mutate(classification = factor(classification, levels = c(paste0(fldr, ">", fldr_ma),  paste0(fldr, "<", fldr_ma)))) %>%
  filter(classification != "Insignificant")

pwc_min <- cluster_wth %>% 
  dunn_test(min_temperature ~ classification, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "classification")
k_min <- cluster_wth %>% kruskal_test(min_temperature ~ classification)
fligner.test(min_temperature ~ classification, data = cluster_wth)

pwc_mean <-cluster_wth %>% 
  dunn_test(mean_temperature ~ classification, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "classification")
k_mean <-cluster_wth %>% kruskal_test(mean_temperature ~ classification)
fligner.test(mean_temperature ~ classification, data = cluster_wth)

pwc_prec <- cluster_wth %>% 
  dunn_test(precipitation ~ classification, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "classification")
k_prec <- cluster_wth %>% kruskal_test(precipitation ~ classification)
fligner.test(precipitation ~ classification, data = cluster_wth)

pwc_hum <-cluster_wth %>% 
  dunn_test(mean_relative_humidity ~ classification, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "classification")
k_hum <-cluster_wth %>% kruskal_test(mean_relative_humidity ~ classification)
fligner.test(mean_relative_humidity ~ classification, data = cluster_wth)

ggpubr::ggarrange(
  ggplot(cluster_wth) +
    geom_boxplot(aes(x = classification, y = min_temperature)) +
    ggpubr::stat_pvalue_manual(pwc_min, hide.ns = TRUE) +
    labs(
      subtitle = get_test_label(k_min, detailed = TRUE),
      caption = get_pwc_label(pwc_min)
    ) +
    theme_classic(),
  ggplot(cluster_wth) +
    geom_boxplot(aes(x = classification, y = mean_temperature)) +
    ggpubr::stat_pvalue_manual(pwc_mean, hide.ns = TRUE) +
    labs(
      subtitle = get_test_label(k_mean, detailed = TRUE),
      caption = get_pwc_label(pwc_mean)
    ) +
    theme_classic(),
  ggplot(cluster_wth) +
    geom_boxplot(aes(x = classification, y = precipitation)) +
    ggpubr::stat_pvalue_manual(pwc_prec, hide.ns = TRUE) +
    labs(
      subtitle = get_test_label(k_prec, detailed = TRUE),
      caption = get_pwc_label(pwc_prec)
    ) +
    theme_classic(),
  ggplot(cluster_wth) +
    geom_boxplot(aes(x = classification, y = mean_relative_humidity)) +
    ggpubr::stat_pvalue_manual(pwc_hum, hide.ns = TRUE) +
    labs(
      subtitle = get_test_label(k_hum, detailed = TRUE),
      caption = get_pwc_label(pwc_hum)
    ) +
    theme_classic()
)
ggsave(paste0(loc.fig, "Spatial_clusters/period_cluster_boxplot_", fldr, "_", fldr_ma, "_", period, ".png"),
       width = 25,  height = 25, units = "cm")

ggpubr::ggarrange(
  ggplot(cluster_wth) +
    geom_density(aes(x = min_temperature, fill = classification), alpha = 0.8) +
    theme_classic(),
  ggplot(cluster_wth) +
    geom_density(aes(x = mean_temperature, fill = classification), alpha = 0.8) +
    theme_classic(),
  ggplot(cluster_wth) +
    geom_density(aes(x = precipitation, fill = classification), alpha = 0.8) +
    theme_classic(),
  ggplot(cluster_wth) +
    geom_density(aes(x =  mean_relative_humidity, fill = classification), alpha = 0.8) +
    theme_classic()
)

ggsave(paste0(loc.fig, "Spatial_clusters/period_cluster_density_", fldr, "_", fldr_ma, "_", period, ".png"),
       width = 25,  height = 25, units = "cm")

# GLM: logistic regression -----------------------------------------------------
cluster_wth <- readRDS(paste0(loc.output, "cluster_", fldr, "_", fldr_ma, "_", period, ".rds")) %>%
  mutate(classification = factor(classification, levels = c(paste0(fldr, ">", fldr_ma),  paste0(fldr, "<", fldr_ma)))) %>%
  filter(classification != "Insignificant")

logistic_model <- glm(classification ~ scale(min_temperature) +  scale(max_temperature) +
                        scale(mean_temperature) + scale(mean_relative_humidity) + 
                        scale(precipitation), 
                      data = cluster_wth, 
                      family = binomial(link = "logit"))
summary(logistic_model)

l_0 = levels(logistic_model$data$classification)[1]
l_1 = levels(logistic_model$data$classification)[2]

boo <- list()
j <- 0

cat("min_temperature -----\n")
var = "min_temperature"
new_data <- data.frame(
  min_temperature = seq(min(cluster_wth$min_temperature), 
                        max(cluster_wth$min_temperature), length.out = 1000),
  max_temperature = mean(cluster_wth$max_temperature),
  mean_temperature = mean(cluster_wth$mean_temperature), 
  mean_relative_humidity = mean(cluster_wth$mean_relative_humidity),
  precipitation = mean(cluster_wth$precipitation)
)
new_data$response_prob <- predict(logistic_model, newdata = new_data, type = "response")

car_tab <- cluster_wth %>%
  dplyr::select(all_of(var), classification) %>%
  rename(value = all_of(var))
a <- ggplot(data = car_tab, aes(x = classification, y = value, colour = classification)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, alpha = 0.06) +
  scale_color_manual(values = c("gray50", "orangered2")) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") 
b <- ggplot(new_data, aes(x = min_temperature, y = response_prob)) +
  geom_line() +
  labs(x = var, y = paste("Probability of", l_1),
       caption = paste0("Estimate (SE): ", round(logistic_model$coefficients[paste0("scale(", var, ")")], 3), 
                        " (", round(summary(logistic_model)$coefficients[8], 3), ")\n",
                        "P-value: ", round(summary(logistic_model)$coefficients[20], 5))
  ) +
  theme_classic()

boo_plot <- (a | b) +
  plot_annotation(title = var, tag_levels = c("A"), tag_suffix = ")") 
j = j+1
boo[[j]] <- boo_plot

cat("max_temperature -----\n")
var = "max_temperature"
new_data <- data.frame(
  max_temperature = seq(min(cluster_wth$max_temperature), 
                        max(cluster_wth$max_temperature), length.out = 1000),
  min_temperature = mean(cluster_wth$min_temperature),
  mean_temperature = mean(cluster_wth$mean_temperature), 
  mean_relative_humidity = mean(cluster_wth$mean_relative_humidity),
  precipitation = mean(cluster_wth$precipitation)
)
new_data$response_prob <- predict(logistic_model, newdata = new_data, type = "response")

car_tab <- cluster_wth %>%
  dplyr::select(all_of(var), classification) %>%
  rename(value = all_of(var))
a <- ggplot(data = car_tab, aes(x = classification, y = value, colour = classification)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, alpha = 0.06) +
  scale_color_manual(values = c("gray50", "orangered2")) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") 
b <- ggplot(new_data, aes(x = max_temperature, y = response_prob)) +
  geom_line() +
  labs(x = var, y = paste("Probability of", l_1),
       caption = paste0("Estimate (SE): ", round(logistic_model$coefficients[paste0("scale(", var, ")")], 3), 
                        " (", round(summary(logistic_model)$coefficients[9], 3), ")\n",
                        "P-value: ", round(summary(logistic_model)$coefficients[21], 5))
  ) +
  theme_classic()

boo_plot <- (a | b) +
  plot_annotation(title = var, tag_levels = c("A"), tag_suffix = ")") 
j = j+1
boo[[j]] <- boo_plot

cat("mean_temperature -----\n")
var = "mean_temperature"
new_data <- data.frame(
  mean_temperature = seq(min(cluster_wth$mean_temperature), 
                         max(cluster_wth$mean_temperature), length.out = 1000),
  min_temperature = mean(cluster_wth$min_temperature), 
  max_temperature = mean(cluster_wth$max_temperature),
  mean_relative_humidity = mean(cluster_wth$mean_relative_humidity),
  precipitation = mean(cluster_wth$precipitation)
)
new_data$response_prob <- predict(logistic_model, newdata = new_data, type = "response")

car_tab <- cluster_wth %>%
  dplyr::select(all_of(var), classification) %>%
  rename(value = all_of(var))
a <- ggplot(data = car_tab, aes(x = classification, y = value, colour = classification)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, alpha = 0.06) +
  scale_color_manual(values = c("gray50", "orangered2")) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") 
b <- ggplot(new_data, aes(x = mean_temperature, y = response_prob)) +
  geom_line() +
  labs(x = var, y = paste("Probability of", l_1),
       caption = paste0("Estimate (SE): ", round(logistic_model$coefficients[paste0("scale(", var, ")")], 3), 
                        " (", round(summary(logistic_model)$coefficients[10], 3), ")\n",
                        "P-value: ", round(summary(logistic_model)$coefficients[22], 5))) +
  theme_classic()

boo_plot <- (a | b) +
  plot_annotation(title = var, tag_levels = c("A"), tag_suffix = ")") 
j = j+1
boo[[j]] <- boo_plot

cat("mean_relative_humidity -----\n")
var = "mean_relative_humidity"
new_data <- data.frame(
  mean_relative_humidity = seq(min(cluster_wth$mean_relative_humidity), 
                               max(cluster_wth$mean_relative_humidity), length.out = 1000),
  min_temperature = mean(cluster_wth$min_temperature), 
  max_temperature = mean(cluster_wth$max_temperature),
  mean_temperature = mean(cluster_wth$mean_temperature),
  precipitation = mean(cluster_wth$precipitation)
)
new_data$response_prob <- predict(logistic_model, newdata = new_data, type = "response")

car_tab <- cluster_wth %>%
  dplyr::select(all_of(var), classification) %>%
  rename(value = all_of(var))
a <- ggplot(data = car_tab, aes(x = classification, y = value, colour = classification)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, alpha = 0.06) +
  scale_color_manual(values = c("gray50", "orangered2")) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") 
b <- ggplot(new_data, aes(x = mean_relative_humidity, y = response_prob)) +
  geom_line() +
  labs(x = var, y = paste("Probability of", l_1),
       caption = paste0("Estimate (SE): ", round(logistic_model$coefficients[paste0("scale(", var, ")")], 3), 
                        " (", round(summary(logistic_model)$coefficients[11], 3), ")\n",
                        "P-value: ", round(summary(logistic_model)$coefficients[23], 5))
  ) +
  theme_classic()

boo_plot <- (a | b) +
  plot_annotation(title = var, tag_levels = c("A"), tag_suffix = ")") 
j = j+1
boo[[j]] <- boo_plot

cat("precipitation -----\n")
var = "precipitation"
new_data <- data.frame(
  precipitation = seq(min(cluster_wth$precipitation), 
                      max(cluster_wth$precipitation), length.out = 1000),
  min_temperature = mean(cluster_wth$min_temperature), 
  mean_temperature = mean(cluster_wth$mean_temperature),
  max_temperature = mean(cluster_wth$max_temperature),
  mean_relative_humidity = mean(cluster_wth$mean_relative_humidity)
)
new_data$response_prob <- predict(logistic_model, newdata = new_data, type = "response")

car_tab <- cluster_wth %>%
  dplyr::select(all_of(var), classification) %>%
  rename(value = all_of(var))
a <- ggplot(data = car_tab, aes(x = classification, y = value, colour = classification)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, alpha = 0.06) +
  scale_color_manual(values = c("gray50", "orangered2")) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") 
b <- ggplot(new_data, aes(x = precipitation, y = response_prob)) +
  geom_line() +
  labs(x = var, y = paste("Probability of", l_1),
       caption = paste0("Estimate (SE): ", round(logistic_model$coefficients[paste0("scale(", var, ")")], 3), 
                        " (", round(summary(logistic_model)$coefficients[12], 3), ")\n",
                        "P-value: ", round(summary(logistic_model)$coefficients[24], 5))
  ) +
  theme_classic()

boo_plot <- (a | b) +
  plot_annotation(title = var, tag_levels = c("A"), tag_suffix = ")") 
j = j+1
boo[[j]] <- boo_plot

ggpubr::ggarrange(boo[[1]], boo[[2]], boo[[3]], boo[[4]], boo[[5]])

ggsave(paste0(loc.fig, "Spatial_clusters/glm/period_boots_", fldr, "_", fldr_ma, "_", period, ".png"),
       width = 35,  height = 22, units = "cm")




