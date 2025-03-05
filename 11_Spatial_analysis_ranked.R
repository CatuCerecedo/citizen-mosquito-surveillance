############################## Spatial analysis ################################
library(ggplot2)
library(tidyverse)
library(sf)
library(parallel)
library(spdep)
library(sfdep)
library(rstatix)
library(patchwork)

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

months = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
y = "2020"

for (m in months) {
  # Load files -------------------------------------------------------------------
  spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
    mutate(
      municipality = if_else(is.na(municipality), "no_name", municipality)
    )
  
  # Load the predictionin which you are interested:
  pred <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                         "/tiger_", m ,"_", y, sub, ".rds")) %>%
    janitor::clean_names() %>%
    st_drop_geometry()
  
  pred$pred_count <- rowMeans(pred[,6:ncol(pred)], na.rm = TRUE)
  
  pred <- pred[c("municipality","id", "prov_name", "lon", "lat", "pred_count")]
  # print(summary(pred$pred_count))
  
  pred <- merge(pred, spain, by = c("municipality", "id", "prov_name"))
  st_geometry(pred) <- "geometry"
  
  a <- ggplot() +
    geom_sf(data = pred, aes(fill = pred_count), color = "transparent",
            size = 0.01, alpha = 0.8, na.rm = TRUE) +
    scale_fill_distiller("", palette = "Spectral", na.value = "transparent") + 
    ggtitle(paste(fldr, m, "-", y)) +
    theme_classic()
  
  pred_ma <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr_ma, mdl_name_ma,
                            "/tiger_", m ,"_", y, sub_ma, ".rds")) %>%
    janitor::clean_names() %>%
    st_drop_geometry()
  
  pred_ma$pred_count <- rowMeans(pred_ma[,6:ncol(pred_ma)], na.rm = TRUE)
  
  pred_ma <- pred_ma[c("municipality","id", "prov_name", "lon", "lat", "pred_count")]
  # print(summary(pred$pred_count))
  
  pred_ma <- merge(pred_ma, spain, by = c("municipality", "id", "prov_name"))
  st_geometry(pred_ma) <- "geometry"
  
  b <- ggplot() +
    geom_sf(data = pred_ma, aes(fill = pred_count), color = "transparent",
            size = 0.01, alpha = 0.8, na.rm = TRUE) +
    scale_fill_distiller("", palette = "Spectral", na.value = "transparent") + 
    ggtitle(paste(fldr_ma, m, "-", y)) +
    theme_classic()
  
  ggpubr::ggarrange(a, b)
  
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
    ggtitle(paste(fldr, m, "-", y)) +
    theme_classic()
  
  d <- ggplot() +
    geom_sf(data = pred_ma, aes(fill = pred_count_rnk), color = "transparent",
            size = 0.01, alpha = 0.8, na.rm = TRUE) +
    scale_fill_distiller("", palette = "Spectral", na.value = "transparent") + 
    ggtitle(paste(fldr_ma, m, "-", y)) +
    theme_classic()
  
  ggpubr::ggarrange(a, b, c, d)
  
  ggsave(paste0(loc.fig, "Spatial_clusters/monthly/ranked_", m, "_", y, "_", fldr, "_", fldr_ma, ".png"),
         width = 25,  height = 25, units = "cm", bg = "white")
  
  pred_ma <- pred_ma %>% 
    rename(pred_ma_rnk = pred_count_rnk) %>%
    select(-pred_count)
  pred <- pred %>% 
    select(-pred_count)
  
  diff_pred <- merge(pred, pred_ma %>% st_drop_geometry(),
                     by = c("municipality", "id", "prov_name", "lon", "lat")
  )
  
  diff_pred$pred_rnk_diff <- diff_pred$pred_count_rnk - diff_pred$pred_ma_rnk
  
  diff_rank_plot <- ggplot() +
    geom_sf(data = diff_pred, aes(fill = pred_rnk_diff), color = "transparent",
            size = 0.01, alpha = 0.8, na.rm = TRUE) +
    scale_fill_distiller("", palette = "Spectral", na.value = "transparent") + 
    ggtitle(paste(fldr, "-", fldr_ma)) +
    theme_classic()
  
  ggpubr::ggarrange(c, d, diff_rank_plot, nrow = 1)
  ggsave(paste0(loc.fig, "Spatial_clusters/monthly/diff_ranked_", m, "_", y, "_", fldr, "_", fldr_ma, ".png"),
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
        gi > 0 & p_folded_sim <= 0.01 ~ paste0(fldr, ">", fldr_ma),
        # gi > 0 & p_folded_sim <= 0.1 ~ "Somewhat hot",
        # gi < 0 & p_folded_sim <= 0.01 ~ "Very cold",
        gi < 0 & p_folded_sim <= 0.01 ~ paste0(fldr, "<", fldr_ma),
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
      title = paste0("Predicted Clusters (", fldr, "-", fldr_ma, ")")
    )
  
  ggpubr::ggarrange(diff_rank_plot, cluster_plot)
  
  ggsave(paste0(loc.fig, "Spatial_clusters/monthly/clusters_", m, "_", y, "_", fldr, "_", fldr_ma, ".png"),
         width = 25,  height = 25, units = "cm", bg = "white")
  
  # Analysing the climatic variables----------------------------------------------
  spain_wth <- mclapply(1:nrow(spain), function(i){
    
    cat(paste0("Number row:", i, "\n"))
    
    data_row <- spain[i, ]
    
    data_point <- st_centroid(data_row)
    
    wth <- readRDS(paste0(loc.output, "monthly_weather_data/prep_", 
                          m, "-", y, ".rds"))[, c(2,9:12)]
    
    data_row <- merge(data_row, wth, by = "id", all.x = TRUE)
  }
  , mc.cores = 8)
  spain_wth <- do.call(rbind, spain_wth)
  
  cluster_wth <- merge(cluster %>% st_drop_geometry(), 
                       spain_wth %>% st_drop_geometry(),
                       by = "id")
  
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
  ggsave(paste0(loc.fig, "Spatial_clusters/monthly/cluster_boxplot_", m, "_", y, "_", fldr, "_", fldr_ma, ".png"),
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
  
  ggsave(paste0(loc.fig, "Spatial_clusters/monthly/cluster_density_", m, "_", y, "_", fldr, "_", fldr_ma, ".png"),
         width = 25,  height = 25, units = "cm")
  
  # Bootstrapping ----------------------------------------------------------------
  # cluster_wth <-  readRDS("p.rds")
  boots_func <- function(variables, R = 1000, mean = TRUE) {
    boo <- list()
    j <- 0
    
    for (var in variables) {
      if (mean == TRUE) {
        # Seleccionar la variable y la clasificación
        car_tab <- cluster_wth %>%
          dplyr::select(all_of(var), classification) %>%
          rename(value = all_of(var))
        
        # Estadísticas descriptivas
        summary_stats <- car_tab %>%
          group_by(classification) %>%
          summarise(
            media = mean(value, na.rm = TRUE),
            mediana = median(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            n = n(),
            .groups = "drop"
          )
        
        # Diferencia observada
        mean1 <- car_tab %>% 
          filter(classification == levels(car_tab$classification)[1]) %>% 
          pull(value) %>% 
          mean(na.rm = TRUE)
        mean2 <- car_tab %>% 
          filter(classification == levels(car_tab$classification)[2]) %>% 
          pull(value) %>% 
          mean(na.rm = TRUE)
        dif_obs <- mean1 - mean2
        
        # Bootstrapping
        values1 <- car_tab %>%
          filter(classification == levels(car_tab$classification)[1])
        values2 <- car_tab %>%
          filter(classification == levels(car_tab$classification)[2])
        
        theta <- numeric(R) # Inicializar vector
        for (i in 1:R) { # Bootstrap resampling
          car_tab$group <- sample(car_tab$classification,
                                  length(car_tab$classification),
                                  replace = TRUE
          )
          Mean = tapply(X = car_tab$value, INDEX = car_tab$group, mean)
          theta[i] = as.numeric(Mean[1] - Mean[2])
          # xx1 <- sample(values1$value, length(values1$value), replace = TRUE)
          # xx2 <- sample(values2$value, length(values2$value), replace = TRUE)
          # theta[i] <- mean(xx1) - mean(xx2)
        }
        
        # Cálculo del p-valor
        p_value <- (sum(abs(theta) >= abs(dif_obs)) + 1) / (R + 1)
        # p_value <- mean(abs(theta) >= abs(dif_obs))
        # p_value_up <- sum(theta >= dif_obs) / R 
        # p_value_down <- (sum(theta <= dif_obs) + 1) / (R + 1)
        
        # Resultados
        cat("---------------------\n")
        cat("Summary of", var, "\n")
        cat("Observed estimate:", dif_obs, "\n")
        cat("Bootstrap mean:", mean(theta), "\n")
        cat("Bootstrap mean 95%CIs:", quantile(theta, probs=c(.025,.975)), "\n")
        cat("P-value:", p_value, "\n")
        
        # cat("P-value_up:", p_value_up, "\n")
        # cat("P-value_down:", p_value_down, "\n")
        
        # Visualización
        a <- ggplot(data = car_tab, aes(x = classification, y = value, colour = classification)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(width = 0.25, alpha = 0.06) +
          scale_color_manual(values = c("gray50", "orangered2")) +
          coord_flip() +
          theme_bw() +
          theme(legend.position = "none") 
        b <- ggplot(data = car_tab, aes(x = value, fill = classification)) +
          geom_density(alpha = 0.5) +
          scale_fill_manual(values = c("gray50", "orangered2")) +
          theme_bw() +
          theme(legend.position = "bottom")
        c <- ggplot(data = as.data.frame(theta), aes(x = theta)) +
          geom_histogram(alpha = 0.5) +
          scale_fill_manual(values = c("gray50")) +
          geom_vline(xintercept = dif_obs, linetype="dotted", 
                     color = "orangered2", size=1.5) +
          labs(caption = paste("Observed mean: ", signif(dif_obs, 4),
                               "\n; Bootstrap mean: ", signif(mean(theta), 4),
                               "\n; Bootstrap mean 95%CIs:",  signif(as.numeric(quantile(theta, probs=c(.025,.975))[1]), 4),
                               signif(as.numeric(quantile(theta, probs=c(.025,.975))[2]), 4),
                               "\n; P-value: ", signif(p_value, 4))) +
          theme_bw() +
          theme(legend.position = "bottom")
        
        boo_plot <- (a | c) +
          plot_annotation(title = var, tag_levels = c("A"), tag_suffix = ")") 
      } else {
        # Seleccionar la variable y la clasificación
        car_tab <- cluster_wth %>%
          dplyr::select(all_of(var), classification) %>%
          rename(value = all_of(var))
        
        # Estadísticas descriptivas
        summary_stats <- car_tab %>%
          group_by(classification) %>%
          summarise(
            media = mean(value, na.rm = TRUE),
            mediana = median(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            n = n(),
            .groups = "drop"
          )
        
        # Diferencia observada
        median1 <- car_tab %>% 
          filter(classification == levels(car_tab$classification)[1]) %>% 
          pull(value) %>% 
          median(na.rm = TRUE)
        median2 <- car_tab %>% 
          filter(classification == levels(car_tab$classification)[2]) %>% 
          pull(value) %>% 
          median(na.rm = TRUE)
        dif_obs <- median1 -median2
        
        # Bootstrapping
        values1 <- car_tab %>%
          filter(classification == levels(car_tab$classification)[1])
        values2 <- car_tab %>%
          filter(classification == levels(car_tab$classification)[2])
        
        theta <- numeric(R) # Inicializar vector
        for (i in 1:R) { # Bootstrap resampling
          car_tab$group <- sample(car_tab$classification,
                                  length(car_tab$classification),
                                  replace = TRUE
          )
          Median = tapply(X = car_tab$value, INDEX = car_tab$group, median)
          theta[i] = as.numeric(Median[1] - Median[2])
          # xx1 <- sample(values1$value, length(values1$value), replace = TRUE)
          # xx2 <- sample(values2$value, length(values2$value), replace = TRUE)
          # theta[i] <- mean(xx1) - mean(xx2)
        }
        
        # Cálculo del p-valor
        p_value <- (sum(abs(theta) >= abs(dif_obs)) + 1) / (R + 1)
        # p_value <- mean(abs(theta) >= abs(dif_obs))
        # p_value_up <- sum(theta >= dif_obs) / R 
        # p_value_down <- (sum(theta <= dif_obs) + 1) / (R + 1)
        
        # Resultados
        cat("---------------------\n")
        cat("Summary of", var, "\n")
        cat("Observed estimate:", dif_obs, "\n")
        cat("Bootstrap median:", median(theta), "\n")
        cat("Bootstrap median 95%CIs:", quantile(theta, probs=c(.025,.975)), "\n")
        cat("P-value:", p_value, "\n")
        
        # cat("P-value_up:", p_value_up, "\n")
        # cat("P-value_down:", p_value_down, "\n")
        
        # Visualización
        a <- ggplot(data = car_tab, aes(x = classification, y = value, colour = classification)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(width = 0.25, alpha = 0.06) +
          scale_color_manual(values = c("gray50", "orangered2")) +
          coord_flip() +
          theme_bw() +
          theme(legend.position = "none") 
        b <- ggplot(data = car_tab, aes(x = value, fill = classification)) +
          geom_density(alpha = 0.5) +
          scale_fill_manual(values = c("gray50", "orangered2")) +
          theme_bw() +
          theme(legend.position = "bottom") 
        c <- ggplot(data = as.data.frame(theta), aes(x = theta)) +
          geom_histogram(alpha = 0.5) +
          scale_fill_manual(values = c("gray50")) +
          geom_vline(xintercept = dif_obs, linetype="dotted", 
                     color = "orangered2", size=1.5) +
          labs(caption = paste("Observed median: ", signif(dif_obs, 4),
                               "\n; Bootstrap median: ", signif(median(theta), 4),
                               "\n; Bootstrap median 95%CIs:",  signif(as.numeric(quantile(theta, probs=c(.025,.975))[1]), 4),
                               signif(as.numeric(quantile(theta, probs=c(.025,.975))[2]), 4),
                               "\n; P-value: ", signif(p_value, 4))) +
          theme_bw() +
          theme(legend.position = "bottom")
        
        boo_plot <- a + c +
          plot_annotation(title = var, tag_levels = c("A"), tag_suffix = ")")
      }
      j = j+1
      boo[[j]] <- boo_plot
    }
    return(boo)
  }
  variables <- c("min_temperature", "mean_temperature", "mean_relative_humidity", "precipitation")
  n <- boots_func(variables, R = 10000, mean = F)
  
  ggpubr::ggarrange(n[[1]], n[[2]], n[[3]], n[[4]])
  
  ggsave(paste0(loc.fig, "Spatial_clusters/monthly/boots_", m, "_", y, "_", fldr, "_", fldr_ma, "_", ".png"),
         width = 30,  height = 20, units = "cm")
}

