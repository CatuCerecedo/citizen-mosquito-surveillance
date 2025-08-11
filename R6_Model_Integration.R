########################## Integration MA on BG model ##########################

#' Devoloping the CITRAP model (using as only predictor CITSCI model outputs)
#' Cross-validation of CITRAP model

################################################################################
library(tidyverse)
library(sf)
library(parallel)
library(dplyr)
library(lubridate)
library(data.table)
library(brms)
library(janitor)
library(cmdstanr)
library(rstanarm)
library(caret)

rm(list = ls())
# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")

# Model names ------------------------------------------------------------------
mdl_ma <- "mtiger7_citsci"
mdl_name_ma <- "/mtiger7_citsci"
fldr_ma <- "COTSCI"
sub_ma <- "_citsci" 

sf::sf_use_s2(FALSE)

# Integrating data -------------------------------------------------------------
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  )

tiger <- readRDS(paste0(loc.data, "trap_surveillance_data.rds")) 


# Creating a monthly average of MA:
months = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
years = c("2020", "2021", "2022")

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
      
      n_col <- names(pred[6:ncol(pred)])
      pred <- pred %>%
        pivot_longer(
        cols = n_col,
        names_to = "date",
        values_to = "ma"
      )
      
      # pred$count <- rowMeans(pred[,6:ncol(pred)], na.rm = TRUE)
      # pred <- pred[c("municipality","id", "prov_name", "lon", "lat", "count")]
    } else {
      p <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr_ma, mdl_name_ma,
                          "/tiger_", months[i], "_", y, sub_ma, ".rds")) %>%
        janitor::clean_names() %>%
        st_drop_geometry()
      
      n_col <- names(p[6:ncol(p)])
      p <- p %>%
        pivot_longer(
          cols = n_col,
          names_to = "date",
          values_to = "ma"
        )
      
       pred <- rbind(pred, p)
    }
  }

pred$date <- as.Date(gsub("x", "", pred$date), format="%Y_%m_%d")

sum_ma <- vector()
prod_ma <- vector()
mean_ma <- vector()
for (i in 1:nrow(tiger)){
  if (i %% 100 == 0) {print(i)}
  sum_ma[i] = pred %>% 
    filter(id == tiger$id[i]) %>%
    filter(date >= tiger$start_date[i] & date <= tiger$end_date[i]) %>%
    pull(ma) %>% 
    sum()
  prod_ma[i] = pred %>% 
    filter(id == tiger$id[i]) %>%
    filter(date >= tiger$start_date[i] & date <= tiger$end_date[i]) %>%
    pull(ma) %>% 
    prod()
  mean_ma[i] = pred %>% 
    filter(id == tiger$id[i]) %>%
    filter(date >= tiger$start_date[i] & date <= tiger$end_date[i]) %>%
    pull(ma) %>% 
    mean()
}

integrated_data <- tiger %>%
  mutate(
    ma = sum_ma
  )

saveRDS(integrated_data, paste0(loc.output, "integrated_data_summa.rds"))

integrated_data <- tiger %>%
  mutate(
    ma = prod_ma
  )

saveRDS(integrated_data, paste0(loc.output, "integrated_data_prodma.rds"))

integrated_data <- tiger %>%
  mutate(
    ma = mean_ma
  )

saveRDS(integrated_data, paste0(loc.output, "integrated_data_meanma.rds"))
# Modeling process -------------------------------------------------------------
nchains = 4
threads_per_chain = 1
iteret = 5000
wup = 2000

mtiger_inte <- brm(females ~ ma + 
                  offset(log(trapping_effort)) +
                  (1 | id) + (1 | year),
                data = integrated_data,
                prior = set_prior("cauchy(0,2.5)", class="b"),
                family = negbinomial(link = "log"),
                iter = iteret,
                chains = nchains,
                cores = nchains,
                backend = "cmdstanr",
                save_pars = save_pars(all = TRUE),
                threads = threading(threads_per_chain),
                control = list(adapt_delta = 0.999))
summary(mtiger_inte)
loo(mtiger_inte)
loo(mtiger_inte, moment_match = TRUE, recompile = TRUE, reloo = TRUE)
bayes_R2(mtiger_inte)
saveRDS(mtiger_inte, file = paste0(loc.output, "mtiger_inte_only_meanma.rds"))

# Predictions ------------------------------------------------------------------
mtiger_inte <- readRDS(file = paste0(loc.output, "mtiger_inte_only_meanma.rds"))
ncores = 8

xy <- st_centroid(spain) %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry() %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
new_points <- xy

fldr = "Monthly_Integration/"
mdl_name = "tiger_inter_monthly_meanma"
sub = ""

# Loading model
model <- mtiger_inte

# Preparing ERA5 folder

for (year in c("2020", "2021", "2022")){
  for (m in c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")){
    
    # Empty box
    pred_points <- xy %>% st_drop_geometry()
    pred_points_sd <- xy %>% st_drop_geometry()
    
    print(paste("month: ", m, "year: ", year))
    
    # Calculate the first date of the month
    first_date <- lubridate::as_date(paste(year, m, "01", sep = "-"))
    
    # Calculate the last date of the month
    last_day <- ceiling_date(lubridate::as_date(paste(year, m, "01", sep = "-")), "month") - 1
    
    # Function to extract daily data -----------------------------------------------
    for (sel_date in seq.Date(first_date, last_day, 1)){
      sel_date <- as_date(sel_date)
      print(sel_date)
      
      data_prep_day <- readRDS(paste0(loc.output, "daily_weather_data/prep_",
                                      sel_date, ".rds")) %>%
        mutate(
          month = month(date),
          year = year(date)
        )
      
      pred_all_ma <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr_ma, mdl_name_ma,
                             "/tiger_", m, "_", year, sub_ma, ".rds")) %>%
        janitor::clean_names() %>%
        st_drop_geometry() %>%
        dplyr::select(municipality, id, prov_name, lon, lat, janitor::make_clean_names(sel_date))
      
      data_prep_day <- merge(data_prep_day, pred_all_ma,
        by = c("municipality", "id", "prov_name", "lon", "lat"),
        all.x = TRUE
      ) %>%
        rename("ma" = janitor::make_clean_names(sel_date))
      
      print("Calculating variables done")
      
      # Function to predict daily data -----------------------------------------------
      
      data_prep_day <- data_prep_day %>% 
        mutate(
          n_traps = 1,
          trapping_effort = 7
        )
      
      nrow_these_pred_points = nrow(data_prep_day)
      max_chunksize = 300000
      chunksize = min(as.integer((nrow_these_pred_points/ncores)), max_chunksize)
      
      pred <- bind_rows(mclapply(seq(1, nrow(data_prep_day), chunksize), function(i){
        print(i)
        
        data_chunk = data_prep_day[i:min(nrow(data_prep_day), (i+(chunksize-1))), ]
        flush.console()
        
        pp <- apply(posterior_predict(model,
                                      newdata = data_chunk,
                                      allow_new_levels = TRUE,
                                      re_formula =  NA,
                                      ndraws = 1000),
                    2, function(x) mean(x)) %>% 
          as.data.frame()
        
        colnames(pp) <- as.Date(sel_date)
        
        data_chunk <- data_chunk %>% dplyr::select(municipality, id, prov_name, lon, lat)
        
        pp <- bind_cols(data_chunk, pp)
        
        return(pp)
        
      }, mc.cores = ncores))
      
      pred_points <- merge(pred_points, pred,  by = c("municipality", "id", "prov_name", "lon", "lat"))
      
      pred_sd <- bind_rows(mclapply(seq(1, nrow(data_prep_day), chunksize), function(i){
        print(i)
        
        data_chunk = data_prep_day[i:min(nrow(data_prep_day), (i+(chunksize-1))), ]
        flush.console()
        
        pp <- apply(posterior_predict(model, 
                                      newdata = data_chunk,
                                      allow_new_levels = TRUE, 
                                      re_formula = NA,
                                      ndraws = 1000), 
                    2, function(x) sd(x)) %>% 
          as.data.frame()
        
        colnames(pp) <- as.Date(sel_date)
        
        data_chunk <- data_chunk %>% dplyr::select(municipality, id, prov_name, lon, lat)
        
        pp <- bind_cols(data_chunk, pp)
        
        return(pp)
        
      }, mc.cores = ncores))
      
      pred_points_sd <- merge(pred_points_sd, pred_sd, by = c("municipality", "id", "prov_name", "lon", "lat"))
      
      
    } # Days
    
    print(paste("Saving pred: ", m))
    saveRDS(pred_points %>% st_drop_geometry(), file = paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                                                              "/tiger_", m, "_", year,  sub, ".rds"))
    
    print(paste("Saving pred sd: ", m))
    saveRDS(pred_points_sd %>% st_drop_geometry(), file = paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                                                                 "/tiger_", m, "_", year, sub, "_sd.rds"))
  } # months
} # years

## Cross-validation ------------------------------------------------------------
# Split the data: train/test 
set.seed(121234)
sample <- sample(c(TRUE, FALSE), nrow(integrated_data), replace=TRUE, prob=c(0.7, 0.3))
train  <- integrated_data[sample, ]
test   <- integrated_data[!sample, ]

nchains = 4
threads_per_chain = 1
iteret = 5000
wup = 2000

mtiger_inte <- brm(females ~ ma + 
                           offset(log(trapping_effort)) +
                           (1 | id) + (1 | year),
                         data = train,
                         prior = set_prior("cauchy(0,2.5)", class="b"),
                         family = negbinomial(link = "log"),
                         iter = iteret,
                         chains = nchains,
                         cores = nchains,
                         backend = "cmdstanr",
                         save_pars = save_pars(all = TRUE),
                         threads = threading(threads_per_chain),
                         control = list(adapt_delta = 0.999))
summary(mtiger_inte)
loo(mtiger_inte_month, moment_match = TRUE, recompile = TRUE, reloo = TRUE)
bayes_R2(mtiger_inte)

# Predictions: 
pp <- apply(posterior_predict(mtiger_inte, 
                              newdata = test,
                              allow_new_levels = TRUE, 
                              # re_formula = NA,
                              ndraws = 1000), 
            2, function(x) mean(x)) %>%
  as.data.frame()

test$predicted <- pp$.

ggplot(test, aes(x = predicted, y = females)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() 
summary(lm(females ~ predicted, data = test))
summary(lm(I(females - 1) ~ predicted - 1, data = test))

library(caret)
cv_metrics <- data.frame(Scor = cor(test$predicted, test$females, method = "spearman"),
                         Pcor = cor(test$predicted, test$females, method = "pearson"),
                         R2 = R2(test$predicted, test$females),
                         RMSE = RMSE(test$predicted, test$females),
                         MAE = MAE(test$predicted, test$females),
                         model = "INTEGRATED model")
cv_metrics

# Working at monthly scale -----------------------------------------------------
## Integrating data ------------------------------------------------------------
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  )

# Creating a monthly average of MA:
months = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
years = c("2020", "2021", "2022")

tiger_integrated <- data.frame()
for(y in years){
  print(y)
  
  for(i in 1:length(months)){
    print(i)
    if (i == 1){
      pred <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr_ma, mdl_name_ma,
                             "/tiger_", months[i], "_", y, sub_ma, ".rds")) %>%
        janitor::clean_names() %>%
        st_drop_geometry()
      
      pred$count <- rowSums(pred[,6:ncol(pred)], na.rm = TRUE)
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
      values_to = "ma"
    )
  pred$year <- y
  pred$month <- as.numeric(pred$month)
  tiger_integrated <- rbind(tiger_integrated, pred)
} 

tiger_integrated <- readRDS(paste0(loc.output, "tiger_integrated_mean_monthly.rds"))
## Modeling process ------------------------------------------------------------
tiger_integrated <- readRDS(paste0(loc.output, "tiger_integrated_mean_monthly.rds"))

nchains = 4
threads_per_chain = 1
iteret = 5000
wup = 2000

mtiger_inte_month <- brm(females ~ ma + 
                     offset(log(n_traps)) + offset(log(n_days)) + 
                     (1 | id) + (1 | year),
                   data = tiger_integrated,
                   prior = set_prior("cauchy(0,2.5)", class="b"),
                   family = negbinomial(link = "log"),
                   iter = iteret,
                   chains = nchains,
                   cores = nchains,
                   backend = "cmdstanr",
                   save_pars = save_pars(all = TRUE),
                   threads = threading(threads_per_chain),
                   control = list(adapt_delta = 0.999))
summary(mtiger_inte_month)
plot(mtiger_inte_month)
loo(mtiger_inte_month, moment_match = TRUE, recompile = TRUE, reloo = TRUE)
bayes_R2(mtiger_inte_month)
saveRDS(mtiger_inte_month, file = paste0(loc.output, "mtiger_inte_monthly_mean_ma.rds"))

## Monthly predictions ---------------------------------------------------------
mtiger_inte <- readRDS(file = paste0(loc.output, "mtiger_inte_monthly_mean_ma.rds"))
ncores = 8

xy <- st_centroid(spain) %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry() %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
new_points <- xy

fldr = "Monthly_Integration/"
mdl_name = "tiger_inter_monthly_meanma"
sub = ""

# Loading model
model <- mtiger_inte

# Preparing ERA5 folder

for (y in c("2020", "2021", "2022")){
  for (m in c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")){
    
    # Empty box
    pred_points <- xy %>% st_drop_geometry()
    pred_points_sd <- xy %>% st_drop_geometry()
    
    print(paste("month: ", m, "year: ", y))
    
    pred_all_ma <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr_ma, mdl_name_ma,
                                  "/tiger_", m, "_", y, sub_ma, ".rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    pred_all_ma$ma <- rowMeans(pred_all_ma[6:ncol(pred_all_ma)])
    pred_all_ma <- pred_all_ma %>% dplyr::select(id, ma)
    
    data_prep_day <- merge(new_points, pred_all_ma, 
                           by = "id", all.x = TRUE)
    
    print("Calculating variables done")
    
    # Function to predict monthly data 
    
    data_prep_day <- data_prep_day %>% 
      mutate(
        attractor = "CO2",
        n_traps = 1,
        n_days = 7
      )
    
    nrow_these_pred_points = nrow(data_prep_day)
    max_chunksize = 300000
    chunksize = min(as.integer((nrow_these_pred_points/ncores)), max_chunksize)
    
    pred_points <- bind_rows(mclapply(seq(1, nrow(data_prep_day), chunksize), function(i){
      print(i)
      
      data_chunk = data_prep_day[i:min(nrow(data_prep_day), (i+(chunksize-1))), ]
      flush.console()
      
      pp <- apply(posterior_predict(model,
                                    newdata = data_chunk,
                                    allow_new_levels = TRUE,
                                    re_formula =  NA,
                                    ndraws = 1000),
                  2, function(x) mean(x)) %>% 
        as.data.frame()
      
      colnames(pp) <- "pred"
      
      data_chunk <- data_chunk %>% dplyr::select(municipality, id, prov_name, lon, lat)
      
      pp <- bind_cols(data_chunk, pp)
      
      return(pp)
      
    }, mc.cores = ncores))
    
    pred_points_sd <- bind_rows(mclapply(seq(1, nrow(data_prep_day), chunksize), function(i){
      print(i)
      
      data_chunk = data_prep_day[i:min(nrow(data_prep_day), (i+(chunksize-1))), ]
      flush.console()
      
      pp <- apply(posterior_predict(model, 
                                    newdata = data_chunk,
                                    allow_new_levels = TRUE, 
                                    re_formula = NA,
                                    ndraws = 1000), 
                  2, function(x) sd(x)) %>% 
        as.data.frame()
      
      colnames(pp) <- "pred"
      
      data_chunk <- data_chunk %>% dplyr::select(municipality, id, prov_name, lon, lat)
      
      pp <- bind_cols(data_chunk, pp)
      
      return(pp)
      
    }, mc.cores = ncores))
    
    print(paste("Saving pred: ", m))
    saveRDS(pred_points %>% st_drop_geometry(), file = paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                                                              "/tiger_", m, "_", y,  sub, "5_traps_7_days.rds"))
    
    print(paste("Saving pred sd: ", m))
    saveRDS(pred_points_sd %>% st_drop_geometry(), file = paste0(loc.output, "PREDICTIONS/", fldr, mdl_name,
                                                                 "/tiger_", m, "_", y, sub, "_sd_5_traps_7_days.rds"))
  } # months
} # years

## Mapping_predictions ---------------------------------------------------------
years <- c("2020", "2021", "2022")
months <- c("01", "02", "03", "04", "05", "06" ,"07", "08", "09", "10", "11", "12")

fldr = "Monthly_Integration/"
mdl_name = "tiger_inter_monthly_meanma"
sub = ""

tiger_maps <- list()
month_values <- list()
iter <- 0

for(y in years){
  pred_month <- list()
  
  for (m in months){
    print(paste0("Plotting: ", m, "-", y))
    iter <- iter + 1
    
    pred <- readRDS(paste0(loc.output, "PREDICTIONS/", fldr, mdl_name, "/tiger_", m ,"_", y, sub, "1_traps_7_days.rds")) %>%
      janitor::clean_names() %>%
      st_drop_geometry()
    
    colnames(pred) <- c("municipality","id", "prov_name", "lon", "lat", "pred_count") 
    
    pred <- merge(pred, spain, by = c("municipality", "id", "prov_name"))
    st_geometry(pred) <- "geometry"
    
    print(max(pred$pred_count, na.rm = TRUE))
    plt <- ggplot() +
      geom_sf(data = pred, aes(fill = pred_count), color = "transparent",
              size = 0.01, alpha = 0.8, na.rm = TRUE) +
      scale_fill_distiller("Predicted\ncounts", palette = "Spectral", 
                           na.value = "transparent", limits = c(0, 23)) + 
      ggtitle(paste0(m, "-", y)) +
      theme_void(base_size = 8, base_family = "Helvetica")
    
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
ggsave(paste0(loc.fig, "monthly_tiger_2020_", mdl_name, ".pdf"), 
       width = 25, height = 18, dpi = 600, units = "cm", device = cairo_pdf)

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
ggsave(paste0(loc.fig, "monthly_tiger_2021_", mdl_name, ".pdf"), 
       width = 25, height = 18, dpi = 600, units = "cm", device = cairo_pdf)

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
ggsave(paste0(loc.fig, "monthly_tiger_2022_", mdl_name, ".pdf"), 
       width = 25, height = 18, dpi = 600, units = "cm", device = cairo_pdf)

rm(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12)

envelop <- fread(paste0(loc.data, "ES_albo_info.csv")) %>% # Private data of known Aedes albopictus distribution
  filter(Year != "año_2023")

natcode <- mapSpain::esp_get_munic() %>%
  mutate(
    natcode = paste0("34", codauto, cpro, LAU_CODE)
  ) %>%
  dplyr::select(natcode) %>%
  filter(
    natcode %in% envelop$NATCODE
  ) 
natcode <- st_make_valid(natcode)
natcode <- st_union(natcode)
natcode <- st_as_sf(natcode)

envelop_integrated <- ggplot() +
  geom_sf(data = pred_month, aes(fill = total), color = "white",
          size = 0.005, na.rm = TRUE) +
  scale_fill_distiller("Predicted\ncounts", palette = "Spectral", 
                       na.value = "transparent", 
                       limits = c(0, 12)
                       ) + 
  geom_sf(data = natcode, fill = "transparent") + 
  # xlim(-13, 5) +
  # ylim(34, 44) +
  theme_void(base_size = 8, base_family = "Helvetica")

envelop_trap <- ggplot() +
  geom_sf(data = pred_month, aes(fill = total), color = "white",
          size = 0.01, na.rm = TRUE) +
  scale_fill_distiller("Predicted\ncounts", palette = "Spectral", 
                       na.value = "transparent", 
                       limits = c(0, 12)
  ) + 
  geom_sf(data = natcode, fill = "transparent") + 
  # xlim(-13, 5) +
  # ylim(34, 44) +
  theme_void(base_size = 8, base_family = "Helvetica")

ggpubr::ggarrange(envelop_trap, envelop_integrated, common.legend = TRUE, 
                  legend = "right", labels = c("a", "b"),
                  label.y = 0.78)

ggsave(file = paste0(loc.fig, "Fig_6.pdf"),
       width = 25, height = 18, dpi = 600, units = "cm", device = cairo_pdf)
## Cross-validation ------------------------------------------------------------
# Split the data: train/test 
set.seed(121234)
sample <- sample(c(TRUE, FALSE), nrow(tiger_integrated), replace=TRUE, prob=c(0.7, 0.3))
train  <- tiger_integrated[sample, ]
test   <- tiger_integrated[!sample, ]

nchains = 4
threads_per_chain = 1
iteret = 5000
wup = 2000

mtiger_inte_month <- brm(females ~ ma + 
                           offset(log(n_traps)) + offset(log(n_days)) + 
                           (1 | id) + (1 | year),
                         data = train,
                         prior = set_prior("cauchy(0,2.5)", class="b"),
                         family = negbinomial(link = "log"),
                         iter = iteret,
                         chains = nchains,
                         cores = nchains,
                         backend = "cmdstanr",
                         save_pars = save_pars(all = TRUE),
                         threads = threading(threads_per_chain),
                         control = list(adapt_delta = 0.999))
summary(mtiger_inte_month)
# loo(mtiger_inte_month)
# loo(mtiger_inte_month, moment_match = TRUE, recompile = TRUE, reloo = TRUE)
bayes_R2(mtiger_inte_month)

# Predictions: 
pp <- apply(posterior_predict(mtiger_inte_month, 
                              newdata = test,
                              allow_new_levels = TRUE, 
                              # re_formula = NA,
                              ndraws = 1000), 
            2, function(x) mean(x)) %>%
  as.data.frame()

test$predicted <- pp$.

ggplot(test, aes(x = predicted, y = females)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() 
summary(lm(females ~ predicted, data = test))
summary(lm(I(females - 1) ~ predicted - 1, data = test))

cv_metrics <- data.frame(Scor = cor(test$predicted, test$females, method = "spearman"),
                         Pcor = cor(test$predicted, test$females, method = "pearson"),
                         R2 = R2(test$predicted, test$females),
                         RMSE = RMSE(test$predicted, test$females),
                         MAE = MAE(test$predicted, test$females),
                         model = "INTEGRATED model")
cv_metrics
