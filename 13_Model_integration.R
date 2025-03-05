########################## Integration MA on BG model ##########################
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

# mdl <- "mtiger16"
# mdl_name <- "/mtiger16"
# fldr <- "Counts"
# sub <- "" # with _

# mdl <- "mbites_7_reduce"
# mdl_name <- "/mbites_7_reduce"
# fldr <- "Bites"
# sub <- "bts" # with _

mdl_ma <- "mtiger7_ma"
mdl_name_ma <- "/mtiger7_ma"
fldr_ma <- "MA"
sub_ma <- "_ma" # with _

sf::sf_use_s2(FALSE)

# Integrating data -------------------------------------------------------------
spain <- readRDS(paste0(loc.output, "spain_mun.rds")) %>%
  mutate(
    municipality = if_else(is.na(municipality), "no_name", municipality)
  )

tiger <- readRDS(paste0(loc.output, "bg_tiger_spain_daily.rds")) %>% 
  filter(y != "2023") %>%
  mutate(
    year = year(end_date),
    month = month(end_date)
  ) 
# %>%
#   group_by(id, year, month) %>%
#   summarise(
#     females = sum(females, na.rm = TRUE),
#     l21precipitation = mean(l21precipitation, na.rm = TRUE),
#     trapping_effort = mean(trapping_effort, na.rm = TRUE)
#   )

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
      # p$count <- rowMeans(p[,6:ncol(p)], na.rm = TRUE)
      # p <- p[c("municipality","id", "prov_name", "lon", "lat", "count")]
      # pred <- cbind(pred, p["count"])
    }
  }
  # colnames(pred) <- c("municipality","id", "prov_name", "lon", "lat", months)
  # pred <- pred %>%
  #   pivot_longer(
  #     cols = "01":"12",
  #     names_to = "month",
  #     values_to = "ma"
  #   )
  # pred$year <- y
  # pred_all_ma <- rbind(pred_all_ma, pred)
} 

# integrated_data <- merge(tiger, pred_all_ma, by = c("id", "month", "year"))

pred$date <- as.Date(gsub("x", "", pred$date), format="%Y_%m_%d")

sum_ma <- vector()
prod_ma <- vector()
for (i in 1:nrow(tiger)){
  if (i %% 100 == 0) {print(i)}
  sum_ma[i] = pred %>% 
    filter(id == tiger$id[i]) %>%
    filter(date >= tiger$start_date[i] & date <= tiger$end_date[i]) %>%
    dplyr::select(ma) %>% 
    sum()
  prod_ma[i] = pred %>% 
    filter(id == tiger$id[i]) %>%
    filter(date >= tiger$start_date[i] & date <= tiger$end_date[i]) %>%
    dplyr::select(ma) %>% 
    prod()
}

integrated_data <- tiger %>%
  mutate(
    ma = sum_ma
  )

saveRDS(integrated_data, paste0(loc.output, "integrated_data_summa.rds"))
# Modeling process -------------------------------------------------------------
nchains = 4
threads_per_chain = 1
iteret = 5000
wup = 2000

mtiger_inte <- brm(females ~ poly(l21mean_temperature, 2) + l21precipitation  + ma + 
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
saveRDS(mtiger_inte, file = paste0(loc.output, "mtiger_inte_daily_sum_with_temperature.rds"))

# Predictions ------------------------------------------------------------------
mtiger_inte <- readRDS(file = paste0(loc.output, "mtiger_inte_daily_sum_with_temperature.rds"))
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
mdl_name = "tiger_inte_with_t"
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
          attractor = "CO2",
          n_traps = 1,
          SE = 1,
          n_total_reporters = 1,
          n_total_reports = 1,
          trapping_effort = 7
          # min_temperature = scale(min_temperature),
          # agricultural = scale(agricultural),
          # forests_scrub = scale(forests_scrub)
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
                    2, function(x) mean(x)) %>% # Or mean
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
