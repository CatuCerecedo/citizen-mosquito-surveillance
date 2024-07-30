########################## Preparing land Cover ################################
#' Percentages by grid

library(tidyverse)
library(sf)
library(janitor)
library(parallel)
library(dplyr)
library(lubridate)
library(terra)
library(data.table)
library(exactextractr)

rm(list = ls())

# Directories ------------------------------------------------------------------

loc.output <- paste0(getwd(), "/EU_Culex/OUTPUT/")
loc.data <- paste0(getwd(), "/EU_Culex/DATA/")

loc.era5 <- paste0(getwd(), "/EU_Culex/ERA5_Download/")

# # In local
# loc.output <- paste0(getwd(), "/OUTPUT/")
# loc.data <- paste0(getwd(), "/DATA/")
# 
# loc.era5 <- paste0(getwd(), "/ERA5_Download/")

sf::sf_use_s2(FALSE)
# Load raster ------------------------------------------------------------------

# Getting a ERA5 template
pname <- "ERA5_EU_hourly_"
tmp_temp <- rast(paste0(loc.era5, pname, "2021.nc"), "t2m") 
template <- tmp_temp[[1]] # Raster to include the mosquito data
# plot(template)
# crs(template)

# Alternative: Martas' script --------------------------------------------------
## from Corine landcover: https://land.copernicus.eu/en/products/corine-land-cover/clc2018/download-by-area
## Load raster ------------------------------------------------------------------
path <- paste0(loc.data, "u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif")
landcover <- rast(path)

# Little test
# esp_can <- st_transform(mapSpain::esp_get_prov("Galicia"), crs(landcover))
# landcover <- crop(landcover, esp_can)
# template <- crop(template, esp_can)
# plot(landcover)
# plot(template)

landcover <- project(landcover, crs(template))

# plot(landcover) # Plot the raster

df_cat <- levels(landcover)[[1]]

## Compute rater percentage for each landcover ---------------------------------
# Load weather raster template (see above)

# # Filter raster Category
# rast_cat <- landcover %in% df_cat[1,2]
# rast_cat <- as.numeric(rast_cat)
# 
# # Project and compute sum of areas of small squares.
# # example in:https://stackoverflow.com/questions/77421977/weighted-sum-in-terraproject-what-are-weights
# rast_per <- resample(rast_cat, template, method= "sum")
# 
# # plot(rast_per)

# Loop through all categories and compute percentage then add as raster layer
# for(i in c(1:nrow(df_cat))){
#   print(paste0("i:",i))
# 
#   # Filter raster Category
#   rast_cat <- landcover %in% df_cat[i,2]
#   rast_cat <- as.numeric(rast_cat)
# 
#   # Project and compute sum of areas of small squares.
#   # example in:https://stackoverflow.com/questions/77421977/weighted-sum-in-terraproject-what-are-weights
#   rast_aux <- resample(rast_cat, template, method= "sum")
# 
#   # plot(it_rast_per)
#   # I can not save in a raster because it is too large --> save each raster
#   saveRDS(rast_aux, file = paste0(loc.output, janitor::make_clean_names(df_cat[1,2]), "_rast.rds"))
#   rm(rast_aux)
#   # rast_per <- list(rast_per, rast_aux)
#   # rast_per <- rast(rast_per)
# }

# Function to process each chunk
process_chunk <- function(chunk_idx, df_cat, landcover, template, loc.output) {
  cat("Processing chunk:", df_cat[chunk_idx, 2], "\n")
  
  # Filter each category
  rast_cat <- landcover %in% df_cat[chunk_idx, 2]
  rast_cat <- as.numeric(rast_cat)
  
  # Projection and weighted sum of pixels
  rast_aux <- terra::resample(rast_cat, template, method= "sum")
  
  # Generating file name
  filename <- paste0(loc.output, make_clean_names(df_cat[chunk_idx, 2]), "_rast.rds")
  
  # Verify no empty name
  if (filename != "") {
    # Guardar el resultado en un archivo .rds
    cat("Saving: ", filename, "\n")
    saveRDS(rast_aux, file = filename)
  } else {
    cat("Error: Filename is empty for category", chunk_idx, "\n")
  }
  
  # Free memory
  rm(rast_aux)
  gc()
}

# Total number of rows in df_cat
n <- nrow(df_cat)

# # Number per each chunk
# chunk_size <- 10

# # Process by chunks
# for (start_idx in seq(1, n, by = chunk_size)) {
#   end_idx <- min(start_idx + chunk_size - 1, n)
#   cat("Number of categories:", end_idx, "\n")
#   
#   chunk_indices <- start_idx:end_idx
#   for (i in chunk_indices) {
#     process_chunk(i, df_cat, landcover, template, loc.output)
#   }
# }

# Process by categories
for (i in 1:n) {
  process_chunk(i, df_cat, landcover, template, loc.output)
}


# Joinning land covers ---------------------------------------------------------
tiger <- readRDS(paste0(loc.output, "tiger_spain_climatic_variables_pixel_monthly.rds"))
clc_surface <- readRDS(paste0(loc.output, "clc_all_pixels.rds"))

tiger <- merge(tiger, clc_surface, by = "pixel_id", all.x = TRUE)
tiger$no_data <- NULL

saveRDS(ma, file = paste0(loc.output, "tiger_spain.rds"))

