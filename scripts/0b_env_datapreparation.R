# load libraries
library(raster)
library(sf)
library(NINA)
library(caret)
library(dplyr)
library(tidyverse)
library(nngeo)
library(ecospat)

# Set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025/"
setwd(wd)

# store rasters in memory
rasterOptions(todisk = FALSE)

# read environmental raster layers
env_files <- list.files("data/environmental_layers/", pattern = "\\.asc$", full.names = TRUE)
all_env_layers <- lapply(env_files, raster)

# create a raster template to homogenize the environmental layers
rasterEx <- raster::extent(c(-180, 180, -90, 90))
rasterDim <- do.call(pmax, as.data.frame(sapply(all_env_layers, dim)))
rRe <- raster::raster(nrow=rasterDim[1],ncol=rasterDim[2], ext = rasterEx, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84")  # Set CRS to WGS84 (EPSG:4326)

# resample layers to template dimensions
all_env_stack <- stack(lapply(all_env_layers, resample, rRe, method = "ngb"))

# shift longitude coordinates from [-180,180] to [0,360]
## Get original extent
ext <- extent(all_env_stack)

## Define two extents
left_ext  <- extent(ext@xmin, 0, ext@ymin, ext@ymax)   # -180 to 0 (to be shifted)
right_ext <- extent(0, ext@xmax, ext@ymin, ext@ymax)   # 0 to 180 (unchanged)

## Crop rasters
left_raster  <- crop(all_env_stack, left_ext)  
right_raster <- crop(all_env_stack, right_ext)

## Shift left part from [-180,0] to [180,360]
left_raster <- shift(left_raster, dx = 360)

# If the rasters have different resolutions, resample them to the same resolution:
if (any(res(left_raster) != res(right_raster))) {
  left_ext <- extent(right_raster)
  right_ext <- extent(left_raster)
  extent(rRe) <- extent(left_ext@xmin, right_ext@xmax, ext@ymin, ext@ymax) # c(0, 360, -90, 90)
  left_raster <- resample(left_raster, rRe, method="bilinear")
  right_raster <- resample(right_raster, rRe, method="bilinear")
}

# Adjust the origin of the shifted raster (left_raster) to match right_raster
if (any(origin(left_raster) != origin(right_raster))) {
  origin(left_raster) <- origin(right_raster)
}

# Merge both rasters
all_env_stack_0_360 <- merge(right_raster, left_raster)

# Add layers' name
names(all_env_stack_0_360) <- names(all_env_stack)

# Simplify names 
names(all_env_stack_0_360) <- sub("Present.Benthic.Mean.Depth.", "", names(all_env_stack_0_360))

# Crop to area of the study
all_env_stack.cropped = crop(all_env_stack_0_360, extent(30, 240, -40, 40))

# save object
saveRDS(all_env_stack.cropped, "./data/all_env_stack.cropped.RDS")

# read object
all_env_stack.cropped <- readRDS("./data/all_env_stack.cropped.RDS")

# Convert to dataframe
all_env_stack.cropped.df <- as.data.frame(all_env_stack.cropped, xy = TRUE)

# remove rwos with complete NA cases
all_env_stack.cropped.df <- all_env_stack.cropped.df %>%
  filter(rowSums(is.na(select(., 3:ncol(.)))) != ncol(.) - 2)

# Select variables with less than 5% missing values
all_env_stack.cropped.df <- all_env_stack.cropped.df %>%
  select(which(colSums(is.na(.)) < 0.05 * nrow(.))) %>%
  filter(rowSums(!is.na(select(., 3:ncol(.)))) > 2)

write.csv(all_env_stack.cropped.df, file = "./data/env_NAfiltered_env.csv", row.names = FALSE)

################################################################################
## ----------------   Filter by physiological info ------------------------ ##
################################################################################
# Environmental filtering
# light > 0 (D Iluz, Z Dubinsky - Zoology, 2015)
# Temperature > 20 degrees Celsius (Methari, V. R., et al. Journal of Coastal Life Medicine, 2015)
# Depth > -200 (Epipelagic Zone)
all_env_stack.cropped.filtered.df <- all_env_stack.cropped.df %>%
  # Apply environmental filters
  filter(Light.bottom.Mean > 0, Temperature.Mean > 20, depth > -200) %>%
  # Remove rows with complete NA cases
  filter(rowSums(is.na(select(., 3:ncol(.)))) < (ncol(.) - 2)) %>%
  # Remove columns with all NAs or no variation
  select(where(~ !all(is.na(.)) && length(unique(na.omit(.))) > 1))

write.csv(all_env_stack.cropped.filtered.df, file = "./data/env_filtered_env.csv", row.names = FALSE)

################################################################################
## ----------------   Filter by coral reef locations ------------------------ ##
################################################################################
# extracting coral reef locations
coral_sf <- st_read("./data/WCMC_CoralReefs2018/01_Data/WCMC008_CoralReef2018_Py_v4.shp")

# shift longitude coordinates from [-180,180] to [0,360]
coral_sf <- st_shift_longitude(coral_sf)

# matching lon-lat coordinates from environmental dataframe with coral reef polygons shape file
## Convert environmental dataframe to an sf object 
all_env_stack.cropped.sf <- all_env_stack.cropped.df %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)  # Ensure it's in WGS 84

## Ensure coral reef data is in the same CRS
coral_sf <- coral_sf %>%
  st_make_valid() %>%  # Fix topology errors
  st_transform(crs = 4326)

## Perform spatial intersection to keep only points inside coral reefs
reef_filtered_env <- st_join(all_env_stack.cropped.sf, coral_sf, left = FALSE)

## Convert back to a dataframe 
reef_filtered_env.df <- reef_filtered_env %>%
  mutate(x = st_coordinates(.)[, 1],  # Extract longitude
         y = st_coordinates(.)[, 2]) %>%
  st_drop_geometry()  # Remove spatial geometry

# extract environmental variable names
env.var <- names(all_env_stack.cropped.filtered.df)[names(all_env_stack.cropped.filtered.df) %in% names(reef_filtered_env.df)]

# subset and order to keep only coordinates and environmental variables
reef_filtered_env.df <- reef_filtered_env.df[, env.var]
all_env_stack.cropped.filtered.df <- all_env_stack.cropped.filtered.df[, env.var]

# write dataset
write.csv(reef_filtered_env.df, file = "./data/reef_filtered_env.csv", row.names = FALSE)


final_all_env_stack.cropped.filtered_reef.df <- unique(rbind(all_env_stack.cropped.filtered.df, reef_filtered_env.df))

write.csv(final_all_env_stack.cropped.filtered_reef.df, file = "./data/env_reef_filtered_env.csv", row.names = FALSE)

