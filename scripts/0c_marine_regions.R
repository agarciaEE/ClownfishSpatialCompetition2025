#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load required libraries
library(sp)
library(plyr)
library(ecospat)
library(ggplot2)
library(sf)
library(dplyr)

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025/"
setwd(wd)

# extract MEOW ECOS shapefile and filter accorginly
regions <- st_read("raw_data/MEOW_ECOS/", "meow_ecos")
regions <- regions[regions$PROV_CODE %in% c(18:41, 58, 55, 9),] # select only regions of interest
regions <- st_transform(regions, crs = 4326) # transform to WGS84

# Dissolve polygons by province
prov_data <- regions %>%
  group_by(PROVINCE) %>%
  dplyr::summarise(
    PROV_CODE = first(PROV_CODE),
    REALM = first(REALM),
    RLM_CODE = first(RLM_CODE),
    Lat_Zone = first(Lat_Zone),
    geometry = st_union(geometry)
  ) %>%
  ungroup()

# shift longitude axis
prov_data <- st_shift_longitude(prov_data)

# save shapefile
st_write(prov_data, "./data/marine_regions.shp", delete_layer = TRUE)

# Keep province ID inside geometry before extracting coordinates
provs_geom <- prov_data %>%
  st_cast("MULTIPOLYGON", warn = FALSE) %>%
  st_cast("POLYGON", warn = FALSE) 

# add geometry id column
provs_geom$id = 1:nrow(provs_geom)

# Extract coordinates while maintaining IDs
provs_coords <- st_coordinates(provs_geom) %>%
  as.data.frame()

# Rename coordinate columns correctly
colnames(provs_coords) <- c("long", "lat", "feature", "id")

# Join polygon attributes back to coordinates using correct ID
provs.points <- left_join(provs_coords, provs_geom %>% st_drop_geometry(), by = "id")

# Save polygon data into csv file
write.csv(provs.points, file = "data/meow_ecos_df.csv", row.names = FALSE)

##### -------------------------------------------------------------------- #####
#####                                                                      #####
#####                After finished 0b_env.var.selection.R                 #####
#####                                                                      #####
##### -------------------------------------------------------------------- #####

# Load datasets
mr_sf <- prov_data["PROVINCE"]  # Marine regions
env <- read.csv("data/selected_environmental_variables.csv") # Environmental dataset

# Convert environmental data to sf spatial points
env_sf <- st_as_sf(env, coords = c("x", "y"), crs = 4326)

mr_sf <- st_make_valid(mr_sf)  # Fix invalid geometries

# Spatial join: Assign marine region to each environmental point
env_sf <- st_join(env_sf, mr_sf, left = TRUE)

# Find rows where marine region is NA
missing_idx <- which(is.na(env_sf$PROVINCE))  

if (length(missing_idx) > 0) {
  # Find the nearest marine region for each missing row
  nearest_idx <- st_nearest_feature(env_sf[missing_idx, ], mr_sf)
  
  # Add the corresponding marine region data
  env_sf$PROVINCE[missing_idx] <- mr_sf$PROVINCE[nearest_idx]
}

# Convert sf object to data frame with coordinates
env_df <- env_sf %>%
  cbind(st_coordinates(.)) %>%
  st_drop_geometry() %>%
  dplyr::select(X, Y, colnames(env)[-c(1:2)])

colnames(env_df) <- tolower(colnames(env_df))

# rewrite env data (in case some env location fell out of the regions)
write.csv(env_df, file = "data/selected_environmental_variables.csv", row.names = FALSE)

marine.regions <- env_sf %>%
  cbind(st_coordinates(.)) %>%
  st_drop_geometry() %>%
  dplyr::select(X, Y, PROVINCE)

colnames(marine.regions) <- tolower(colnames(marine.regions))

# Add correspondent realms
mr <- read.csv("data/meow_ecos_df.csv") # Marine regions dataset
for (i in unique(marine.regions$province)) {
  marine.regions$realm[marine.regions$province == i] = mr[mr$PROVINCE == i, "REALM"][1]
}

# modify province and realms names to be more user-friendly (replace spaces & special characters for underscore)
marine.regions$province = gsub("[^[:alnum:]]", "_", marine.regions$province)
marine.regions$realm = gsub("[^[:alnum:]]", "_", marine.regions$realm)

write.csv(marine.regions, file = "data/marine_regions.csv", row.names = FALSE)
