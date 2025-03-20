# Load required packages
library(dplyr)
library(stringr)
library(sp)
library(sf)
library(reshape)
library(rgbif)
library(ecospat)
library(raster)
library(robis)
library(readr)
library(ggplot2)
library(geosphere)

# Set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025/"
setwd(wd)

###################################################################################
########################### Anemones Occurrences Data #############################
###################################################################################
# Define target anemone species
anemones_sps <- c("Heteractis_magnifica", "Heteractis_malu", "Heteractis_crispa",  
                  "Cryptodendrum_adhaesivum", "Macrodactyla_doreensis", "Stichodactyla_mertensii",
                  "Stichodactyla_gigantea", "Heteractis_aurora", "Stichodactyla_haddoni", "Entacmaea_quadricolor")

anem_selected_species <- str_replace_all(anemones_sps, "_", " ")

anem_keys <- sapply(anem_selected_species, function(sp) name_suggest(q=sp)$data$key[1])

# Download occurrences from GBIF
anemones_GBIF <- do.call(rbind,
                         lapply(anem_keys, function(key) occ_search(taxonKey = key, 
                                                                    fields = c("taxonKey", "scientificName", "decimalLongitude", "decimalLatitude"))$data))

# Filter occurrences for target species
anem_keys <- setNames(names(anem_keys), as.character(anem_keys))
anemones_GBIF <- anemones_GBIF %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  mutate(species = anem_keys[as.character(taxonKey)]) %>%
  select(species, decimalLongitude, decimalLatitude)

anemones_GBIF <- na.exclude(anemones_GBIF)

# Rename columns
colnames(anemones_GBIF) <- c("species", "x", "y")

# Save dataset
write.csv(anemones_GBIF, "./data/raw/GBIF_anemones_dataset.csv", row.names = FALSE)

### MERGING DATASETS

# Read additional datasets
anemones_HEX <- read.csv("./data/raw/HEX_anemones_dataset.csv") 
colnames(anemones_HEX) <- c("species", "x", "y")

anemones_KGS <- read.csv("./data/raw/KGS_anemones_dataset.csv") %>%
  select(Species, decimalLongitude, decimalLatitude)
colnames(anemones_KGS) <- c("species", "x", "y")

# Merge datasets
anem.all.occs <- bind_rows(anemones_GBIF, anemones_HEX, anemones_KGS) %>%
  drop_na() %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  mutate(species = str_replace_all(species, "\\.", "_")) %>%
  distinct(species, x, y, .keep_all = TRUE)  %>% # Remove duplicates
  select(x, y, species)

# Save merged dataset
write.csv(anem.all.occs, "./data/raw/allDB.anem.occ.dataset.csv", row.names = FALSE)

###################################################################################
########################### Amphiprion occurrences data ###########################
###################################################################################

amphiprion_sps = c("Amphiprion_akallopisos"  , 
                   "Amphiprion_akindynos"  ,   "Amphiprion_allardi" ,
                   "Amphiprion_barberi"     ,  "Amphiprion_bicinctus" ,
                   "Amphiprion_chagosensis"  , "Amphiprion_chrysogaster",
                   "Amphiprion_chrysopterus",  "Amphiprion_clarkii" ,
                   "Amphiprion_ephippium"   ,  "Amphiprion_frenatus" ,
                   "Amphiprion_fuscocaudatus", "Amphiprion_latezonatus" ,
                   "Amphiprion_latifasciatus", "Amphiprion_leucokranos" ,
                   "Amphiprion_mccullochi"  ,  "Amphiprion_melanopus" ,
                   "Amphiprion_nigripes"    ,  "Amphiprion_ocellaris" ,
                   "Amphiprion_omanensis"   ,  "Amphiprion_pacificus" ,
                   "Amphiprion_percula"     ,  "Amphiprion_perideraion" ,
                   "Amphiprion_polymnus"    ,  "Amphiprion_rubrocinctus" ,
                   "Amphiprion_sandaracinos",  "Amphiprion_sebae"   ,
                   "Amphiprion_tricinctus", "Premnas_biaculeatus")

amph_selected_species <- str_replace_all(amphiprion_sps, "_", " ")

amph_keys <- sapply(amph_selected_species, function(sp) name_suggest(q=sp)$data$key[1])

# Download occurrences from GBIF
amphiprion_GBIF <- do.call(rbind,
                         lapply(amph_keys, function(key) occ_search(taxonKey = key, 
                                                                    fields = c("taxonKey", "scientificName", "decimalLongitude", "decimalLatitude"))$data))

# Filter occurrences for target species
amph_keys <- setNames(names(amph_keys), as.character(amph_keys))
amphiprion_GBIF <- amphiprion_GBIF %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  mutate(species = amph_keys[as.character(taxonKey)]) %>%
  select(species, decimalLongitude, decimalLatitude)

amphiprion_GBIF <- na.exclude(amphiprion_GBIF)

# Rename columns
colnames(amphiprion_GBIF) <- c("species", "x", "y")

# Save dataset
write.csv(amphiprion_GBIF, "./data/raw/GBIF_amphiprion_dataset.csv", row.names = FALSE)

# Read additional datasets
## OBIS
Amphiprion.IOBIS <- occurrence("Amphiprion")
Premnas.IOBIS <- occurrence("Premnas biaculeatus")

Amphiprion.IOBIS = Amphiprion.IOBIS %>%
  select(species, decimalLongitude, decimalLatitude)
Premnas.IOBIS = Premnas.IOBIS %>%
  select(species, decimalLongitude, decimalLatitude)

Amphiprion.IOBIS <- rbind(Amphiprion.IOBIS, Premnas.IOBIS)
colnames(Amphiprion.IOBIS) <- c("species", "x", "y")

# Save dataset
write.csv(Amphiprion.IOBIS, "./data/raw/IOBIS_amphiprion_dataset.csv", row.names = FALSE)

## FISHBASE
Amphiprion.FISHBASE <- read.csv("./data/raw/FISHBASE_amphiprion_dataset.csv")

Amphiprion.FISHBASE = Amphiprion.FISHBASE %>%
  select(species, lon, lat)
colnames(Amphiprion.FISHBASE) <- c("species", "x", "y")

## RLS
Amphiprion.RLS<- read.csv("./data/raw/RLS_amphiprion_dataset.csv")

Amphiprion.RLS = Amphiprion.RLS %>%
  select(species, lon, lat)
colnames(Amphiprion.RLS) <- c("species", "x", "y")

### MERGING DATASETS
amph.all.occs <- bind_rows(amphiprion_GBIF, Amphiprion.IOBIS, Amphiprion.FISHBASE, Amphiprion.RLS) %>%
  drop_na() %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  mutate(species = str_replace_all(species, "\\.", "_")) %>%
  distinct(species, x, y, .keep_all = TRUE)  # Remove duplicates

# correct names
unique(amph.all.occs$species)[!unique(amph.all.occs$species) %in% amphiprion_sps]
amph.all.occs$species <- stringr::str_replace_all(amph.all.occs$species, "Amphiprion_biaculeatus", "Premnas_biaculeatus")
amph.all.occs$species <- stringr::str_replace_all(amph.all.occs$species, "AMPHIPRION_PERCULA", "Amphiprion_percula")
amph.all.occs$species <- stringr::str_replace_all(amph.all.occs$species, "Amphiprion_sebae_", "Amphiprion_sebae")

# drop unknown species
amph.all.occs = amph.all.occs %>%
  filter(species %in% amphiprion.sps)  %>%
  distinct(species, x, y, .keep_all = TRUE)  %>% # Remove duplicates
  select(x, y, species)

# remove entries wirh NA values
amph.all.occs = na.exclude(amph.all.occs)

table(amph.all.occs$species)
nrow(amph.all.occs)

write.csv(amph.all.occs, file =  "./data/raw/allDB.amph.occ.dataset.csv", row.names = FALSE)

##### -------------------------------------------------------------------- #####
#####                                                                      #####
#####             After first step on 0b_env.var.selection.R               #####
#####                                                                      #####
##### -------------------------------------------------------------------- #####

################################################################################
## ---------------------   Check occurrences data --------------------------- ##
################################################################################

### ------------------------------------------------------------------------ ###
##      Check if occurrences fall inside the environmental data set
### ------------------------------------------------------------------------ ###
### read occurrence datasets
amph.occ = read.csv("./data/raw/allDB.amph.occ.dataset.csv")
anem.occ = read.csv("./data/raw/allDB.anem.occ.dataset.csv")

amphiprion.sps <- unique(amph.occ$species)
anemone.sps <- unique(anem.occ$species)

### combine  occs  datasets
all.occ <- rbind(amph.occ, anem.occ)

### shift longitude to 0-360
all.occ$x[all.occ$x < 0] = all.occ$x[all.occ$x < 0] + 360

### read environmental dataset
env.df <- read.csv("./data/selected_environmental_variables.csv")

# convert environmental dataset to raster stack
env.stack <- rast(env.df)

# get maps resolution
res <- max(res(env.stack))

# match occurrence locations with environmental locations with a margin of two times the resolution
all.occ.env.inn <- ecospat.sample.envar(dfsp = all.occ, colspxy = 1:2, colspkept = 1:3,
                                        dfvar = env.df, colvarxy = 1:2, colvar = 1:ncol(env.df), 
                                        resolution = res * 2)

# remove not matching rows
all.occ.env.inn <- all.occ.env.inn[!is.na(all.occ.env.inn[,4]),]

# replace xy occurrence coordinates by matched environmental xy coordinates for precision
all.occ.env.inn <- cbind(all.occ.env.inn[,4:5], 
                         species = all.occ.env.inn[,3], all.occ.env.inn[,6:ncol(all.occ.env.inn)])

# get number of observations per species
all.occ.env.inn.tab <- table(all.occ.env.inn$species)

# print
print(all.occ.env.inn.tab[all.occ.env.inn.tab > 4])

# split datasets into clownfish and sea anemones
amph_env_filtered.df <- all.occ.env.inn %>%
  dplyr::filter(species %in% amphiprion.sps) %>%
  dplyr::select(x, y, species)

anem_env_filtered.df <- all.occ.env.inn %>%
  dplyr::filter(species %in% anemone.sps) %>%
  dplyr::select(x, y, species)

# write filtered datasets
write.csv(amph_env_filtered.df, file = "./data/amph_env_filtered.csv", row.names = FALSE)
write.csv(anem_env_filtered.df, file = "./data/anem_env_filtered.csv", row.names = FALSE)

### ------------------------------------------------------------------------ ###
##    Filtering occurrences species per species based on known distributions
### ------------------------------------------------------------------------ ###

amph_env_filtered.df = read.csv("./data/amph_env_filtered.csv")
anem_env_filtered.df = read.csv("./data/anem_env_filtered.csv")

# Sea anemones
# Most with widespread distributions
# Filter species absent in the West Indian Ocean
invalid_species <- c("Macrodactyla_doreensis", "Heteractis_malu")

anem_env_filtered.df <- anem_env_filtered.df %>%
  filter(!(species %in% invalid_species & x < 100)) %>%
  drop_na()

# Remove duplicate records
anem_env_filtered.df <- anem_env_filtered.df %>%
  dplyr::distinct()

nrow(anem_env_filtered.df)
anem_tab <- table(anem_env_filtered.df$species)
summary(as.vector(anem_tab))

write.csv(anem_env_filtered.df, file = "./data/anem_occ_env_final_dataset.csv", row.names = FALSE)

# Clownfishes
### Sources: fishbase.org and amphiprionology blog 
### Using marine regions MEOUW_ECOS

prov_data <- st_read("./data/marine_regions.shp")
provs.points <- read.csv("data/meow_ecos_df.csv")

prov_legend <- provs.points %>%
  arrange(long) %>% # Order by longitude
  dplyr::select(PROVINCE, PROV_CODE) %>%
  distinct(.keep_all = TRUE) 

# Plot marine regions with names
ggplot(data = prov_data) +
  geom_sf(fill = "lightblue", color = "black", alpha = 0.5) +  # Plot polygons
  geom_sf_text(aes(label = PROVINCE), size = 3, color = "black") +  # Add labels
  theme_minimal() +
  labs(title = "Marine Regions", x = "Longitude", y = "Latitude")

species_regions <- list(
   Amphiprion_akallopisos    = c("Western Indian Ocean", "West and South Indian Shelf", 
                                 "Central Indian Ocean Islands", "Andaman", "Java Transitional"),
   
   Amphiprion_clarkii        = c("Somali/Arabian", "West and South Indian Shelf", "Central Indian Ocean Islands",
                                 "Andaman", "Java Transitional", "Sunda Shelf", "South China Sea", 
                                 "West Central Australian Shelf", "Northwest Australian Shelf", "Western Coral Triangle",
                                 "Warm Temperate Northwest Pacific", "South Kuroshio", "Sahul Shelf", 
                                 "Tropical Northwestern Pacific", "Northeast Australian Shelf",  
                                 "Eastern Coral Triangle"),
   
   Amphiprion_ephippium      = c("Andaman"),
   
   Amphiprion_frenatus       = c("Andaman", "Java Transitional", "Western Coral Triangle", "Sunda Shelf",
                                 "South China Sea", "South Kuroshio", "Warm Temperate Northwest Pacific",
                                 "Tropical Northwestern Pacific"),
   
   Amphiprion_ocellaris      = c("Andaman", "Java Transitional",
                                 "Northwest Australian Shelf",  "Western Coral Triangle", "Sunda Shelf",
                                 "South China Sea", "South Kuroshio", "Warm Temperate Northwest Pacific", 
                                 "Tropical Northwestern Pacific", "Sahul Shelf"),
   
   Amphiprion_sebae          = c("Somali/Arabian", "West and South Indian Shelf", "South China Sea", 
                                 "Andaman", "Central Indian Ocean Islands", 
                                 "Sunda Shelf", "Java Transitional", "Western Coral Triangle"),
   
   Amphiprion_bicinctus      = c("Somali/Arabian", "Red Sea and Gulf of Aden"),
   
   Amphiprion_chagosensis    = c("Central Indian Ocean Islands"),
   
   Amphiprion_nigripes       = c("Central Indian Ocean Islands", "West and South Indian Shelf"),
   
   Amphiprion_chrysopterus   = c("Western Coral Triangle", "Marshall, Gilbert and Ellis Islands", "Tropical Northwestern Pacific", 
                                 "Eastern Coral Triangle", "Sahul Shelf", "Northeast Australian Shelf", 
                                 "Tropical Southwestern Pacific", "Central Polynesia", "Southeast Polynesia"),
   
   Amphiprion_akindynos      = c("East Central Australian Shelf", "Northeast Australian Shelf", 
                                "Sahul Shelf", "Tropical Southwestern Pacific"),
   
   Amphiprion_latezonatus    = c("East Central Australian Shelf", "Tropical Southwestern Pacific", "Lord Howe and Norfolk Islands"),
   
   Amphiprion_melanopus      = c("Java Transitional", "Western Coral Triangle", "Sunda Shelf",
                                 "Marshall, Gilbert and Ellis Islands", "Tropical Northwestern Pacific", 
                                 "Eastern Coral Triangle", "Sahul Shelf", "Northeast Australian Shelf", 
                                 "East Central Australian Shelf",  
                                 "Tropical Southwestern Pacific", "Central Polynesia"),
   
   Amphiprion_leucokranos    = c("Eastern Coral Triangle"),
   
   Amphiprion_percula        = c("Western Coral Triangle", "Marshall, Gilbert and Ellis Islands", "Eastern Coral Triangle",
                                 "Sahul Shelf", "Northeast Australian Shelf", "East Central Australian Shelf",
                                 "Southeast Polynesia", "Tropical Southwestern Pacific"),
   
   Amphiprion_perideraion    = c("Andaman", "Java Transitional", "West Central Australian Shelf",  
                                 "Northwest Australian Shelf",  "Western Coral Triangle", "Sunda Shelf",
                                 "South China Sea", "South Kuroshio", "Warm Temperate Northwest Pacific", 
                                 "Marshall, Gilbert and Ellis Islands", "Tropical Northwestern Pacific",
                                 "Eastern Coral Triangle", "Sahul Shelf",  "Northeast Australian Shelf", 
                                 "East Central Australian Shelf", "Tropical Southwestern Pacific", 
                                 "Central Polynesia"),
   
   Amphiprion_polymnus       = c("Andaman", "Java Transitional", "Northwest Australian Shelf", "Northeast Australian Shelf", "Sahul Shelf",
                                 "Eastern Coral Triangle", "South China Sea", "Sunda Shelf", "Western Coral Triangle",
                                 "Tropical Northwestern Pacific", "South Kuroshio", "Warm Temperate Northwest Pacific"),
   
   Amphiprion_sandaracinos   = c("Eastern Coral Triangle", "South China Sea", "Sunda Shelf", "Western Coral Triangle",
                                 "Java Transitional", "South Kuroshio", "Northwest Australian Shelf", "Sahul Shelf"),
   
   Premnas_biaculeatus       = c("Eastern Coral Triangle", "Java Transitional", "Northeast Australian Shelf", 
                                "Sahul Shelf", "Sunda Shelf", "Western Coral Triangle", "Tropical Southwestern Pacific",
                                "Tropical Northwestern Pacific", "South Kuroshio", "Andaman", "Northwest Australian Shelf"),
   
   Amphiprion_mccullochi     = c("Lord Howe and Norfolk Islands"),
   
   Amphiprion_rubrocinctus   = c("Northwest Australian Shelf", "Sahul Shelf", "West Central Australian Shelf"),
   
   Amphiprion_omanensis      = c("Somali/Arabian"),
   
   Amphiprion_barberi        = c("Tropical Southwestern Pacific", "Central Polynesia", "Southeast Polynesia"),
   
   Amphiprion_tricinctus     = c("Tropical Northwestern Pacific", "Marshall, Gilbert and Ellis Islands"),
   
   Amphiprion_allardi        = c("Western Indian Ocean"),
   
   Amphiprion_chrysogaster   = c("Western Indian Ocean"),
   
   Amphiprion_fuscocaudatus  = c("Western Indian Ocean"),
   
   Amphiprion_pacificus      = c("Eastern Coral Triangle", "Tropical Southwestern Pacific"),
   
   Amphiprion_latifasciatus  = c("Western Indian Ocean")
)

# Plot marine regions with names
sp = names(species_regions)[29]
ggplot(data = prov_data[prov_data$PROVINCE %in% species_regions[[sp]],]) +
  geom_sf(fill = "lightblue", color = "black", alpha = 0.5) +  # Plot polygons
  geom_sf_text(aes(label = PROVINCE), size = 3, color = "black") +  # Add labels
  theme_minimal() +
  labs(title = sp, x = "Longitude", y = "Latitude")

# Merge species occurrences with their respective regions
amph_env_filtered.sf <- amph_env_filtered.df %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)  # Ensure it's in WGS 84

filtered_occurrences <- amph_env_filtered.sf %>%
  rowwise() %>%
  filter({
    # Get species name
    species_name <- species
    
    # Get the regions for the given species
    allowed_regions <- species_regions[[species_name]]
    
    if (is.null(allowed_regions)) return(FALSE)  # Skip if species has no defined regions
    
    # Subset the marine regions to only those corresponding to this species
    species_region_polygons <- prov_data %>%
      filter(PROVINCE %in% allowed_regions)
    
    # Check if occurrence falls within the filtered region polygons
    any(st_contains(species_region_polygons, geometry, sparse = FALSE))
  }) %>%
  ungroup()

# Convert sf object to data frame with coordinates
amph_env_filtered.df <- filtered_occurrences %>%
  cbind(st_coordinates(.)) %>%
  st_drop_geometry() %>%
  dplyr::select(X, Y, species)

# rename columns
colnames(amph_env_filtered.df) <- c("x", "y", "species")

# View the first few rows
head(amph_env_filtered.df)

# Remove duplicate records
amph_env_filtered.df <- amph_env_filtered.df %>%
  dplyr::distinct()

amph_tab <- table(amph_env_filtered.df$species)
nrow(amph_env_filtered.df)
summary(as.vector(amph_tab))

write.csv(amph_env_filtered.df, "./data/amph_occ_env_final_dataset.csv", row.names = FALSE)

