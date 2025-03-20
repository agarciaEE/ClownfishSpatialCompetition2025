library(NINA)
library(plyr)
library(raster)
library(tidyr)
library(dplyr)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

# load custom functions
source("./scripts/functions.R")

# store in memory
rasterOptions(todisk = FALSE)

ROU_data <- read.csv("./results/ROU_data.csv")
D_df <- read.csv("results/niche_overlaps.csv")

amph_EN <- readRDS("./Rdata/amphENMs.RDS")
anem_EN <- readRDS("./Rdata/anemENMs.RDS")
amph_BC <- readRDS("./Rdata/amphEBMs.RDS")

# association matrix
int.mat = read.csv(paste0("data/interaction_matrix.csv"), row.names = 1)
int.mat[int.mat > 0] = 1

# PA maps
amph_EN.PA <- prob_to_PA(amph_EN$maps, th = setNames(pmax(amph_EN$eval$threshold$threshold, 0.00001), rownames(amph_EN$eval$threshold)))
amph_BC.PA <- prob_to_PA(amph_BC$maps, th = setNames(pmax(amph_BC$eval$threshold$threshold, 0.00001), rownames(amph_BC$eval$threshold)))
anem_EN.PA <- prob_to_PA(anem_EN$maps, th = setNames(pmax(anem_EN$eval$threshold$threshold, 0.00001), rownames(anem_EN$eval$threshold)))

# Helper function to extract raster data
extract_raster <- function(raster_data, coords) {
  extracted <- raster::extract(raster_data, coords)
  extracted[is.na(extracted)] <- 0
  extracted
}

# Pre-allocate INT.pred
INT.pred <- data.frame(
  region = amph_EN$clus$province, 
  amph_EN$clus[,1:2],
  amph.richness = 0, anem.richness = 0,
  nG = 0, nS = 0,
  Region.interactions = 0,
  Local.interactions = 0,
  GG = 0, GS = 0, SS = 0, 
  HGG = 0, HGS = 0, HSS = 0,
  ENGG = 0, ENGS = 0, ENSS = 0,
  ROGG = 0, ROGS = 0, ROSS = 0
  
)

INT.data <- list()
# Loop through regions
for (region in names(amph_EN$z.mod)) {

  region_mask <- INT.pred$region == region
  coords <- INT.pred[region_mask, 2:3]
  
  clown_dis <- extract_raster(amph_BC.PA, coords)
  anem_dis <- extract_raster(anem_EN.PA, coords)
  
  INT.pred$amph.richness[region_mask] <- rowSums(clown_dis, na.rm = TRUE)
  INT.pred$anem.richness[region_mask] <- rowSums(anem_dis, na.rm = TRUE)
  
  Xvar <- names(which(colSums(clown_dis) != 0))
  Xint <- names(which(colSums(anem_dis) != 0))

  if (length(Xvar) == 0) { next }
  
  names(Xvar) <- sapply(Xvar, function(x) ROU_data$behavior_reg[match(paste(region, x), paste(ROU_data$region, ROU_data$species))])

  generalist_species <- Xvar[names(Xvar) == "generalist"]
  specialist_species <- Xvar[names(Xvar) == "specialist"]
  
  # Count generalists and specialists
  INT.pred[region_mask, "nG"] <- ifelse(length(generalist_species) > 0, rowSums(clown_dis[, generalist_species, drop = FALSE]), 0)
  INT.pred[region_mask, "nS"] <- ifelse(length(specialist_species) > 0, rowSums(clown_dis[, specialist_species, drop = FALSE]), 0)

  if (length(Xvar) > 1 && length(Xint) > 1) {
   
    INT.data[[region]] <- data.frame()
    
    # Pairwise interactions
    comb_sps <- combn(Xvar, 2)

    INT.pred[region_mask, "Region.interactions"] <- ncol(comb_sps)
    INT.pred[region_mask, "Local.interactions"] <- choose(INT.pred$amph.richness[region_mask], 2)
    
    pairwise_interaction <- function(spa, spb, correct = FALSE) {
      interaction <- (clown_dis[, spa] + clown_dis[, spb]) == 2
      if (correct) {
        interaction <- interaction * aint
      }
      as.integer(interaction)
    }
    
    for (i in seq_len(ncol(comb_sps))) {
      spa <- comb_sps[1, i]
      spb <- comb_sps[2, i]

      # Find shared interactions
      Xintab <- names(which(colSums(int.mat[c(spa, spb), Xint, drop = FALSE]) == 2))
      aint <- if (length(Xintab) == 0) 0 else rowSums(anem_dis[, Xintab, drop = FALSE]) > 0
      
      # Interaction type classification
      behavior_spa <- toupper(substr(names(Xvar[Xvar == spa]), 1, 1))
      behavior_spb <- toupper(substr(names(Xvar[Xvar == spb]), 1, 1))
      type <- paste(sort(c(behavior_spa, behavior_spb)), collapse = "")

      # number of interactions based on distribution
      INT.pred[region_mask, type] <- INT.pred[region_mask, type] +  pairwise_interaction(spa, spb)
      # number of interactions corrected by common hosts presence
      h_overlap <- pairwise_interaction(spa, spb, TRUE)
      INT.pred[region_mask, paste0("H", type)] <- INT.pred[region_mask, paste0("H", type)] + h_overlap
      
      # add environmental & resource overlap
      D_row <- unique(which(D_df$region == region & 
                     ((D_df$spa == spa & D_df$spb == spb) | (D_df$spa == spb & D_df$spb == spa))))
      
      if (length(D_row) == 1) {
        
        Denv <- h_overlap * ifelse(is.na(D_df[D_row, "Denv"]), 0, D_df[D_row, "Dmut"])
        Rove <- h_overlap * ifelse(is.na(D_df[D_row, "Roverlap"]), 0, D_df[D_row, "Roverlap"])
        
        INT.pred[region_mask, paste0("EN", type)] <- INT.pred[region_mask, paste0("EN", type)] + Denv
        INT.pred[region_mask, paste0("RO", type)] <- INT.pred[region_mask, paste0("RO", type)] + Rove
        
        INT.data[[region]] <- rbind(INT.data[[region]], data.frame(sp1 = spa, sp2 = spb, interaction = type, Denv, Rove))
        
      }
      
    }
    
  } 
  
}

# Total number of interactions
INT.pred$DTotal <- rowSums(INT.pred[, c("GG", "GS", "SS")], na.rm = TRUE) # based on distribution
INT.pred$HTotal <- rowSums(INT.pred[, c("HGG", "HGS", "HSS")], na.rm = TRUE) # based on habitat overlap
INT.pred$Ptotal <- with(INT.pred, ifelse(Local.interactions > 0, HTotal / Local.interactions, 0)) # proportion of habitat interactions relative to the total local interactions based on species richness

# Average environmental overlap per location
INT.pred$ENTotal <- with(INT.pred, ifelse(Local.interactions > 0,
                                          rowSums(INT.pred[, c("ENGG", "ENGS", "ENSS")], na.rm = TRUE) / Local.interactions,
                                          0))

INT.pred$ROTotal <- with(INT.pred, ifelse(Local.interactions > 0,
                                          rowSums(INT.pred[, c("ROGG", "ROGS", "ROSS")], na.rm = TRUE) / Local.interactions,
                                          0))

# divide each overlap by the number of interactions to obtain the average overlap
INT.pred$ENGG <- with(INT.pred, ifelse(HGG > 0, ENGG / HGG, 0))
INT.pred$ENGS <- with(INT.pred, ifelse(HGS > 0, ENGS / HGS, 0))
INT.pred$ENSS <- with(INT.pred, ifelse(HSS > 0, ENSS / HSS, 0))

INT.pred$ROGG <- with(INT.pred, ifelse(HGG > 0, ROGG / HGG, 0))
INT.pred$ROGS <- with(INT.pred, ifelse(HGS > 0, ROGS / HGS, 0))
INT.pred$ROSS <- with(INT.pred, ifelse(HSS > 0, ROSS / HSS, 0))

# Replace NaN, Inf, and NA with 0
INT.pred <- INT.pred %>%
  mutate(across(everything(), ~ ifelse(is.na(.) | is.nan(.) | is.infinite(.), 0, .)))

write.csv(INT.pred, file = "results/spatial_results.csv", row.names = F)

########

summary_INT_reg <- as.data.frame(INT.pred %>%
                                   group_by(region) %>%
                                   summarise_all(mean, na.rm = T)) 

summary_INT_reg[,-1] <- round(summary_INT_reg[,-1], 3)

write.csv(summary_INT_reg, file = "results/summary_spatial.csv", row.names = F)

cat("\nclownfish richness sorted:\n")
summary_INT_reg[order(summary_INT_reg$amph.richness, decreasing = T),] # clownfish richness
cat("\nenvironmental overlap sorted:\n")
summary_INT_reg[order(summary_INT_reg$ENTotal, decreasing = T),] # environmental overlap
cat("\nhost-specific overlap sorted:\n")
summary_INT_reg[order(summary_INT_reg$ROTotal, decreasing = T),] # host-specific overlap
cat("\ndifference between host-specific and ecological overlap sorted:\n")
summary_INT_reg[order(summary_INT_reg$ROTotal-summary_INT_reg$ENTotal, decreasing = F),] # difference between host-specific and ecological overlap

#correlation test
cat("\nPearson's correlation clownfish richness vs Host-specific overlap:\n")
cor.test(INT.pred$amph.richness, INT.pred$ROTotal)
cat("\nPearson's correlation clownfish richness vs environmental overlap:\n")
cor.test(INT.pred$amph.richness, INT.pred$ENTotal)
cat("\nPearson's correlation Host-specific overlap vs environmental overlap:\n")
cor.test(INT.pred$ROTotal, INT.pred$ENTotal)
