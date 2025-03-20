# load libraries
library(raster)
library(NINA)
library(stringr)
library(ggplot2)
library(ggpubr)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

# load custom functions
source("./scripts/functions.R")

# Metric functions
##################

estimate_niche_ROU <- function(en_niche, bc_niche, w_niche) {
  
  # niche occupancy
  en_o = en_niche$w
  bc_o = bc_niche$w
  w_o = w_niche$w
  
  if (raster::cellStats(bc_o, "max") == 0) {
    w_o = bc_o
  }
  
  # compare rasters
  if (!raster::compareRaster(en_o, bc_o, w_o, stopiffalse = FALSE)) {
    stop("Rasters of different extent.")
  }
  
  # compute ROU metrics
  R <- raster::cellStats(en_o > 0 & bc_o == 0, sum, na.rm = TRUE) / raster::cellStats(en_o > 0, sum, na.rm = TRUE) # Restricted niche
  O <- raster::cellStats(en_o > 0 & bc_o > 0, sum, na.rm = TRUE) / raster::cellStats(en_o > 0, sum, na.rm = TRUE) # Occupied niche
  U <- 1 - raster::cellStats(bc_o > 0, sum, na.rm = TRUE) / raster::cellStats(w_o > 0, sum, na.rm = TRUE) # Unexploited niche
  
  return(c("R" = R, 
           "O" = O, 
           "U" = U))
}

niche_position <- function(z, cor = FALSE, quantile = 0.5) {
  
  R <- length(z$x)
  v <- raster::as.data.frame(if (cor) z$z.cor else z$z.uncor, xy = TRUE)
  v[is.na(v)] <- 0
  v <- v[v[, 3] != 0, ]
  
  if (nrow(v) == 0) {
    return(NA)
  }
  
  qt <- quantile(v[, 3], quantile, na.rm = TRUE)
  opt <- v[v[, 3] >= qt, ]
  
  ctr <- stats::cov.wt(opt[, 1:2], wt = opt[, 3])$center
  
  return(ctr)
}

centroid_shift <- function(x, y) {
  
  Cx = niche_position(x, cor = FALSE, quantile = 0.5)
  Cy = niche_position(y, cor = FALSE, quantile = 0.5)
  
  CS = as.numeric(sqrt((Cx[1] - Cy[1])^2 + (Cx[2] - Cy[2])^2))
  
  return(CS)
}

environmental_shift <- function(x, y) {
  
  Cx = niche_position(x, cor = FALSE, quantile = 0.5)
  Cy = niche_position(y, cor = FALSE, quantile = 0.5)
  
  alpha = Morpho::angle.calc(Cx, Cy)/pi
  
  return(alpha)
  
}

niche_dissimilarity <- function(x, y, cor = FALSE) {
  
  if (!raster::compareRaster(x$Z, y$Z, stopiffalse = FALSE)) {
    stop("Rasters of different extent.")
  }
  
  if (cor) {
    p1 <- t(raster::as.matrix(x$z.cor)/sum(raster::as.matrix(x$z.cor)))
    p2 <- t(raster::as.matrix(y$z.cor)/sum(raster::as.matrix(y$z.cor)))
  }
  else {
    p1 <- t(raster::as.matrix(x$z.uncor)/sum(raster::as.matrix(x$z.uncor)))
    p2 <- t(raster::as.matrix(y$z.uncor)/sum(raster::as.matrix(y$z.uncor)))
  }
  
  D <- 1 - (0.5 * (sum(abs(p1 - p2))))
  
  return( 1 - D)
}

estimate_spatial_ROU <- function(en_map, bc_map, w_map) {
  
  # compare rasters
  if (!raster::compareRaster(en_map, bc_map, w_map, stopiffalse = FALSE)) {
    stop("Rasters of different extent.")
  }
  
  # compute ROU metrics
  R <- raster::cellStats(en_map > 0 & bc_map == 0, sum, na.rm = TRUE) / raster::cellStats(en_map > 0, sum, na.rm = TRUE) # Restricted distribution
  O <- raster::cellStats(en_map > 0 & bc_map > 0, sum, na.rm = TRUE) / raster::cellStats(en_map > 0, sum, na.rm = TRUE) # Occupied distribution
  U <- 1 - raster::cellStats(bc_map > 0, sum, na.rm = TRUE) / raster::cellStats(w_map > 0, sum, na.rm = TRUE) # Unexploited distribution
  
  return(c("R" = R, 
           "O" = O, 
           "U" = U))
}

# store in memory
rasterOptions(todisk = FALSE)

amph_EN <- readRDS("./Rdata/amphENMs.RDS")
anem_EN <- readRDS("./Rdata/anemENMs.RDS")
amph_BC <- readRDS("./Rdata/amphEBMs.RDS")

# association matrix
int.mat = read.csv("data/interaction_matrix.csv", row.names = 1)
int.mat[int.mat > 0] = 1
##################

#  merge spatial distributions dataset and marine regions
amph_ENdis <- merge(amph_EN$pred.dis, amph_EN$clus[,1:3], by = c("x", "y"))
amph_BCdis <- merge(amph_BC$pred.dis, amph_BC$clus[,1:3], by = c("x", "y"))
anem_ENdis <- merge(anem_EN$pred.dis, anem_EN$clus[,1:3], by = c("x", "y"))

ROU_data <- data.frame()
for (reg in names(amph_BC$z.mod)) {
  
  for (sp in names(amph_BC$z.mod[[reg]])) {

    # specific associations
    Xvar = colnames(int.mat)[int.mat[sp,] == 1]
    # sea anemone species in region
    hosts = names(anem_EN$z.mod[[reg]])
    # number of hosts in region
    hosts_reg <- sum(Xvar %in% hosts) 
    # proportion of hosts over total of sea anemone species in region
    host_use = hosts_reg/length(hosts)
    # region specific behavior
    behavior_reg <- ifelse(hosts_reg > 2, "generalist", "specialist")
    # global behavior
    behavior_glob <- ifelse(sum(int.mat[sp,]) > 2, "generalist", "specialist")
    
    en_niche <- amph_EN$z.mod[[reg]][[sp]]
    bc_niche <- amph_BC$z.mod[[reg]][[sp]]
    w_niche <- amph_BC$w[[reg]][[sp]]
    
    ## Niche
    ROU_niche <- estimate_niche_ROU(en_niche, bc_niche, w_niche)
    
    Cen_shift <- centroid_shift(en_niche, bc_niche)
    
    Env_shift <- environmental_shift(en_niche, bc_niche)
    
    Niche_diss <-  niche_dissimilarity(en_niche, bc_niche)
    
    ## Spatial

    sp_EN_th <- amph_EN$eval$threshold[sp,1]
    sp_BC_th <- amph_BC$eval$threshold[sp,1]
    w_EN_th <- setNames(anem_EN$eval$threshold$threshold,
                        rownames(anem_EN$eval$threshold))
    
    # get regional distributions
    sp_ENdis <- amph_ENdis %>%
      dplyr::filter(province == reg) %>%
      dplyr::select(x, y, all_of(sp))
  
    sp_BCdis <- amph_BCdis %>%
      dplyr::filter(province == reg) %>%
      dplyr::select(x, y, all_of(sp))
    
    w_ENdis <- anem_ENdis %>%
      dplyr::filter(province == reg) %>%
      dplyr::select(-province)
    
    # convert to rasters
    sp_ENras <- raster_projection(sp_ENdis, crs = 4326)
    sp_BCras <- raster_projection(sp_BCdis, crs = 4326)
    w_ENras <- raster_projection(w_ENdis,  crs = 4326)

    # convert probabilities to Presence/Absences
    sp_ENras_PA <- prob_to_PA(sp_ENras, sp_EN_th)
    sp_BCras_PA <- prob_to_PA(sp_BCras, sp_BC_th)
    
    
    w_ENras_PA <- sum(prob_to_PA(w_ENras, w_EN_th))
    w_ENras_PA[w_ENras_PA > 1] = 1 # Presence/Absence of any host 

    # Spatial ROU
    ROU_spatial <- estimate_spatial_ROU(sp_ENras_PA, sp_BCras_PA, w_ENras_PA)
    
    ROU_data <- rbind(ROU_data, data.frame(species = sp, 
                                           behavior = behavior_glob,
                                           region = reg,
                                           num_hosts = hosts_reg,
                                           p_host_shared = host_use,
                                           behavior_reg = behavior_reg,
                                           Restricted_niche = ROU_niche["R"], 
                                           Occupied_niche = ROU_niche["O"], 
                                           Unexploited_niche = ROU_niche["U"],
                                           Centroid_shift = Cen_shift,
                                           Environmental_shift = Env_shift,
                                           Niche_dissimilarity = Niche_diss,
                                           Restricted_spatial = ROU_spatial["R"], 
                                           Occupied_spatial = ROU_spatial["O"], 
                                           Unexploited_spatial = ROU_spatial["U"]
                                           ))
  }
}

dir.create("results", showWarnings = FALSE)
write.csv(ROU_data, file = "./results/ROU_data.csv", row.names = FALSE)

