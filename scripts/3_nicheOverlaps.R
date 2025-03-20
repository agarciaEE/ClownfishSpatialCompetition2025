# load libraries
library(NINA)
library(plyr)
library(raster)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

# load custom functions
source("./scripts/functions.R")

# niche overlap functions
#########################

niche_overlap <- function(x, y, cor = FALSE) {
  
  if (is.null(x) | is.null(y)) {
    return(NA)
  }
  
  if (raster::cellStats(x$z, "max") == 0 | raster::cellStats(y$z, "max") == 0) {
    return(NA)
  }
  
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
  
  return(D)
  
}

resource_overlap <- function(x, y, w.list, Xvar, Yvar, cor = FALSE) {
  
  if (is.null(x) | is.null(y)) {
    return(NA)
  }
  
  if (raster::cellStats(x$z, "max") == 0 | raster::cellStats(y$z, "max") == 0) {
    return(NA)
  }
  
  if (!raster::compareRaster(c(x$Z, y$Z, 
                               lapply(w.list, function(w) w$Z)), stopiffalse = FALSE)) {
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
  
  # get environmental niche overlap
  Denv <- 1 - (0.5 * (sum(abs(p1 - p2))))
  
  if ( Denv == 0 ) { return(Denv) }

  Deco <- c()
  for (var in Yvar){
    if (var %in% Xvar){
      
      # Refine niche layers  
      if (cor) {
        z1 <- x$z.cor * w.list[[var]]$z.cor
        z2 <- y$z.cor * w.list[[var]]$z.cor
      }
      else {
        z1 <- x$z.uncor * w.list[[var]]$z.uncor
        z2 <- y$z.uncor * w.list[[var]]$z.uncor
      }
      
      # integrate values to 1
      p1 <- as.matrix(z1)/sum(as.matrix(z1))
      p2 <- as.matrix(z2)/sum(as.matrix(z2))
      
      # compute D
      Deco <- c(Deco, 1 - (0.5 * (sum(abs(p1 - p2)))))
    }
    else{
      Deco <- c(Deco, 0)
    }
  }
  Deco <- mean(Deco, na.rm = TRUE) / max(mean(Deco, na.rm = TRUE), Denv) # normalize by the environmental overlap (or if it is higher by itself as it cannot be higher than 1)
  
  return(Deco)
}

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

# store in memory
rasterOptions(todisk = FALSE)

ROU_data <- read.csv("./results/ROU_data.csv")

amph_EN <- readRDS("./Rdata/amphENMs.RDS")
anem_EN <- readRDS("./Rdata/anemENMs.RDS")
amph_BC <- readRDS("./Rdata/amphEBMs.RDS")

# association matrix
int.mat = read.csv(paste0("data/interaction_matrix.csv"), row.names = 1)
int.mat[int.mat > 0] = 1
#########################

## Species pairwise niche overlaps
#####

pairs_list <- plyr::ldply(
                    sapply(amph_EN$z.mod, function(reg) {
                    if (length(reg) > 1) {
                      t(data.frame(combn(names(reg),2), 
                                   row.names = c("spa", "spb")))
                    } else {
                      NULL
                    }
                    }),
                    .id = "region"
                  )

for (i in 1:nrow(pairs_list)) {
  
  reg <- as.character(pairs_list$region[i])
  spa <- as.character(pairs_list$spa[i])
  spb <- as.character(pairs_list$spb[i])
  
  # species behavior
  spa_behavior <- ROU_data$behavior_reg[ROU_data$region == reg &
                                        ROU_data$species == spa]
  spb_behavior <- ROU_data$behavior_reg[ROU_data$region == reg &
                                          ROU_data$species == spb]
  
  # environmental niches
  spa_en <- amph_EN$z.mod[[reg]][[spa]]
  spb_en <- amph_EN$z.mod[[reg]][[spb]]
  
  # mutualism refined niches
  spa_bc <- amph_BC$z.mod[[reg]][[spa]]
  spb_bc <- amph_BC$z.mod[[reg]][[spb]]

  # host niche lists
  w_list <- anem_EN$z.mod[[reg]]

  # Shared hosts between spA and spB
  Xvar <- colnames(int.mat)[colSums(int.mat[c(spa,spb),]) == 2]
  Xvar <- Xvar[Xvar %in% names(anem_EN$z.mod[[reg]])]
  
  # Hosts used by either spa or spb
  Yvar <- colnames(int.mat)[colSums(int.mat[c(spa,spb),]) > 0]
  Yvar <- Yvar[Yvar %in% names(anem_EN$z.mod[[reg]])]
  
  # proportion of shared hosts over total hosts available
  p_host_shared = length(Xvar)/(length(Yvar))

  Denv <- niche_overlap(spa_en, spb_en)
  Dmut <- niche_overlap(spa_bc, spb_bc)

  if (!is.na(Dmut)) {
    Roverlap <- resource_overlap(spa_bc, spb_bc, 
                                 w.list = w_list, 
                                 Xvar = Xvar, 
                                 Yvar = Yvar, cor = FALSE)
  } else {
    Roverlap <- NA
  }
 
  # overlaps
  pairs_list$Denv[i] <- Denv
  pairs_list$Dmut[i] <- ifelse(is.na(Dmut), 0, Dmut)
  pairs_list$Roverlap[i] <- ifelse(is.na(Roverlap), 0, Roverlap)

  # host and interaction info
  pairs_list$shared_hosts[i] <- length(Xvar)
  pairs_list$p_host_shared[i] <- p_host_shared
  pairs_list$interaction[i] <- paste(sort(c(spa_behavior, spb_behavior)), collapse = "-")
  
}

pairs_list$interaction <- factor(pairs_list$interaction, levels = c("specialist-specialist", "generalist-specialist", "generalist-generalist"))
pairs_list$c_host_shared = factor(sapply(pairs_list$p_host_shared, function(i) if(i == 0){"none"} else if(i == 1){"all"} else {"some"}), levels = c("none", "some", "all"))

write.csv(pairs_list, file = "results/niche_overlaps.csv", row.names = FALSE)

## Compare estimates of niche and resource overlap
median(pairs_list$Denv, na.rm = T)
IQR(pairs_list$Denv, na.rm = T)
median(pairs_list$Dmut, na.rm = T)
IQR(pairs_list$Dmut, na.rm = T)
median(pairs_list$Roverlap, na.rm = T)
IQR(pairs_list$Roverlap, na.rm = T)

wilcox.test(pairs_list$Denv, pairs_list$Dmut, paired = T)
wilcox.test(pairs_list$Denv, pairs_list$Roverlap, paired = T)
wilcox.test(pairs_list$Dmut, pairs_list$Roverlap, paired = T)

## Check niche overlap and interaction

table(pairs_list$c_host_shared, pairs_list$interaction)
pairs_list[pairs_list$c_host_shared == "all",]

kruskal.test(pairs_list$Denv, pairs_list$interaction)
pairwise.wilcox.test(pairs_list$Denv, pairs_list$interaction)
# SS
median(pairs_list[pairs_list$interaction == "specialist-specialist", "Denv"], na.rm = T)
IQR(pairs_list[pairs_list$interaction == "specialist-specialist", "Denv"], na.rm = T)
#GG
median(pairs_list[pairs_list$interaction == "generalist-generalist", "Denv"], na.rm = T)
IQR(pairs_list[pairs_list$interaction == "generalist-generalist", "Denv"], na.rm = T)
#GS
median(pairs_list[pairs_list$interaction == "generalist-specialist", "Denv"], na.rm = T)
IQR(pairs_list[pairs_list$interaction == "generalist-specialist", "Denv"], na.rm = T)

kruskal.test(pairs_list$Dmut, pairs_list$interaction) # no differences

kruskal.test(pairs_list$Roverlap, pairs_list$interaction)
pairwise.wilcox.test(pairs_list$Roverlap, pairs_list$interaction)
# SS
median(pairs_list[pairs_list$interaction == "specialist-specialist", "Roverlap"], na.rm = T)
IQR(pairs_list[pairs_list$interaction == "specialist-specialist", "Roverlap"], na.rm = T)
#GG
median(pairs_list[pairs_list$interaction == "generalist-generalist", "Roverlap"], na.rm = T)
IQR(pairs_list[pairs_list$interaction == "generalist-generalist", "Roverlap"], na.rm = T)
#GS
median(pairs_list[pairs_list$interaction == "generalist-specialist", "Roverlap"], na.rm = T)
IQR(pairs_list[pairs_list$interaction == "generalist-specialist", "Roverlap"], na.rm = T)

wilcox.test(pairs_list[pairs_list$interaction == "generalist-generalist" & pairs_list$c_host_shared != "none", "Roverlap"],
            pairs_list[pairs_list$interaction == "specialist-specialist" & pairs_list$c_host_shared != "none", "Roverlap"])

wilcox.test(pairs_list[pairs_list$interaction == "generalist-specialist" & pairs_list$c_host_shared != "none", "Roverlap"],
            pairs_list[pairs_list$interaction == "specialist-specialist" & pairs_list$c_host_shared != "none", "Roverlap"])

wilcox.test(pairs_list[pairs_list$interaction == "generalist-specialist" & pairs_list$c_host_shared != "none", "Roverlap"],
            pairs_list[pairs_list$interaction == "generalist-generalist" & pairs_list$c_host_shared != "none", "Roverlap"])
pairwise.wilcox.test(pairs_list$Roverlap[pairs_list$c_host_shared != "none"], pairs_list$interaction[pairs_list$c_host_shared != "none"])

