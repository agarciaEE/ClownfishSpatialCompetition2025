library(raster)
library(NINA)
library(stringr)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(tidyr)
library(dplyr)

# Perform clownfish models with inferred sea anemone distribution as predictors
setwd("~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025/")

source("./scripts/functions.R")
rasterOptions(todisk = FALSE)

env <- read.csv("./data/selected_environmental_variables.csv")
clus <- read.csv("./data/marine_regions.csv")
anem_occ <- read.csv("./data/anem_occ_env_final_dataset.csv")
amph_occ <- read.csv("./data/amph_occ_env_final_dataset.csv")

anemENM <- readRDS("./Rdata/anemENMs.RDS")

# perform clownfish models with biotic as predictors
int.mat <- read.csv("./data/interaction_matrix.csv", row.names = 1)
int.mat[int.mat > 0] = 1

preds <- merge(env, anemENM$pred.dis, by = c("x", "y"))
amphENM_Host = ENmodel(amph_occ, preds, clus = clus, R = 100, extrapolate = TRUE, extrapolateIf = FALSE, FORCE = TRUE, eval = FALSE, 
                  type.pred = "probability", crs = 4326, res = 1)

amphENM_Host$pred.dis[is.na(amphENM_Host$pred.dis)] = 0
amphENM_Host$eval <- models_evaluation(amphENM_Host, sample.pseudoabsences = TRUE, 
                                  alpha = 0.5, transformation = "none", trans.param = 20, 
                                  plot = FALSE, rep = 100, th = NULL,  best.th = "TSS")

saveRDS(amphENM_Host, "./Rdata/amphENMs_envhostpredictors.RDS")

amphENM <- readRDS("./Rdata/amphENMs.RDS")
amphEBM <- readRDS("./Rdata/amphEBMs.RDS")

# Prepare data
df <- data.frame(
  species = rownames(amphENM_Host$eval$tab),
  ENM      = amphENM$eval$tab[rownames(amphENM_Host$eval$tab), "AUC"],
  ENM_host = amphENM_Host$eval$tab[rownames(amphENM_Host$eval$tab), "AUC"],
  EBM      = amphEBM$eval$tab[rownames(amphENM_Host$eval$tab), "AUC"]
)

# Long format
df_long <- df %>%
  pivot_longer(cols = c(ENM, ENM_host, EBM),
               names_to = "model",
               values_to = "AUC") %>%
  mutate(model = factor(recode(model,
                        "ENM" = "ENM",
                        "ENM_host" = "ENM_Host",
                        "EBM" = "MR_ENM"), levels = c("ENM", "ENM_Host", "MR_ENM")))

# Plot
ggplot(df_long, aes(x = model, y = AUC, fill = model, color = model)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "white", alpha = 0.9,
               outlier.shape = NA, fatten = 3, lwd = 1) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  #scale_y_continuous(limits = c(0.5, 1.01)) +
  scale_x_discrete(labels = c("Environmental predictors \n ENM",
                              "Env + Host predictors \n ENM",
                              "Mutualism-refined \n ENM")) +
  labs(x = "", y = "AUC", title = "Performance of Ecological Models") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 16),
        axis.text = element_text(color = "black")) +
  stat_compare_means(
    comparisons = list(
      c("ENM", "ENM_Host"),
      c("ENM_Host", "MR_ENM"),
      c("ENM", "MR_ENM")
    ),
    paired = TRUE, step.increase = 0.075,
    label = "p.signif", hide.ns = FALSE, size= 10,
    method.args = list(exact = FALSE)
  )
ggsave("./figures/AUCmodels_comparison.pdf", width = 10, height = 9)

pdf("./figures/host_pred_contrbib.pdf", width = 10, height = 7)
ecospat.plot.contrib(amphENM_Host$pca$co, amphENM_Host$pca$eig)
dev.off()

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
  
  D <- 1 - (0.5 * (sum(abs(p1 - p2), na.rm = TRUE)))
  
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
amph_BC <- readRDS("./Rdata/amphENMs_hostpredictors.RDS")
amph_BC <- amphENM_Host
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
   
    if (length(hosts[hosts %in% Xvar]) == 0) next
    
    en_niche <- amph_EN$z.mod[[reg]][[sp]]
    bc_niche <- amph_BC$z.mod[[reg]][[sp]]
    bc_niche$w <- raster::rasterize(amph_EN$env.scores[,3:4], en_niche$w, 
                      field = raster::extract(bc_niche$w, amph_BC$env.scores[,3:4]), 
                      fun = mean, na.rm = TRUE)
    bc_niche$w[bc_niche$w > 0] = 1
    bc_niche$z.uncor <- raster::rasterize(amph_EN$env.scores[,3:4], en_niche$w, 
                                    field = raster::extract(bc_niche$z.uncor, amph_BC$env.scores[,3:4]), 
                                    fun = mean, na.rm = TRUE)
    bc_niche$Z <- en_niche$Z
    w_niche_all <- anem_EN$z.mod[[reg]][hosts[hosts %in% Xvar]]

    w_niche <- w_niche_all[[1]]
    w_niche$w <- sum(stack(sapply(w_niche_all, function(i) i$w)), na.rm = TRUE)
    w_niche$w[w_niche$w > 0] = 1

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
sum(apply(amphENM_Host$pred.dis, 2, function(i) sum(i)) == 0)

dir.create("results", showWarnings = FALSE)
ROU_data$Unexploited_spatial[!is.finite(ROU_data$Unexploited_spatial)] = 0
write.csv(ROU_data, file = "./results/ROU_data_envhost_predictors.csv", row.names = FALSE)

# Load ROU estimates
data <- read.csv("./results/ROU_data_envhost_predictors.csv")
data <- data[, c("species", "Restricted_niche", "Unexploited_niche", "Occupied_niche", "behavior_reg")]
colnames(data) <- c("species", "Restricted", "Unexploited", "Occupied", "behavior")

# Load interactions matrix to order species in based to number of interactions
int.mat = read.csv("data/interaction_matrix.csv", row.names = 1)
int.mat[int.mat > 0] = 1

mean(ROU_data$Occupied_niche)
mean(ROU_data$Restricted_niche)
mean(ROU_data$Unexploited_niche)

mean(ROU_data$Occupied_spatial)
mean(ROU_data$Restricted_spatial)
mean(ROU_data$Unexploited_spatial, na.rm = TRUE)

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# compute relative proportions to integrate to 1
data$Restricted <- as.numeric(data$Restricted)
data$Occupied <- as.numeric(data$Occupied)
data$Unexploited <- as.numeric(data$Unexploited)

rw = data$Restricted + data$Unexploited + data$Occupied
data$Occupied = data$Occupied/rw
data$Unexploited = data$Unexploited/rw
data$Restricted = data$Restricted/rw
data$Occupied = data$Restricted + data$Occupied
data$Unexploited = data$Unexploited + data$Occupied

speciesLevels <- names(sort(rowSums(int.mat > 0 )))
speciesLevels <- speciesLevels[speciesLevels %in% unique(data$species)]

# Convert to table of data parameters and values
data <- gather(data, parameter, value, 2:4)
data$parameter <- factor(data$parameter, levels = rev(c("Restricted", "Occupied", "Unexploited")))
data$species = factor(gsub("Premnas", "Amphiprion", data$species), levels = gsub("Premnas", "Amphiprion", speciesLevels))

# plot
ROUcols <- c("#E69F00", "#2171B5", "#27AD81FF")
behavior_cols <- c("generalist" = "#1F3B73",  # Deep Blue  
                   "specialist" = "#F28E80")  # Soft Coral  

# Summarize the most predominant behavior for each species
behavior_data <- data %>%
  group_by(species, behavior) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  group_by(species) %>%
  slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
  dplyr::select(species, behavior)

behavior_data$nhost <- rowSums(int.mat[gsub("Amphiprion_biaculeatus", "Premnas_biaculeatus", as.character(behavior_data$species)),])

# (a) species average proportions across regions
fig2a <- ggplot(data, aes(x = species, y = value, fill = parameter)) +
  stat_summary(aes(fill = parameter), fun.data = "mean_sdl", fun.args = list(mult = 1), 
               geom = "bar", position = "identity") +
  stat_summary(data = data[data$parameter == "Restricted", ], fun.data = "mean_sdl", 
               col = alpha("white", 0.5), fun.args = list(mult = 1), geom = "pointrange", 
               position = "identity", size = 1, shape = 45) +
  stat_summary(data = data[data$parameter == "Occupied", ], fun.data = "mean_sdl", 
               col = alpha("white", 0.5), fun.args = list(mult = 1), geom = "pointrange", 
               position = "identity", size = 1, shape = 45) +
  scale_fill_manual(
    name = "",
    labels = c("Restricted", "Occupied", "Unexploited"),
    values = ROUcols, 
    limits = c("Restricted", "Occupied", "Unexploited") 
  ) +
  scale_x_discrete(labels = gsub("_", " ", levels(data$species)), expand = c(0, 0)) +
  scale_y_continuous("ROU proportions", breaks = seq(0, 1, 0.25), limits = c(0, 1.1), expand = c(0, 0)) +
  theme_classic() +
  coord_flip() +
  theme(axis.title.x = element_text(color = "black", size = 20, vjust = 0),
        axis.title.y = element_text(color = "black", size = 20, vjust = 2),
        axis.text.x = element_text(color = "black", size = 18),
        axis.text.y = element_text(color = "black", size = 18, face = "italic"),
        axis.line.x = element_blank(),
        plot.margin = unit(c(0, 2, 1, 2), "lines"),
        legend.position = "top",
        legend.text = element_text(size = 16),
        axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid")) +
  geom_text(data = behavior_data[, c("species", "nhost")], aes(x = species, y = 1.05, label = nhost, fill = NULL), 
            hjust = 0, vjust = 0.5, size = 6) +
  geom_point(data = behavior_data, aes(x = species, y = 1.025, color = behavior), 
             shape = 19, size = 8, fill = "black") +
  scale_color_manual(values = behavior_cols,
                     guide = "none")

# (b) generalist versus specialist comparisons
dodge <- position_dodge(width = 0.4)
symnum.args <- list(cutpoints = c(0,0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "n.s."))
data <- ROU_data
data <- read.csv("./results/ROU_data_host_predictors.csv")
data <- data[, c("species", "behavior_reg", "behavior", "Restricted_niche", "Unexploited_niche", "Occupied_niche")]
colnames(data) <- c("species", "behavior_reg", "behavior", "Restricted", "Unexploited", "Occupied")

## Check G/S differences
#####
## plot

RN_GS <- ggplot(data, aes(x = behavior_reg, y = Restricted, fill = behavior_reg, col = behavior_reg)) +
  scale_y_continuous(name = "Restricted", expand = c(0,0), breaks = seq(0,1,0.25), limits = c(-0.1,1.2))+
  scale_x_discrete(expand = c(0,0)) +
  geom_violin() +
  geom_boxplot(width = 0.02, outlier.colour = NA, col = darken(behavior_cols,0.5), fill = NA) +
  stat_summary(fun=median, geom="point", size=2, shape = 19, position = dodge, col = "white") +
  scale_fill_manual(values= behavior_cols) +
  scale_color_manual(values= behavior_cols) +
  theme_classic() +
  coord_cartesian(xlim=c(0.5,2.5), ylim=c(-0.1,1.2)) +
  annotate(x=0.5, xend=0.5, y=0, yend=1, lwd = 1, colour="black", geom="segment") +
  annotate(x=1, xend=2, y=-0.1, yend=-0.1, lwd = 1, colour="black", geom="segment") +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks = element_line(size = 0.75),
        axis.ticks.length=unit(.2, "cm"),
        plot.margin = unit(c(2,2,2,2), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 14, vjust = 5, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12,  vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none",
        legend.text = element_text(color = "black", size = 12, family = "mono" ),
        axis.line = element_blank()) +
  stat_compare_means(comparisons = list(c("generalist", "specialist")), hide.ns = T, size = 8,  symnum.args = symnum.args)

ON_GS <- ggplot(data, aes(x = behavior_reg, y = Occupied, fill = behavior_reg, col = behavior_reg)) +
  scale_y_continuous(name = "Occupied", expand = c(0,0), breaks = seq(0,1,0.25), limits = c(-0.1,1.2))+
  scale_x_discrete(expand = c(0,0)) +
  geom_violin() +
  geom_boxplot(width = 0.02, outlier.colour = NA, col = darken(behavior_cols,0.5), fill = NA) +
  stat_summary(fun=median, geom="point", size=2, shape = 19, position = dodge, col = "white") +
  scale_fill_manual(values= behavior_cols) +
  scale_color_manual(values= behavior_cols) +
  theme_classic() +
  coord_cartesian(xlim=c(0.5,2.5), ylim=c(-0.1,1.2)) +
  annotate(x=0.5, xend=0.5, y=0, yend=1, lwd = 1, colour="black", geom="segment") +
  annotate(x=1, xend=2, y=-0.1, yend=-0.1, lwd = 1, colour="black", geom="segment") +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks = element_line(size = 0.75),
        axis.ticks.length=unit(.2, "cm"),
        plot.margin = unit(c(2,2,2,2), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 14, vjust = 5, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12,  vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none",
        legend.text = element_text(color = "black", size = 12, family = "mono" ),
        axis.line = element_blank()) +
  stat_compare_means(comparisons = list(c("generalist", "specialist")), hide.ns = T, size = 8,  symnum.args = symnum.args)

UN_GS <- ggplot(data, aes(x = behavior_reg, y = Unexploited, fill = behavior_reg, col = behavior_reg)) +
  scale_y_continuous(name = "Unexploited", expand = c(0,0), breaks = seq(0,1,0.25), limits = c(-0.1,1.2))+
  scale_x_discrete(expand = c(0,0)) +
  geom_violin() +
  geom_boxplot(width = 0.02, outlier.colour = NA, col = darken(behavior_cols,0.5), fill = NA) +
  stat_summary(fun=median, geom="point", size=2, shape = 19, position = dodge, col = "white") +
  scale_fill_manual(values= behavior_cols) +
  scale_color_manual(values= behavior_cols) +
  theme_classic() +
  coord_cartesian(xlim=c(0.5,2.5), ylim=c(-0.1,1.2)) +
  annotate(x=0.5, xend=0.5, y=0, yend=1, lwd = 1, colour="black", geom="segment") +
  annotate(x=1, xend=2, y=-0.1, yend=-0.1, lwd = 1, colour="black", geom="segment") +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks = element_line(size = 0.75),
        axis.ticks.length=unit(.2, "cm"),
        plot.background = element_rect(fill = "white", color = "white"),
        plot.margin = unit(c(2,2,2,2), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 14, vjust = 5, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12,  vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none",
        legend.text = element_text(color = "black", size = 12, family = "mono" ),
        axis.line = element_blank()) +
  stat_compare_means(comparisons = list(c("generalist", "specialist")), hide.ns = T, size = 8,  symnum.args = symnum.args)

fig2bd <- ggarrange(RN_GS, ON_GS, UN_GS,
                    ncol = 1, nrow = 3, labels = paste0("(", letters[2:4], ")"), font.label = list(size = 20))

fig2 <- ggarrange(fig2a, fig2bd, ncol = 2, labels = paste0("(", letters[1], ")"), widths = c(1.75,1), font.label = list(size = 20))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/SuppFig2.speciesROU_and_ROU_vs_behavior_envhostpredictors.pdf", plot = fig2, width=18, height=12, units="in")


