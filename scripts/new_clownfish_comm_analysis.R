# Load required packages
library(ape)
library(phytools)
library(geiger)
library(vegan)
library(bipartite)
library(BioGeoBEARS)
library(ENMeval)
library(tidyr)
library(dplyr)
library(vegan)
library(terra)
library(rJava)
library(dismo)
library(ade4)
library(ecospat)
library(nlme)
library(ggpubr)
library(cowplot)
library(ggplot2)

setwd("~/Unil/Research/PNC_clownfish/")

# Data preparation ############################################################

# environmental data
env_data <- read.csv("~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025/data/env_selected_dataset.csv")
env_stack <- rast(env_data)
coords <- env_data[,1:2]

# marine regions
marine_regions <- read.csv("~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025/data/marine_regions.csv")

# occurrence data
clownfish_occ <- read.csv("~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025/data/amph_occ_env_final_dataset.csv")
anemone_occ <- read.csv("~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025/data/anem_occ_env_final_dataset.csv")

# interaction matrix
int.mat <- read.csv("~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025/data/interaction_matrix_new_prop.csv", row.names = 1)
int.mat_bin <- int.mat
int.mat_bin[int.mat_bin > 0] = 1

# behavior dataset
behavior_df <- read.csv("~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025/data/amph_behavior_new.csv")

# host categories
host_categories <- setNames(behavior_df$category, behavior_df$species)
host_categories <- as.factor(host_categories)

# clownfish phylogeny
clownfish_tree <- read.tree("~/Unil/Research/1_NicheCompetition/calibrated_tree.tre")
clownfish_tree$tip.label <- gsub("A.", "Amphiprion", clownfish_tree$tip.label)
clownfish_tree$tip.label[clownfish_tree$tip.label == "Amphiprion_biaculeatus"] = "Premnas_biaculeatus"

# add environmental variables to species occurrences
clownfish_occ <- na.exclude(cbind(clownfish_occ, terra::extract(env_stack, clownfish_occ[,1:2])[,-1]))
anemone_occ <- na.exclude(cbind(anemone_occ, terra::extract(env_stack, anemone_occ[,1:2])[,-1]))

# get unique species
clownfish_species <- unique(clownfish_occ$species)
anemone_species <- unique(anemone_occ$species)

# Species Distribution Modeling ###############################################



# Environmental Niche Modeling ################################################

# perform PCA
env_pca <- dudi.pca(env_data[,-c(1:2)], center = TRUE, scale = TRUE, nf = 2, scannf = FALSE)
env_scores <- cbind(coords, env_pca$li)

# get species PCA scores
clownfish_scores <- cbind(suprow(env_pca, data.frame(clownfish_occ)[, colnames(env_data[,-c(1:2)])])$li, species = clownfish_occ$species) 
anemone_scores <- cbind(suprow(env_pca, data.frame(anemone_occ)[, colnames(env_data[,-c(1:2)])])$li, species = anemone_occ$species) 

plot(env_scores[,3:4])
points(clownfish_scores[,1:2], pch = 19, col = "lightblue")
points(anemone_scores[,1:2], pch = 19, col = "pink")

global_species_scores_ext <- c(range(c(clownfish_scores$Axis1, anemone_scores$Axis1)),
                               range(c(clownfish_scores$Axis2, anemone_scores$Axis2)))
# subset env_scores to remove scores far from the occurrences and narrow down the environmental space to their occupancy (increase resolution)
kepp_rows <- which(env_scores$Axis1 >= (global_species_scores_ext[1] * 1.5) & 
                   env_scores$Axis1 <= (global_species_scores_ext[2] * 1.5) & 
                   env_scores$Axis2 >= (global_species_scores_ext[3] * 1.5) & 
                   env_scores$Axis2 <= (global_species_scores_ext[4] * 1.5))

env_scores <- env_scores[kepp_rows, ] 
env_data <- env_data[kepp_rows, ] 
coords <- coords[kepp_rows, ] 
env_scores <- merge(env_scores, marine_regions, by = c("x", "y"), all.x = TRUE)

plot(env_scores[,3:4])
points(clownfish_scores[,1:2], pch = 19, col = "lightblue")
points(anemone_scores[,1:2], pch = 19, col = "pink")

## GLOBAL ####################

# combine clownfish and sea anemone datasets
all_occ <- rbind(clownfish_occ, anemone_occ)
all_scores <- rbind(clownfish_scores, anemone_scores)
all_species <- c(clownfish_species, anemone_species)

R = 100
ENMs <- list()
for (sp in all_species){ 
  
  sp_scores <- na.exclude(all_scores[all_scores$species == sp,1:2])
  
  if (nrow(sp_scores) > 4){
    ENMs[[sp]] <- ecospat.grid.clim.dyn(glob = env_scores[,3:4], 
                                        glob1 = env_scores[,3:4],
                                        sp = sp_scores[,1:2],
                                        R = R)
  } else {
    message(sp, " failed. At least 5 locations are required to fit an home range.
            Since failed species are highly endemic with very narrow ranges in relation to the geographical resolution,
            locations will be duplicated to be able to fit the kernel distribution. Porbabilities will be normalized to [0,1],
            so this tweak would not have a strong effect on the inference.")
    
    ENMs[[sp]] <- ecospat.grid.clim.dyn(glob = env_scores[,3:4], 
                                        glob1 = env_scores[,3:4],
                                        sp = rbind(sp_scores[,1:2],sp_scores[,1:2]),
                                        R = R)
  }
}

clownfish_ENMs <- ENMs[names(ENMs) %in% clownfish_species]
anemone_ENMs <- ENMs[names(ENMs) %in% anemone_species]

clownfish_species <- names(clownfish_ENMs)
anemone_species <- names(anemone_ENMs)

## REGIONAL MODELS ####################
R = 100
subset_environments = TRUE 
expand_factor = 0.2
# FALSE Considering all species inhabited environments for each region even if they are not present to estimate the species niche usage in the region,
# TRUE would consider the specific niche of the species for the subset of the environments present in the region disregarding the overall environmental suitability
regions <- unique(env_scores$province)
pred.dis <- data.frame()
regENMs <- list()
for (reg in regions) {
  
  regENMs[[reg]] <- list()

  reg_scores <- na.exclude(env_scores[env_scores$province == reg,1:5])
  reg_scores_ext <- c(range(reg_scores$Axis1), range(reg_scores$Axis2))
  reg_coords_ext <- c(range(reg_scores$x), range(reg_scores$y))
  
  # Expand the the niche extent
  x_range <- (reg_scores_ext[2] - reg_scores_ext[1]) * expand_factor
  y_range <- (reg_scores_ext[4] - reg_scores_ext[3]) * expand_factor
  expanded_scores_ext <- c(reg_scores_ext[1] - x_range,
                           reg_scores_ext[2] + x_range,
                           reg_scores_ext[3] - y_range,
                           reg_scores_ext[4] + y_range)
  
  reg_species <- unique(all_occ %>%
    dplyr::filter(x >= reg_coords_ext[1], 
                  x <= reg_coords_ext[2],
                  y >= reg_coords_ext[3],
                  y <= reg_coords_ext[4]) %>%
    pull(species))
  
  for (sp in reg_species){ 
    
    sp_scores <- na.exclude(all_scores[all_scores$species == sp,1:2])
   
    # subset species environments to those of the region?
    if (subset_environments & sp != "Amphiprion_mccullochi") {
      sp_scores <- sp_scores %>%
        filter(Axis1 >= expanded_scores_ext[1], 
               Axis1 <= expanded_scores_ext[2],
               Axis2 >= expanded_scores_ext[3],
               Axis2 <= expanded_scores_ext[4])
    }
  
    if (nrow(sp_scores) > 4){
      regENMs[[reg]][[sp]] <- ecospat.grid.clim.dyn(glob = reg_scores[,3:4], 
                                                    glob1 = reg_scores[,3:4],
                                                    sp = sp_scores[,1:2],
                                                    R = R, th.sp = NULL, th.env = NULL,
                                                    extend.extent = expanded_scores_ext - reg_scores_ext)
      ecospat.plot.niche(regENMs[[reg]][[sp]])
      
    } else if (nrow(sp_scores) > 2) {
      message(sp, " failed in ", reg, ". At least 5 locations are required to fit an home range.
            Since failed species are highly endemic with very narrow ranges in relation to the geographical resolution,
            locations will be duplicated to be able to fit the kernel distribution. Porbabilities will be normalized to [0,1],
            so this tweak would not have a strong effect on the inference.")
      
      z <- ecospat.grid.clim.dyn(glob = reg_scores[,3:4], 
                                 glob1 = reg_scores[,3:4],
                                 sp = rbind(sp_scores[,1:2],sp_scores[,1:2]),
                                 R = R, th.sp = NULL, th.env = NULL,
                                 extend.extent =  expanded_scores_ext - reg_scores_ext)
      
      z$sp <- sp_scores[,1:2]
      z$z <- z$z / 2
      z$z.cor <- z$z / z$Z
      regENMs[[reg]][[sp]] <- z
      
      ecospat.plot.niche(regENMs[[reg]][[sp]])
      
    } else {
      message(sp, " failed in ", reg, ". At least 5 locations are required to fit an home range.")
    }
    
  }
  
  # get projected distributions of environmental niches
  distributions <- cbind(reg_scores[,1:2], sapply(all_species, function(sp) {
    
    if (sp %in% names(regENMs[[reg]])) {
      # extract niche densities from env space
      raster::extract(regENMs[[reg]][[sp]]$z, 
                      reg_scores[,3:4])[,2]
    } else {
      0
    }
  }, simplify = FALSE))
  
  pred.dis <- rbind(pred.dis, distributions)
}

# split niche models between clownfish and anemones
clownfish_regENMs <- sapply(regENMs, function(reg) reg[names(reg) %in% clownfish_species])
clownfish_regENMs <- clownfish_regENMs[sapply(clownfish_regENMs, length) > 0]

anemone_regENMs <- sapply(regENMs, function(reg) reg[names(reg) %in% anemone_species])
anemone_regENMs <- anemone_regENMs[sapply(anemone_regENMs, length) > 0]

## Species pairwise niche overlaps #############################################

# niche overlap functions #########################

niche_overlap <- function(x, y, cor = FALSE) {
  
  if (is.null(x) | is.null(y)) {
    return(NA)
  }

  # Check if maximum values are zero
  if (terra::global(x$z, "max", na.rm = TRUE)$max == 0 | terra::global(y$z, "max", na.rm = TRUE)$max == 0) {
    return(NA)
  }
  
  # Check if rasters have the same properties
  if (!identical(x$Z, y$Z)) {
    stop("Rasters of different extent.")
  }
  
  
  # If 'cor' is TRUE, process the 'z.cor' data
  if (cor) {
    p1 <- t(terra::as.matrix(x$z.cor) / sum(terra::as.matrix(x$z.cor), na.rm = TRUE))
    p2 <- t(terra::as.matrix(y$z.cor) / sum(terra::as.matrix(y$z.cor), na.rm = TRUE))
  } else {
    # If 'cor' is FALSE, process the 'z.uncor' data
    p1 <- t(terra::as.matrix(x$z.uncor) / sum(terra::as.matrix(x$z.uncor), na.rm = TRUE))
    p2 <- t(terra::as.matrix(y$z.uncor) / sum(terra::as.matrix(y$z.uncor), na.rm = TRUE))
  }
  
  # Calculate D
  D <- 1 - (0.5 * (sum(abs(p1 - p2))))
  
  return(D)
  
}

resource_overlap <- function(x, y, w.list, Xvar, Yvar) {
  
  if (is.null(x) | is.null(y)) {
    return(NA)
  }

  # Check if maximum values are zero
  if (terra::global(x$z, "max", na.rm = TRUE)$max == 0 | terra::global(y$z, "max", na.rm = TRUE)$max == 0) {
    return(NA)
  }
  
  # Check if rasters have the same properties
  if (!identical(x$Z, y$Z)) {
    stop("Rasters of different extent.")
  }
  
  ecoD <- c()
  for (var in Yvar){
    if (var %in% Xvar){
      
      # Refine niche layers      
      z1 <- x$z.uncor * w.list[[var]]$z.uncor
      z2 <- y$z.uncor * w.list[[var]]$z.uncor
      
      # integrate values to 1
      p1 <- terra::as.matrix(z1)/sum(terra::as.matrix(z1), na.rm = TRUE)
      p2 <- terra::as.matrix(z2)/sum(terra::as.matrix(z2), na.rm = TRUE)
      
      # compute D
      D <- 1 - (0.5 * (sum(abs(p1 - p2))))
      ecoD <- c(ecoD, D)
    }
    else{
      ecoD <- c(ecoD, 0)
    }
  }

  if (length(ecoD) > 0 & any(!is.na(ecoD))) {
    ecoD <- mean(ecoD, na.rm = TRUE)
  }
  else {
    ecoD <- NA
  }
  
  return(ecoD)
}

# -------------------------------------------------

## GLOBAL ####################

pairs_list <- expand.grid(spa = clownfish_species, 
                          spb = clownfish_species)
pairs_list <- pairs_list[pairs_list$spa != pairs_list$spb,]

for (i in 1:nrow(pairs_list)) {
  
  spa <- as.character(pairs_list$spa[i])
  spb <- as.character(pairs_list$spb[i])
  
  # species behavior
  spa_behavior <- behavior_df$category[behavior_df$species == spa]
  spb_behavior <- behavior_df$category[behavior_df$species == spb]
  
  # environmental niches
  spa_en <- clownfish_ENMs[[spa]]
  spb_en <- clownfish_ENMs[[spb]]
  
  # host niche lists
  w_list <- anemone_ENMs
  
  # Shared hosts between spA and spB
  Xvar <- colnames(int.mat_bin)[colSums(int.mat_bin[c(spa,spb),]) == 2]
  Xvar <- Xvar[Xvar %in% names(w_list)]
  
  # Hosts used by either spa or spb
  Yvar <- colnames(int.mat_bin)[colSums(int.mat_bin[c(spa,spb),]) > 0]
  Yvar <- Yvar[Yvar %in% names(w_list)]
  
  # proportion of shared hosts over total hosts available
  p_host_shared = length(Xvar)/(length(Yvar))
  
  Denv <- niche_overlap(spa_en, spb_en)

  Roverlap <- resource_overlap(spa_en, spb_en, 
                                 w.list = w_list, 
                                 Xvar = Xvar, 
                                 Yvar = Yvar)
  # overlaps
  pairs_list$Denv[i] <- Denv
  pairs_list$Roverlap[i] <- Roverlap
  
  # host and interaction info
  pairs_list$shared_hosts[i] <- length(Xvar)
  pairs_list$p_host_shared[i] <- p_host_shared
  pairs_list$c_host_shared[i] <- ifelse(p_host_shared == 0, "none", ifelse(i == 1, "all", "some"))
  pairs_list$interaction[i] <- paste(sort(c(spa_behavior, spb_behavior)), collapse = "-")
  
}

pairs_list$interaction <- factor(pairs_list$interaction, levels = unique(pairs_list$interaction))
pairs_list$c_host_shared <- factor(pairs_list$c_host_shared, levels = c("none", "some", "all"))

## PLOT

pairs_list$pair = 1:nrow(pairs_list)

median(pairs_list$Denv, na.rm = T)
IQR(pairs_list$Denv, na.rm = T)

median(pairs_list$Roverlap, na.rm = T)
IQR(pairs_list$Roverlap, na.rm = T)

wilcox.test(pairs_list$Denv, pairs_list$Roverlap, paired = T)

dodge <- position_dodge(width = 0.4)
symnum.args <- list(cutpoints = c(0,0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "n.s."))
colors <- c("#E69F00",  "#2171B5", "#27AD81FF")

## select cuts
bins <- 20
cols <- c("darkred", "red", "white", "blue",  "darkblue")
colGradient <- colorRampPalette(cols)
cut.cols <- colGradient(bins)
cuts = cut(seq(-1,1,0.1),20)[-1]
names(cuts) <- sapply(cuts,function(t) cut.cols[which(as.character(t) == levels(cuts))])
breaks = seq(-1,1,0.1)

## env vs resource-use changes
env_res.hist = data.frame(diff = rowSums(cbind(pairs_list$Denv,-pairs_list$Roverlap), na.rm = T), pair = pairs_list$pair)
env_res.hist[env_res.hist$diff>0, "group"] = "ENV"
env_res.hist[env_res.hist$diff<0, "group"] = "RES"
group_tags <- cut(env_res.hist$diff,
                  breaks=breaks,
                  include.lowest=TRUE,
                  right=FALSE,
                  labels=cuts)
env_res.hist$tag = group_tags
env_res.hist$tag = factor(env_res.hist$tag, levels = sort(unique(env_res.hist$tag)))

env_res_hist <- ggplot(na.exclude(env_res.hist), aes(x=diff, fill=tag, col = tag)) +
  geom_histogram(breaks=seq(-1, 1, by=0.1),
                 alpha=0.8,
                 position="identity") +
  annotate(x=0, xend=0, y=0, yend=150, lwd = 0.5, linetype = 2,  colour="black", geom="segment") +
  geom_histogram(data = env_res.hist[is.na(env_res.hist$group),], breaks=seq(-1.05, 1.05, by=0.1),
                 alpha=0.8,
                 position="identity", col = "black", fill = "black") +
  geom_hline(yintercept = 0, linetype = 1, linewidth = 2, col = "white") +
  annotate(x=median(env_res.hist$diff), xend=median(env_res.hist$diff), y=0, yend=150, lwd = 1,linetype = 1,  colour=if(median(env_res.hist$diff)<0){"red"}else{"blue"}, geom="segment") +
  scale_y_continuous("count", limits = c(-20,150) , expand = c(0,0)) +
  scale_x_continuous(name = "difference", labels = c(1,0.5,0,0.5,1), breaks = seq(-1,1,0.5), expand = c(0,0), limits = c(-1.4,1.4)) +
  scale_fill_manual(values = names(cuts)[cuts %in% env_res.hist$tag]) +
  scale_color_manual(values =  names(cuts)[cuts %in% env_res.hist$tag]) +
  theme_classic() +
  coord_cartesian(xlim=c(-1.4,1.4), ylim=c(-20,170)) +
  annotate(x=-1.4, xend=-1.4, y=0, yend=150, lwd = 1, colour="black", geom="segment") +
  annotate(x=-1, xend=1, y=-20, yend=-20, lwd = 1, colour="black",  geom="segment") +
  theme(panel.border=element_blank(),
        panel.background = element_rect(fill ="transparent", color = NA),
        plot.background = element_rect(fill ="transparent", color = NA),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(3,1,3,1), "cm"),
        axis.title.x = element_text(color = "black", size = 20, vjust = 0),
        legend.title = element_text(color = "black", face = "bold"),
        legend.text = element_text(color = "black"),
        axis.title.y = element_text(color = "black",  size = 20, vjust = 1, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 18,  vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 18 ),
        axis.ticks = element_line(linewidth = 0.75),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "none",
        axis.line = element_blank())

D_gather <- gather(pairs_list, type, D, 3:4)
D_gather$type <- factor(D_gather$type, levels = c("Denv", "Roverlap"))

## ecological to resource use
niche_overlap <- ggplot(D_gather, aes(x = type, y = D, fill = type, col = type)) +
  geom_violin(col = NA) +
  geom_rect(data=NULL,aes(xmin=1,xmax=2,ymin=-Inf,ymax=Inf),
            fill="white" , col = NA, alpha = 0.5) +
  geom_line(data = D_gather[pairs_list$pair %in% env_res.hist[which(env_res.hist$diff >0), "pair"],],
            aes(group = pair), alpha = 0.2, col = "blue") +
  geom_line(data = D_gather[pairs_list$pair %in% env_res.hist[which(env_res.hist$diff <0), "pair"],],
            aes(group = pair), alpha = 0.6, col = "red") +
  geom_line(data = D_gather[pairs_list$pair %in% env_res.hist[which(env_res.hist$diff ==0), "pair"],],
            aes(group = pair), alpha = 0.2, col = "grey40") +
  geom_boxplot(width = 0.05, outlier.colour = NA, outlier.size = 0.8,
               col = "white") +
  stat_summary(fun=median, geom="point", size=3, col = "white", shape = 19, position = dodge) +
  theme_classic() +
  scale_fill_manual("", labels = c("Ecological", "Resource-use"), values= rev(as.vector(colors))[c(1,3)]) +
  scale_color_manual("", labels = c("Ecological", "Resource-use"), values= rev(as.vector(colors))[c(1,3)]) +
  scale_y_continuous(name = "Schoener's D", expand = c(0,0),  seq(0,1,0.25), limits = c(-0.1,1.2)) +
  scale_x_discrete(labels = c("Ecological\nniche overlap", "Resource-use \noverlap"), expand = c(0,0)) +
  coord_cartesian(xlim=c(0.5,2.5), ylim=c(-0.1,1.2)) +
  annotate(x=0.5, xend=0.5, y=0, yend=1, lwd = 1, colour="black", geom="segment") +
  annotate(x=1, xend=2, y=-0.1, yend=-0.1, lwd = 1, colour="black",  geom="segment") +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks = element_line(linewidth = 0.75),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "none",
        axis.line = element_blank(),
        plot.margin = margin(1, 6, 1, 1, "cm"),
        title = element_blank(),
        legend.text = element_text( color = "black", size=18),
        legend.title = element_text(color = "black", size = 18, vjust = 3),
        axis.title.x = element_blank(),
        axis.title.y = element_text( color = "black", size=24, vjust = 5),
        axis.text.x = element_text( color = "black", size=20, vjust = 0),
        axis.text.y = element_text( color = "black", size=20),
        ) +
  stat_compare_means(comparisons = list(c("Denv", "Roverlap")), paired = T, size  = 12, method = "wilcox.test", label = "..p.signif..", label.y =0.95, symnum.args = symnum.args) 

# Create the combined plot with an offsite, larger inset
fig_overlap <- ggdraw() +
  draw_plot(niche_overlap, x = 0, y = 0, width = 1, height = 1) +  # Main plot
  draw_plot(env_res_hist, x = 0.625, y = 0.55, width = 0.4, height = 0.4)  # Larger and shifted inset

dir.create("figures")
ggsave(paste0("figures/species_niche_overlaps_global.pdf"), plot = fig_overlap, width=12, height= 10, units="in")

## REGIONAL ####################

reg_pairs_list <- plyr::ldply(
  sapply(clownfish_regENMs, function(reg) {
    if (length(reg) > 1) {
      t(data.frame(combn(names(reg),2), 
                   row.names = c("spa", "spb")))
    } else {
      NULL
    }
  }),
  .id = "region"
)

for (i in 1:nrow(reg_pairs_list)) {
  
  reg <- as.character(reg_pairs_list$region[i])
  spa <- as.character(reg_pairs_list$spa[i])
  spb <- as.character(reg_pairs_list$spb[i])
  
  # species behavior
  spa_behavior <- behavior_df$category[behavior_df$species == spa]
  spb_behavior <- behavior_df$category[behavior_df$species == spb]
  
  # environmental niches
  spa_en <- clownfish_regENMs[[reg]][[spa]]
  spb_en <- clownfish_regENMs[[reg]][[spb]]
  
  # host niche lists
  w_list <- anemone_regENMs[[reg]]

  # Shared hosts between spA and spB
  Xvar <- colnames(int.mat_bin)[colSums(int.mat_bin[c(spa,spb),]) == 2]
  Xvar <- Xvar[Xvar %in% names(w_list)]
  
  # Hosts used by either spa or spb
  Yvar <- colnames(int.mat_bin)[colSums(int.mat_bin[c(spa,spb),]) > 0]
  Yvar <- Yvar[Yvar %in% names(w_list)]
  
  # proportion of shared hosts over total hosts available
  p_host_shared = length(Xvar)/(length(Yvar))
  
  Denv <- niche_overlap(spa_en, spb_en)
  
  Roverlap <- resource_overlap(spa_en, spb_en, 
                               w.list = w_list, 
                               Xvar = Xvar, 
                               Yvar = Yvar)
  # overlaps
  reg_pairs_list$Denv[i] <- Denv
  reg_pairs_list$Roverlap[i] <- Roverlap
  
  # host and interaction info
  reg_pairs_list$shared_hosts[i] <- length(Xvar)
  reg_pairs_list$p_host_shared[i] <- p_host_shared
  reg_pairs_list$c_host_shared[i] <- ifelse(p_host_shared == 0, "none", ifelse(i == 1, "all", "some"))
  reg_pairs_list$interaction[i] <- paste(sort(c(spa_behavior, spb_behavior)), collapse = "-")
  
}

reg_pairs_list$interaction <- factor(reg_pairs_list$interaction, levels = unique(pairs_list$interaction))
reg_pairs_list$c_host_shared <- factor(reg_pairs_list$c_host_shared, levels = c("none", "some", "all"))

## PLOTS

reg_pairs_list$pair = 1:nrow(reg_pairs_list)

median(reg_pairs_list$Denv, na.rm = T)
IQR(reg_pairs_list$Denv, na.rm = T)

median(reg_pairs_list$Roverlap, na.rm = T)
IQR(reg_pairs_list$Roverlap, na.rm = T)

wilcox.test(reg_pairs_list$Denv, reg_pairs_list$Roverlap, paired = T)

dodge <- position_dodge(width = 0.4)
symnum.args <- list(cutpoints = c(0,0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "n.s."))
colors <- c("#E69F00",  "#2171B5", "#27AD81FF")

## select cuts
bins <- 20
cols <- c("darkred", "red", "white", "blue",  "darkblue")
colGradient <- colorRampPalette(cols)
cut.cols <- colGradient(bins)
cuts = cut(seq(-1,1,0.1),20)[-1]
names(cuts) <- sapply(cuts,function(t) cut.cols[which(as.character(t) == levels(cuts))])
breaks = seq(-1,1,0.1)

## env vs resource-use changes
env_res.hist = data.frame(diff = rowSums(cbind(reg_pairs_list$Denv,-reg_pairs_list$Roverlap), na.rm = T), pair = reg_pairs_list$pair)
env_res.hist[env_res.hist$diff>0, "group"] = "ENV"
env_res.hist[env_res.hist$diff<0, "group"] = "RES"
group_tags <- cut(env_res.hist$diff,
                  breaks=breaks,
                  include.lowest=TRUE,
                  right=FALSE,
                  labels=cuts)
env_res.hist$tag = group_tags
env_res.hist$tag = factor(env_res.hist$tag, levels = sort(unique(env_res.hist$tag)))

env_res_hist <- ggplot(na.exclude(env_res.hist), aes(x=diff, fill=tag, col = tag)) +
  geom_histogram(breaks=seq(-1, 1, by=0.1),
                 alpha=0.8,
                 position="identity") +
  annotate(x=0, xend=0, y=0, yend=150, lwd = 0.5, linetype = 2,  colour="black", geom="segment") +
  geom_histogram(data = env_res.hist[is.na(env_res.hist$group),], breaks=seq(-1.05, 1.05, by=0.1),
                 alpha=0.8,
                 position="identity", col = "black", fill = "black") +
  geom_hline(yintercept = 0, linetype = 1, linewidth = 2, col = "white") +
  annotate(x=median(env_res.hist$diff), xend=median(env_res.hist$diff), y=0, yend=150, lwd = 1,linetype = 1,  colour=if(median(env_res.hist$diff)<0){"red"}else{"blue"}, geom="segment") +
  scale_y_continuous("count", limits = c(-20,150) , expand = c(0,0)) +
  scale_x_continuous(name = "difference", labels = c(1,0.5,0,0.5,1), breaks = seq(-1,1,0.5), expand = c(0,0), limits = c(-1.4,1.4)) +
  scale_fill_manual(values = names(cuts)[cuts %in% env_res.hist$tag]) +
  scale_color_manual(values =  names(cuts)[cuts %in% env_res.hist$tag]) +
  theme_classic() +
  coord_cartesian(xlim=c(-1.4,1.4), ylim=c(-20,170)) +
  annotate(x=-1.4, xend=-1.4, y=0, yend=150, lwd = 1, colour="black", geom="segment") +
  annotate(x=-1, xend=1, y=-20, yend=-20, lwd = 1, colour="black",  geom="segment") +
  theme(panel.border=element_blank(),
        panel.background = element_rect(fill ="transparent", color = NA),
        plot.background = element_rect(fill ="transparent", color = NA),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(3,1,3,1), "cm"),
        axis.title.x = element_text(color = "black", size = 20, vjust = 0),
        legend.title = element_text(color = "black", face = "bold"),
        legend.text = element_text(color = "black"),
        axis.title.y = element_text(color = "black",  size = 20, vjust = 1, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 18,  vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 18 ),
        axis.ticks = element_line(linewidth = 0.75),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "none",
        axis.line = element_blank())

D_gather <- gather(reg_pairs_list, type, D, 4:5)
D_gather$type <- factor(D_gather$type, levels = c("Denv", "Roverlap"))

## ecological to resource use
niche_overlap <- ggplot(D_gather, aes(x = type, y = D, fill = type, col = type)) +
  geom_violin(col = NA) +
  geom_rect(data=NULL,aes(xmin=1,xmax=2,ymin=-Inf,ymax=Inf),
            fill="white" , col = NA, alpha = 0.5) +
  geom_line(data = D_gather[reg_pairs_list$pair %in% env_res.hist[which(env_res.hist$diff >0), "pair"],],
            aes(group = pair), alpha = 0.2, col = "blue") +
  geom_line(data = D_gather[reg_pairs_list$pair %in% env_res.hist[which(env_res.hist$diff <0), "pair"],],
            aes(group = pair), alpha = 0.6, col = "red") +
  geom_line(data = D_gather[reg_pairs_list$pair %in% env_res.hist[which(env_res.hist$diff ==0), "pair"],],
            aes(group = pair), alpha = 0.2, col = "grey40") +
  geom_boxplot(width = 0.05, outlier.colour = NA, outlier.size = 0.8,
               col = "white") +
  stat_summary(fun=median, geom="point", size=3, col = "white", shape = 19, position = dodge) +
  theme_classic() +
  scale_fill_manual("", labels = c("Ecological", "Resource-use"), values= rev(as.vector(colors))[c(1,3)]) +
  scale_color_manual("", labels = c("Ecological", "Resource-use"), values= rev(as.vector(colors))[c(1,3)]) +
  scale_y_continuous(name = "Schoener's D", expand = c(0,0),  seq(0,1,0.25), limits = c(-0.1,1.2)) +
  scale_x_discrete(labels = c("Ecological\nniche overlap", "Resource-use \noverlap"), expand = c(0,0)) +
  coord_cartesian(xlim=c(0.5,2.5), ylim=c(-0.1,1.2)) +
  annotate(x=0.5, xend=0.5, y=0, yend=1, lwd = 1, colour="black", geom="segment") +
  annotate(x=1, xend=2, y=-0.1, yend=-0.1, lwd = 1, colour="black",  geom="segment") +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks = element_line(linewidth = 0.75),
        axis.ticks.length=unit(.2, "cm"),
        legend.position = "none",
        axis.line = element_blank(),
        plot.margin = margin(1, 6, 1, 1, "cm"),
        title = element_blank(),
        legend.text = element_text( color = "black", size=18),
        legend.title = element_text(color = "black", size = 18, vjust = 3),
        axis.title.x = element_blank(),
        axis.title.y = element_text( color = "black", size=24, vjust = 5),
        axis.text.x = element_text( color = "black", size=20, vjust = 0),
        axis.text.y = element_text( color = "black", size=20),
  ) +
  stat_compare_means(comparisons = list(c("Denv", "Roverlap")), paired = T, size  = 12, method = "wilcox.test", label = "..p.signif..", label.y =0.95, symnum.args = symnum.args) 

# Create the combined plot with an offsite, larger inset
fig_overlap <- ggdraw() +
  draw_plot(niche_overlap, x = 0, y = 0, width = 1, height = 1) +  # Main plot
  draw_plot(env_res_hist, x = 0.625, y = 0.55, width = 0.4, height = 0.4)  # Larger and shifted inset

dir.create("figures")
ggsave(paste0("figures/species_niche_overlaps_regional.pdf"), plot = fig_overlap, width=12, height=10, units="in")

 # Species Distributions #######################################################

# get species distributions
clownfish_distributions <- pred.dis[, c("x", "y", clownfish_species)]
# convert occurrence densities to probability of occurrence
clownfish_distributions <- cbind(clownfish_distributions[,1:2], 
                                 apply(clownfish_distributions[,-c(1:2)], 2, function(i) i / max(i)))

# convert to Presence/Absences (Presence > 0) (th > 0)
clownfish_PA <- cbind(clownfish_distributions[,1:2], 
                      apply(clownfish_distributions[,-c(1:2)], 2, function(i) as.numeric(i > quantile(i, 0.25))))

# get maps
clownfish_maps <- rast(clownfish_distributions)
clownfish_PAmaps <- rast(clownfish_PA)

# COMMUNITY ECOLOGY ANALYSES ##################################################

# get community matrix based on abundance (number of locations present) per region
community_matrix <- merge(clownfish_PA, marine_regions, by = c("x", "y"), all.x = TRUE) %>%
  dplyr::group_by(province) %>%
  dplyr::summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
  as.data.frame()
rownames(community_matrix) <- community_matrix$province
community_matrix <- community_matrix[,-c(1:3)]
community_matrix <- community_matrix[rowSums(community_matrix) > 0,]

# get by region environmental data (average values)
community_env_data <- merge(env_data, marine_regions, by = c("x", "y")) %>%
  group_by(province) %>%
  summarise(across(where(is.numeric), list(mean = ~mean(.x, na.rm = TRUE), 
                                           sd = ~sd(.x, na.rm = TRUE)))) %>%
  as.data.frame()

rownames(community_env_data) <- community_env_data$province
community_env_data <- community_env_data[,-c(1:5)]
community_env_data <- community_env_data[rownames(community_matrix),]

nmds <- metaMDS(community_matrix, distance = "bray")
stressplot(nmds)  # Shepard plot to check stress vs. distances
print(nmds$stress)  # Print the stress value

# Plot NMDS results
plot(nmds, type = "t", main = "Community Composition")
realms <- sapply(rownames(community_matrix), function(i) marine_regions$realm[marine_regions$province == i][1])
ordihull(nmds, groups = realms, draw = "polygon", col = rainbow(5))

ordiplot(nmds, display = "species", type = "n")
text(nmds, display = "species", col = "red")

adonis2(community_matrix ~ ., data = community_env_data, method = "bray")

envfit(nmds, community_env_data, permutations = 999)

db_rda <- dbrda(community_matrix ~ ., data = community_env_data, dist="bray",
                sqrt.dist = TRUE)
plot(db_rda)
sppscores(db_rda) <- wisconsin(community_matrix)
plot(db_rda)

common_species <- intersect(clownfish_tree$tip.label, colnames(community_matrix))
clownfish_tree_subset <- keep.tip(clownfish_tree, common_species)

pd_values <- picante::pd(community_matrix, clownfish_tree_subset)$PD  # Calculate PD for each region
names(pd_values) <- rownames(community_matrix)
bp <- barplot(sort(pd_values), 
              main = "Phylogenetic Diversity by Region", 
              xlab = "", 
              ylab = "PD", 
              xaxt = "n")  # Remove default x-axis labels
text(x = bp, 
     y = par("usr")[3] - 0.1,  # Adjust position
     labels = gsub("_", " ", names(sort(pd_values))), 
     srt = 45,  # Rotate text 45 degrees
     adj = 1, 
     xpd = TRUE, 
     cex = 0.7)  # Adjust text size


# 1. Species Richness Patterns ----------------------------------------------

# Calculate observed richness per region
richness <- specnumber(community_matrix)

# Rarefaction (if sample effort varies, assuming community_region values are counts)
min_sample <- 3  # or another value if justified
rarefied_richness <- rarefy(community_matrix, sample = min_sample)

# Diversity indices: Margalef's and Pielou's evenness, Shannon, Simpson etc.
# Margalef's index: (S - 1)/log(N)
margalef <- (richness - 1) / log(rowSums(community_matrix))
# Shannon diversity
shannon <- vegan::diversity(community_matrix, index = "shannon")
# Simpson diversity
simpson <- vegan::diversity(community_matrix, index = "simpson")
# Pielou's evenness: Shannon / log(richness)
pielou <- shannon / log(richness)

diversity_indices <- data.frame(
  Richness = richness,
  Rarefied = rarefied_richness,
  Margalef = margalef,
  Shannon = shannon,
  Simpson = simpson,
  Pielou = pielou
)

#################################

library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggrepel)

compute_niche_mds <- function(reg_pairs_list, variable = "Denv") {
  
  # Step 1: Compute Global Dissimilarity Matrix -----------------------------
  
  full_species <- unique(c(reg_pairs_list$spa, reg_pairs_list$spb))  # All species globally
  
  global_matrix <- matrix(NA, nrow = length(full_species), ncol = length(full_species),
                          dimnames = list(full_species, full_species))
  
  # Compute global dissimilarity per species pair
  global_dissimilarity <- reg_pairs_list %>%
    dplyr::group_by(spa, spb) %>%
    dplyr::summarise(dissimilarity = 1 - mean(.data[[variable]], na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = spb, values_from = dissimilarity)
  
  # Fill global matrix
  for (i in 1:nrow(global_dissimilarity)) {
    sp1 <- global_dissimilarity$spa[i]
    for (j in 2:ncol(global_dissimilarity)) {
      sp2 <- colnames(global_dissimilarity)[j]
      if (!is.na(global_dissimilarity[i, j])) {
        global_matrix[sp1, sp2] <- as.numeric(global_dissimilarity[i, j])
        global_matrix[sp2, sp1] <- as.numeric(global_dissimilarity[i, j]) 
      }
    }
  }
  
  diag(global_matrix) <- 0
  global_matrix[is.na(global_matrix)] <- 1 
  
  # Compute global MDS
  global_mds <- cmdscale(as.dist(global_matrix), k = 2)
  global_mds_df <- data.frame(species = rownames(global_mds), X = global_mds[,1], Y = global_mds[,2])
  
  # Step 2: Compute Regional MDS and Align with Procrustes -------------------
  
  regional_matrix_list <- list()
  regional_mds_list <- list()
  aligned_mds_list <- list()
  
  for(region in unique(reg_pairs_list$region)) {
    
    region_data <- reg_pairs_list %>%
      dplyr::filter(region == !!region) %>%
      dplyr::select(spa, spb, all_of(variable)) %>%
      tidyr::pivot_wider(names_from = spb, values_from = variable)
    
    species_in_region <- unique(c(region_data$spa, colnames(region_data)[-1]))
    
    if (length(species_in_region) > 2) {
      
      region_matrix <- matrix(NA, nrow = length(species_in_region), ncol = length(species_in_region),
                              dimnames = list(species_in_region, species_in_region))
      
      for (i in 1:nrow(region_data)) {
        sp1 <- region_data$spa[i]
        for (j in 2:ncol(region_data)) {
          sp2 <- colnames(region_data)[j]
          if (!is.na(region_data[i, j])) {
            region_matrix[sp1, sp2] <- 1 - as.numeric(region_data[i, j])  
            region_matrix[sp2, sp1] <- 1 - as.numeric(region_data[i, j])  
          }
        }
      }
      
      diag(region_matrix) <- 0
      region_matrix[is.na(region_matrix)] <- 1  
      
      regional_matrix_list[[region]] <- region_matrix
      
      regional_mds <- cmdscale(as.dist(region_matrix), k = 2)
      regional_mds_list[[region]] <- regional_mds
      
      common_species <- intersect(rownames(regional_mds), global_mds_df$species)
      
      if (length(common_species) > 2) {
        proc <- procrustes(
          global_mds_df %>% dplyr::filter(species %in% common_species) %>% dplyr::select(X, Y),
          regional_mds[common_species, ]
        )
        
        aligned_mds <- predict(proc, regional_mds)
        aligned_mds_list[[region]] <- data.frame(species = rownames(aligned_mds), X = aligned_mds[,1], Y = aligned_mds[,2])
      }
    }
  }
  
  # Step 3: Compute Species Niche Centroids and Niche Breadth --------------
  
  all_regions_df <- dplyr::bind_rows(aligned_mds_list, .id = "region")
  
  species_niche <- all_regions_df %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(niche_X = mean(X, na.rm = TRUE),
                     niche_Y = mean(Y, na.rm = TRUE),
                     breadth_X = sd(X, na.rm = TRUE),
                     breadth_Y = sd(Y, na.rm = TRUE),
                     .groups = "drop")
  
  
  return(list(global_matrix = global_matrix,
              global_mds = global_mds_df,
              regional_matrix = regional_matrix_list,
              regional_mds = regional_mds_list,
              aligned_mds = aligned_mds_list,
              species_niche = species_niche))
}

Denv_niche_mds <- compute_niche_mds(reg_pairs_list, "Denv")
Roverlap_niche_mds <- compute_niche_mds(reg_pairs_list, "Roverlap")

Denv_niche_mds$species_niche$behavior <- Roverlap_niche_mds$species_niche$behavior <- factor(gsub("_", " ", sapply(Denv_niche_mds$species_niche$species, function(sp) behavior_df$category[behavior_df$species == sp])),
                                                                                             levels = c("generalist", "EQ specialist", "HM specialist", "SD specialist"))

host_cols <- setNames(c( "black", "red", "orange", "lightblue"), c("generalist", "EQ specialist", "HM specialist", "SD specialist"))

ggplot(Denv_niche_mds$species_niche, aes(x = niche_X, y = niche_Y, label = species, color = behavior)) +
  geom_errorbar(aes(ymin = niche_Y - breadth_Y, ymax = niche_Y + breadth_Y), 
                width = 0.2, linewidth = 0.8, alpha = 0.25) +  
  geom_errorbarh(aes(xmin = niche_X - breadth_X, xmax = niche_X + breadth_X), 
                 height = 0.2, linewidth = 0.8, alpha = 0.25) +  
  geom_point(size = 6) +  
  geom_text_repel(aes(label = paste0("italic('", sub("Premnas_", "P. ", sub("Amphiprion_", "A. ", species)), "')")), 
                  size = 4, max.overlaps = 15, box.padding = 0.5, color = "black", parse = TRUE) +  
  scale_color_manual(values = host_cols) +  # Custom colors for behavior
  theme_minimal(base_size = 16) +  
  labs(
    x = "Principal Coordinate 1",
    y = "Principal Coordinate 2",
    title = "Species Niche Positions based on environmental dissimilarity",
    color = "Behavior Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),  
    panel.grid = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  
  )
ggsave(paste0("figures/species_environmental_niche_pos_from_Denv.pdf"), width=12, height=12, units="in")


ggplot(Roverlap_niche_mds$species_niche, aes(x = niche_X, y = niche_Y, label = species, color = behavior)) +
  geom_errorbar(aes(ymin = niche_Y - breadth_Y, ymax = niche_Y + breadth_Y), 
                width = 0.2, linewidth = 0.8, alpha = 0.25) +  
  geom_errorbarh(aes(xmin = niche_X - breadth_X, xmax = niche_X + breadth_X), 
                 height = 0.2, linewidth = 0.8, alpha = 0.25) +  
  geom_point(size = 6) +  
  geom_text_repel(aes(label = paste0("italic('", sub("Premnas_", "P. ", sub("Amphiprion_", "A. ", species)), "')")), 
                  size = 4, max.overlaps = 15, box.padding = 0.5, color = "black", parse = TRUE) +  
  scale_color_manual(values = host_cols) +  # Custom colors for behavior
  theme_minimal(base_size = 16) +  
  labs(
    x = "Principal Coordinate 1",
    y = "Principal Coordinate 2",
    title = "Species Niche Positions based on resrouce dissimilarity",
    color = "Behavior Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),  
    panel.grid = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  
  )
ggsave(paste0("figures/species_resource_niche_pos_from_Roverlap.pdf"), width=12, height=12, units="in")


NO_diss_matrix <- Denv_niche_mds$global_matrix
RO_diss_matrix <- Roverlap_niche_mds$global_matrix

common_species <- intersect(clownfish_tree$tip.label, colnames(NO_diss_matrix))
clownfish_tree_subset <- keep.tip(clownfish_tree, common_species)

NO_diss_matrix <- NO_diss_matrix[common_species, common_species]
RO_diss_matrix <- RO_diss_matrix[common_species, common_species]

# Compute pairwise phylogenetic distances
phylo_dist <- cophenetic(clownfish_tree_subset)

# Mantel test: Correlation between phylogenetic distance and niche overlap
mantel_test_NO <- mantel(phylo_dist, NO_diss_matrix, method="pearson", permutations=999)
mantel_test_RO <- mantel(phylo_dist, RO_diss_matrix, method="pearson", permutations=999)

print(mantel_test_NO)
print(mantel_test_RO)

# Perform Principal Coordinates Analysis (PCoA)
pcoa_NO <- cmdscale(NO_diss_matrix, k=4, eig = TRUE)  
pcoa_RO <- cmdscale(RO_diss_matrix, k=4, eig = TRUE)  

NO_var_perc <- round(pcoa_NO$eig^2 / sum(pcoa_NO$eig^2) * 100,2)[1:4]
RO_var_perc <- round(pcoa_RO$eig^2 / sum(pcoa_RO$eig^2) * 100,2)[1:4]

sum(NO_var_perc)
sum(RO_var_perc)

# Convert PCoA scores to named vectors
RO_niche_df <- data.frame(pcoa_RO$points, species = rownames(pcoa_RO$points))
colnames(RO_niche_df)[1:4] <- paste0("PCo", 1:4)

RO_niche_df$behavior <- factor(gsub("_", " ", sapply(RO_niche_df$species, function(i) behavior_df$category[behavior_df$species==i])), 
                               levels = names(host_cols))

RO_pco12 <- ggplot(RO_niche_df, aes(x = PCo1, y = PCo2, label = species, fill = behavior)) +
                  geom_point(size = 5, shape = 21) +  
                  geom_text_repel(aes(label = paste0("italic('", sub("Premnas_", "P. ", sub("Amphiprion_", "A. ", species)), "')")), 
                                  size = 4, max.overlaps = 15, box.padding = 0.5, color = "black", parse = TRUE) +                    
                  scale_fill_manual(values = host_cols) +  
                  theme_minimal(base_size = 16) +
                  labs(
                    x = paste0("PCo 1 (", RO_var_perc[1], "%)"),
                    y = paste0("PCo 2 (", RO_var_perc[2], "%)"),
                    fill = "Behavior Type"
                  ) +
                  theme(
                    plot.margin = margin(1, 1, 1, 1, "cm"),
                    plot.title = element_text(hjust = 0.5),  
                    panel.grid = element_blank(),  
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  
                  )

RO_pco34 <- ggplot(RO_niche_df, aes(x = PCo3, y = PCo4, label = species, fill = behavior)) +
                  geom_point(size = 5, shape = 21) +  
                  geom_text_repel(aes(label = paste0("italic('", sub("Premnas_", "P. ", sub("Amphiprion_", "A. ", species)), "')")), 
                                  size = 4, max.overlaps = 15, box.padding = 0.5, color = "black", parse = TRUE) +  
                  scale_fill_manual(values = host_cols) +  
                  theme_minimal(base_size = 16) +
                  labs(
                    x = paste0("PCo 3 (", RO_var_perc[3], "%)"),
                    y = paste0("PCo 4 (", RO_var_perc[4], "%)"),
                    fill = "Behavior Type"
                  ) +
                  theme(
                    plot.margin = margin(1, 1, 1, 1, "cm"),
                    plot.title = element_text(hjust = 0.5),  
                    panel.grid = element_blank(),  
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  
                  )

ggarrange(RO_pco12, RO_pco34, ncol = 2, common.legend = TRUE, labels = letters[1:2], legend = "bottom",  font.label = list(size = 35))
ggsave(paste0("figures/species_resource_niche_pos_from_Roverlap_global.pdf"), width=18, height=12, units="in")

NO_niche_df <- data.frame(pcoa_NO$points, species = rownames(pcoa_NO$points))
colnames(NO_niche_df)[1:4] <- paste0("PCo", 1:4)

NO_niche_df$behavior <- factor(gsub("_", " ", sapply(NO_niche_df$species, function(i) behavior_df$category[behavior_df$species==i])), 
                               levels = names(host_cols))

NO_pco12 <- ggplot(NO_niche_df, aes(x = PCo1, y = PCo2, label = species, fill = behavior)) +
  geom_point(size = 5, shape = 21) +  
  geom_text_repel(aes(label = paste0("italic('", sub("Premnas_", "P. ", sub("Amphiprion_", "A. ", species)), "')")), 
                  size = 4, max.overlaps = 15, box.padding = 0.5, color = "black", parse = TRUE) +                    
  scale_fill_manual(values = host_cols) +  
  theme_minimal(base_size = 16) +
  labs(
    x = paste0("PCo 1 (", NO_var_perc[1], "%)"),
    y = paste0("PCo 2 (", NO_var_perc[2], "%)"),
    fill = "Behavior Type"
  ) +
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"),
    plot.title = element_text(hjust = 0.5),  
    panel.grid = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  
  )

NO_pco34 <- ggplot(NO_niche_df, aes(x = PCo3, y = PCo4, label = species, fill = behavior)) +
  geom_point(size = 5, shape = 21) +  
  geom_text_repel(aes(label = paste0("italic('", sub("Premnas_", "P. ", sub("Amphiprion_", "A. ", species)), "')")), 
                  size = 4, max.overlaps = 15, box.padding = 0.5, color = "black", parse = TRUE) +  
  scale_fill_manual(values = host_cols) +  
  theme_minimal(base_size = 16) +
  labs(
    x = paste0("PCo 3 (", NO_var_perc[3], "%)"),
    y = paste0("PCo 4 (", NO_var_perc[4], "%)"),
    fill = "Behavior Type"
  ) +
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"),
    plot.title = element_text(hjust = 0.5),  
    panel.grid = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  
  )

ggarrange(NO_pco12, NO_pco34, ncol = 2, common.legend = TRUE, labels = letters[1:2], legend = "bottom",  font.label = list(size = 35))
ggsave(paste0("figures/species_environmental_niche_pos_from_Denv_global.pdf"), width=18, height=12, units="in")

# Blomberg's K and Pagel's λ for PCoA axis
PCo_cols <- grep("PCo", colnames(NO_niche_df), value = TRUE)

get_phylosig <- function(df) {
  
  phylosig_df <- data.frame(matrix(NA, nrow = length(PCo_cols), ncol = 4, dimnames = list(PCo_cols, c("K", "K.p", "lambda", "lambda.p"))))

  for (pco in 1:nrow(phylosig_df)) {
    pco_vec <- setNames(df[,pco], df$species)
    
    K <- phylosig(clownfish_tree_subset, pco_vec, method = "K", test = 999)
    lambda <- phylosig(clownfish_tree_subset, pco_vec, method = "lambda", test = 999)
    
    phylosig_df[pco, ] <- c(K$K, K$P, lambda$lambda, lambda$P)
  }
  
  return(phylosig_df)
}

NO_phylosig <- get_phylosig(NO_niche_df)
RO_phylosig <- get_phylosig(RO_niche_df)

# Create a phylogenetic variance-covariance matrix
vcv_phylo <- vcv(clownfish_tree_subset, corr=TRUE)

clownfish_tree_subset$root.edge <- NULL

library(phylolm)
model <- phylolm(PCo3 ~ behavior, 
                 data = RO_niche_df[clownfish_tree_subset$tip.label,], 
                 phy = clownfish_tree_subset, 
                 model = "BM")  
summary(model)


# Perform PGLS regression
RO_dat <- data.frame(species = names(RO_niche_trait1), niche_trait1 = RO_niche_trait1, niche_trait2 = RO_niche_trait2)
RO_pgls_model <- gls(niche_trait1 ~ niche_trait2, correlation=corBrownian(1, clownfish_tree_subset, form = ~ species), method="ML",
                     data = RO_niche_df)
summary(RO_pgls_model)

NO_dat <- data.frame(species = names(NO_niche_trait1), niche_trait1 = NO_niche_trait1, niche_trait2 = NO_niche_trait2)
NO_pgls_model <- gls(niche_trait1 ~ niche_trait2, correlation=corBrownian(1, clownfish_tree_subset, form = ~ species), method="ML",
                     data = NO_dat)
summary(NO_pgls_model)

library(betapart)
library(picante)
library(FD)
devtools::install_github("heibl/phyloclim")
library(phyloclim)

# Compute biogeographic dissimilarity
comm_matrix <- community_matrix
comm_matrix[comm_matrix>0] = 1
biogeo_dist <- beta.pair(t(comm_matrix))$beta.sor  # Sørensen dissimilarity
biogeo_dist <- as.dist(as.matrix(biogeo_dist)[common_species,common_species])

# Test correlation with phylogenetic distances
mantel_biogeo <- mantel(phylo_dist, biogeo_dist, method="pearson", permutations=999)
print(mantel_biogeo)

# Calculate Net Relatedness Index (NRI) and Nearest Taxon Index (NTI)
phylo_structure <- ses.mpd(comm_matrix[,clownfish_tree_subset$tip.label], cophenetic(clownfish_tree_subset), null.model="taxa.labels", runs=999)
print(phylo_structure)
plot(phylo_structure$mpd.rand.mean, phylo_structure$mpd.obs)
abline(0, 1)

Niche_centroid_global <- t(sapply(ENMs, function(sp) niche_position(sp$z.uncor)))

Niche_centroid_reg <- sapply(names(clownfish_regENMs), function(reg) t(sapply(clownfish_regENMs[[reg]], function(sp) niche_position(sp$z.uncor))))
Niche_centroid_reg <- Niche_centroid_reg[lapply(Niche_centroid_reg, ncol) > 0]
Niche_centroid_reg <- sapply(Niche_centroid_reg, function(i) data.frame(species = rownames(i), i), simplify = FALSE)
Niche_centroid_reg <- plyr::ldply(Niche_centroid_reg, .id = "region")
Niche_centroid_reg$behavior <- sub("_", " ", host_categories[Niche_centroid_reg$species])

env_envelope <- env_scores %>% 
  slice(chull(Axis1, Axis2))  

species_niche_centroid <- Niche_centroid_reg  %>% 
  dplyr::group_by(species) %>%
  dplyr::summarise(Axis1 = mean(x, na.rm = TRUE),
                   Axis2 = mean(y, na.rm = TRUE),
                   behavior = behavior[1])

species_niche_hull <- Niche_centroid_reg %>% 
  dplyr::group_by(species) %>%
  slice(chull(x, y))

ggplot(data = env_envelope, aes(x = Axis1, y = Axis2),) +
  #geom_polygon(color = "grey", alpha = 0.5) + 
  geom_polygon(data = species_niche_hull, aes(x = x, y = y, group = species, fill = behavior), alpha = 0.5) +  
  geom_point(data = species_niche_centroid, aes(fill = behavior), size = 5, shape = 21) +  
  geom_text_repel(data = species_niche_centroid, aes(label = paste0("italic('", sub("Premnas_", "P. ", sub("Amphiprion_", "A. ", species)), "')")), 
                  size = 4, max.overlaps = 15, box.padding = 0.5, color = "black", parse = TRUE) +  
  scale_fill_manual(values = host_cols) +  
  theme_minimal(base_size = 16) +
  #coord_cartesian(xlim = c(-1,1.5), ylim = c(-0.7, 0)) +
  coord_cartesian(xlim = range(Niche_centroid_reg$x), ylim = range(Niche_centroid_reg$y)) +
  labs(
    x = "x",
    y = "y",
    fill = "Behavior Type"
  ) +
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"),
    plot.title = element_text(hjust = 0.5),  
    panel.grid = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  
  )


niche_data <- data.frame(centroid.x = Niche_centroid_global[common_species, 1],
                         centroid.y = Niche_centroid_global[common_species, 2],
                         NO_niche_df[common_species, 1:4],
                         RO_niche_df[common_species, -5])
colnames(niche_data)[3:10] <- paste0( rep(c("env", "res"), each = 4), gsub("\\.1", "", colnames(niche_data)[3:10] ))

# Calculate functional diversity indices (e.g., Rao's Q)
functional_diversity <- dbFD(
  x = niche_data, a = community_matrix[,common_species],
  calc.FRic = TRUE,
  calc.FDiv = TRUE,
  clust.type = "ward"
)

print(functional_diversity)

# Convert results into a data frame
FD_results <- data.frame(
  region = rownames(functional_diversity$CWM),  # Species names
  functional_diversity$CWM,
  RaoQ = functional_diversity$RaoQ,  # Rao's Quadratic Entropy
  FRic = functional_diversity$FRic,  # Functional Richness
  FDiv = functional_diversity$FDiv   # Functional Divergence
)

FD_results$behavior <- factor(FD_results$behavior, levels = names(host_cols))

# Functional Richness vs. Rao’s Q**
ggplot(FD_results, aes(x = FRic, y = RaoQ, fill = behavior)) +
  geom_point(size = 5, shape = 21, color = "black", alpha = 0.7) +  
  geom_text_repel(aes(label = gsub("_", " ",  region)), size = 4,
                  max.overlaps = 15, box.padding = 0.5, color = "black") +  
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = host_cols) +  
  labs(x = "Functional Richness", 
       y = "Rao's Quadratic Entropy",
       fill = "Dominant behavior type") +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.margin = margin(1, 1, 1, 1, "cm"),
        panel.grid = element_blank(),  
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "bottom")  
ggsave(paste0("figures/community_FunctionalDiversity.pdf"), width=10, height=10, units="in")

#  Comparing Functional Divergence (FDiv) Across Behaviors**
Fdiv <- ggplot(FD_results, aes(x = behavior, y = FDiv, fill = behavior)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(shape = 21, size = 3, color = "black", width = 0.2) +
  scale_fill_manual(values = host_cols) +
  theme_minimal(base_size = 16) +
  labs(x = "Dominant behavior type", y = "Functional Divergence") +
  theme(plot.margin = margin(1, 1, 1, 1.5, "cm"),
        axis.title.x = element_text(vjust = -1))

# Comparing Functional Richness (FRic) Across Behaviors**
Fric <- ggplot(FD_results, aes(x = behavior, y = FRic, fill = behavior)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(shape = 21, size = 3, color = "black", width = 0.2) +
  scale_fill_manual(values = host_cols) +
  theme_minimal(base_size = 16) +
  labs(x = "Dominant behavior type", y = "Functional Richness") +
  theme(plot.margin = margin(1, 1, 1, 1.5, "cm"),
        axis.title.x = element_text(vjust = -1))


# Comparing Functional Entropy (RaoQ) Across Behaviors**
RaoQ <- ggplot(FD_results, aes(x = behavior, y = RaoQ, fill = behavior)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(shape = 21, size = 3, color = "black", width = 0.2) +
  scale_fill_manual(values = host_cols) +
  theme_minimal(base_size = 16) +
  labs(x = "Dominant behavior type", y = "Rao's Quadratic Entropy") +
  theme(plot.margin = margin(1, 1, 1, 1.5, "cm"),
        axis.title.x = element_text(vjust = -1))

ggarrange(Fdiv, Fric, RaoQ, ncol = 1, common.legend = TRUE, labels = letters[1:3], legend = "none",  font.label = list(size = 35))
ggsave(paste0("figures/community_FD_boxplots.pdf"), width=12, height=18, units="in")


# Graph analysis
interaction_network <- graph_from_biadjacency_matrix(community_matrix)
modularity <- computeModules(interaction_network)
nestedness <- nestedness(interaction_network)
specialization <- H2fun(interaction_network)

# Null biogeographic test
null_test <- custom_null_biogeographic_test(niche_overlap, phylogeny)

## from heibl github phyloclim package
nbConnectingNodes <- function(phy, npair){
  ntips <- length(phy$tip.label)
  nds <- getMRCA(phy, npair)
  nds <- descendants(phy, nds, internal = TRUE)
  if (identical(sort(nds), sort(npair)))
    nb <- 1										else {
      nds <- nds[nds > ntips]
      check <- function(x, npair)
        any(npair %in% descendants(phy, x))
      id <- sapply(nds, check, npair = npair)
      nds <- nds[id] 
      nb <- length(nds) + 1 
    }
  nb
}

descendants <-
  function(tree, node, internal = FALSE, string = FALSE){
    
    tips <- seq(along = tree$tip.label)
    x <- tree$edge[,2][tree$edge[,1] == node]
    repeat{
      xx <- x
      x <- sort(unique(c(x, tree$edge[,2][tree$edge[,1] %in% x])))
      if (identical(x, xx)) break
    }
    # return tip number if input is tip number:
    # -----------------------------------------
    if (length(x) == 0) x <- node
    if (!internal)
      x <- x[x %in% tips]
    if (string)
      x <- tree$tip.label[x]
    x
  }

nested.mean.overlap <- function(phy, node, olap){
  
  # match ordering of phy and olap
  # ------------------------------
  id <- match(phy$tip.label, rownames(olap))
  olap <- olap[id, id]
  
  # get daughter nodes
  # ------------------
  d2 <- phy$edge[phy$edge[, 1] == node, 2]
  
  # get descendents of both daughter nodes
  # --------------------------------------
  C1 <- descendants(phy, d2[1])
  C2 <- descendants(phy, d2[2])
  
  # calculate mean overlap
  # ----------------------
  o <- 0
  for (j in C1){
    for (k in C2){
      n <- nbConnectingNodes(phy, c(j, k))
      o <- o + 0.5 ^ (n - 1) * olap[j, k]
    }
  }
  o
}

age.range.correlation <- function(phy, overlap, tri = "upper", 
                                  n = 10000){
  
  # check input
  # -----------
  if (!inherits(phy, "phylo")) 
    stop("object 'phy' is not of class 'phylo'")
  if (!is.ultrametric(phy))
    stop("object 'phy' must be ultrametric")
  
  # ages:
  # -----
  age <- branching.times(phy)	
  
  # make matrix symmetrical
  # -----------------------
  ovlap <- overlap
  if (tri == "upper")
    ovlap[lower.tri(ovlap)] <- t(ovlap)[lower.tri(ovlap)]
  if (tri == "lower")
    ovlap[upper.tri(ovlap)] <- t(ovlap)[upper.tri(ovlap)]
  
  # match matrix to tree
  # --------------------
  id <- match(phy$tip.label, rownames(ovlap))
  ovlap <- ovlap[id, id]
  
  # calculate 'nested mean overlap'
  # -------------------------------
  overlap <- sapply(names(age), nested.mean.overlap, phy = phy, 		
                    olap = ovlap)
  
  x <- cbind(age, overlap)	
  x.lm <- lm(overlap ~ age)
  
  # randomization:
  # --------------
  randomization <- function(phy, o, n, age) {
    id <- sample(seq_len(nrow(o)))
    colnames(o) <- rownames(o) <- colnames(o)[id]
    
    o_values <- sapply(names(age), nested.mean.overlap, phy = phy, olap = o)
    
    lm_fit <- lm(o_values ~ age)
    lm_fit$coefficients
  }
  
  random.x <- t(replicate(n, randomization(phy, ovlap, n, age)))
  
  ## fraction of intercepts and slopes from the randomization,
  ## which are greater than the observed values
  ## ------------------------------------------
  f.intercept <- sum(random.x[, "(Intercept)" ] > x.lm$coefficients["(Intercept)"]) / n
  f.slope <- sum(random.x[, "age" ] > x.lm$coefficients["age"]) / n

  ## 2-sided p-values
  ## ----------------
  f <- c(f.intercept, f.slope)
  p <- sapply(f, function(x) 2 * min(x, 1 - x))
  sig <- setNames(c(f, p), c("(intercept)", "age"))
  
  list(age.range.correlation = x, 
       linear.regression = x.lm, 
       sig = sig,
       MonteCarlo.replicates = random.x	
  )
}

# Age-range correlation test
NO_arc_test <- age.range.correlation(clownfish_tree_subset, NO_overlap_matrix, tri = "upper", n = 1000)
NO_arc_test$sig

RO_arc_test <- age.range.correlation(clownfish_tree_subset, RO_overlap_matrix, tri = "upper", n = 1000)
RO_arc_test$sig

 # Phylogenetic Comparative Methods
bm_model <- fitContinuous(phylogeny, trait_data, model = "BM")
ou_model <- fitContinuous(phylogeny, trait_data, model = "OU")
model_selection <- aic(bm_model, ou_model)

# Ancestral state reconstruction
anc_states <- ace(trait_data, phylogeny, type = "continuous")

# Biogeographic reconstruction
bio_reconstruction <- BioGeoBEARS::runBSM(phylogeny, range_data)

# Phylogenetic diversity patterns
pd <- pd(community_matrix, clownfish_tree_subset)
mpd <- mpd(community_matrix, cophenetic(clownfish_tree_subset))
mntd <- mntd(community_matrix, cophenetic(clownfish_tree_subset))

ses_pd <- ses.pd(community_matrix, clownfish_tree_subset)
ses_mpd <- ses.mpd(community_matrix, cophenetic(clownfish_tree_subset))
ses_mntd <- ses.mntd(community_matrix, cophenetic(clownfish_tree_subset))

# Adjust p-values for multiple comparisons
p_adjusted <- p.adjust(c(K$P, lambda$P, null_test$p, arc_test$p), method = "BH")

# Account for spatial autocorrelation
mems <- MEM.autocor(coordinates(occurrence_data))
