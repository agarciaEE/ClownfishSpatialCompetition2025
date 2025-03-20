# Load dependencies
library(tidyverse)
library(colorspace)
library(ggplot2)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

# load custom functions
source("./scripts/functions.R")

# Load ROU estimates
data <- read.csv("./results/ROU_data.csv")
data <- data[, c("species", "Restricted_niche", "Unexploited_niche", "Occupied_niche", "behavior_reg")]
colnames(data) <- c("species", "Restricted", "Unexploited", "Occupied", "behavior")

# Load interactions matrix to order species in based to number of interactions
int.mat = read.csv("data/interaction_matrix.csv", row.names = 1)
int.mat[int.mat > 0] = 1

data[order(data$Restricted_niche, decreasing = TRUE),]
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

data <- read.csv("./results/ROU_data.csv")
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
ggsave("figures/Fig2.speciesROU_and_ROU_vs_behavior.pdf", plot = fig2, width=18, height=12, units="in")
