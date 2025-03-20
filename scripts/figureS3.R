# load libraries
library(raster)
library(NINA)
library(stringr)
library(ggplot2)
library(ggpubr)
library(colorspace)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

# load custom functions
source("./scripts/functions.R")

# store in memory
rasterOptions(todisk = FALSE)

amph_EN <- readRDS("./Rdata/amphENMs.RDS")
anem_EN <- readRDS("./Rdata/anemENMs.RDS")
amph_BC <- readRDS("./Rdata/amphEBMs.RDS")

# association matrix
int.mat = read.csv("data/interaction_matrix.csv", row.names = 1)
int.mat[int.mat > 0] = 1

# marine regions
marine_regions <- read.csv("./data/marine_regions.csv")
##################

# species regions
comm_mat <- amph_EN$tab

ROUcols <- c("#E69F00", "#2171B5", "#27AD81FF")

# ROU geographical projections
##########################

# Convert to P/A distributions
amph_EN.PA <- prob_to_PA(amph_EN$maps, th = setNames(pmax(amph_EN$eval$threshold$threshold, 0.00001), rownames(amph_EN$eval$threshold)))
amph_BC.PA <- prob_to_PA(amph_BC$maps, th = setNames(pmax(amph_BC$eval$threshold$threshold, 0.00001), rownames(amph_BC$eval$threshold)))
anem_EN.PA <- prob_to_PA(anem_EN$maps, th = setNames(pmax(anem_EN$eval$threshold$threshold, 0.00001), rownames(anem_EN$eval$threshold)))

amph_species <- names(amph_BC.PA)

R <- O <- U <- stack()
for (sp in amph_species) {
  
  # specific associations
  Xvar = colnames(int.mat)[int.mat[sp,] == 1]
  
  # get ecological and mutualism-refined maps
  en_rasPA <- amph_EN.PA[[sp]]
  bc_rasPA <- amph_BC.PA[[sp]]

  # get map of host availability
  w_ENras_PA <- sum(anem_EN.PA[[Xvar]])
  w_ENras_PA[w_ENras_PA > 1] = 1 # Presence/Absence of any host 
  
  # mask inhabited marine regions only
  regions <- names(comm_mat[sp,])[comm_mat[sp,] == 1]
  xy.reg <- marine_regions[!marine_regions$province %in% regions,1:2]
  cellIDs <- cellFromXY(w_ENras_PA, xy.reg)
  
  w_ENras_PA[cellIDs] = NA
  en_rasPA[cellIDs] = NA
  bc_rasPA[cellIDs] = NA
  
  R <- stack(R, en_rasPA - en_rasPA * bc_rasPA)
  O <- stack(O, en_rasPA * bc_rasPA)
  U <- stack(U, w_ENras_PA - w_ENras_PA * bc_rasPA)
  
}
R <- rast(mean(R, na.rm = TRUE))
O <- rast(mean(O, na.rm = TRUE))
U <- rast(mean(U, na.rm = TRUE))

# Define color palettes for each map
R_colors <- colorRampPalette(c(lighten(ROUcols[1], 0.75), darken(ROUcols[1], 0.25)))(100)
O_colors <- colorRampPalette(c(lighten(ROUcols[2], 0.75), darken(ROUcols[2], 0.25)))(100)
U_colors <- colorRampPalette(c(lighten(ROUcols[3], 0.75), darken(ROUcols[3], 0.25)))(100)

# Convert raster to dataframe for ggplot
R_df <- na.exclude(as.data.frame(R, xy = TRUE))
O_df <- na.exclude(as.data.frame(O, xy = TRUE))
U_df <- na.exclude(as.data.frame(U, xy = TRUE))

# Load world map for reference
world <- ne_countries(scale = "medium", returnclass = "sf")

# Function to create ggplot for each dataset
plot_map <- function(data, fill_colors, title) {
  ggplot() +
    geom_sf(data = world, fill = "gray90", color = "grey60") +
    geom_tile(data = data, aes(x = x, y = y, fill = layer)) +
    scale_fill_gradientn(colors = fill_colors, name = title) +  
    scale_y_continuous(breaks = seq(-30, 30, 10), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(50, 200, 50), expand = c(0,0)) +
    coord_sf(xlim = range(data$x, na.rm = TRUE), ylim = range(data$y, na.rm = TRUE)) +
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude") +
    theme(plot.margin = margin(20, 20, 20, 20),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white", color = NA),
          legend.key.size = unit(1.5, "cm"),
          legend.ticks = element_line(linewidth = c(0, 1, 1, 1, 0), color ="white"),
          panel.grid = element_line(color = "grey90"),
          text = element_text(size = 24))
}

# Generate plots
R_plot <- plot_map(R_df, R_colors, "Restricted\n")
O_plot <- plot_map(O_df, O_colors, "Occupied\n")
U_plot <- plot_map(U_df, U_colors, "Unexploited\n")

figS3 <- ggarrange(R_plot, O_plot, U_plot, nrow = 3, labels = c("(a)", "(b)", "(c)"), font.label = list(size = 30))

ggsave("figures/FigS3.ROU_geographical_distribution.pdf", plot = figS3, width=20, height=19, units="in")

