# Load dependencies
library(spdep)
library(spatialreg)
library(reshape2)
library(spaMM)
library(NINA)
library(ggplot2)
library(ggtext)  
library(metR)
library(MetBrewer)
library(dplyr)
library(DHARMa)
library(ggpubr)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

INT.pred <- read.csv("results/spatial_results.csv")

### Too many data points to process the spatial models. Aggregate raster to reduce data points
INT.ras <- NINA::raster_projection(INT.pred[,-1])
INT.ras10 <- raster::aggregate(INT.ras, fact = 10, fun = mean, na.rm = T)
#### Spatial Analyses on spatial point interactions dataset
spat_int_df <- raster::as.data.frame(INT.ras10, xy = T)

# remove NA's and points of no interaction
spat_int_df <- spat_int_df[!is.na(spat_int_df$amph.richness),]
spat_int_df <- spat_int_df[spat_int_df$amph.richness > 0,]
spat_int_df[is.na(spat_int_df)] = 0

fit.glm.spatial <- readRDS("Rdata/all_spatialGLMM.RDS")
fit.glm.spatialGG <- readRDS("Rdata/GG_spatialGLMM.RDS")
fit.glm.spatialGS <- readRDS("Rdata/GS_spatialGLMM.RDS")
fit.glm.spatialSS <- readRDS("Rdata/SS_spatialGLMM.RDS")

# Create a grid of predictor values, ensuring a proper range
sp.pred_RO <- expand.grid(
  amph.richness = seq(2, max(spat_int_df$amph.richness, na.rm = TRUE), length.out = 50),
  ENTotal = seq(0, 1, length.out = 50),
  x = mean(spat_int_df$x, na.rm = TRUE),  # Keeping x constant
  y = mean(spat_int_df$y, na.rm = TRUE)   # Keeping y constant
)
sp.pred_ROGG <- expand.grid(
  nG = seq(2, max(spat_int_df$amph.richness, na.rm = TRUE), length.out = 50),
  ENGG = seq(0, 1, length.out = 50),
  x = mean(spat_int_df$x, na.rm = TRUE),  # Keeping x constant
  y = mean(spat_int_df$y, na.rm = TRUE)   # Keeping y constant
)
sp.pred_ROGS <- expand.grid(
  amph.richness = seq(2, max(spat_int_df$amph.richness, na.rm = TRUE), length.out = 50),
  ENGS = seq(0, 1, length.out = 50),
  x = mean(spat_int_df$x, na.rm = TRUE),  # Keeping x constant
  y = mean(spat_int_df$y, na.rm = TRUE)   # Keeping y constant
)
sp.pred_ROSS <- expand.grid(
  nS = seq(2, max(spat_int_df$amph.richness, na.rm = TRUE), length.out = 50),
  ENSS = seq(0, 1, length.out = 50),
  x = mean(spat_int_df$x, na.rm = TRUE),  # Keeping x constant
  y = mean(spat_int_df$y, na.rm = TRUE)   # Keeping y constant
)
# Predict ROGG values using the fitted model
sp.pred_RO$predicted_RO <- predict(fit.glm.spatial, newdata = sp.pred_RO)
sp.pred_RO$predicted_RO[sp.pred_RO$predicted_RO < 0] = 0
sp.pred_ROGG$predicted_ROGG <- predict(fit.glm.spatialGG, newdata = sp.pred_ROGG)
sp.pred_ROGG$predicted_ROGG[sp.pred_ROGG$predicted_ROGG < 0] = 0
sp.pred_ROGS$predicted_ROGS <- predict(fit.glm.spatialGS, newdata = sp.pred_ROGS)
sp.pred_ROGS$predicted_ROGS[sp.pred_ROGS$predicted_ROGS < 0] = 0
sp.pred_ROSS$predicted_ROSS <- predict(fit.glm.spatialSS, newdata = sp.pred_ROSS)
sp.pred_ROSS$predicted_ROSS [sp.pred_ROSS$predicted_ROSS  < 0] = 0

# Define custom breaks for contours
breaks_seq <- seq(0, 1, by = 0.1)
colnames(sp.pred_RO) <- colnames(sp.pred_ROGG) <- colnames(sp.pred_ROGS) <- colnames(sp.pred_ROSS) <- c("ns", "EN", "x", "y", "predRO")
sp.pred_RO_df <- rbind(cbind(sp.pred_RO, model = "all"),
                       cbind(sp.pred_ROGG, model = "GG"),
                       cbind(sp.pred_ROGS, model = "GS"),
                       cbind(sp.pred_ROSS, model = "SS"))

labeller_fn <- as_labeller(c(
  "all" = "All Interactions",
  "GG" = "Generalist-Generalist",
  "GS" = "Generalist-Specialist",
  "SS" = "Specialist-Specialist"
))

zlim <- ceiling(max(sp.pred_RO_df$predRO) * 10) / 10

# Create contour plot
RO_all <- ggplot(sp.pred_RO_df[sp.pred_RO_df$model == "all",], aes(x = ns, y = EN, z = predRO)) +
  geom_raster(aes(fill = predRO), interpolate = TRUE) +
  #scale_fill_viridis_c(name = "Predicted \nResource Overlap", limits = c(0,0.6)) + 
  scale_fill_gradientn(
    colors = MetBrewer::met.brewer("Hokusai3"),  # Try "Hokusai3" or "VanGogh1"
    name = "Resource-use\nOverlap\n",
    transform = "identity",
    limits = c(0, zlim), 
    na.value = MetBrewer::met.brewer("Hokusai3")[6]
  ) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_contour(breaks = breaks_seq, color = "white") +
  geom_text_contour(breaks = breaks_seq, stroke = 0.2) +
  labs(
    x = "Number of species",
    y = "Ecological Niche Overlap"
  ) +
  theme_minimal(base_size = 14) +  # Clean theme with larger text
  theme(
    legend.position = "right",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 60),
    legend.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 16)
  )

RO_int <- ggplot(sp.pred_RO_df[sp.pred_RO_df$model != "all",], aes(x = ns, y = EN, z = predRO)) +
  geom_raster(aes(fill = predRO), interpolate = TRUE) +
  #scale_fill_viridis_c(name = "Predicted \nResource Overlap", limits = c(0,0.6)) + 
  scale_fill_gradientn(
    colors = MetBrewer::met.brewer("Hokusai3"),  # Try "Hokusai3" or "VanGogh1"
    name = "Resource\nOverlap\n",
    transform = "identity",
    limits = c(0, zlim), 
    na.value = MetBrewer::met.brewer("Hokusai3")[6],
    guide = "none"
  ) +
  facet_wrap(.~model, ncol = 3, labeller = labeller_fn) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_contour(breaks = breaks_seq, color = "white") +
  geom_text_contour(breaks = breaks_seq, stroke = 0.2) +
  labs(
    x = "Number of species",
    y = "Ecological Niche Overlap",
  ) +
  theme_minimal(base_size = 14) +  # Clean theme with larger text
  theme(
    legend.position = "none",
    plot.margin = margin(t = 10, r = 100, b = 10, l = 40),
    strip.text = element_text(size = 16),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 16)
  )

fig5 <- ggarrange(RO_all, RO_int, nrow = 2, heights = c(1.75, 1), common.legend = FALSE, legend = "right",
                   labels = paste0("(", letters[1:2], ")"), font.label = list(size = 20))

ggsave("figures/Fig5.spatialGLMM_ROpred_interactions.pdf", plot = fig5, width = 10, height = 10)
