# Load dependencies
library(ggplot2)
library(ggpubr)
library(colorspace)
library(ggthemes)
library(sf)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggnewscale)
library(ggExtra)
library(patchwork)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

INT.pred <- read.csv("results/spatial_results.csv")

# 1. BASELINE MAPS
# Load world map for reference
world <- ne_countries(scale = "medium", returnclass = "sf")

richness_palette <- function(a, b, max_a = 1, max_b = 1) {
  # Calculate RGB components using bilinear interpolation formulas
  red <- ((max_a * max_b) - max_a*b + a*b) / (max_a * max_b)
  green <- ((max_a * max_b) - max_a*b - max_b*a + 2*a*b) / (max_a * max_b)
  blue <- (max_a - a) / max_a
  
  # Clamp values between 0 and 1
  red <- pmin(pmax(red, 0), 1)
  green <- pmin(pmax(green, 0), 1)
  blue <- pmin(pmax(blue, 0), 1)
  
  # Convert to RGB color code
  rgb(red, green, blue)
}

# Define richness ranges
max_amph <- max(INT.pred$amph.richness)
max_anem <- max(INT.pred$anem.richness)

# Create a grid of values for the color scale
legend_data <- expand.grid(amph_richness = seq(0,max_amph, 0.1), anem_richness = seq(0, max_anem, 0.1))

# Apply the function to generate colors
legend_data <- legend_data %>%
  mutate(color = richness_palette(amph_richness, anem_richness, max_amph, max_anem))

INT.pred <- INT.pred %>%
  mutate(color = richness_palette(amph.richness, anem.richness, max_amph, max_anem))

# Create the richness map (Main Plot)
p_combined <- ggplot() +
  geom_sf(data = world, fill = "gray90", color = "grey60") +
  geom_tile(data = INT.pred, aes(x = x, y = y, fill = color), show.legend = FALSE) +
  scale_fill_identity() +  # Use precomputed colors
  coord_sf(xlim = range(INT.pred$x), ylim = range(INT.pred$y)) +
  scale_y_continuous(breaks = seq(-30, 30, 10), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(50, 200, 50), expand = c(0,0)) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude") +
  theme(plot.margin = margin(0, 0, 0, 0),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_line(color = "grey90"),
    text = element_text(size = 24)
  )

# Create aggregated data for margins
longitude_agg <- INT.pred %>%
  group_by(bin = round(x)) %>%
  summarise(amph.richness = mean(amph.richness, na.rm = TRUE),
            anem.richness = mean(anem.richness, na.rm = TRUE))

latitude_agg <- INT.pred %>%
  group_by(bin = round(y)) %>%
  summarise(amph.richness = mean(amph.richness, na.rm = TRUE),
            anem.richness = mean(anem.richness, na.rm = TRUE))

# Top margin plot: Density of species richness along longitude
top_margin <- ggplot(longitude_agg, aes(x = bin)) +
  geom_area(aes(y = amph.richness), fill = "#D96666", alpha = 0.5) +
  geom_area(aes(y = anem.richness), fill = "#6699D9", alpha = 0.5) +
  scale_x_continuous(breaks = seq(50, 200, 50), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,10,2), limits = c(0,max_amph), expand = c(0,0)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 240),
        panel.grid.major = element_line(color = "grey90"),
        panel.background = element_rect(fill ="white", color = NA))


# Right margin plot: Density of species richness along latitude
right_margin <- ggplot(latitude_agg, aes(x = bin)) +
  geom_area(aes(y = amph.richness), fill = "#D96666", alpha = 0.5) +
  geom_area(aes(y = anem.richness), fill = "#6699D9", alpha = 0.5) +
  theme_void() +
  scale_x_continuous(breaks = seq(-30, 30, 10), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,10,2), limits = c(0,10), expand = c(0,0)) +
  theme(plot.margin = margin(0, 0, 0, 0),
        panel.grid.major = element_line(color = "grey90"),
        panel.background = element_rect(fill ="white", color = NA)) +
  coord_flip()


# Create 2D Color Legend (Richness Gradient)
p_legend <- ggplot(legend_data, aes(x = anem_richness, y = amph_richness, fill = color)) +
  geom_tile() +
  geom_abline(color = "white", linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = seq(2,10,2), color = "white", linewidth = 0.25) +
  geom_vline(xintercept = seq(2,10,2), color = "white", linewidth = 0.25) +
  scale_fill_identity() + 
  scale_x_continuous(
    expand = c(0,0), 
    breaks = seq(0, 10, 2),
    position = "top"  
  ) +
  scale_y_continuous(
    expand = c(0,0), 
    breaks = seq(0, 10, 2),
    position = "right"  
  ) +
  theme_minimal() +
  labs( x = "Sea Anemone Richness", y = "Clownfish Richness") +
  theme(
    panel.background = element_rect(color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.text = element_text(size = 18),
    axis.text.y.right = element_text(margin = margin(l = 10)),  
    axis.title = element_text(size = 18),
    axis.title.y = element_text(
      angle = 90,  
      vjust = 0.5,  
      hjust = 0.5,  
      margin = margin(l = 10, r = 10)  
    ),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid = element_blank()
  )


# Combine plots using patchwork
SRplot <- wrap_plots(
  top_margin, 
  p_legend, 
  p_combined, 
  right_margin,
  ncol = 2, nrow = 2,
  widths = c(4, 1),
  heights = c(1.5, 3.5)
) 
SRplot


### (b)
generate_2D_gradient <- function(a, b) {
  # Ensure inputs are within the range [0, 1]
  if (a < 0 || a > 1 || b < 0 || b > 1) {
    stop("Values of 'a' and 'b' must be between 0 and 1.")
  }
  
  # Define the corner colors
  white <- c(1, 1, 1)   # RGB for white
  red <- c(1, 0, 0)     # RGB for red
  blue <- c(0, 0, 1)    # RGB for blue
  yellow <- c(1, 1, 0)  # RGB for yellow
  
  # Compute the interpolated color
  color <- (1 - a) * (1 - b) * white + 
    a * (1 - b) * red + 
    (1 - a) * b * blue + 
    a * b * yellow
  
  # Convert to hexadecimal format
  rgb(color[1], color[2], color[3], maxColorValue = 1)
}

INT.pred <- INT.pred %>%
  mutate(color = mapply(generate_2D_gradient, ENTotal, ROTotal))

max_a <- 1
max_b <- 1

# Create a grid of values for the color scale
legend_data <- expand.grid(ENTotal = seq(0, max_a, 0.01), ROTotal = seq(0, max_b, 0.01))

# Apply the function to generate colors
legend_data <- legend_data %>%
  dplyr::mutate(color = mapply(generate_2D_gradient, ENTotal, ROTotal))

# Create the richness map (Main Plot)
p_combined <- ggplot() +
  geom_sf(data = world, fill = "gray90", color = "grey60") +
  geom_tile(data = INT.pred, aes(x = x, y = y, fill = color), show.legend = FALSE) +
  scale_fill_identity() +  # Use precomputed colors
  scale_y_continuous(breaks = seq(-30, 30, 10), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(50, 200, 50), expand = c(0,0)) +
  coord_sf(xlim = range(INT.pred$x), ylim = range(INT.pred$y)) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude") +
  theme(plot.margin = margin(0, 0, 0, 0),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_line(color = "grey90"),
    text = element_text(size = 24)
  )

# Create aggregated data for margins
longitude_agg <- INT.pred %>%
  dplyr::group_by(bin = round(x)) %>%
  dplyr::summarise(ENTotal = mean(ENTotal, na.rm = TRUE),
            ROTotal = mean(ROTotal, na.rm = TRUE))

latitude_agg <- INT.pred %>%
  dplyr::group_by(bin = round(y)) %>%
  dplyr::summarise(ENTotal = mean(ENTotal, na.rm = TRUE),
            ROTotal = mean(ROTotal, na.rm = TRUE))

# Top margin plot: Density of species richness along longitude
top_margin <- ggplot(longitude_agg, aes(x = bin)) +
  geom_area(aes(y = ENTotal), fill = "#D96666", alpha = 0.5) +
  geom_area(aes(y = ROTotal), fill = "#6699D9", alpha = 0.5) +
  scale_x_continuous(breaks = seq(50, 200, 50), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1), expand = c(0,0)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 240),
        panel.grid.major = element_line(color = "grey90"),
        panel.background = element_rect(fill ="white", color = NA))


# Right margin plot: Density of species richness along latitude
right_margin <- ggplot(latitude_agg, aes(x = bin)) +
  geom_area(aes(y = ENTotal), fill = "#D96666", alpha = 0.5) +
  geom_area(aes(y = ROTotal), fill = "#6699D9", alpha = 0.5) +
  theme_void() +
  scale_x_continuous(breaks = seq(-30, 30, 10), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1), expand = c(0,0)) +
  theme(plot.margin = margin(0, 0, 0, 0),
        panel.grid.major = element_line(color = "grey90"),
        panel.background = element_rect(fill ="white", color = NA)) +
  coord_flip()

# Create 2D Color Legend (Overlap Gradient)
p_legend <- ggplot(legend_data, aes(x = ROTotal, y = ENTotal, fill = color)) +
  geom_tile() +
  geom_abline(color = "white", linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = seq(0.2,max_a,0.2), color = "white", linewidth = 0.25) +
  geom_vline(xintercept = seq(0.2,max_b,0.2), color = "white", linewidth = 0.25) +
  scale_fill_identity() + 
  scale_x_continuous(
    expand = c(0,0), 
    breaks = seq(0, max_a, 0.2),
    position = "top"  
  ) +
  scale_y_continuous(
    expand = c(0,0), 
    breaks = seq(0, max_b, 0.2),
    position = "right"  
  ) +
  theme_minimal() +
  labs(x = "Resource-use Overlap", y = "Ecological Niche Overlap") +
  theme(
    panel.background = element_rect(color = "black"),
    axis.ticks = element_line(linewidth = 1),
    axis.text = element_text(size = 18),
    axis.text.y.right = element_text(margin = margin(l = 10)),  
    axis.title = element_text(size = 18),
    axis.title.y = element_text(
      angle = 90,  
      vjust = 0.5,  
      hjust = 0.5,  
      margin = margin(l = 10, r = 10)  
    ),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid = element_blank()
  )

# Combine plots using patchwork
Overlap_plot <- wrap_plots(
  top_margin, 
  p_legend, 
  p_combined, 
  right_margin,
  ncol = 2, nrow = 2,
  widths = c(4, 1),
  heights = c(1.5, 3.5)
) 
Overlap_plot

fig4 <- ggarrange(SRplot, Overlap_plot, nrow = 2, labels = paste0("(", letters[1:2], ")"), font.label = list(size = 30), label.x = 0.1)

ggsave("figures/Fig4.geographical_richness_overlap.pdf", plot = fig4, width=24, height=19, units="in")

