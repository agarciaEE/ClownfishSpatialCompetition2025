# Load dependencies
library(ggplot2)
library(ggpubr)
library(colorspace)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

dodge <- position_dodge(width = 0.4)
symnum.args <- list(cutpoints = c(0,0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "n.s."))

ROU_data <- read.csv("results/ROU_data.csv")

behavior_cols <- c("generalist" = "#1F3B73",  # Deep Blue  
                   "specialist" = "#F28E80")  # Soft Coral  

## Check G/S differences
#####
# tests on niche metrics
BE_GS <- ggplot(ROU_data, aes(x = behavior_reg, y = Niche_dissimilarity, col = behavior_reg, fill = behavior_reg)) +
  scale_y_continuous(name = "Niche dissimilarity", expand = c(0,0), breaks = seq(0,1,0.25), limits = c(-0.1,1.2))+
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
        plot.margin = unit(c(2,2,2,3), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 14, vjust = 5, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12,  vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none",
        legend.text = element_text(color = "black", size = 12, family = "mono" ),
        axis.line = element_blank()) +
  stat_compare_means( comparisons = list(c("generalist", "specialist")), hide.ns = T, size = 8,  symnum.args = symnum.args)


CS_GS <- ggplot(ROU_data, aes(x = behavior_reg, y = Centroid_shift, fill = behavior_reg, col = behavior_reg)) +
  scale_y_continuous(name = "Centroid Shift", expand = c(0,0), breaks = seq(0,ceiling(max(ROU_data$Centroid_shift, na.rm = T)), 
                                                                            ceiling(max(ROU_data$Centroid_shift, na.rm = T))/4), limits = c(-0.1*ceiling(max(ROU_data$Centroid_shift, na.rm = T)),
                                                                                                                                ceiling(max(ROU_data$Centroid_shift, na.rm = T))*1.2)) +
  scale_x_discrete(expand = c(0,0)) +
  geom_violin() +
  geom_boxplot(width = 0.02, outlier.colour = NA, col = darken(behavior_cols,0.5), fill = NA) +
  stat_summary(fun=median, geom="point", size=2, shape = 19, position = dodge, col = "white") +
  scale_fill_manual(values= behavior_cols) +
  scale_color_manual(values= behavior_cols) +
  theme_classic() +
  coord_cartesian(xlim=c(0.5,2.5), ylim=c(-0.1*ceiling(max(ROU_data$Centroid_shift, na.rm = T)),ceiling(max(ROU_data$Centroid_shift, na.rm = T))*1.2)) +
  annotate(x=0.5, xend=0.5, y=0, yend=ceiling(max(ROU_data$Centroid_shift, na.rm = T)), lwd = 1, colour="black", geom="segment") +
  annotate(x=1, xend=2, y=-0.1*ceiling(max(ROU_data$Centroid_shift, na.rm = T)), yend=-0.1*ceiling(max(ROU_data$Centroid_shift, na.rm = T)), lwd = 1, colour="black", geom="segment") +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks = element_line(size = 0.75),
        axis.ticks.length=unit(.2, "cm"),
        plot.margin = unit(c(2,2,2,3), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 14, vjust = 5, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12,  vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none",
        legend.text = element_text(color = "black", size = 12, family = "mono" ),
        axis.line = element_blank()) +
  stat_compare_means( comparisons = list(c("generalist", "specialist")), hide.ns = T, size = 8,  symnum.args = symnum.args)

ER_GS <- ggplot(ROU_data, aes(x = behavior_reg, y = Environmental_shift, fill = behavior_reg, col = behavior_reg)) +
  scale_y_continuous(name = "Environmental shift", expand = c(0,0), breaks = seq(0,1,0.25), limits = c(-0.1,1.3), labels = c(0, expression(pi/4), expression(pi/2), expression(pi*3/4), expression(pi))) +
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
        plot.margin = unit(c(2,2,2,3), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 14, vjust = 5, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12,  vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none",
        legend.text = element_text(color = "black", size = 12, family = "mono" ),
        axis.line = element_blank()) +
  stat_compare_means( comparisons = list(c("generalist", "specialist")), hide.ns = T, size = 8,  symnum.args = symnum.args)

figS5a <- ggarrange(CS_GS,  ER_GS, BE_GS,
                   ncol = 1, nrow = 3, labels = paste0("(", letters[1:3], ")"), font.label = list(size = 20))

SRP_GS <- ggplot(ROU_data, aes(x = behavior_reg, y = Restricted_spatial, fill = behavior_reg, col = behavior_reg)) +
  scale_y_continuous(name = "Restricted distribution", expand = c(0,0), breaks = seq(0,1,0.25), limits = c(-0.1,1.2))+
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
        plot.margin = unit(c(2,2,2,3), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 14, vjust = 5, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12,  vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none",
        legend.text = element_text(color = "black", size = 12, family = "mono" ),
        axis.line = element_blank()) +
  stat_compare_means(comparisons = list(c("generalist", "specialist")), hide.ns = T, size = 8,  symnum.args = symnum.args)

SOP_GS <- ggplot(ROU_data, aes(x = behavior_reg, y = Occupied_spatial, fill = behavior_reg, col = behavior_reg)) +
  scale_y_continuous(name = "Occupied distribution", expand = c(0,0), breaks = seq(0,1,0.25), limits = c(-0.1,1.2))+
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
        plot.margin = unit(c(2,2,2,3), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 14, vjust = 5, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12,  vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none",
        legend.text = element_text(color = "black", size = 12, family = "mono" ),
        axis.line = element_blank()) +
  stat_compare_means(comparisons = list(c("generalist", "specialist")), hide.ns = T, size = 8,  symnum.args = symnum.args)

SUP_GS <- ggplot(ROU_data, aes(x = behavior_reg, y = Unexploited_spatial, fill = behavior_reg, col = behavior_reg)) +
  scale_y_continuous(name = "Unexploited distribution", expand = c(0,0), breaks = seq(0,1,0.25), limits = c(-0.1,1.2))+
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
        plot.margin = unit(c(2,2,2,3), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 14, vjust = 5, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12,  vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none",
        legend.text = element_text(color = "black", size = 12, family = "mono" ),
        axis.line = element_blank()) +
  stat_compare_means(comparisons = list(c("generalist", "specialist")), hide.ns = T, size = 8,  symnum.args = symnum.args)

figS5b <- ggarrange(SRP_GS, SOP_GS, SUP_GS,
                   ncol = 1, nrow = 3, labels = paste0("(", letters[4:6], ")"), font.label = list(size = 20))

figS5 <- ggarrange(figS5a, figS5b) 

ggsave("figures/FigS5.extra_niche parameters_vs_behavior.pdf", plot = figS5, width=12, height=12, units="in")
