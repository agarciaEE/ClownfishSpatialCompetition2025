# Load dependencieslibrary(ggplot2)
library(ggpubr)
library(ggplot2)
library(colorspace)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

# Load 3U estimates
ROU_data <- read.csv("results/ROU_data.csv")
symnum.args <- list(cutpoints = c(0,0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "n.s."))


SvsN_UP <- ggplot(ROU_data, aes(x = Restricted_niche, y = Restricted_spatial)) +
  scale_x_continuous(name = "Niche proportions",limits = c(0,1), expand = c(0.005,0.005))+
  scale_y_continuous(name = "Spatial proportions", limits = c(0,1), expand = c(0.005,0.005))+
  geom_point(col = "#E69F00", shape = 19) +
  geom_abline(intercept = 0, linetype =2, col =  "#E69F00", lwd = 1) +
  theme_classic() +
  border(color = "black", size = 1, linetype = NULL) +
  theme(axis.title.x = element_text(color = "black", size = 12, vjust = 1, hjust = 0.5),
        axis.title.y = element_text(color = "black", size = 12, vjust = 1, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        plot.margin = unit(c(2,2,2,2), "lines"),
        axis.ticks = element_blank(),
        axis.line.y = element_line(color = "black"),
        legend.position = "none",
        legend.text = element_text(color = "black", size = 14),
        legend.title = element_text(color = "black", size = 14),
        axis.line = element_line(colour = "black",
                                 size = 0.5, linetype = "solid"))
SvsN_SP <- ggplot(ROU_data, aes(x = Occupied_niche, y = Occupied_spatial)) +
  scale_x_continuous(name = "Niche proportions",limits = c(0,1), expand = c(0.005,0.005))+
  scale_y_continuous(name = "Spatial proportions", limits = c(0,1), expand = c(0.005,0.005))+
  geom_point(col = "#2171B5", shape = 19) +
  geom_abline(intercept = 0, linetype =2, col =  "#2171B5", lwd = 1) +
  theme_classic() +
  border(color = "black", size = 1, linetype = NULL) +
  theme(axis.title.x = element_text(color = "black", size = 12, vjust = 1, hjust = 0.5),
        axis.title.y = element_text(color = "black", size = 12, vjust = 1, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        plot.margin = unit(c(2,2,2,2), "lines"),
        axis.ticks = element_blank(),
        axis.line.y = element_line(color = "black"),
        legend.position = "none",
        legend.text = element_text(color = "black", size = 14),
        legend.title = element_text(color = "black", size = 14),
        axis.line = element_line(colour = "black",
                                 size = 0.5, linetype = "solid"))
SvsN_EP <- ggplot(ROU_data, aes(x = Unexploited_niche, y = Unexploited_spatial)) +
  scale_x_continuous(name = "Niche proportions",limits = c(0,1), expand = c(0.005,0.005))+
  scale_y_continuous(name = "Spatial proportions", limits = c(0,1), expand = c(0.005,0.005))+
  geom_point(col = "#27AD81FF", shape = 19) +
  geom_abline(intercept = 0, linetype =2, col =  "#27AD81FF", lwd = 1) +
  theme_classic() +
  border(color = "black", size = 1, linetype = NULL) +
  theme(axis.title.x = element_text(color = "black", size = 12, vjust = 1, hjust = 0.5),
        axis.title.y = element_text(color = "black", size = 12, vjust = 1, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        plot.margin = unit(c(2,2,2,2), "lines"),
        axis.ticks = element_blank(),
        axis.line.y = element_line(color = "black"),
        legend.position = "none",
        legend.text = element_text(color = "black", size = 14),
        legend.title = element_text(color = "black", size = 14),
        axis.line = element_line(colour = "black",
                                 size = 0.5, linetype = "solid"))


# test for difference from equal proportion of ROU_data in niche and spatial projections
median(ROU_data$Restricted_spatial-ROU_data$Restricted_niche)
GRp <- wilcox.test(ROU_data$Restricted_spatial-ROU_data$Restricted_niche, na.rm = T)
median(ROU_data$Occupied_spatial-ROU_data$Occupied_niche)
GOp <- wilcox.test(ROU_data$Occupied_spatial-ROU_data$Occupied_niche, na.rm = T)
median(ROU_data$Unexploited_spatial-ROU_data$Unexploited_niche, na.rm = TRUE)
GUp <- wilcox.test(ROU_data$Unexploited_spatial-ROU_data$Unexploited_niche, na.rm = T)


SvsN_USEdiff <- ggplot(ROU_data) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(aes(x = 1, y = Restricted_spatial-Restricted_niche), fill = "#E69F00", outlier.fill= "#E69F00", outlier.shape = 21, outlier.size = 2)  +
  geom_boxplot(aes(x = 2, y = Occupied_spatial-Occupied_niche), fill = "#2171B5", outlier.fill= "#2171B5", outlier.shape = 21, outlier.size = 2) +
  geom_boxplot(aes(x = 3, y = Unexploited_spatial-Unexploited_niche), fill = "#27AD81FF", outlier.fill= "#27AD81FF", outlier.shape = 21, outlier.size = 2) +
  scale_x_continuous("Proportions", breaks =  c(1,2,3), labels = c("Restricted", "Occupied", "Unexploited")) +
  scale_y_continuous("Geographically\n\nUnderrepresented                   Overrepresented", breaks = seq(-0.5,0.5,0.5), limits = c(-0.6,0.6)) +
  theme_classic() +
  annotate("text", x = 1, y = 0.55, label = symnum.args$symbols[which(GRp$p.value < symnum.args$cutpoints)[1]], size = 8) +
  annotate("text", x = 2, y = 0.55, label = symnum.args$symbols[which(GOp$p.value < symnum.args$cutpoints)[1]], size = 8) +
  annotate("text", x  = 3, y = 0.55, label = symnum.args$symbols[which(GUp$p.value < symnum.args$cutpoints)[1]], size = 8) +
  border(color = "black", size = 1, linetype = NULL) +
  theme(axis.title.x = element_text(color = "black", size = 16, vjust = 0, hjust = 0.5),
        axis.title.y = element_text(color = "black", size = 16, vjust = 3, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 16,  hjust = 0.5),
        axis.text.y = element_text(color = "black", size = 16,  hjust = 0.1),
        plot.margin = unit(c(2,2,2,2), "lines"),
        axis.ticks = element_blank(),
        axis.line.y = element_line(color = "black"),
        legend.position = "none",
        legend.text = element_text(color = "black", size = 14),
        legend.title = element_text(color = "black", size = 14),
        axis.line = element_line(colour = "black",
                                 size = 0.5, linetype = "solid"))

#plot
figS4 <- ggarrange(ggarrange(SvsN_UP, SvsN_SP, SvsN_EP, nrow = 3, labels = paste0("(", letters[1:3], ")"),
                                font.label = list(size = 20)), ggarrange(ggplot() + theme_void(), SvsN_USEdiff, ggplot() + theme_void(), nrow = 3, heights = c(0.2, 1, 0.2), labels = "(d)",  font.label = list(size = 20)),
                      ncol = 2, widths = c(1,2))

ggsave("figures/FigS4.ROU_nicheVSspatial.pdf", plot = figS4, width=10, height=10, units="in")
