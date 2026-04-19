# Load dependencies
library(ggpubr)
library(ggplot2)
library(colorspace)

# Set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025/"
setwd(wd)

# Load 3U estimates
ROU_data <- read.csv("results/ROU_data.csv")
marine_regions <- read.csv("data/marine_regions.csv")

## Regional effect on ROU_data parameters
#set realm and provs colors
col.realm = c("#21908CFF", "#FDE725FF", "#5DC863FF","#2C728EFF","#BB3754FF")
names(col.realm) = c("Western Indo Pacific" , "Central Indo Pacific", "Eastern Indo Pacific", "Temperate Northern Pacific", "Temperate Australasia" )
ROU_data$realm <- factor(sapply(ROU_data$region, function(r) gsub("_", " ", marine_regions$realm[marine_regions$province == r][1])), levels = names(col.realm))
ROU_data$region <- factor(ROU_data$region, levels = unique(ROU_data[order(ROU_data$realm), "region"]))
cluster_realms <-data.frame(x1=c(0.5,6.5,18.5,21.5,22.5), x2=c(6.5,18.5,21.5,22.5, 24.5), y1=rep(0,5), y2=rep(1,5), t= levels(ROU_data$realm), r=c(1,2,3,4,5))

RP_regions <- ggplot(ROU_data) +
  scale_x_discrete(labels = gsub("\\_", " ",levels(ROU_data$region)), expand = c(0,0)) +
  scale_y_continuous(name = "Proportion\n Restricted niche", breaks = seq(0,1,0.25), limits = c(-0.1,1.1), expand = c(0,0))+
  geom_rect(data=cluster_realms, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=col.realm, color=NA, alpha=0.3) +
  annotate("text", (cluster_realms$x1 + cluster_realms$x2)/2, 0.95, label = c("WIP", "CIP", "EIP", "TNP", "TA"), col = "grey30")  +
  geom_boxplot(aes(x = region, y = Restricted_niche), width = 0.5, outlier.colour = NA, fill = "#E69F00") +
  theme_classic() +
  coord_cartesian(xlim=c(0,24.5), ylim=c(-0.05,1.1)) +
  annotate(x=0, xend=0, y=0, yend=1, lwd = 1, colour="black", geom="segment") +
  annotate(x=1, xend=24, y=-0.05, yend=-0.05, lwd = 1, colour="black", geom="segment") +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks = element_line(size = 0.75),
        axis.ticks.length=unit(.2, "cm"),
        plot.margin = unit(c(3,2,0,2), "lines"),
        plot.title = element_text(size  = 20, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12, vjust = 2, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_blank()) +
  stat_compare_means(aes(x = region, y = Restricted_niche
  ), size = 4, label.y = 1.05, label.x = 2)

# stat_compare_means(comparisons = lapply(1:ncol(combn(levels(ROU_data$region),2)), function(i) combn(levels(ROU_data$region),2)[,i]), step.increase = 0.1, hide.ns = T, size = 6,  symnum.args = symnum.args,)
UP_regions <-  ggplot(ROU_data) +
  scale_x_discrete(labels = gsub("\\_", " ",levels(ROU_data$region)), expand = c(0,0)) +
  scale_y_continuous(name = "Proportion\n Unexploited niche", breaks = seq(0,1,0.25), limits = c(-0.1,1.1), expand = c(0,0))+
  geom_rect(data=cluster_realms, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=col.realm, color=NA, alpha=0.3) +
  annotate("text", (cluster_realms$x1 + cluster_realms$x2)/2, 0.95, label = c("WIP", "CIP", "EIP", "TNP", "TA"), col = "grey30")  +
  geom_boxplot(aes(x = region, y = Unexploited_niche), width = 0.5, outlier.colour = NA, fill = "#27AD81FF") +
  theme_classic() +
  coord_cartesian(xlim=c(0,24.5), ylim=c(-0.05,1.1)) +
  annotate(x=0, xend=0, y=0, yend=1, lwd = 1, colour="black", geom="segment") +
  annotate(x=1, xend=24, y=-0.05, yend=-0.05, lwd = 1, colour="black", geom="segment") +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks = element_line(size = 0.75),
        axis.ticks.length=unit(.2, "cm"),
        plot.margin = unit(c(3,2,0,2), "lines"),
        plot.title = element_text(size  = 20, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12, vjust = 2, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12,  vjust = 1, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_blank()) +
  stat_compare_means(aes(x = region, y = Unexploited_niche), size = 4, label.y = 1.05, label.x = 2)
#stat_compare_means(comparisons = lapply(levels(ROU_data$region)[-4], function(i) c(levels(ROU_data$region)[4], i)), label.y = seq(0.7,1,0.1), hide.ns = T, size = 9,  symnum.args = symnum.args)
OP_regions <- ggplot(ROU_data) +
  scale_x_discrete(labels = gsub("\\_", " ",levels(ROU_data$region)), expand = c(0,0)) +
  scale_y_continuous(name = "Proportion\n Occupied niche", breaks = seq(0,1,0.25), limits = c(-0.1,1.1), expand = c(0,0))+
  geom_rect(data=cluster_realms, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=col.realm, color=NA, alpha=0.3) +
  annotate("text", (cluster_realms$x1 + cluster_realms$x2)/2, 0.05, label = c("WIP", "CIP", "EIP", "TNP", "TA"), col = "grey30")  +
  geom_boxplot(aes(x = region, y = Occupied_niche), width = 0.5, outlier.colour = NA, fill = "#2171B5") +
  theme_classic() +
  coord_cartesian(xlim=c(0,24.5), ylim=c(-0.05,1.1)) +
  annotate(x=0, xend=0, y=0, yend=1, lwd = 1, colour="black", geom="segment") +
  annotate(x=1, xend=24, y=-0.05, yend=-0.05, lwd = 1, colour="black", geom="segment") +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks = element_line(size = 0.75),
        axis.ticks.length=unit(.2, "cm"),
        plot.margin = unit(c(3,2,0,2), "lines"),
        plot.title = element_text(size  = 20, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12, vjust = 2, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_blank()) +
  stat_compare_means(aes(x = region, y = Occupied_niche), size = 4, label.y = 1.05, label.x = 2)

SRP_regions <- ggplot(ROU_data) +
  scale_x_discrete(labels = gsub("\\_", " ",levels(ROU_data$region)), expand = c(0,0)) +
  scale_y_continuous(name = "Proportion\n Restricted distribution", breaks = seq(0,1,0.25), limits = c(-0.1,1.1), expand = c(0,0))+
  geom_rect(data=cluster_realms, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=col.realm, color=NA, alpha=0.3) +
  annotate("text", (cluster_realms$x1 + cluster_realms$x2)/2, 0.95, label = c("WIP", "CIP", "EIP", "TNP", "TA"), col = "grey30")  +
  geom_boxplot(aes(x = region, y = Restricted_spatial), width = 0.5, outlier.colour = NA, fill = "#E69F00") +
  theme_classic() +
  coord_cartesian(xlim=c(0,24.5), ylim=c(-0.05,1.1)) +
  annotate(x=0, xend=0, y=0, yend=1, lwd = 1, colour="black", geom="segment") +
  annotate(x=1, xend=24, y=-0.05, yend=-0.05, lwd = 1, colour="black", geom="segment") +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks = element_line(size = 0.75),
        axis.ticks.length=unit(.2, "cm"),
        plot.margin = unit(c(3,2,0,2), "lines"),
        plot.title = element_text(size  = 20, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12, vjust = 2, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_blank()) +
  stat_compare_means(aes(x = region, y = Restricted_spatial), size = 4, label.y = 1.05, label.x = 2)

SUP_regions <-  ggplot(ROU_data) +
  scale_x_discrete(labels = gsub("\\_", " ",levels(ROU_data$region)), expand = c(0,0)) +
  scale_y_continuous(name = "Proportion\n Unexploited distribution", breaks = seq(0,1,0.25), limits = c(-0.1,1.1), expand = c(0,0))+
  geom_rect(data=cluster_realms, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=col.realm, color=NA, alpha=0.3) +
  annotate("text", (cluster_realms$x1 + cluster_realms$x2)/2, 0.95, label = c("WIP", "CIP", "EIP", "TNP", "TA"), col = "grey30")  +
  geom_boxplot(aes(x = region, y = Unexploited_spatial), width = 0.5, outlier.colour = NA, fill = "#27AD81FF") +
  theme_classic() +
  coord_cartesian(xlim=c(0,24.5), ylim=c(-0.05,1.1)) +
  annotate(x=0, xend=0, y=0, yend=1, lwd = 1, colour="black", geom="segment") +
  annotate(x=1, xend=24, y=-0.05, yend=-0.05, lwd = 1, colour="black", geom="segment") +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks = element_line(size = 0.75),
        axis.ticks.length=unit(.2, "cm"),
        plot.margin = unit(c(3,2,0,2), "lines"),
        plot.title = element_text(size  = 20, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12, vjust = 2, hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 12,  vjust = 1, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_blank()) +
  stat_compare_means(aes(x = region, y = Unexploited_spatial), size = 4, label.y = 1.05, label.x = 2)

SOP_regions <- ggplot(ROU_data) +
  scale_x_discrete(labels = gsub("\\_", " ",levels(ROU_data$region)), expand = c(0,0)) +
  scale_y_continuous(name = "Proportion\n Occupied distribution", breaks = seq(0,1,0.25), limits = c(-0.1,1.1), expand = c(0,0))+
  geom_rect(data=cluster_realms, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=col.realm, color=NA, alpha=0.3) +
  annotate("text", (cluster_realms$x1 + cluster_realms$x2)/2, 0.05, label = c("WIP", "CIP", "EIP", "TNP", "TA"), col = "grey30")  +
  geom_boxplot(aes(x = region, y = Occupied_spatial), width = 0.5, outlier.colour = NA, fill = "#2171B5") +
  theme_classic() +
  coord_cartesian(xlim=c(0,24.5), ylim=c(-0.05,1.1)) +
  annotate(x=0, xend=0, y=0, yend=1, lwd = 1, colour="black", geom="segment") +
  annotate(x=1, xend=24, y=-0.05, yend=-0.05, lwd = 1, colour="black", geom="segment") +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks = element_line(size = 0.75),
        axis.ticks.length=unit(.2, "cm"),
        plot.margin = unit(c(3,2,0,2), "lines"),
        plot.title = element_text(size  = 20, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12, vjust = 2, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_blank()) +
  stat_compare_means(aes(x = region, y = Occupied_spatial), size = 4, label.y = 1.05, label.x = 2)

figS2 <- ggarrange(ggarrange(RP_regions, OP_regions, UP_regions, nrow = 3, heights = c(0.66, 0.66, 1)),
                                    ggarrange(SRP_regions, SOP_regions, SUP_regions, nrow = 3, heights = c(0.66, 0.66, 1)),
                                    ncol = 2, labels = paste0("(", letters[1:2], ")"), font.label = list(size = 20))

ggsave("figures/FigS2.ROUproportions_by_province.pdf", plot = figS3, width=15, height=12, units="in")
