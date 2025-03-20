# Load dependencies
library(ggpubr)
library(ggplot2)
library(colorspace)
library(tidyr)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

# load custom functions
source("./scripts/functions.R")

dodge <- position_dodge(width = 0.4)
symnum.args <- list(cutpoints = c(0,0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "n.s."))
colors <- c("#E69F00",  "#2171B5", "#27AD81FF")

D_df <- read.csv("results/niche_overlaps.csv")
D_df$pair = 1:nrow(D_df)

## select cuts
bins <- 20
cols <- c("darkred", "red", "white", "blue",  "darkblue")
colGradient <- colorRampPalette(cols)
cut.cols <- colGradient(bins)
cuts = cut(seq(-1,1,0.1),20)[-1]
names(cuts) <- sapply(cuts,function(t) cut.cols[which(as.character(t) == levels(cuts))])
breaks = seq(-1,1,0.1)

## env vs eff changes
env_eff.hist = data.frame(diff = rowSums(cbind(D_df$Denv,-D_df$Dmut), na.rm = T), pair = D_df$pair)
env_eff.hist[env_eff.hist$diff>0, "group"] = "ENV"
env_eff.hist[env_eff.hist$diff<0, "group"] = "EFF"
group_tags <- cut(env_eff.hist$diff,
                  breaks=breaks,
                  include.lowest=TRUE,
                  right=FALSE,
                  labels=cuts)
env_eff.hist$tag = group_tags
env_eff.hist$tag = factor(env_eff.hist$tag, levels = sort(unique(env_eff.hist$tag)))

env_eff_hist <- ggplot(na.exclude(env_eff.hist), aes(x=diff, fill=tag, col = tag)) +
  geom_histogram(breaks=seq(-1, 1, by=0.1),
                 alpha=0.8,
                 position="identity") +
  geom_histogram(data = env_eff.hist[is.na(env_eff.hist$group),], breaks=seq(-1.05, 1.05, by=0.1),
                 alpha=0.8,
                 position="identity", col = "grey", fill = "grey") +
  geom_hline(yintercept = 0, linetype = 1, linewidth = 1, col = "white") +
  annotate(x=median(env_eff.hist$diff), xend=median(env_eff.hist$diff), y=0, yend=150, lwd = 1,linetype = 1,  colour=if(median(env_eff.hist$diff)<0){"red"}else{"blue"}, geom="segment") +
  annotate(x=0, xend=0, y=0, yend=150, lwd = 0.5,linetype = 2,  colour="black", geom="segment") +
  scale_y_continuous("count", limits = c(-20,150) , expand = c(0,0)) +
  scale_x_continuous(name = "difference", labels = c(1,0.5,0,0.5,1), breaks = seq(-1,1,0.5), expand = c(0,0), limits = c(-1.4,1.4)) +
  scale_fill_manual(values = names(cuts)[cuts %in% env_eff.hist$tag]) +
  scale_color_manual(values =  names(cuts)[cuts %in% env_eff.hist$tag]) +
  theme_classic() +
  coord_cartesian(xlim=c(-1.4,1.4), ylim=c(-20,170)) +
  annotate(x=-1.4, xend=-1.4, y=0, yend=150, lwd = 1, colour="black", geom="segment") +
  annotate(x=-1, xend=1, y=-20, yend=-20, lwd = 1, colour="black",  geom="segment") +
  theme(panel.border=element_blank(),
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

## env vs resource-use changes
env_res.hist = data.frame(diff = rowSums(cbind(D_df$Denv,-D_df$Roverlap), na.rm = T), pair = D_df$pair)
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
  geom_histogram(data = env_res.hist[is.na(env_res.hist$group),], breaks=seq(-1.05, 1.05, by=0.1),
                 alpha=0.8,
                 position="identity", col = "grey", fill = "grey") +
  geom_hline(yintercept = 0, linetype = 1, linewidth = 1, col = "white") +
  annotate(x=median(env_res.hist$diff), xend=median(env_res.hist$diff), y=0, yend=150, lwd = 1,linetype = 1,  colour=if(median(env_res.hist$diff)<0){"red"}else{"blue"}, geom="segment") +
  annotate(x=0, xend=0, y=0, yend=150, lwd = 0.5,linetype = 2,  colour="black", geom="segment") +
  scale_y_continuous("count", limits = c(-20,150) , expand = c(0,0)) +
  scale_x_continuous(name = "difference", labels = c(1,0.5,0,0.5,1), breaks = seq(-1,1,0.5), expand = c(0,0), limits = c(-1.4,1.4)) +
  scale_fill_manual(values = names(cuts)[cuts %in% env_res.hist$tag]) +
  scale_color_manual(values =  names(cuts)[cuts %in% env_res.hist$tag]) +
  theme_classic() +
  coord_cartesian(xlim=c(-1.4,1.4), ylim=c(-20,170)) +
  annotate(x=-1.4, xend=-1.4, y=0, yend=150, lwd = 1, colour="black", geom="segment") +
  annotate(x=-1, xend=1, y=-20, yend=-20, lwd = 1, colour="black",  geom="segment") +
  theme(panel.border=element_blank(),
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

D_gather <- gather(D_df, type, D, 4:6)
D_gather$type <- factor(D_gather$type, levels = c("Denv", "Dmut", "Roverlap"))
D_gather$interaction <-   factor(D_gather$interaction, levels = c("generalist-generalist",  "generalist-specialist", "specialist-generalist", "specialist-specialist"))
D_gather$interaction[D_gather$interaction == "generalist-specialist"] = "specialist-generalist"
D_gather$interaction <-   factor(D_gather$interaction, levels = c("generalist-generalist", "specialist-generalist", "specialist-specialist"))

D_enveff <- gather(D_df, type, D, 4:5)
D_enveff$type <- factor(D_enveff$type, levels = c("Denv", "Dmut"))

D_effeco <- gather(D_df, type, D, c(4,6))
D_effeco$type <- factor(D_effeco$type, levels = c("Denv", "Roverlap"))

## ecological to mutualsm-refined
niche_overlap <- ggplot(D_gather[D_gather$type != "Roverlap",], aes(x = type, y = D, fill = type, col = type)) +
  geom_violin(col = NA) +
  geom_rect(data=NULL,aes(xmin=1,xmax=2,ymin=-Inf,ymax=Inf),
            fill="white" , col = NA, alpha = 0.5) +
  geom_line(data = D_enveff[D_df$pair %in% env_eff.hist[which(env_eff.hist$diff >0), "pair"],],
            aes(group = pair), alpha = 0.2, col = "blue") +
  geom_line(data = D_enveff[D_df$pair %in% env_eff.hist[which(env_eff.hist$diff <0), "pair"],],
            aes(group = pair), alpha = 0.6, col = "red") +
  geom_line(data = D_enveff[D_df$pair %in% env_eff.hist[which(env_eff.hist$diff ==0), "pair"],],
            aes(group = pair), alpha = 0.2, col = "grey40") +
  geom_boxplot(width = 0.05, outlier.colour = NA, outlier.size = 0.8,
               col = "white") +
  stat_summary(fun=median, geom="point", size=3, col = "white", shape = 19, position = dodge) +
  theme_classic() +
  scale_fill_manual("", labels = c("Ecological", "Mutualism-refined"), values= rev(as.vector(colors))[1:2]) +
  scale_color_manual("", labels = c("Ecological", "Mutualism-refined"), values= rev(as.vector(colors))[1:2]) +
  scale_y_continuous(name = "Schoener's D", expand = c(0,0),  seq(0,1,0.25), limits = c(-0.1,1.2)) +
  scale_x_discrete(labels = c("Ecological \nniche overlap", "Mutualism-refined\nniche overlap"), expand = c(0,0)) +
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
        plot.margin = margin(1, 1, 1, 1, "cm"),
        title = element_blank(),
        legend.text = element_text( color = "black", size=18),
        legend.title = element_text(color = "black", size = 18, vjust = 3),
        axis.title.x = element_blank(),
        axis.title.y = element_text( color = "black", size=24, vjust = 5),
        axis.text.x = element_text( color = "black", size=20, vjust = 0),
        axis.text.y = element_text( color = "black", size=20)) +
  stat_compare_means(comparisons = list(c("Denv", "Dmut")), paired = T, size  = 8, method = "wilcox.test", label = "..p.signif..", label.y =0.95, symnum.args = symnum.args)
fig4a <- ggarrange(niche_overlap, env_eff_hist,
                   ncol = 2, common.legend = F, widths = c(1,0.5))

## ecological to resource use
niche_overlap <- ggplot(D_gather[D_gather$type != "Dmut",], aes(x = type, y = D, fill = type, col = type)) +
  geom_violin(col = NA) +
  geom_rect(data=NULL,aes(xmin=1,xmax=2,ymin=-Inf,ymax=Inf),
            fill="white" , col = NA, alpha = 0.5) +
  geom_line(data = D_effeco[D_df$pair %in% env_res.hist[which(env_res.hist$diff >0), "pair"],],
            aes(group = pair), alpha = 0.2, col = "blue") +
  geom_line(data = D_effeco[D_df$pair %in% env_res.hist[which(env_res.hist$diff <0), "pair"],],
            aes(group = pair), alpha = 0.6, col = "red") +
  geom_line(data = D_effeco[D_df$pair %in% env_res.hist[which(env_res.hist$diff ==0), "pair"],],
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
        plot.margin = margin(1, 1, 1, 1, "cm"),
        title = element_blank(),
        legend.text = element_text( color = "black", size=18),
        legend.title = element_text(color = "black", size = 18, vjust = 3),
        axis.title.x = element_blank(),
        axis.title.y = element_text( color = "black", size=24, vjust = 5),
        axis.text.x = element_text( color = "black", size=20, vjust = 0),
        axis.text.y = element_text( color = "black", size=20)) +
  stat_compare_means(comparisons = list(c("Denv", "Roverlap")), paired = T, size  = 8, method = "wilcox.test", label = "..p.signif..", label.y =0.95, symnum.args = symnum.args) 
fig4b <- ggarrange(niche_overlap, env_res_hist,
                   ncol = 2, common.legend = F, widths = c(1,0.5))

fig4 <- ggarrange(fig4a, fig4b, labels = paste0("(", letters[1:2], ")"), font.label = list(size = 20),
                  nrow = 2, common.legend = F, heights = c(1,1))

ggsave("figures/Fig3.species_niche_overlaps.pdf", plot = fig4, width=15, height=12, units="in")

