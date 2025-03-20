# Load dependencies
library(ggplot2)
library(ggpubr)
library(colorspace)
library(tidyr)
library(ggbeeswarm)
library(dplyr)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

dodge <- position_dodge(width = 0.5)
symnum.args <- list(cutpoints = c(0,0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "n.s."))
colors <- c("#E69F00",  "#2171B5", "#27AD81FF")
# Define color palette
colors <- c("#E63946", "#F4A261", "#457B9D")  # Example colors: Red, Orange, Blue

D_df <- read.csv("results/niche_overlaps.csv")

D_gather <- gather(D_df, type, D, 4:6)

D_gather$interaction[D_gather$interaction == "specialist-specialist"] = "specialist-specialist"
D_gather$interaction <-   factor(D_gather$interaction, levels = c("generalist-generalist", "generalist-specialist", "specialist-specialist"))

D_gather$c_host_shared = sapply(D_gather$p_host_shared, function(i) if(i == 0){"non-sharing"} else if(i == 1){"sharing-all"} else {"partial-sharing"})
D_gather$c_host_shared <- factor(D_gather$c_host_shared, levels = c("non-sharing", "partial-sharing", "sharing-all"))

D_gather$type = sapply(D_gather$type, function(i) if(i == "Denv"){"Ecological"} else if(i == "Dmut"){"Mutualism-refined"} else {"Resource-use"})
D_gather$type <- factor(D_gather$type, levels = c("Ecological", "Mutualism-refined", "Resource-use"))

# Compute number of observations per group
n_labels_int <- D_gather %>%
  dplyr::group_by(type, interaction) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>%
  dplyr::select(interaction, n) %>%
  distinct()

n_labels_host <- D_gather %>%
  dplyr::group_by(type, c_host_shared) %>%
  dplyr::summarize(n = n(), .groups = "drop")%>%
  dplyr::select(c_host_shared, n) %>%
  distinct()

# Function to add number of observations
add_n_labels <- function(data) {
  geom_text(data = data, aes(x = 0.75, y = 1.2, label = paste0("n=", n)), 
            color = "grey40", size = 5, fontface = "bold")
}
# Define pairwise comparisons per facet
comparisons_list <- list(c("Ecological", "Mutualism-refined"), c("Ecological", "Resource-use"), c("Mutualism-refined", "Resource-use"))

# Combined plot for interaction types
D_int <- ggplot(D_gather, aes(x = type, y = D, fill = interaction)) +
  geom_violin(alpha = 0.6, color = NA, scale = "width") +  # Violin plot for distribution
  geom_quasirandom(size = 1.5, alpha = 0.6, shape = 21) +  
  geom_boxplot(width = 0.15, position = position_dodge(0.9), outlier.shape = NA, fill = "white", color = "black", alpha = 0.5) +  
  stat_summary(fun = median, geom = "point", size = 3, color = "black", shape = 18, position = position_dodge(0.9)) +  
  add_n_labels(n_labels_int) +
  theme_classic() +
  scale_fill_manual(name = "Interaction Type", values = colors) +
  scale_y_continuous(name = "Schoener's D", breaks = seq(0,1, 0.25), expand = c(0,0)) +
  coord_cartesian(ylim = c(0,1.25)) +
  facet_wrap(~ interaction, ncol = 3) +
  theme(plot.margin = margin(20, 20, 20, 60),
    legend.position = "right",
    legend.text = element_text(size = 16),  
    legend.title = element_text(size = 18, face = "bold"), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.text = element_blank()
  ) +
  stat_compare_means(data = D_gather %>% filter(interaction == "generalist-generalist"), 
                     method = "wilcox.test", comparisons = comparisons_list, 
                     label = "p.signif", label.y = c(0.9, 1, 1.1), size = 5, symnum.args = symnum.args) +
  stat_compare_means(data = D_gather %>% filter(interaction == "generalist-specialist"), 
                     method = "wilcox.test", comparisons = comparisons_list, 
                     label = "p.signif", label.y = c(0.9, 1, 1.1), size = 5, symnum.args = symnum.args) +
  stat_compare_means(data = D_gather %>% filter(interaction == "specialist-specialist"), 
                     method = "wilcox.test", comparisons = comparisons_list, 
                     label = "p.signif", label.y = c(0.9, 1, 1.1), size = 5, symnum.args = symnum.args) 
  
# Combined plot for host-sharing categories
D_host_shared <- ggplot(D_gather, aes(x = type, y = D, fill = c_host_shared)) +
  geom_violin(alpha = 0.6, color = NA, scale = "width") +
  geom_quasirandom(size = 1.5, alpha = 0.6, shape = 21) +
  geom_boxplot(width = 0.15, position = position_dodge(0.9), outlier.shape = NA, fill = "white", color = "black", alpha = 0.5) +
  stat_summary(fun = median, geom = "point", size = 3, color = "black", shape = 18, position = position_dodge(0.9)) +
  add_n_labels(n_labels_host) +
  scale_fill_manual(name = "Host Sharing", values = colors) +
  scale_y_continuous(name = "Schoener's D", breaks = seq(0,1, 0.25), expand = c(0,0)) +
  facet_wrap(~ c_host_shared, ncol = 3) +
  coord_cartesian(ylim = c(0,1.25)) +
  theme_classic() +
  theme(plot.margin = margin(20, 20, 20, 60),
    legend.position = "right",
    legend.text = element_text(size = 16),  
    legend.title = element_text(size = 18, face = "bold"), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.text = element_blank()
  )  +
  stat_compare_means(data = D_gather %>% filter(c_host_shared == "non-sharing"), 
                     method = "wilcox.test", comparisons = comparisons_list, 
                     label = "p.signif", label.y = c(0.9, 1, 1.1), size = 5, symnum.args = symnum.args) +
  stat_compare_means(data = D_gather %>% filter(c_host_shared == "partial-sharing"), 
                     method = "wilcox.test", comparisons = comparisons_list, 
                     label = "p.signif", label.y = c(0.9, 1, 1.1), size = 5, symnum.args = symnum.args) +
  stat_compare_means(data = D_gather %>% filter(c_host_shared == "sharing-all"), 
                     method = "wilcox.test", comparisons = comparisons_list, 
                     label = "p.signif", label.y = c(0.9, 1, 1.1), size = 5, symnum.args = symnum.args)
  
# Arrange plots
figS6 <- ggarrange(D_int, D_host_shared, nrow = 2, labels = c("(a)", "(b)"), font.label = list(size = 20))

# Save figure
ggsave("figures/FigS6.niche_overlap_vs_interaction_and_hostshared.pdf", figS6, width = 20, height = 12)
  