library(NINA)
library(kableExtra)
library(knitr)
library(raster)
library(tidyr)
library(ggplot2)
library(ggpubr)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

dodge <- position_dodge(width = 0.4)
symnum.args <- list(cutpoints = c(0,0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "n.s."))

int.mat = read.csv("data/interaction_matrix.csv", row.names = 1)
int.mat[int.mat > 0] = 1

# store in memory
rasterOptions(todisk = FALSE)

amph_EN <- readRDS("./Rdata/amphENMs.RDS")
anem_EN <- readRDS("./Rdata/anemENMs.RDS")
amph_BC <- readRDS("./Rdata/amphEBMs.RDS")

# PA maps
amph_EN.PA <- raster_projection(prob_to_PA(amph_EN$maps, th = setNames(pmax(amph_EN$eval$threshold$threshold, 0.00001), 
                                                                       rownames(amph_EN$eval$threshold))))
amph_BC.PA <- raster_projection(prob_to_PA(amph_BC$maps, th = setNames(pmax(amph_BC$eval$threshold$threshold, 0.00001), 
                                                                       rownames(amph_BC$eval$threshold))))
anem_EN.PA <- raster_projection(prob_to_PA(anem_EN$maps, th = setNames(pmax(anem_EN$eval$threshold$threshold, 0.00001), 
                                                                       rownames(anem_EN$eval$threshold))))

EN_congruence <- data.frame()
BC_congruence <- data.frame()
env.scores <- cbind(region = amph_EN$clus[,3], amph_EN$env.scores)
for (s in names(amph_EN.PA)){
  # get associations
  Xvar = colnames(int.mat)[which(int.mat[s,] == 1)]
  
  ### Ecological niches
  # extract values from coordinates
  en_c <- raster::extract(amph_EN.PA[[s]], env.scores[, 2:3])
  en_c[is.na(en_c)] = 0
  
  bc_c <- raster::extract(amph_BC.PA[[s]], env.scores[, 2:3])
  bc_c[is.na(bc_c)] = 0
  
  a <- sapply(Xvar, function(i) raster::extract(anem_EN.PA[[i]], env.scores[, 2:3]))
  a[is.na(a)] = 0
  
  # sum host presences
  if (dim(a)[2] > 1){
    a <- rowSums(a, na.rm = T)
    a[a > 1] = 1
  }
  a <- as.numeric(a)
  
  # combine clownfish and host presences
  en_ca <- en_c+a  
  bc_ca <- bc_c+a  
  
  # compute metrics and add to data frame
  EN_congruence <- rbind(EN_congruence, data.frame(species = s,
                                                   PP = sum(en_ca == 2), 
                                                   PA = sum(en_c == 1 & a == 0), 
                                                   AP = sum(en_c == 0 & a == 1), 
                                                   AA =  sum(en_ca == 0), 
                                                   jaccard = prabclus::jaccard(cbind(en_c,a))[1,2]))
  
  BC_congruence <- rbind(BC_congruence, data.frame(species = s, 
                                                   PP = sum(bc_ca == 2), 
                                                   PA = sum(bc_c == 1 & a == 0), 
                                                   AP = sum(bc_c == 0 & a == 1), 
                                                   AA =  sum(bc_ca == 0), 
                                                   jaccard = prabclus::jaccard(cbind(bc_c,a))[1,2]))
}

eval_df <- rbind(cbind(species = rownames(amph_EN$eval$tab),
                       amph_EN$eval$tab, 
                       TOR = EN_congruence$PP/(EN_congruence$PP + EN_congruence$PA), 
                       TAR = EN_congruence$AA/(EN_congruence$AA + EN_congruence$AP),
                       model = "Ecological Niche"), 
                 cbind(species = rownames(amph_BC$eval$tab),
                       amph_BC$eval$tab, 
                       TOR = BC_congruence$PP/(BC_congruence$PP + BC_congruence$PA),
                       TAR = BC_congruence$AP/(BC_congruence$AA + BC_congruence$AP),
                       model = "Mutualism-refined Niche"))

eval_df$species <- factor(eval_df$species, levels = rownames(amph_EN$eval$tab))
eval_df$model <- factor(eval_df$model, levels = c("Ecological Niche", "Mutualism-refined Niche"))

eval_gather <- gather(eval_df[,c("species", "model", "TPR", "TNR", "TSS", "AUC", "TOR", "TAR")], metric, value, 3:8)
eval_gather$metric <- factor(eval_gather$metric, levels = c("TPR", "TNR", "TSS", "AUC", "TOR", "TAR"))

metric_names <- setNames(c("specificity", "sensitivity", 
                           "True Skill Statistic", "Area Under the ROC Curve", 
                           "True Co-occurrence Rate", "True Co-absence Rate"),
                  levels(eval_gather$metric))

me_all <- ggplot(eval_gather, aes(x = model, y = value, fill = model)) +  
  geom_boxplot() + 
  scale_fill_manual(values = c("#27AD81FF","#2171B5")) +
  scale_y_continuous(limits = c(0,1.2), breaks = seq(0,1, 0.25)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 14),
        axis.title = element_blank(),
        strip.text = element_text(size = 16)) +
  stat_compare_means(method = "t.test", comparisons = list(c("Ecological Niche", "Mutualism-refined Niche")), hide.ns = T, size = 6, paired = T,
                      symnum.args = symnum.args) +
  facet_wrap(. ~ metric, ncol = 2, labeller = labeller(metric = metric_names))

behavior = setNames(ifelse(rowSums(int.mat) > 2, "generalist", "specialist"), 
                    rownames(int.mat))

behavior_cols <- c("generalist" = "#1F3B73",  # Deep Blue  
                   "specialist" = "#F28E80")  # Soft Coral  

eval_gather$behavior <- sapply(eval_gather$species, function(i) behavior[i])
me_behavior <- ggplot(eval_gather, aes(x = behavior, y = value, fill = behavior)) +  
  geom_boxplot(color = rep(c("grey60", "grey20"), 12)) + 
  scale_fill_manual(values = behavior_cols) +
  scale_y_continuous(limits = c(0,1.1), breaks = seq(0,1, 0.25)) +
  theme_pubr() +
  theme(legend.position = "none") +
  stat_compare_means(method = "t.test", comparisons = list(c("generalist", "specialist")), hide.ns = T, size = 6,
                      symnum.args = symnum.args) +
  facet_wrap(model ~ metric, ncol = 6, labeller = labeller(metric = metric_names))

figS12 <- ggarrange(me_all, me_behavior, ncol = 2, widths = c(0.4, 1.2), labels = paste0("(", letters[1:2], ")"), font.label = list(size = 20))
ggsave("figures/FigS12.Models evaluation.pdf", plot = figS12, width = 25, height = 15, units = "in")

EN_congruence_region <- data.frame()
env.scores <- cbind(region = amph_EN$clus[,3], amph_EN$env.scores)
for (s in names(amph_EN.PA)){
  Xvar = colnames(int.mat)[which(int.mat[s,] == 1)]
  for (e in names(reverse_list(amph_EN$z.mod)[[s]])){
    c <- raster::extract(amph_EN.PA[[s]], env.scores[env.scores$region == e, 2:3])
    a <- sapply(Xvar, function(i) raster::extract(anem_EN.PA[[i]], env.scores[env.scores$region == e, 2:3]))
    c[is.na(c)] = 0
    a[is.na(a)] = 0
    if (dim(a)[2] > 1){
      a <- rowSums(a)
      a[a>1] = 1
    }
    a <- as.numeric(a)
    ca <- c+a
    PP <- sum(ca == 2)
    PA <- sum(c == 1 & a == 0)
    AP <- sum(c == 0 & a == 1)
    AA <- sum(ca == 0)
    jacc <- prabclus::jaccard(cbind(c,a))[1,2]
    EN_congruence_region <- rbind(EN_congruence_region, data.frame(species = s, region = e, PP = PP, PA = PA, AP = AP, AA = AA, jaccard = jacc))
  }
}

BC_congruence_region <- data.frame()
env.scores <- cbind(region = amph_BC$clus[,3], amph_BC$env.scores)
for (s in names(amph_BC.PA)){
  Xvar = colnames(int.mat)[which(int.mat[s,] == 1)]
  for (e in names(reverse_list(amph_BC$z.mod)[[s]])){
    c <- raster::extract(amph_BC.PA[[s]], env.scores[env.scores$region == e, 2:3])
    a <- sapply(Xvar, function(i) raster::extract(anem_EN.PA[[i]], env.scores[env.scores$region == e, 2:3]))
    c[is.na(c)] = 0
    a[is.na(a)] = 0
    if (dim(a)[2] > 1){
      a <- rowSums(a)
      a[a>1] = 1
    }
    a <- as.numeric(a)
    ca <- c+a
    PP <- sum(ca == 2)
    PA <- sum(c == 1 & a == 0)
    AP <- sum(c == 0 & a == 1)
    AA <- sum(ca == 0)
    jacc <- prabclus::jaccard(cbind(c,a))[1,2]
    BC_congruence_region <- rbind(BC_congruence_region, data.frame(species = s, region = e, PP = PP, PA = PA, AP = AP, AA = AA, jaccard = jacc))
  }
}

marine.regions <- read.csv("data/marine_regions.csv")
marine.regions$realm <- factor(marine.regions$realm, levels =  c("Western_Indo_Pacific" , "Central_Indo_Pacific", "Eastern_Indo_Pacific", "Temperate_Northern_Pacific", "Temperate_Australasia" ))
marine.regions$province <- factor(marine.regions$province, levels = unique(marine.regions[order(marine.regions$realm), "province"]))

EN_congruence_region$region <- factor(EN_congruence_region$region, levels = levels(marine.regions$province))
BC_congruence_region$region <- factor(BC_congruence_region$region, levels = levels(marine.regions$province))
EN_congruence_region$species <- as.factor(EN_congruence_region$species)
BC_congruence_region$species <- as.factor(BC_congruence_region$species)
EN_congruence_region$behavior <- sapply(EN_congruence_region$species, function(i) behavior[behavior$species == i, 2])
BC_congruence_region$behavior <- sapply(BC_congruence_region$species, function(i) behavior[behavior$species == i, 2])

theme <- theme_bw() +
  theme(axis.text.y = element_text(size=12, color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(size=10, angle = 45, hjust = 1, color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
  
data <- rbind(cbind(EN_congruence_region, model = "EN"), cbind(BC_congruence_region, model = "BC"))
data$model <- factor(data$model, levels = c("EN", "BC"))
data$species <- gsub("Amphiprion_", "A. ", data$species)
data$species <- gsub("Premnas_", "P. ", data$species)
data$region <- factor(data$region, levels = levels(marine.regions$province))
levels(data$region) <- gsub("_", " ", levels(data$region))

TORreg <- ggplot(data, aes(x = region, y = PP/(PP+PA))) + 
  geom_boxplot(fill = "grey80") +  
  labs(x = "", y = "True Co-occurrence Rate") +
  theme + 
  facet_wrap(model ~ ., ncol = 1, labeller = labeller(model = setNames(c("Ecological niche model", "Corrected niche model"),
                                                      c("EN", "BC")))) +
  theme(plot.margin=unit(c(1,1,1,2), 'cm'),
        strip.text = element_text(size=15))

TORsp <- ggplot(data, aes(x = species, y = PP/(PP+PA))) + 
  geom_boxplot(fill = "grey80") + 
  labs(x = "", y = "True Co-occurrence Rate") +
  theme +
  theme(plot.margin=unit(c(1,2,1,1), 'cm'),
        axis.text.x = element_text(size=10, angle = 45, hjust = 1, color = "black", face = "italic")) +
  facet_wrap(model ~ ., ncol = 1, labeller = labeller(model = setNames(c("Ecological niche model", "Corrected niche model"),
                                                      c("EN", "BC")))) +
  theme(strip.text = element_text(size=15))


ggarrange(TORreg, TORsp, labels = letters[1:2])
ggsave("suppfigures/FigS13.True Cooccurrence Rate per province and species.png", width = 3740, height = 2560, units = "px")
