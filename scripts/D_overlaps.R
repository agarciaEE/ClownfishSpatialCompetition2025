# Load required libraries
library(dplyr)
library(lme4)
library(dunn.test)

D_df <- read.csv("results/niche_overlaps.csv")
D_df$pair = 1:nrow(D_df)
colnames(D_df) <- c("region", "spa", "spb", 
                    'overlap_ecological', 'overlap_mutualism', 'overlap_resource_use',
                    "shared_hosts", "p_host_shared", "specialization", 
                    "c_host_shared", "pair")
D_df$species_pair <- paste0(D_df$spa, "-", D_df$spb)

# 1. Pearson correlation between niche overlap measures
cor_matrix <- cor(na.exclude(D_df[, c("overlap_ecological", "overlap_mutualism", "overlap_resource_use")]), 
                  method = "pearson")

# 2. Kruskal-Wallis test for each overlap measure
kruskal_ecological <- kruskal.test(overlap_ecological ~ specialization, data = D_df)
kruskal_mutualism <- kruskal.test(overlap_mutualism ~ specialization, data = D_df)
kruskal_resource <- kruskal.test(overlap_resource_use ~ specialization, data = D_df)

# 3. Post-hoc Dunn test with Bonferroni correction
dunn_ecological <- dunn.test(D_df$overlap_ecological, D_df$specialization, 
                             method = "bonferroni")
dunn_mutualism <- dunn.test(D_df$overlap_mutualism, D_df$specialization, 
                            method = "bonferroni")
dunn_resource <- dunn.test(D_df$overlap_resource_use, D_df$specialization, 
                           method = "bonferroni")

# 4. Linear mixed-effects models
lmm_ecological <- lmer(overlap_ecological ~ specialization + (1|region), data = na.exclude(D_df))
lmm_mutualism <- lmer(overlap_mutualism ~ specialization + (1|region), data = na.exclude(D_df))
lmm_resource <- lmer(overlap_resource_use ~ specialization + (1|region), data = na.exclude(D_df))

# Summary of results
summary(lmm_ecological)
summary(lmm_mutualism)
summary(lmm_resource)

# ANOVA for fixed effects
anova(lmm_ecological)
anova(lmm_mutualism)
anova(lmm_resource)

# Load required libraries
library(rstatix)
library(ggplot2)
library(dplyr)
library(corrplot)
library(lme4)
library(ggeffects)

# Assuming your models are named lmm_ecological, lmm_mutualism, and lmm_resource

# 1. Visualize fixed effects
plot_fixed_effects <- function(model, title) {
  fixef_data <- as.data.frame(summary(model)$coefficients)
  fixef_data$term <- rownames(fixef_data)
  
  ggplot(fixef_data, aes(x = term, y = Estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = Estimate - 1.96 * `Std. Error`, 
                      ymax = Estimate + 1.96 * `Std. Error`), 
                  width = 0.2) +
    coord_flip() +
    labs(title = title, x = "Term", y = "Estimate") +
    theme_minimal()
}


plot_fixed_effects(lmm_ecological, "Fixed Effects - Ecological Overlap")
plot_fixed_effects(lmm_mutualism, "Fixed Effects - Mutualism Overlap")
plot_fixed_effects(lmm_resource, "Fixed Effects - Resource Use Overlap")

# 2. Visualize random effects
plot_random_effects <- function(model, title) {
  ranef_data <- as.data.frame(ranef(model)$region)
  ranef_data$region <- rownames(ranef_data)
  
  # Extract standard errors of random effects
  re_summary <- as.data.frame(summary(model)$varcor$region)
  vc <- VarCorr(model)
  se <- sqrt(attr(vc, "sc")^2)

  ggplot(ranef_data, aes(x = reorder(region, `(Intercept)`), y = `(Intercept)`)) +
    geom_point() +
    geom_errorbar(aes(ymin = `(Intercept)` - 1.96 * se, 
                      ymax = `(Intercept)` + 1.96 * se), 
                  width = 0.2) +
    coord_flip() +
    labs(title = title, x = "Region", y = "Random Intercept") +
    theme_minimal()
}

plot_random_effects(lmm_ecological, "Random Effects by Region - Ecological Overlap")
plot_random_effects(lmm_mutualism, "Random Effects by Region - Mutualism Overlap")
plot_random_effects(lmm_resource, "Random Effects by Region - Resource Use Overlap")

# 3. Visualize predicted values
plot_predicted_values <- function(model, data, title) {
  # Get predictions and standard errors
  preds <- predict(model, newdata = data, re.form = NA, se.fit = TRUE)
  data$predicted <- preds$fit
  data$se <- preds$se.fit
  
  ggplot(data, aes(x = specialization, y = predicted)) +
    geom_boxplot() +
    geom_errorbar(aes(ymin = predicted - se, ymax = predicted + se), width = 0.2) +
    scale_y_continuous(limits = c(0,1)) +
    labs(title = title, x = "Specialization", y = "Predicted Overlap") +
    theme_minimal()
}

plot_predicted_values(lmm_ecological, niche_data, "Predicted Ecological Overlap by Specialization")
plot_predicted_values(lmm_mutualism, niche_data, "Predicted Mutualism Overlap by Specialization")
plot_predicted_values(lmm_resource, niche_data, "Predicted Resource Use Overlap by Specialization")
