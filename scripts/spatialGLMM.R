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
library(pgirmess)

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
nrow(spat_int_df)

nb_spat_int <- nb2listw(knn2nb(knearneigh(spat_int_df[, c("x", "y")], k = 8)))

avail_thr <- parallel::detectCores(logical=FALSE) - 1L

fit.glm.spatial <- spaMM::fitme(ROTotal ~ amph.richness * ENTotal + Matern(1|x+y), data = spat_int_df, family = gaussian(),
                                control.HLfit=list(NbThreads=max(avail_thr, 1L)))
summary(fit.glm.spatial)

coefs <- as.data.frame(summary(fit.glm.spatial)$beta_table)
lower <- coefs[,'Estimate'] - 1.96*coefs[, 'Cond. SE']
upper <- coefs[,'Estimate'] + 1.96*coefs[, 'Cond. SE']
dfs <- length(fit.glm.spatial$y) - fit.glm.spatial$dfs$pforpv
p_values <- 2 * pt(-abs(coefs$`t-value`), dfs, lower.tail = TRUE) 
CI <- data.frame(lower, upper, row.names = row.names(coefs))
CI$sig <- sapply(1:nrow(CI), function(i) if (0 >= CI$lower[i] & 0 <= CI$upper[i]){ "n.s."} else { "sig" })
coefs <- cbind(coefs, dfs, p_values, CI)

# Check spatial autocorrelation
# Extract residuals from the model
residuals_glm.spatial <- resid(fit.glm.spatial)
# Compute Moran's I to check remaining spatial autocorrelation
moran.test(residuals_glm.spatial, listw = nb_spat_int)

# Compute correlogram of the residuals
nbc <- 50
# correlogram
cor_r <- pgirmess::correlog(coords=spat_int_df[,c("x", "y")],
                            z=residuals_glm.spatial,
                            method="Moran", nbclass=nbc)

correlograms <- as.data.frame(cor_r)
correlograms$variable <- "residuals_glmm" 
# Plot correlogram
p_cor1 <- ggplot(subset(correlograms, variable=="residuals_glmm"), aes(dist.class, coef)) + 
  geom_hline(yintercept = 0, col="grey") +
  geom_line(col="steelblue") + 
  geom_point(data = correlograms[correlograms$p.value >= 0.05,], col="steelblue", size = 2, shape = 4) +
  geom_point(data = correlograms[correlograms$p.value < 0.05,], col="steelblue", size = 3, shape = 18) +
  xlab("distance [in decimal degrees]") + 
  ylab("Moran's coefficient") +
  labs(title = "Moran's I correlogram") +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 16))

nu <- fit.glm.spatial$ranFix$corrPars$`1`$nu
rho <- fit.glm.spatial$ranFix$corrPars$`1`$rho
distances <- dist(spat_int_df[,c("x","y")])
MaternCorrelation <- MaternCorr(distances, nu = nu, rho = rho)
p_cor2 <- ggplot(data = data.frame(Distance = as.numeric(distances), Correlation = as.numeric(MaternCorrelation)), 
       aes(x = Distance, y = Correlation)) +
  geom_point(color = "red", size = 2, alpha = 0.8) +  # Red points with slight transparency
  geom_smooth(method = "loess", color = "black", linetype = "dashed", linewidth = 1, se = FALSE) +  # Add smooth trend line
  labs(
    x = "Distance between pairs of locations [decimal degrees]",
    y = "Estimated correlation",
    title = "Spatial Correlation Decay"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

figS9 <- ggarrange(p_cor1, p_cor2, ncol = 2, labels = c("(a)", "(c)"), font.label = list(size = 20))
ggsave("figures/FigS9.glmm_residuals_Morans_I_correlogram.pdf", plot = figS9, width = 15, height = 10, units = "in")

correlograms$p.value <- format(correlograms$p.value, digits = 3, scientific = T)
kable(correlograms[,1:4], align=c(rep('c',times=4)), row.names = F, digits = 3, format = "latex") %>% 
  kable_styling(bootstrap_options = c("condensed"),
                full_width = F,
                position = "center") %>%
  kable_classic() %>%
  save_kable("tables/TableS3.Morans_I_spatial_correlogram.tex")

# Model validation
##################
# take 20% to act as validation set
set.seed(1)
validation_rows <- sample(1:nrow(spat_int_df), floor(nrow(spat_int_df)*0.2))
spat_int_df_train <- spat_int_df[-validation_rows,]
spat_int_df_valid <- spat_int_df[validation_rows,]

# Fit model using 80%
fit.glm.spatial_validation <-  spaMM::fitme(ROTotal ~ amph.richness * ENTotal + Matern(1|x+y), data = spat_int_df_train, family = gaussian(),
                                            control.HLfit=list(NbThreads=max(avail_thr, 1L)))

# Evaluation
predictions_validation <- predict(fit.glm.spatial_validation, spat_int_df_valid)
ggplot() + geom_point(aes(as.vector(predictions_validation), spat_int_df_valid$ROTotal)) +
  geom_abline(col = "red", linewidth = 1.25)
Metrics::mse(predictions_validation, spat_int_df_valid$ROTotal)

# GG
fit.glm.spatialGG <- spaMM::fitme(ROGG ~ nG * ENGG + Matern(1|x+y), data = spat_int_df[spat_int_df$nG > 1,], family = gaussian(),
                                  control.HLfit=list(NbThreads=max(avail_thr, 1L)))
summary(fit.glm.spatialGG)
# Extract residuals from the model
residuals_glm.spatialGG<- resid(fit.glm.spatialGG)
# Compute Moran's I to check remaining spatial autocorrelation
nb_spat_intGG <- nb2listw(knn2nb(knearneigh(spat_int_df[spat_int_df$nG > 1, c("x", "y")], k = 8)))
moran.test(residuals_glm.spatialGG, listw = nb_spat_intGG)

sims <-  simulateResiduals(fit.glm.spatialGG)
plot(sims)

coefs <- as.data.frame(summary(fit.glm.spatialGG)$beta_table)
lower <- coefs[,'Estimate'] - 1.96*coefs[, 'Cond. SE']
upper <- coefs[,'Estimate'] + 1.96*coefs[, 'Cond. SE']
dfs <- length(fit.glm.spatialGG$y) - fit.glm.spatialGG$dfs$pforpv
p_values <- 2 * pt(-abs(coefs$`t-value`), dfs, lower.tail = TRUE) 
CI <- data.frame(lower, upper, row.names = row.names(coefs))
CI$sig <- sapply(1:nrow(CI), function(i) if (0 >= CI$lower[i] & 0 <= CI$upper[i]){ "n.s."} else { "sig" })
coefs <- cbind(coefs, dfs, p_values, CI)

# Model independently by interation type
########################################
# GS
fit.glm.spatialGS <- spaMM::fitme(ROGS ~ amph.richness * ENGS + Matern(1|x+y), data = spat_int_df, family = gaussian(),
                                  control.HLfit=list(NbThreads=max(avail_thr, 1L)))
summary(fit.glm.spatialGS)
# Extract residuals from the model
residuals_glm.spatialGS<- resid(fit.glm.spatialGS)
# Compute Moran's I to check remaining spatial autocorrelation
moran.test(residuals_glm.spatialGS, listw = nb_spat_int)
sims <- simulateResiduals(fit.glm.spatialGS)
plot(sims)

coefs <- as.data.frame(summary(fit.glm.spatialGS)$beta_table)
lower <- coefs[,'Estimate'] - 1.96*coefs[, 'Cond. SE']
upper <- coefs[,'Estimate'] + 1.96*coefs[, 'Cond. SE']
dfs <- length(fit.glm.spatialGS$y) - fit.glm.spatialGS$dfs$pforpv
p_values <- 2 * pt(-abs(coefs$`t-value`), dfs, lower.tail = TRUE) 
CI <- data.frame(lower, upper, row.names = row.names(coefs))
CI$sig <- sapply(1:nrow(CI), function(i) if (0 >= CI$lower[i] & 0 <= CI$upper[i]){ "n.s."} else { "sig" })
coefs <- cbind(coefs, dfs, p_values, CI)

# SS
fit.glm.spatialSS <- spaMM::fitme(ROSS ~ nS * ENSS + Matern(1|x+y), data = spat_int_df[spat_int_df$nS > 1,], family = gaussian(),
                                  control.HLfit=list(NbThreads=max(avail_thr, 1L)))
summary(fit.glm.spatialSS)
# Extract residuals from the model
residuals_glm.spatialSS<- resid(fit.glm.spatialSS)
# Compute Moran's I to check remaining spatial autocorrelation
nb_spat_intSS <- nb2listw(knn2nb(knearneigh(spat_int_df[spat_int_df$nS > 1, c("x", "y")], k = 8)))
moran.test(residuals_glm.spatialSS, listw = nb_spat_intSS)
sims <- simulateResiduals(fit.glm.spatialSS)
plot(sims)

coefs <- as.data.frame(summary(fit.glm.spatialSS)$beta_table)
lower <- coefs[,'Estimate'] - 1.96*coefs[, 'Cond. SE']
upper <- coefs[,'Estimate'] + 1.96*coefs[, 'Cond. SE']
dfs <- length((fit.glm.spatialSS)$y) - (fit.glm.spatialSS)$dfs$pforpv
p_values <- 2 * pt(-abs(coefs$`t-value`), dfs, lower.tail = TRUE) 
CI <- data.frame(lower, upper, row.names = row.names(coefs))
CI$sig <- sapply(1:nrow(CI), function(i) if (0 >= CI$lower[i] & 0 <= CI$upper[i]){ "n.s."} else { "sig" })
coefs <- cbind(coefs, dfs, p_values, CI)

saveRDS(fit.glm.spatial, "Rdata/all_spatialGLMM.RDS")
saveRDS(fit.glm.spatialGG, "Rdata/GG_spatialGLMM.RDS")
saveRDS(fit.glm.spatialGS, "Rdata/GS_spatialGLMM.RDS")
saveRDS(fit.glm.spatialSS, "Rdata/SS_spatialGLMM.RDS")

as.table_GLMM <- function(x, row.names = NULL){
  sum <- invisible(summary(x))
  
  tab <- rbind(data.frame(round(sum$beta_table, 3)),
               c(round(x$ranFix$corrPars$`1`$nu, 3), "", ""),
               c(round(x$ranFix$corrPars$`1`$rho, 3), "", ""),
               c(round(x$ranFix$lambda, 3), "", ""),
               c(round(x$ranFix$phi,6), "", ""),
               c(round(sum$likelihoods, 3), "", ""))
  if (is.null(row.names)){
    rownames(tab)[5:9] <- c("nu", "rho", "lambda (x+y)", "phi", "logLik")
  } else {
    rownames(tab)[-1]<- c(row.names, "nu", "rho", "lambda (x+y)", "phi", "logLik")
  }
  return(tab)
}

tab_glm.spatial <- as.table_GLMM(fit.glm.spatial, row.names = c("num. species", "Ecological overlap", "num. species x Ecological overlap"))
tab_glm.spatialGG <- as.table_GLMM(fit.glm.spatialGG, row.names = c("num. generalists", "Ecological overlap", "num. species x Ecological overlap"))
tab_glm.spatialGS <- as.table_GLMM(fit.glm.spatialGS, row.names = c("num. species", "Ecological overlap", "num. species x Ecological overlap"))
tab_glm.spatialSS <- as.table_GLMM(fit.glm.spatialSS, row.names = c("num. specialists", "Ecological overlap", "num. species x Ecological overlap"))

write.csv(tab_glm.spatial, "results/overall_spatGLMmodel.csv")
write.csv(tab_glm.spatialGG, "results/GG_spatGLMmodel.csv")
write.csv(tab_glm.spatialGS, "results/GS_spatGLMmodel.csv")
write.csv(tab_glm.spatialSS, "results/SS_spatGLMmodel.csv")
