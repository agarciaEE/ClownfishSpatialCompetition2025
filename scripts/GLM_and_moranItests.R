# load libraries
library(pgirmess)
library(NINA)
library(geoR)
library(ggplot2)
library(knitr)
library(kableExtra)
library(flextable)
library(sp)
library(spdep)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

# load custom functions
source("./scripts/functions.R")

# load dataset
INT.pred <- read.csv("./results/spatial_results.csv")

cor.test(INT.pred$ROTotal, INT.pred$amph.richness)
cor.test(INT.pred$ROTotal, INT.pred$ENTotal)
cor.test(INT.pred$amph.richness, INT.pred$ENTotal)

#quasi-poisson
fit.glm.qp <- glm(ROTotal ~ amph.richness * ENTotal, data = INT.pred, family = quasipoisson())
summary(fit.glm.qp)
ggplot() + geom_point(aes(fit.glm.qp$fitted, INT.pred$ROTotal))
# gaussian
fit.glm.ga <- glm(ROTotal ~ amph.richness * ENTotal, data = INT.pred, family = gaussian())
summary(fit.glm.ga)
ggplot() + geom_point(aes(fit.glm.ga$fitted, INT.pred$ROTotal))

# Compute correlogram of the residuals
nbc <- 20
# sample points for the correlogram as data set is too large
nrws <- sample(1:nrow(INT.pred), nrow(INT.pred) * 0.05, replace = F)
# correlogram
cor_r <- pgirmess::correlog(coords=INT.pred[nrws,c("x", "y")],
                            z=fit.glm.ga$residuals[nrws],
                            method="Moran", nbclass=nbc)

correlograms <- as.data.frame(cor_r)
correlograms$variable <- "residuals_glm" 

# Plot correlogram
ggplot(subset(correlograms, variable=="residuals_glm"), aes(dist.class, coef)) + 
  geom_hline(yintercept = 0, col="grey") +
  geom_line(col="steelblue") + 
  geom_point(data = correlograms[correlograms$p.value >= 0.05,], col="steelblue", size = 2, shape = 4) +
  geom_point(data = correlograms[correlograms$p.value < 0.05,], col="steelblue", size = 3, shape = 18) +
  xlab("distance") + 
  ylab("Moran's coefficient") +
  labs(title = "Moran's I correlogram") +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 9, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 7), 
        axis.title = element_text(size = 7))
ggsave("figures/FigS11.Morans_I_correlogram.pdf", width = 7, height = 5, units = "in")

# compute nearing neighbour data points
sp.coords <- coordinates(INT.pred[,c("x", "y")])
sp.sea.nb <- knn2nb(knearneigh(sp.coords, k = 8, longlat = T))
sp.sea.w <-nb2listw(sp.sea.nb, style="W", zero.policy = TRUE)

sp.moran.glm <- lm.morantest(fit.glm.ga, sp.sea.w)
print(sp.moran.glm)

# compute mean residuals of adjacent data points
sp.res_nb <- sapply(sp.sea.nb, function(x) mean(fit.glm.ga$residuals[x]))
sp.corr <- cor.test(fit.glm.ga$residuals, sp.res_nb, na.rm = T)

## moran's I tests
moran.mc(fit.glm.ga$residuals, sp.sea.w, nsim=999) # positive spatial autocorrelation
moran.mc(INT.pred$ROTotal, sp.sea.w, nsim=999) # positive spatial autocorrelation
## significant spatial Autocorrelation
## explore through tests of spatial regression models to detect which ones may be the more explanatory
sp.LM<-lm.LMtests(fit.glm.ga, sp.sea.w, test="all")
print(sp.LM)

## Variogram
### Too many data points to process the variogram Aggregate raster to reduce data points
INT.ras <- NINA::raster_projection(INT.pred[,-1])
INT.ras10 <- raster::aggregate(INT.ras, fact = 10, fun = mean, na.rm = T)
#### Spatial Analyses on sppatial point interactions dataset
spat_int_df <- raster::as.data.frame(INT.ras10, xy = T)

# remove NA's and points of no interaction
spat_int_df <- spat_int_df[!is.na(spat_int_df$amph.richness),]
spat_int_df <- spat_int_df[spat_int_df$amph.richness > 0,]
spat_int_df[is.na(spat_int_df)] = 0
nrow(spat_int_df)

sp10.EC_geo <- geoR::as.geodata(spat_int_df[,c("x","y","ROTotal")])
#plot(sp10.EC_geo, lowes=T)
MaxDist <- max(dist(spat_int_df[,c("x","y")])) / 2

sp10.Vario <- variog(sp10.EC_geo, max.dist = MaxDist, uvec = seq(0.01,MaxDist,0.2))
sp10.VarioMod_mtn<-variofit(sp10.Vario, cov.model = "matern")

pdf("figures/FigS10.Spatial_dependencies_and_variogram.pdf", width = 10, height = 5 )
par(mfrow = c(1,2))
plot(fit.glm.ga$residuals, sp.res_nb, xlab='Residuals', ylab='Mean adjacent residuals', main = "Spatial dependencies in the residuals")
abline(0,1, col = "red", lty = 2)
text(0.1, -0.3, paste("Pearson's correlation =", round(sp.corr$estimate, 3), "\n p-value < 2.2e-16"),
     cex= 0.6)
plot(sp10.Vario,pch=16, main = "Variogram")
lines(sp10.VarioMod_mtn,col="blue",lwd=2)
legend("bottomright", legend = "Matern model", col="blue", lwd=2, box.lty = 0, bg = "transparent")
dev.off()

correlograms$p.value <- format(correlograms$p.value, digits = 3, scientific = T)
kable(correlograms[,1:4], align=c(rep('c',times=4)), row.names = F, digits = 3, format = "latex") %>% 
  kable_styling(bootstrap_options = c("condensed"),
                full_width = F,
                position = "center") %>%
  kable_classic() %>%
  save_kable("tables/TableS3.Morans_I_spatial_correlogram.tex")
