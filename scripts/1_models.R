
setwd("~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025/")

source("./scripts/functions.R")
rasterOptions(todisk = FALSE)

env <- read.csv("./data/selected_environmental_variables.csv")
clus <- read.csv("./data/marine_regions.csv")
anem_occ <- read.csv("./data/anem_occ_env_final_dataset.csv")
amph_occ <- read.csv("./data/amph_occ_env_final_dataset.csv")

# Sea anemones
anemENM = ENmodel(anem_occ, env, clus = clus, R = 100, extrapolate = TRUE, extrapolateIf = FALSE, FORCE = FALSE, eval = FALSE, 
                  type.pred = "probability", crs = 4326, res = 1)

anemENM$eval <- models_evaluation(anemENM, sample.pseudoabsences = TRUE, 
                                  alpha = 0.5, transformation = "none", trans.param = 20, 
                                  plot = FALSE, rep = 100, th = NULL,  best.th = "TSS")

saveRDS(anemENM, "./Rdata/anemENMs.RDS")

# Clownfish
amphENM = ENmodel(amph_occ, env, clus = clus, R = 100, extrapolate = TRUE, extrapolateIf = FALSE, FORCE = TRUE, eval = FALSE, 
                  type.pred = "probability", crs = 4326, res = 1)


amphENM$eval <- models_evaluation(amphENM, sample.pseudoabsences = TRUE, 
                                  alpha = 0.5, transformation = "none", trans.param = 20, 
                                  plot = FALSE, rep = 100, th = NULL,  best.th = "TSS")

saveRDS(amphENM, "./Rdata/amphENMs.RDS")

# perform biotic corrections
int.mat <- read.csv("./data/interaction_matrix.csv", row.names = 1)
int.mat[int.mat > 0] = 1

amphEBM <- BCmodel(amphENM, anemENM,  D = 1, A.matrix = int.mat, method = "densities",
                        type.pred = "probability", relative.niche = FALSE, cor = FALSE, eval = FALSE)

amphEBM$eval <- models_evaluation(amph_EBM_new, sample.pseudoabsences = TRUE, 
                                       alpha = 0.5, transformation = "none", trans.param = 20, 
                                       plot = FALSE, rep = 100, th = NULL,  best.th = "TSS")

saveRDS(amphEBM, "./Rdata/amphEBMs.RDS")

