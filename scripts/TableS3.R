library(knitr)
library(kableExtra)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

tab_glm.spatial <- read.csv("results/overall_spatGLMmodel.csv", row.names = 1)
tab_glm.spatialGG <- read.csv("results/GG_spatGLMmodel.csv", row.names = 1)
tab_glm.spatialGS <- read.csv("results/GS_spatGLMmodel.csv", row.names = 1)
tab_glm.spatialSS <- read.csv("results/SS_spatGLMmodel.csv", row.names = 1)

tab_glm.spatial[is.na(tab_glm.spatial)] = ""
tab_glm.spatialGG[is.na(tab_glm.spatialGG)] = ""
tab_glm.spatialGS[is.na(tab_glm.spatialGS)] = ""
tab_glm.spatialSS[is.na(tab_glm.spatialSS)] = ""

tabA <- kbl(tab_glm.spatial, align=c(rep('c',times=3)), row.names = T, digits = 3, format = "latex") %>% 
  kable_styling(bootstrap_options = c("condensed"),
                full_width = F,
                position = "center") %>%
  kable_classic() %>%
  add_header_above(c("All interactions" = 4)) %>%
  pack_rows("Fixed effects (beta)", 1, 4) %>%
  pack_rows("Random effects", 5, 7) %>%
  pack_rows("Residual variance", 8, 8) %>%
  pack_rows("Likelihood values", 9, 9) 

tabB <- kbl(tab_glm.spatialGG, align=c(rep('c',times=3)), row.names = T, digits = 3, format = "latex") %>% 
  kable_styling(bootstrap_options = c("condensed"),
                full_width = F,
                position = "center") %>%
  kable_classic() %>%
  add_header_above(c("Generalist-Generalist" = 4)) %>%
  pack_rows("Fixed effects (beta)", 1, 4) %>%
  pack_rows("Random effects", 5, 7) %>%
  pack_rows("Residual variance", 8, 8) %>%
  pack_rows("Likelihood values", 9, 9) 

tabC <- kbl(tab_glm.spatialGS, align=c(rep('c',times=3)), row.names = T, digits = 3, format = "latex") %>% 
  kable_styling(bootstrap_options = c("condensed"),
                full_width = F,
                position = "center") %>%
  kable_classic() %>%
  add_header_above(c("Generalist-Specialist" = 4)) %>%
  pack_rows("Fixed effects (beta)", 1, 4) %>%
  pack_rows("Random effects", 5, 7) %>%
  pack_rows("Residual variance", 8, 8) %>%
  pack_rows("Likelihood values", 9, 9) 

tabD <- kbl(tab_glm.spatialSS, align=c(rep('c',times=3)), row.names = T, digits = 3, format = "latex") %>% 
  kable_styling(bootstrap_options = c("condensed"),
                full_width = F,
                position = "center") %>%
  kable_classic() %>%
  add_header_above(c("Specialist-Specialist" = 4)) %>%
  pack_rows("Fixed effects (beta)", 1, 4) %>%
  pack_rows("Random effects", 5, 7) %>%
  pack_rows("Residual variance", 8, 8) %>%
  pack_rows("Likelihood values", 9, 9) 

kables(list(tabA, tabB, tabC, tabD), "latex") %>% kable_classic() %>%
  save_kable("tables/TableS3.Spatial GLMMs.tex")
