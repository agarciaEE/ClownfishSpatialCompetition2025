library(cramer)
library(kableExtra)
library(knitr)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

dir.create("tables", showWarnings = FALSE)

ROU_data <- read.csv("results/ROU_data.csv")

colnames(ROU_data)
parameters <- c("Restricted_niche", "Occupied_niche", "Unexploited_niche",
                "Restricted_spatial", "Occupied_spatial", "Unexploited_spatial",
                "Centroid_shift", "Environmental_shift", "Niche_dissimilarity")

gQ2 <- sQ2 <- gIQR <- sIQR <- statistic <- pval <- c()
for (p in parameters){
  gQ2 <- c(gQ2, median(ROU_data[ROU_data$behavior_reg == "generalist",p], na.rm = T))
  gIQR <- c(gIQR, IQR(ROU_data[ROU_data$behavior_reg == "generalist",p], na.rm = T))
  
  sQ2 <- c(sQ2, median(ROU_data[ROU_data$behavior_reg == "specialist",p], na.rm = T))
  sIQR <- c(sIQR, IQR(ROU_data[ROU_data$behavior_reg == "specialist",p], na.rm = T))
  test <- kruskal.test(ROU_data[,p], ROU_data$behavior_reg)
  
  statistic <- c(statistic, test$statistic)
  pval <- c(pval, test$p.value)
}

data <- data.frame(gQ2, gIQR, sQ2, sIQR, 
                   statistic, format(pval, digits = 3, scientific = T))

rownames(data) = gsub("_", " ", parameters)

latex_table <- kable(data, align=c(rep('c',times=6)),  digits = 3, 
      col.names =  c("Q2", "IQR", "Q2", "IQR", "Statistic", "p-value"),
      format = "latex") %>% 
  kable_styling(bootstrap_options = c("striped", "condensed"),
                full_width = F,
                position = "center") %>%
  kable_classic() %>% add_header_above(c(" " = 1, "Generalist" = 2, "Specialist" = 2, " " = 2)) 

cat(latex_table, file = "tables/TableS2_ROU_GvsS_summary.tex")
