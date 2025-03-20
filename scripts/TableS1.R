library(cramer)
library(kableExtra)
library(knitr)

# set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025"
setwd(wd)

dir.create("tables", showWarnings = FALSE)

# store in memory
rasterOptions(todisk = FALSE)

amph_EN <- readRDS("./Rdata/amphENMs.RDS")
anem_EN <- readRDS("./Rdata/anemENMs.RDS")
amph_BC <- readRDS("./Rdata/amphEBMs.RDS")

ROU_data <- read.csv("./results/ROU_data.csv")

for (i in 1:nrow(ROU_data)){
  reg <- ROU_data$region[i]
  sp <- ROU_data$species[i]
  data1 <- as.matrix(amph_EN$z.mod[[reg]][[sp]]$z.uncor)
  if (is.null(amph_BC$z.mod[[reg]][[sp]]$z.uncor)){
    data2 <- data1
    data2[] = 0
  } else{
    data2 <- as.matrix(amph_BC$z.mod[[reg]][[sp]]$z.uncor)
  }
  # Cramer test specific for nonparametric tests on multivariate data
  cr <- cramer.test(data1, data2)
  ROU_data[i, "Cramer_statistic"] = cr$statistic
  ROU_data[i, "Cramer_pvalue"] = format(cr$p.value, digits = 3, scientific = T)
  # Kolmogorov-Smirnov test to compare univariate distributions
  #ks <- ks.test(data1, data2)
  #ROU_data[i, "KS.statistic"] = ks$statistic
  #ROU_data[i, "KS.pvalue"] = ks$p.value
}

sum(as.numeric(ROU_data$Cramer_pvalue) < 0.05, na.rm = T) # number of significant tests
length(na.exclude(ROU_data$Cramer_pvalue)) # number of tests performed
sum(as.numeric(ROU_data$Cramer_pvalue) < 0.05, na.rm = T) / length(na.exclude(ROU_data$Cramer_pvalue))

data <- na.exclude(ROU_data[,c("region", "species", "Cramer_statistic", "Cramer_pvalue")])
data$region <- gsub("_", " ", data$region)
data$species <- gsub("Amphiprion_", "A. ", data$species)
data$species <- gsub("Premnas_", "A. ", data$species)

# Create LaTeX table output and store it in a variable
latex_table <- kable(
  list(data),
  align = rep("c", 4),
  row.names = FALSE,
  digits = 3,
  col.names = c("Province", "Species", "Statistic", "p-value"),
  format = "latex"
) %>%
  kable_styling(
    latex_options = c("condensed", "hold_position"),
    full_width = FALSE,
    position = "center"
  ) %>%
  kable_classic()

# Write the LaTeX table output to a .tex file
cat(latex_table, file = "tables/TableS1_Cramers_test.tex")
write.csv(data, "results/Cramers_test.csv", row.names = F)
 
