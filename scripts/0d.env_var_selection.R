# Load required libraries
library(FactoMineR)  
library(factoextra)  
library(caret)       
library(car)         
library(missMDA)     
library(ecospat)
library(dplyr)
library(tidyr)
library(vegan)
library(missMDA)
library(car)
library(usdm)

# Set working directory
wd = "~/Unil/Research/1_NicheCompetition/ClownfishSpatialCompetition2025/"
setwd(wd)

# Load and preprocess environmental data
env.vars_df <- read.csv("./data/env_reef_filtered_env.csv")

# convert to raster object
env_rast <- rast(env.vars_df)

# Remove variables with >5% missing values (adjust threshold as needed)
missing_threshold <- 0.05
missing_prop <- colMeans(is.na(env.vars_df))
env.vars_df <- env.vars_df[, missing_prop <= missing_threshold]

# Impute missing values using PCA-based imputation
ncp_opt <- estim_ncpPCA(env.vars_df, method = "Regularized")$ncp
env.vars_imputed <- imputePCA(env.vars_df, ncp = ncp_opt)$completeObs

# Perform PCA on standardized environmental dataset
pca_res <- PCA(scale(env.vars_imputed[,-c(1:2)]), graph = FALSE, ncp = ncol(env.vars_imputed) - 2)

# Compute the eigenvalues for all PCs
eig_values <- as.data.frame(get_eigenvalue(pca_res))

# Calculate the broken-stick distribution
stick_values <- rev(cumsum(rev(1/(1:ncol(env.vars_imputed)))))

# Compare PCA variance explained with broken-stick values
selected_ncp <- sum(eig_values$eigenvalue >= stick_values)  # Find the best PC cutoff

# Select variables with highest contributions to retained PCs
loading_matrix <- abs(pca_res$var$coord[, 1:selected_ncp])
selected_vars <- apply(loading_matrix, 2, function(x) names(sort(x, decreasing = TRUE)[1:3]))
selected_vars <- unique(as.vector(selected_vars))  # Unique top variables

# Subset environmental data with selected variables
env_selected <- env.vars_imputed[, colnames(env.vars_imputed) %in% selected_vars]

# function to find the variables set that maximizes PC1 and PC2 variance explained
find_best_pca_subset <- function(env_data, min_vars = 2, max_vars = NULL, max_combinations = 1000) {
  # Select only numeric columns and handle missing values
  env_data <- env_data %>%
    dplyr::select(where(is.numeric)) %>%
    mutate_all(~ ifelse(is.na(.), mean(., na.rm = TRUE), .))  # Impute missing values with column mean
  
  # Set max_vars if NULL (default to total number of variables)
  if (is.null(max_vars)) {
    max_vars <- ncol(env_data)
  }
  
  best_result <- list()
  
  for (n_vars in seq(min_vars, max_vars)) {
    best_result[[as.character(n_vars)]] <- list(variables = NULL, variance_explained = 0)
    
    # Generate variable combinations
    var_combinations <- combn(names(env_data), n_vars, simplify = FALSE)
  
    # Reduce the number of combinations if too many
    if (length(var_combinations) > max_combinations) {
      var_combinations <- sample(var_combinations, max_combinations)
    }
    
    best_vars <- NULL
    best_variance <- 0
    
    for (vars in var_combinations) {
      pca_result <- PCA(env_data[, vars], scale.unit = TRUE, graph = FALSE)
      variance_pc1_pc2 <- sum(pca_result$eig[1:2, 2])  # Variance explained by PC1 & PC2
      
      if (variance_pc1_pc2 > best_variance) {
        best_variance <- variance_pc1_pc2
        best_vars <- vars
      }
    }
    
    # Store best result
    best_result[[as.character(n_vars)]] <- list(
      variables = best_vars,
      variance_explained = best_variance
    )
  }
  
  return(best_result)
}

# Run the function to find the best variable subset
best_pca_subset <- find_best_pca_subset(as.data.frame(env_selected), min_vars = 6, max_vars = NULL)

summary <- data.frame(nvar = sapply(best_pca_subset, function(i) length(i$variables)),
                      variance_explained = sapply(best_pca_subset, function(i) i$variance_explained))
plot(summary$nvar, summary$variance_explained)

# Compute fist and second derivative 
second_derivative <- diff(diff(summary$variance_explained))

# Find the elbow point 
elbow_index <- which.min(second_derivative) + 1  # Adding 1 to match it with nvar

# Best index (number of variables)
best_n_vars <- summary$nvar[elbow_index]

var.sel <- best_pca_subset[[elbow_index]]$variables

# Final PCA visualization
pca_final <- PCA(scale(env_selected[,var.sel]), graph = FALSE)
fviz_pca_var(pca_final, col.var = "contrib", gradient.cols = c("blue", "red"), repel = TRUE)

final_var.sel <- as.data.frame(env.vars_imputed[, colnames(env.vars_imputed) %in% c("x", "y", var.sel)])
final_env_rast <- rast(final_var.sel)

# Save the final dataset
write.csv(final_var.sel, "./data/selected_environmental_variables.csv", row.names = FALSE)

# CHECK OCCURRENCES ON ENVIRONMENTAL PCA
### read occurrence datasets
amph.occ = read.csv("./data/amph_occ_env_final_dataset.csv")
anem.occ = read.csv("./data/anem_occ_env_final_dataset.csv")

### combine  occs  datasets
all.occ <- rbind(amph.occ, anem.occ)

# get environmental variables linked to species occurrences
all.occ_env = cbind(all.occ, terra::extract(final_env_rast, all.occ[,1:2])[,-1])

env_pca <- ade4::dudi.pca(final_var.sel[,-c(1:2)], scannf = FALSE, nf = 2)
occ_pca <- ade4::dudi.pca(all.occ_env[,-c(1:3)], scannf = FALSE, nf = 2)

pdf("./figures_new//pca_env_contributions.pdf")
par(mfrow = c(2,1))
ecospat.plot.contrib(env_pca$co, env_pca$eig)
ecospat.plot.contrib(occ_pca$co, occ_pca$eig)
dev.off()

occ.env.inn<-na.exclude(ecospat.sample.envar(dfsp=final_var.sel,colspxy=1:2,colspkept=1:2,dfvar=all.occ,colvarxy=1:2,colvar=3,resolution=max(res(final_env_rast))))

pdf("./figures_new/pca_occ_distribution.pdf")
par(mfrow = c(1,1))
ade4::s.class(env_pca$li[,1:2],
        fac= factor(rownames(env_pca$li) %in% rownames(occ.env.inn), 
                    levels = c("FALSE", "TRUE" ),
                    labels = c("background", "occs")), 
        col=c("red","blue"), 
        csta = 0,
        cellipse = 2,
        cpoint = .5,
        pch = 16)
dev.off()


