# Clownfish-Anemone Interaction Analysis

## Overview
This repository contains data, scripts, figures, and tables for analyzing ecological niche overlap (EN) vs. resource overlap (RO) in clownfish-anemone interactions across the Indo-Pacific. It includes statistical models, geographical analyses, and niche comparisons across interaction types.

## Repository Structure

### 1. `data/`
Contains the processed datasets used in analyses:
- **`amph_occ_env_final_dataset.csv`**: Processed occurrence and environmental data for clownfish.
- **`anem_occ_env_final_dataset.csv`**: Processed occurrence and environmental data for anemones.
- **`interaction_matrix.csv`**: Interaction matrix of clownfish and anemones.
- **`marine_regions.*`**: Shapefile components for marine region data.
- **`meow_ecos_df.csv`**: Processed marine ecoregions data.
- **`selected_environmental_variables.csv`**: Selected environmental variables for analyses.

### 2. `figures/`
Contains figures generated from the analyses:
- **`Fig1.framework_scheme.pdf`**: Conceptual framework of the study.
- **`Fig2.speciesROU_and_ROU_vs_behavior.pdf`**: Resource overlap distributions by species and behavior.
- **`Fig3.species_niche_overlaps.pdf`**: Niche overlap comparisons.
- **`Fig4.geographical_richness_overlap.pdf`**: Geographic distribution of richness overlap.
- **`Fig5.spatialGLMM_ROpred_interactions.pdf`**: Spatial GLMM predictions.
- **Supplementary Figures (`FigS1` to `FigS9`)**: Additional supporting figures.

### 3. `scripts/`
R scripts for data processing, modeling, and visualization:
- **Data Preparation:**
  - `0a_occ_datapreparation.R`, `0b_env_datapreparation.R`, `0c_marine_regions.R`, `0d.env_var_selection.R`
- **Modeling and Analysis:**
  - `1_reg_models.R`, `2_ROU_analysis.R`, `3_nicheOverlaps.R`, `4_spatialCompetition.R`, `spatialGLMM.R`
- **Statistical Tests and Model Evaluation:**
  - `D_overlaps.R`, `GLM_and_moranItests.R`, `Models_evaluation.R`
- **Visualization and Table Generation:**
  - `figure2.R`, `figure3.R`, `figure4.R`, `figure5.R`, `figureS1.R`, `figureS2.R`, `figureS3.R`, `figureS4.R`, `figureS5.R`, `figureS6.R`, `figureS7.R`
  - `TableS1.R`, `TableS2.R`, `TableS3.R`
- **Functions:**
  - `functions.R`, `spatialGLMM.R`

### 4. `tables/`
Contains tables summarizing statistical results:
- **`TableS1_Cramers_test.tex`**: Cramér's test results.
- **`TableS2_ROU_GvsS_summary.tex`**: Summary of ROU differences between generalists and specialists.
- **`TableS3.Spatial GLMMs.tex`**: Summary of spatial GLMM results.

## Usage
1. Clone the repository:
   ```sh
   git clone <repository_url>
   cd <repository>
   ```
2. Prepare the environment:
   - Ensure R and required packages are installed.
   - Use the scripts in `scripts/` for data processing and analysis.
3. Data Preparation (in case of running from the raw data (available soon in Dryad repository):
   - Run `0a_occ_datapreparation.R`, `0b_env_datapreparation.R`, `0c_marine_regions.R`, `0d.env_var_selection.R` in order.
3. Run models:
   - Execute `1_reg_models.R` to fit the models.
4. Carry out analysis:
   - Execute `2_ROU_analysis.R`, `3_nicheOverlaps.R`, `4_spatialCompetition.R` `spatialGLMM.R` to analyse the models.
   - Execute `D_overlaps.R`, `GLM_and_moranItests.R`, `Models_evaluation.R` for extra statistical tests.
4. Generate figures and tables:
   - Run `figureX.R` scripts to generate figures.
4. Check tables and results:
   - Run `tableX.R` scripts to generate tables.

## License

ClownfishSpatialCompetition2025 © 2025 by Alberto Garcia Jimenez is licensed under Creative Commons Attribution-ShareAlike 4.0 International.
To view a copy of this license, visit https://creativecommons.org/licenses/by-sa/4.0/ or see the LICENSE file.

## Citation
If using this dataset or analysis, please cite accordingly (citation details to be added).

## Contact
For inquiries, please contact agarcia26286@gmail.com
