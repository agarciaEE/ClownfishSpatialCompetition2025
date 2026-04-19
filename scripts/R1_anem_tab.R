# assumes: anem_EN$tab is a matrix/data.frame with species as rownames
# and marine regions as columns containing 0/1 presence-absence.

library(dplyr)
library(tibble)
library(knitr)
library(kableExtra)

library(knitr)
library(kableExtra)

# Convert table object to data.frame
anem_tab <- as.data.frame.matrix(anem_EN$tab)

# Add species column (rownames) and format for LaTeX
anem_tab <- data.frame(
  species = rownames(anem_tab),
  anem_tab,
  row.names = NULL
)

# Italicize species names (LaTeX / PDF output)
anem_tab$species <- gsub("_", " ", anem_tab$species)
colnames(anem_tab) <- gsub("_", " ", colnames(anem_tab))
colnames(anem_tab)[1] <- ""

anem_tab_t <- t(anem_tab) 

# Create kable table
kbl(
  anem_tab_t,
  booktabs = TRUE,
  escape = FALSE,
  caption = "Presence (1) / absence (0) of clownfish-hosting sea anemones across marine regions."
) %>%
  kable_styling(full_width = FALSE, font_size = 10) %>%
  row_spec(1, italic = TRUE) %>%
  column_spec(1, bold = TRUE)

