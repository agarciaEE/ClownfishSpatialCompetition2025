install.packages(c("igraph", "ggraph", "tidygraph", "ggplot2", "dplyr"))

library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(dplyr)

# Create an edge list for species interactions
INT.data_reg <- INT.data %>%
  filter(region == "Sunda_Shelf") %>%
  select(sp1, sp2, region, interaction, Denv, Rove) %>%
  group_by(sp1, sp2, region, interaction) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")


# Create a node list with unique species
nodes <- unique(c(INT.data_reg$sp1, INT.data_reg$sp2))
nodes_df <- data.frame(name = nodes)

# Convert edges to an igraph object
graph <- graph_from_data_frame(INT.data_reg, vertices = nodes_df, directed = FALSE)

ggraph(graph, layout = "fr") +  
  geom_edge_link(aes(width = Denv, color = interaction), alpha = 0.6) +  
  geom_node_point(size = 5, color = "black") +  
  geom_node_text(aes(label = name), repel = TRUE) +  
  scale_edge_width(range = c(0.5, 3)) +  
  scale_edge_color_manual(values = c("GG" = "red", "GS" = "green", "SS" = "blue")) +  
  theme_void() +  
  ggtitle("Species Interaction Network")


# Assign behavior categories
nodes_df$behavior <- ifelse(nodes_df$name %in% behavior$species[behavior$mutualism == "generalist"], "Generalist", "Specialist")

# Create a behavior-based graph
graph <- as_tbl_graph(graph) %>%
  mutate(behavior = nodes_df$behavior[match(name, nodes_df$name)])

ggraph(graph, layout = "fr") +  
  geom_edge_link(aes(width = Denv, color = interaction), alpha = 0.6) +  
  geom_node_point(aes(color = behavior), size = 5) +  
  geom_node_text(aes(label = name), repel = TRUE) +  
  scale_edge_width(range = c(0.5, 3)) +  
  scale_edge_color_manual(values = c("GG" = "red", "GS" = "green", "SS" = "blue")) +  
  scale_color_manual(values = c("Generalist" = "orange", "Specialist" = "cyan")) +  
  theme_void() +  
  ggtitle("Behavioral Categories in Species Interactions")


ggraph(graph, layout = "kk") +  
  geom_edge_link(aes(width = Denv, color = Rove), alpha = 0.7) +  
  geom_node_point(aes(color = behavior), size = 5) +  
  geom_node_text(aes(label = name), repel = TRUE) +  
  scale_edge_width(range = c(0.5, 4)) +  
  scale_edge_color_gradient(low = "yellow", high = "red") +  
  scale_color_manual(values = c("Generalist" = "orange", "Specialist" = "cyan")) +  
  theme_void() +  
  ggtitle("Environmental and Resource Overlap in Species Interactions")


regions <- unique(edges$region)

for (reg in regions) {
  sub_edges <- edges %>% filter(region == reg)
  sub_nodes <- unique(c(sub_edges$species1, sub_edges$species2))
  sub_graph <- graph_from_data_frame(sub_edges, vertices = data.frame(name = sub_nodes), directed = FALSE)
  
  print(
    ggraph(sub_graph, layout = "fr") +  
      geom_edge_link(aes(width = ENTotal, color = ROTotal), alpha = 0.7) +  
      geom_node_point(size = 5, color = "darkblue") +  
      geom_node_text(aes(label = name), repel = TRUE) +  
      scale_edge_width(range = c(0.5, 3)) +  
      scale_edge_color_gradient(low = "yellow", high = "red") +  
      theme_void() +  
      ggtitle(paste("Species Interaction Network in", reg))
  )
}

library(igraph)

library(igraph)
library(sf)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(geosphere)  # For distance calculations

# Function to find the most central coordinate within the species' distribution
find_center <- function(df) {
  if (nrow(df) == 1) return(df)  # If only one point, return it
  x <- df$x
  y <- df$y
  # Compute sum of distances to all other points
  idx <- which.min(sqrt((x - mean(x, na.rm = TRUE))^2 + (y - mean(y, na.rm = TRUE))^2))

  df[idx,]

}

species_coords <- na.exclude(as.data.frame(amph_BC.PA, xy = TRUE)) %>%
  pivot_longer(cols = -c(x,y), names_to = "species", values_to = "presence") %>%
  filter(presence == 1) %>%
  group_by(species) %>%
  group_modify(~ find_center(.x)) %>%
  ungroup()

# Convert to spatial points
species_sf <- st_as_sf(species_coords, coords = c("x", "y"), crs = 4326) %>%
  st_shift_longitude()

# Extract edges from graph
edges_df <- as_data_frame(graph, what = "edges")

# Merge with coordinates for start and end points
edges_mapped <- edges_df %>%
  left_join(species_coords, by = c("from" = "species")) %>%
  rename(lat_from = y, lon_from = x) %>%
  left_join(species_coords, by = c("to" = "species")) %>%
  rename(lat_to = y, lon_to = x)

# Plot world map with network
map.world <- map_data("world") %>%
  mutate(long = ifelse(long < 0, long + 360, long)) %>%
  filter(long > 20, long < 240, lat > -40, lat < 40)

# Plot
ggplot() +
  geom_polygon(data = map.world, aes(x = long, y = lat, group = group),
               fill = "gray90", color = "gray70") +  # Base map
  geom_sf(data = species_sf, color = "red", size = 3) +  # Plot species locations
  geom_curve(data = edges_mapped, aes(x = lon_from, y = lat_from,  
                                       xend = lon_to, yend = lat_to, 
                                       color = interaction),  # Color by interaction type
              curvature = 0.2, arrow = arrow(length = unit(0.02, "inches")),
              size = 0.7, alpha = 0.5) +  # Add interaction lines
  scale_color_manual(values = c("GG" = "red", "GS" = "green", "SS" = "blue")) +  
  geom_text_repel(data = species_coords, aes(x = x, y = y, label = species), size = 3) +  # Label species
  theme_minimal() +
  ggtitle("Geographical Network of Species Interactions")



region_levels <- unique(INT.pred$region[order(INT.pred$x)])

# Create an edge list for species interactions
INT.data_all <- INT.data %>%
  select(sp1, sp2, region, interaction, Denv, Rove) %>%
  group_by(sp1, sp2, region, interaction) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

INT.data_all$region <- factor(INT.data_all$region, levels = region_levels)

# Create a node list with unique species
nodes <- unique(c(INT.data_all$sp1, INT.data_all$sp2))
nodes_df <- data.frame(name = nodes)

# Assign behavior categories
nodes_df$behavior <- ifelse(nodes_df$name %in% behavior$species[behavior$mutualism == "generalist"], "Generalist", "Specialist")

# Convert edges to an igraph object
graph <- graph_from_data_frame(INT.data_all, vertices = nodes_df, directed = FALSE)

ggraph(graph, layout = "fr") +  
  geom_edge_link(aes(width = Denv, color = interaction), alpha = 0.6) +  
  geom_node_point(aes(color = behavior), size = 5) +  
  geom_node_text(aes(label = name), repel = TRUE) +  
  scale_edge_width(range = c(0.5, 3)) +  
  scale_edge_color_manual(values = c("GG" = "red", "GS" = "green", "SS" = "blue")) +  
  scale_color_manual(values = c("Generalist" = "orange", "Specialist" = "cyan")) +  
  theme_void() +  
  ggtitle("Behavioral Categories in Species Interactions")


ggraph(graph, layout = "kk") +  
  geom_edge_link(aes(width = Denv, color = Rove), alpha = 0.7) +  
  geom_node_point(aes(color = behavior), size = 5) +  
  geom_node_text(aes(label = name), repel = TRUE) +  
  scale_edge_width(range = c(0.5, 4)) +  
  scale_edge_color_gradient(low = "yellow", high = "red") +  
  scale_color_manual(values = c("Generalist" = "orange", "Specialist" = "cyan")) +  
  theme_void() +  
  ggtitle("Environmental and Resource Overlap in Species Interactions")

# Network Structure & Complexity
# Degree Centrality: Find the most connected species (hubs) by counting their interactions.
V(graph)$degree <- degree(graph)

# Betweenness Centrality:  Identify species that act as bridges between different groups.
V(graph)$betweenness <- betweenness(graph, normalized = TRUE)

# Clustering Coefficient: Measures how interconnected species are.
clustering_coeff <- transitivity(graph, type = "global")

# Network Density: Proportion of actual connections relative to possible connections.
network_density <- edge_density(graph)

print(clustering_coeff)
print(network_density)

# Community Detection
#  Identify if species segregate into different functional or ecological guilds.
communities <- cluster_louvain(graph)
plot(graph, vertex.color = communities$membership)

# Interaction Asymmetry & Specialization
# Interaction Strength: Measure how frequently species interact.
V(graph)$generalist_score <- strength(graph)  # Sum of edge weights per node

# Plot Degree Centrality
ggraph(graph, layout = "circle") +
  geom_edge_link(alpha = 0.5) +
  geom_node_text(aes(label = name), repel = TRUE) +  
  geom_node_point(aes(size = degree, color =  as.factor(V(graph)$behavior))) +
  scale_color_manual(values = c("Generalist" = "orange", "Specialist" = "cyan")) +  
  theme_void() + ggtitle("Degree Centrality (Species Connectivity)")

# Betweenness Centrality Histogram
ggplot(data.frame(Betweenness = V(graph)$betweenness), aes(x = Betweenness)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue") +
  theme_minimal() + ggtitle("Betweenness Centrality Distribution")

# Print Network-Wide Metrics
print(paste("Clustering Coefficient:", round(clustering_coeff, 3)))
print(paste("Network Density:", round(network_density, 3)))

key_species <- data.frame(
  Species = V(graph)$name,
  Degree = V(graph)$degree,
  Betweenness = V(graph)$betweenness
) %>% arrange(desc(Degree))

ggplot(key_species, aes(x = reorder(Species, -Degree), y = Degree)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Top Species by Degree Centrality")


# Compare Networks Across Species
# Network Similarity: Compare species interactions between different regions (Jaccard similarity).
# Compute Jaccard similarity matrix
similarity_matrix <- similarity(graph, method = "jaccard")

# Convert matrix to dataframe with row and column names as region names
rownames(similarity_matrix) <- V(graph)$name
colnames(similarity_matrix) <- V(graph)$name

library(ape)
library(vegan)

phylo_tree <- ape::read.tree("../../3_EcoPopGen/data/calibrated_tree.tre")

phylo_tree$tip.label <- gsub("A._", "Amphiprion_", phylo_tree$tip.label)
phylo_tree$tip.label[phylo_tree$tip.label == "Amphiprion_biaculeatus"] = "Premnas_biaculeatus"

# Extract species names from both datasets
species_in_matrix <- rownames(similarity_matrix)  # Species from similarity matrix
species_in_tree <- phylo_tree$tip.label            # Species from phylogenetic tree

# Find common species
common_species <- intersect(species_in_matrix, species_in_tree)

# Subset the similarity matrix
species_similarity <- similarity_matrix[common_species, common_species]

# Prune the phylogenetic tree to keep only common species
phylo_tree <- drop.tip(phylo_tree, setdiff(species_in_tree, common_species))

# Check if everything matches
all(rownames(species_similarity) %in% phylo_tree$tip.label)  # Should be TRUE
all(phylo_tree$tip.label %in% rownames(species_similarity))  # Should be TRUE

# Compute pairwise phylogenetic distances
phylo_dist <- cophenetic(phylo_tree) 

# Ensure species names in both matrices match
phylo_dist <- phylo_dist[common_species, common_species]

# Perform Mantel test to check correlation
mantel_test <- mantel(species_similarity, phylo_dist, method = "spearman", permutations = 999)
print(mantel_test)


library(ggplot2)
library(ggtree)
library(ape)
library(dplyr)
library(reshape2)
library(patchwork)

species_levels <- phylo_tree$tip.label

# Convert similarity matrix to long format for ggplot2
species_similarity_reordered <- species_similarity[species_levels, species_levels]

similarity_long <- melt(as.matrix(species_similarity_reordered))
colnames(similarity_long) <- c("Species1", "Species2", "Similarity")

similarity_long$Species1 = factor(similarity_long$Species1, levels = rev(species_levels))
similarity_long$Species2 = factor(similarity_long$Species2, levels = rev(species_levels))


# Ensure species order matches between tree and heatmap
similarity_long <- similarity_long %>%
  filter(Species1 %in% tree_species_order, Species2 %in% tree_species_order)

# Create the heatmap using ggplot2
heatmap_plot <- ggplot(similarity_long, aes(x = Species1, y = Species2, fill = Similarity)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Jaccard Similarity Between Species")

# Create the tree plot
tree_plot <- ggtree(phylo_tree) + 
  theme_tree2()  # Cleaner tree visualization

# Combine the tree and heatmap
tree_plot | heatmap_plot


# Convert to long format for visualization
similarity_df <- as.data.frame(as.table(species_similarity))
colnames(similarity_df) <- c("Species1", "Species2", "Similarity")


# Visualize with a heatmap
ggplot(similarity_df, aes(x = Species1, y = Species2, fill = Similarity)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Jaccard Similarity Between Regions (Ordered)")


hc <- hclust(as.dist(1 - species_similarity_reordered))
plot(hc, main = "Hierarchical Clustering of Species")


library(igraph)
library(ggplot2)
library(ggraph)

# Add environmental and resource overlap as edge attributes
E(graph)$Denv <- edge_attr(graph, "Denv")  # Environmental Overlap
E(graph)$Rove <- edge_attr(graph, "Rove")  # Resource Overlap

# Normalize values for visualization
E(graph)$Denv_scaled <- scales::rescale(E(graph)$Denv, to = c(1, 5)) # Edge width
E(graph)$Rove_scaled <- scales::rescale(E(graph)$Rove, to = c(1, 5))

# Define color scale based on Denv or Rove
edge_colors <- ifelse(E(graph)$Rove > 0.5, "blue", "red") # Blue = high env overlap, Red = low

# Plot
ggraph(graph, layout = "fr") + 
  geom_edge_link(aes(width = Rove, color = edge_colors), alpha = 0.8) + 
  geom_node_point(aes(size = degree(graph), color = as.factor(V(graph)$behavior)), alpha = 0.9) +
  scale_color_manual(values = c("generalist" = "green", "specialist" = "purple")) +
  theme_void() + ggtitle("Species Interaction Network with Environmental Overlap")
