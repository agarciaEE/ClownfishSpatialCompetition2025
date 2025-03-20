args = commandArgs(TRUE)
width <- if(is.na(args[1])){472}else{args[1]}
height <- if(is.na(args[2])){231}else{args[2]}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load libraries
library(ggplot2)
library(rgdal)
library(colorspace)
library(ggthemes)
library(dplyr)
library(sp)
library(maptools)

# load marine regions dataset
mr  = read.csv("data/marine_regions.csv") # environmental locations with designated marine regions

provs  = read.csv("data/meow_ecos_df.csv") # marine regions polygon data set
provs$REALM = factor(provs$REALM, levels = c("Western Indo-Pacific" , "Central Indo-Pacific", "Eastern Indo-Pacific", "Temperate Northern Pacific", "Temperate Australasia" )) # order realms
provs$PROVINCE = factor(provs$PROVINCE, levels = unique(provs[order(provs$REALM),"PROVINCE"])) # order provs in base to realms

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#set realm and provs colors
col.realm = c("#21908CFF", "#FDE725FF", "#5DC863FF","#2C728EFF","#BB3754FF")
names(col.realm) = c("Western Indo-Pacific" , "Central Indo-Pacific", "Eastern Indo-Pacific", "Temperate Northern Pacific", "Temperate Australasia" )

col.prov = rep(col.realm)[as.factor(unlist(lapply(levels(provs$PROVINCE), function(i) provs[provs$PROVINCE == as.character(i), "REALM"][1])))]
names(col.prov) = levels(provs$PROVINCE)

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# get provs indexes and key legend
ctr = data.frame(matrix(NA, length(levels(provs$id)),2), row.names = levels(provs$id))
for(i in levels(provs$PROVINCE)){
  c <- cov.wt(provs[provs$PROVINCE == i,1:2])$center
  r = gsub("[^[:alnum:]]", "_", as.character(i))
  ctr[i,] <- mr[mr$province == r,1:2][which.min(sqrt((mr[mr$province == r,1] - c[1])^2 + (mr[mr$province == r,2] - c[2])^2)),]
}
ctr$n <- 1:length(levels(provs$PROVINCE))
ctr$id <- rownames(ctr)
ctr <- dplyr::mutate(ctr, field_key = paste0(c(rep("   ", 9), rep("  ", 18)),
                                             paste(n, id, sep = "     ")),
                     field_key = factor(field_key, levels = field_key))
levels(provs$id) = levels(ctr$field_key)

ctr[5,c("X1", "X2")] <- c(60,20) # relocate number of Somali/Arabian Peninsula to avoid confusion with Red Sea and Gulf of Aden

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# get wprld map
wmp <- map_data("world")
wmp[which(wmp$lon < 0),"long"] = wmp[which(wmp$lon < 0),"long"] + 360
wmp <- wmp[wmp$long > 20 & wmp$long < 240,]
wmp <- wmp[wmp$lat < 40 & wmp$lat > -40,]

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# plot
mp <- ggplot(provs,aes(long, lat, group=id,  col = PROVINCE)) + theme_bw() +
  geom_polygon(aes(fill = REALM)) +
  annotate(x=180, xend=180, y=-13.75, yend=-4.4, lwd = 2, colour= col.realm[3], geom="segment", lineend = "round") +
  annotate(x=180, xend=180, y=26.3, yend=30, lwd = 1, colour= col.realm[3], geom="segment", lineend = "round") +
  annotate(x=179.9, xend=179.9, y=-4, yend=-2, lwd = 0.8, colour= col.realm[3], geom="segment", lineend = "round") +
  annotate(x=180, xend=180, y=-23.2, yend=-14.5, lwd = 1.5, colour= col.realm[2], geom="segment", lineend = "round") +
  coord_equal() +
  geom_polygon(data = wmp,
               aes(x = long, y = lat, group = group), fill = "grey85", col = "grey75") +
  geom_text(data = ctr, aes(x = X1, y = X2, label = n, group = NA),size = 6, col = "black", fontface = 'bold') +
  scale_y_continuous("latitude", breaks = seq(-30, 30, 10), expand = c(0,0)) +
  scale_x_continuous("longitude", labels = c(0,50,100,150,-160,-110), expand = c(0,0)) +
  scale_fill_manual(name = "REALMS", values = col.realm, labels = paste0("            ", levels(provs$REALM))) +
  scale_color_manual("PROVINCE", values = lighten(col.prov,0.25), labels = paste0("  ", ctr[names(col.prov), "field_key"])) +
  guides(color = guide_legend(override.aes = list(name = "regions", fill = col.prov, col = NA, size = 3),
                              ncol=4, keywidth = 1.1, keyheight = 0.4, label.vjust = 2.25)) +
  theme_map() +
  theme(plot.title = element_text(size = 20, hjust = 0.5, vjust = -1),
        plot.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.text = element_text(margin = margin(5, 0, 0, -21)),
        axis.title.x = element_text(color = "black", size = 16, vjust = -1 ),
        axis.title.y = element_text(color = "black", size = 16, vjust = 2),
        axis.text.x = element_text(color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 16))

ggsave("figures/FigS1.Selected_MEOW_provs_realms.pdf", plot = mp, width=15, height=12, units="in")
