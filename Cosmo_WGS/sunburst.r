rm(list=ls())
library(plotly)
library(treeio)
library(randomcoloR)
library(gridExtra)
if (!require("processx")) install.packages("processx")
library(processx)
library(TreeTools)
library(ggtree)
library(RColorBrewer)
library(rgdal)
library(ggplot2)
library(rworldmap)
library(cleangeo)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(mapproj)
library("rnaturalearth")
library(rnaturalearthdata)
library(scatterpie)

Sys.setenv("PATH" = "/Users/kathryncampbell/miniconda3/bin")

Sys.getenv("PATH")

lineage_info<-read.csv("Cosmo_WGS/Outputs/Cosmo_WGS_lineage_info.csv")
node_data<-read.csv("Cosmo_WGS/Outputs/Cosmo_WGS_node_data.csv")
tree<-read.tree("Cosmo_WGS/Trees/Cosmo_WGS_aligned.fasta.contree")
metadata<-read.csv("Cosmo_WGS/Cosmo_WGS_metadata.csv")
sequence_data<-read.csv("Cosmo_WGS/Outputs/Cosmo_WGS_sequence_data.csv")

previous<-data.frame(assignment = unique(metadata$alignment.name), parent = "Cosmopolitan", n_seqs = NA)
previous$parent[1]<-""
for (i in 1:length(previous$assignment)) {
  previous$n_seqs[i]<-length(which(metadata$alignment.name == previous$assignment[i]))
}
previous$assignment<-gsub("AL_", "", previous$assignment)
previous$assignment<-gsub("Cosmopolitan_", "", previous$assignment)

node_data$parent<-NA
node_data$parent[1]<-""

for (i in 2:length(node_data$Node)) {
  parent<-node_data$cluster[which(node_data$Node %in% ancestor(tree, node_data$Node[i]))]
  node_data$parent[i]<-parent[length(parent)]
}

lineage_info$parent<-NA

for (i in 1:length(lineage_info$cluster)) {
  lineage_info$parent[i]<-node_data$parent[which(node_data$cluster == lineage_info$cluster[i])]
  
}

lineage_info$colour<-NA

Colours<-c("Reds","Purples","YlOrRd","Greens","Greys",
           "BuGn","RdPu","Oranges","BuPu","YlGn","PuRd","YlGnBu")

clades<-c("AF1a","AF1b","AM2a","CA1","CA2","CE","Cosmopolitan EE",
          "ME1a","ME2","NEE","WE","Cosmopolitan_")

lineage<-lineage_info$cluster[-c(grep("_", lineage_info$cluster))]
cols<-brewer.pal(9, "Blues")
pal<-colorRampPalette(c(cols))
pal<-rev(pal(length(lineage)))
lineage_info$colour[-c(grep("_", lineage_info$cluster))]<-pal

for (i in 1:length(clades)) {
  lineage<-grep(clades[i], lineage_info$cluster)
  cols<-brewer.pal(3, Colours[i])
  pal<-colorRampPalette(c(cols))
  pal<-rev(pal(length(lineage)))
  lineage_info$colour[(grep(clades[i], lineage_info$cluster))]<-pal
}

lineage_info$cluster<-gsub("Cosmopolitan ", "", lineage_info$cluster)
lineage_info$parent<-gsub("Cosmopolitan ", "", lineage_info$parent)

dput(lineage_info$colour)

new<-plot_ly(
  labels = c(lineage_info$cluster),
  parents = c(lineage_info$parent),
  values = c(lineage_info$n_seqs),
  type = "sunburst",
  width = 1000,
  height = 850,
  marker = list(colors = list("#08306B", "#08468B", "#105BA4", "#2171B5", "#3787C0", "#4F9BCB", 
                              "#6BAED6", "#8DC0DD", "#DE2D26", "#FEE0D2", "#756BB1", "#7C73B5", 
                              "#837CBA", "#8B84BE", "#928DC3", "#9A96C7", "#A19ECC", "#A9A7D0", 
                              "#B0B0D5", "#B8B8D9", "#BEBFDD", "#C4C4DF", "#C9C9E2", "#CECEE5", 
                              "#D4D3E7", "#D9D8EA", "#DEDDED", "#E4E2EF", "#E9E7F2", "#EFEDF5", 
                              "#FFEDA0", "#31A354", "#A1D99B", "#E5F5E0", "#636363", "#BDBDBD", 
                              "#F0F0F0", "#E5F5F9", "#C51B8A", "#FDE0DD", "#E6550D", "#E9631C", 
                              "#ED722C", "#F1813C", "#F5904B", "#F99F5B", "#FDAE6B", "#FDB77B", 
                              "#FDC08C", "#FDCA9C", "#FDD3AD", "#FDDCBD", "#FEE6CE", "#8856A7", 
                              "#907EBB", "#99A7CF", "#ABC5DF", "#C5D8E9", "#E0ECF4", "#31A354", 
                              "#5AB667", "#83C97A", "#ADDD8E", "#C5E79C", "#DEF1AA", "#F7FCB9", 
                              "#E7E1EF", "#EDF8B1", "#ABCFE5", "#C6DBEF", "#D6E5F4", "#E6F0F9", 
                              "#F7FBFF"))
)

new

#previous_sun plot now done in N_all

# orca(new, "Cosmo_WGS/Figures/sunburst.png", width = 800, height = 700)
# Better resolution to just save from plot window

# Plot a nice figure to save
plot_tree<-ggtree(tree, colour = "grey50", ladderize = T) %<+% sequence_data +
  geom_tippoint(aes(color=cluster), size=3)  +
  ggtitle("Cosmo WGS Lineage Tree")+
  theme(plot.title = element_text(size = 40, face = "bold"))+ 
  scale_color_manual(values=c(lineage_info$colour))

plot_tree

lineage_info$node<-NA

node_data$cluster<-gsub("Cosmopolitan ", "", node_data$cluster)

for (i in 1:length(lineage_info$cluster)) {
  lineage_info$node[i]<-node_data$Node[which(node_data$cluster == lineage_info$cluster[i])]
  lineage_info$group[i]<-node_data$overlaps[which(node_data$cluster == lineage_info$cluster[i])]
}

group0<-which(lineage_info$group == 0)

for (i in 1:length(group0)) {
  plot_tree<-
    collapse(plot_tree, lineage_info$node[group0[i]], 'max', fill=lineage_info$colour[group0[i]], alpha=1)
}

plot_tree<-plot_tree + theme(legend.position = "none")

ggsave("Cosmo_WGS/Figures/figure_lineage_tree.png",
       plot = plot_tree,
       height = 25, width = 40)

world<-ne_countries()

#' **Cleaning the data**
#' Find which country names do not match between data and map file
map_countries <- as.character(world$admin)
countries <- unique(metadata$country)
no_match <- setdiff(countries, map_countries); message(length(no_match), " countries are mis-matched: \n", paste0(no_match, collapse="\n"))

metadata$country[grepl("USA|United States", metadata$country)] <- "United States of America"
metadata$country[grepl("Tanzania", metadata$country)] <- "United Republic of Tanzania"
metadata$country[grepl("Serbia", metadata$country)] <- "Republic of Serbia"
metadata$country[grepl("Czechia", metadata$country)] <- "Czech Republic"
metadata$country[grepl("Grenada", metadata$country)] <- "Grenada"

#' How many entries are there for the remaining mis-matches?
table(metadata$country[which(metadata$country %in% no_match)])

countries <- unique(metadata$country)
#' ** Tabulating the data**
country_table<-table(metadata$country); country_table
length(country_table)

sequence_data<-sequence_data[-c(which(is.na(sequence_data$cluster))),]

lineage_table <- data.frame(matrix(ncol = (length(unique(sequence_data$cluster))+1), nrow = length(countries)))
x <- c("country", unique(sequence_data$cluster))
colnames(lineage_table) <- x

lineage_table$country<-countries

for (i in 1:length(countries)) {
  country_place<-which(lineage_table$country == countries[i])
  lineages<-unique(sequence_data$cluster[which(sequence_data$ID %in% metadata$ID[which(metadata$country == countries[i])])])
  
  for (j in 1:length(lineages)) {
    country_lineage<-which(colnames(lineage_table) == lineages[j])
    lineage_table[country_place, country_lineage]<-length(which((sequence_data$cluster[which(sequence_data$ID %in% 
                                                                                          metadata$ID[which(metadata$country == countries[i])])]) == lineages[j]))
  }
}

lineage_table[is.na(lineage_table)] <- 0

lineage_table$LAT<-NA
lineage_table$LON<-NA

lineage_table<-lineage_table[-c(which(lineage_table$country == "-")),]
lineage_table<-lineage_table[-c(which(lineage_table$country == "Grenada")),]

for (i in 1:length(lineage_table$country)) {
  lineage_table$LAT[i]<-world@polygons[which(world$admin == lineage_table$country[i])][[1]]@labpt[2]
  lineage_table$LON[i]<-world@polygons[which(world$admin == lineage_table$country[i])][[1]]@labpt[1]
}

world <- ne_countries(scale = "medium", returnclass = "sf")
plot<-
  ggplot(data = world) +
  geom_sf() 

lineage_table <- with(lineage_table, lineage_table[abs(LON) < 150 & abs(LAT) < 70,])
n <- nrow(lineage_table)
lineage_table$region <- factor(1:n)
lineage_table$radius<-NA
for (i in 1:length(lineage_table$country)) {
  lineage_table$radius[i]<-sum(lineage_table[i,2:74])
}

lineage_table$radius[which(lineage_table$radius %in% 6:10)]<-6
lineage_table$radius[which(lineage_table$radius %in% 11:20)]<-7
lineage_table$radius[which(lineage_table$radius %in% 21:50)]<-8
lineage_table$radius[which(lineage_table$radius %in% 51:100)]<-9
lineage_table$radius[which(lineage_table$radius %in% 101:200)]<-10
lineage_table$radius[which(lineage_table$radius > 200)]<-11

head(lineage_table)

lineage_table$radius<-lineage_table$radius/1.5

colour_table<-data.frame(lineage = colnames(lineage_table[2:74]), colour = NA)
colour_table$lineage<-gsub("Cosmopolitan ", "", colour_table$lineage)

for (i in 1:length(colour_table$lineage)) {
  colour_table$colour[i]<-lineage_info$colour[which(lineage_info$cluster == colour_table$lineage[i])]
}

plot_world<-plot + geom_scatterpie(aes(x=LON, y=LAT, group=region, r=radius),
                       data=lineage_table, cols=c(colnames(lineage_table)[2:74]), color=NA, alpha=.8)+ 
  theme(legend.position = "none") + scale_fill_manual(values = c(colour_table$colour))+
  coord_sf(xlim = c(-155, 155), ylim = c(-50, 70))

plot_zoom<-plot + geom_scatterpie(aes(x=LON, y=LAT, group=region, r=(radius/2)),
                             data=lineage_table, cols=c(colnames(lineage_table)[2:74]), color=NA, alpha=.8)+ 
  theme(legend.position = "none") + scale_fill_manual(values = c(colour_table$colour))+
  coord_sf(xlim = c(-15,55), ylim = c(-30,62))

plot_world
plot_zoom
ggsave("Cosmo_WGS/Figures/pie_map.png", plot = plot_world, width = 45, height = 20)
ggsave("Cosmo_WGS/Figures/pie_map_zoom.png", plot = plot_zoom)
