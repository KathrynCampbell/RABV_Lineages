rm(list=ls())
library(plotly)
library(treeio)
library(randomcoloR)
library(gridExtra)
if (!require("processx")) install.packages("processx")
library(processx)
library(TreeTools)
library(ggtree)
library(scales)
library(RColorBrewer)
library(rgdal)
library(ggplot2)
library(scatterpie)
library(phytools)
library(phangorn)

Sys.setenv("PATH" = "/Users/kathryncampbell/miniconda3/bin")

Sys.getenv("PATH")

lineage_info<-read.csv("Tanz_WGS/Outputs/Tanz_WGS_lineage_info.csv")
node_data<-read.csv("Tanz_WGS/Outputs/Tanz_WGS_node_data.csv")
outgroup<-read.tree("Tanz_WGS/Trees/Tanz_outgroup_aligned.fasta.contree")
tree<-read.tree("Tanz_WGS/Trees/Tanz_WGS_aligned.fasta.contree")
metadata<-read.csv("Tanz_WGS/Tanz_WGS_metadata_outgroup.csv")
sequence_data<-read.csv("Tanz_WGS/Outputs/Tanz_WGS_sequence_data.csv")

lineage_info$parent<-NA

for (i in 1:length(lineage_info$cluster)) {
  lineage_info$parent[i]<-all_lineage$parent[which(all_lineage$lineage == lineage_info$cluster[i])]
}


names(node_data)<-c("Node", "n_tips", "cluster")

shape_region<-readOGR("Tanz_WGS/Shapefiles/gadm40_TZA_1.shp")
shape_district<-readOGR("Tanz_WGS/Shapefiles/gadm40_TZA_2.shp")

outgroup<-reroot(outgroup, getMRCA(
  outgroup, tip = which(outgroup$tip.label %in% metadata$sequence.sequenceID[(
    grep("Asian", metadata$assignment))])))

outgroup<-tree_subset(outgroup, "KR906752", levels_back = 2)

metadata<-metadata[-c(grep("Asian", metadata$assignment)),]

lineage_info$colour<-NA

Colours<-c("Purples","Greens","Reds","Greys","Oranges", "BuGn")

clades<-c("AF1b_A","AF1b_B","AF1b_C")

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

new<-plot_ly(
  labels = c(lineage_info$cluster),
  parents = c(lineage_info$parent),
  values = c(lineage_info$n_seqs),
  type = "sunburst",
  width = 1000,
  height = 850,
  marker = list(colors = lineage_info$colour)
)

htmlwidgets::saveWidget(plotly::as_widget(new), "Tanz_WGS/Figures/sunburst.html")

# Plot a nice figure to save
plot_tree<-ggtree(outgroup, colour = "grey50", ladderize = T) %<+% sequence_data +
  geom_tippoint(color="grey50", size=4)+
  geom_tippoint(aes(color=cluster), size=3)  +
  theme(plot.title = element_text(size = 40, face = "bold"))+ 
  scale_color_manual(values=c(lineage_info$colour)) + 
  theme(legend.position = "none")+
  ggplot2::ylim(0, 200) +
  geom_treescale()

genotype<-data.frame(lineage = sequence_data$cluster)
rownames(genotype)<-sequence_data$ID

plot_tree<-gheatmap(plot_tree, genotype, offset=-0, width=.1, font.size=3, color = NA, 
                    colnames_angle=-45, hjust=0) +
  scale_fill_manual(values=c(lineage_info$colour), name = "")+ 
  theme(legend.position = "none")

plot_tree
ggsave("Tanz_WGS/Figures/figure_lineage_tree.png",
       plot = plot_tree)


#' **Cleaning the data**
#' Find which country names do not match between data and map file
map_regions <- as.character(shape_region@data$NAME_1)
map_districts <- as.character(shape_district@data$NAME_2)
areas <- unique(metadata$sequence.gb_place_sampled)
no_match <- setdiff(areas, map_regions); message(length(no_match), " countries are mis-matched: \n", paste0(no_match, collapse="\n"))

metadata$sequence.gb_place_sampled<-gsub("region|Region|District|district|island", "", metadata$sequence.gb_place_sampled)
metadata$sequence.gb_place_sampled<-gsub(" ", "", metadata$sequence.gb_place_sampled)
metadata$sequence.gb_place_sampled[grepl("DaresSalaam", metadata$sequence.gb_place_sampled)] <- "Dar es Salaam"
metadata$sequence.gb_place_sampled[grepl("Chakechake", metadata$sequence.gb_place_sampled)] <- "Pemba"
metadata$sequence.gb_place_sampled[grepl("KusiniUnguja", metadata$sequence.gb_place_sampled)] <- "Kusini"

#' How many entries are there for the remaining mis-matches?
areas <- unique(metadata$sequence.gb_place_sampled)
no_match <- setdiff(areas, map_regions); message(length(no_match), " countries are mis-matched: \n", paste0(no_match, collapse="\n"))

no_match <- setdiff(no_match, map_districts); message(length(no_match), " countries are mis-matched: \n", paste0(no_match, collapse="\n"))
#Serengeti at district level
#Kusini Unguja at district level, just "Kusini"

shape_region@data[grep("Pemba", shape_region@data$NAME_1),]
#Pemba 2 regions; maybe average to get middle?

areas <- unique(metadata$sequence.gb_place_sampled)
#' ** Tabulating the data**
areas_table<-table(metadata$sequence.gb_place_sampled); areas_table
length(areas_table)

lineage_table <- data.frame(matrix(ncol = (length(unique(sequence_data$cluster))+1), nrow = length(areas)))
x <- c("area", unique(sequence_data$cluster))
colnames(lineage_table) <- x

lineage_table$area<-areas

for (i in 1:length(areas)) {
  area_place<-which(lineage_table$area == areas[i])
  lineages<-unique(sequence_data$cluster[which(sequence_data$ID %in% metadata$sequence.sequenceID[which(metadata$sequence.gb_place_sampled == areas[i])])])
  for (j in 1:length(lineages)) {
    area_lineage<-which(colnames(lineage_table) == lineages[j])
    lineage_table[area_place, area_lineage]<-length(which((sequence_data$cluster[which(sequence_data$ID %in% 
                                                                                         metadata$sequence.sequenceID[which(metadata$sequence.gb_place_sampled == areas[i])])]) == lineages[j]))
  }
}

lineage_table[is.na(lineage_table)] <- 0

lineage_table$LAT<-NA
lineage_table$LON<-NA

lineage_table<-lineage_table[-c(which(lineage_table$area == "-")),]

lineage_table_problems<-lineage_table[grep("Serengeti|Pemba|Kusini", lineage_table$area),]
lineage_table<-lineage_table[-c(grep("Serengeti|Pemba|Kusini", lineage_table$area)),]

for (i in 1:length(lineage_table$area)) {
  lineage_table$LAT[i]<-shape_region@polygons[which(shape_region@data$NAME_1 == lineage_table$area[i])][[1]]@labpt[2]
  lineage_table$LON[i]<-shape_region@polygons[which(shape_region@data$NAME_1 == lineage_table$area[i])][[1]]@labpt[1]
}

for (i in 1) {
  lineage_table_problems$LAT[i]<-shape_district@polygons[which(shape_district@data$NAME_2 == lineage_table_problems$area[i])][[1]]@labpt[2]
  lineage_table_problems$LON[i]<-shape_district@polygons[which(shape_district@data$NAME_2 == lineage_table_problems$area[i])][[1]]@labpt[1]
}

lineage_table_problems$LAT[3]<-shape_district@polygons[which(shape_district@data$NAME_2 == lineage_table_problems$area[3])][[1]]@labpt[2]
lineage_table_problems$LON[3]<-shape_district@polygons[which(shape_district@data$NAME_2 == lineage_table_problems$area[3])][[1]]@labpt[1]

lineage_table_problems$LAT[2]<-mean(c(
  shape_region@polygons[[grep(lineage_table_problems$area[2], shape_region@data$NAME_1)[1]]]@labpt[2], 
  shape_region@polygons[[grep(lineage_table_problems$area[2], shape_region@data$NAME_1)[2]]]@labpt[2]
))

lineage_table_problems$LON[2]<-mean(c(
  shape_region@polygons[[grep(lineage_table_problems$area[2], shape_region@data$NAME_1)[1]]]@labpt[1], 
  shape_region@polygons[[grep(lineage_table_problems$area[2], shape_region@data$NAME_1)[2]]]@labpt[1]
))

lineage_table<-rbind(lineage_table, lineage_table_problems)

plot<-
  ggplot() +
  geom_polygon(data = shape_region, aes( x = long, y = lat, group = group), fill="grey80", color="white") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

lineage_table <- with(lineage_table, lineage_table[abs(LON) < 150 & abs(LAT) < 70,])
n <- nrow(lineage_table)
lineage_table$region <- factor(1:n)
lineage_table$radius<-NA
for (i in 1:length(lineage_table$area)) {
  lineage_table$radius[i]<-sum(lineage_table[i,2:15])
}

lineage_table$radius
lineage_table$radius[which(lineage_table$radius %in% 1:5)]<-2
lineage_table$radius[which(lineage_table$radius %in% 6:10)]<-3
lineage_table$radius[which(lineage_table$radius > 10)]<-4

lineage_table$radius<-lineage_table$radius/6

colour_table<-data.frame(lineage = colnames(lineage_table[2:15]), colour = NA)
colour_table$lineage<-gsub("Cosmopolitan ", "", colour_table$lineage)

for (i in 1:length(colour_table$lineage)) {
  colour_table$colour[i]<-lineage_info$colour[which(lineage_info$cluster == colour_table$lineage[i])]
}

plot_world<-plot + geom_scatterpie(aes(x=LON, y=LAT, group=region, r=radius),
                                   data=lineage_table, cols=c(colnames(lineage_table)[2:15]), color=NA)+ 
  theme(legend.position = "none") + scale_fill_manual(values = c(colour_table$colour))+ coord_equal(ratio=1)+
  theme(panel.grid.major = element_blank())+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())

plot_world

ggsave("Tanz_WGS/Figures/pie_map.png", plot = plot_world)

