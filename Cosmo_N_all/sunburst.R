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
library(rworldmap)
library(cleangeo)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(mapproj)
library("rnaturalearth")
library(rnaturalearthdata)
library(scatterpie)
library(phytools)
library(phangorn)

Sys.setenv("PATH" = "/Users/kathryncampbell/miniconda3/bin")

Sys.getenv("PATH")

lineage_info<-read.csv("Cosmo_N_all/Outputs/Cosmo_N_all_lineage_info.csv")
node_data<-read.csv("Cosmo_N_all/Outputs/Cosmo_N_all_node_data.csv")
tree<-read.tree("Cosmo_N_all/Trees/Cosmo_N_all_aligned.fasta.contree")
outgroup<-read.tree("Cosmo_N_all/Trees/Cosmo_N_all_outgroup_aligned.fasta.contree")
outgroup_meta<-read.csv("Cosmo_WGS/Cosmo_WGS_metadata_outgroup.csv")
metadata<-read.csv("Cosmo_N_all/Cosmo_N_all_metadata.csv")
sequence_data<-read.csv("Cosmo_N_all/Outputs/Cosmo_N_all_sequence_data.csv")

metadata<-metadata[which(metadata$ID %in% sequence_data$ID),]

metadata<-rbind(metadata, outgroup_meta[grep("Asian", outgroup_meta$assignment),])

outgroup<-reroot(outgroup, getMRCA(
  outgroup, tip = which(outgroup$tip.label %in% metadata$ID[(
    grep("Asian", metadata$assignment))])))

metadata<-metadata[-c(grep("Asian", metadata$assignment)),]

outgroup<-tree_subset(outgroup, 2138, levels_back = 0)

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

Colours<-c("Reds","Purples","YlOrBr","PuBuGn","YlOrRd","OrRd","PuBu","Pastel1","Greens","Greys",
           "GnBu","BuGn","RdPu","Oranges","BuPu","YlGn","PuRd","YlGnBu")

clades<-c("AF1a","AF1b","AF4","AM1","AM2a","AM3a","AM3b","AM4","CA1","CA2","CA3","CE","Cosmopolitan EE",
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
  width = 700,
  height = 700,
  marker = list(colors = list("#08306B", "#083675", "#083D7F", "#084489", "#084B94", "#09529D", 
                              "#0E59A2", "#1360A7", "#1966AD", "#1E6DB2", "#2474B6", "#2B7BBA", 
                              "#3282BD", "#3989C1", "#4090C5", "#4896C8", "#519CCB", "#59A2CF", 
                              "#62A8D2", "#6BAED6", "#75B3D8", "#80B9DA", "#8BBFDC", "#95C5DF", 
                              "#A0CAE1", "#A8CEE4", "#B0D2E7", "#B9D5EA", "#DE2D26", "#E23C31", 
                              "#E74C3D", "#EB5B49", "#F06B54", "#F57A60", "#F98A6C", "#FC9779", 
                              "#FCA388", "#FCB096", "#FDBBA5", "#FDC8B4", "#FDD3C3", "#FEE0D2", 
                              "#756BB1", "#7F77B7", "#8A84BE", "#9590C4", "#A09DCB", "#ABAAD2", 
                              "#B6B6D8", "#BFC0DD", "#C7C8E1", "#CFCFE5", "#D7D6E9", "#DFDEED", 
                              "#E7E5F1", "#EFEDF5", "#FFF7BC", "#1C9099", "#A6BDDB", "#ECE2F0", 
                              "#F03B20", "#F98A3D", "#FEC568", "#FFEDA0", "#E34A33", "#E86043", 
                              "#ED7753", "#F28D63", "#F7A473", "#FDBB84", "#FDC391", "#FDCC9F", 
                              "#FDD6AC", "#FDDFBA", "#FEE8C8", "#2B8CBE", "#A6BDDB", "#ECE7F2", 
                              "#FBB4AE", "#31A354", "#42AB5E", "#53B369", "#64BB74", "#75C47F", 
                              "#87CC8A", "#98D495", "#A6DBA0", "#B0DFAA", "#BBE3B5", "#C5E8C0", 
                              "#D0ECCA", "#DAF0D5", "#E5F5E0", "#636363", "#9F9F9F", "#CECECE", 
                              "#F0F0F0", "#43A2CA", "#75BFBF", "#A8DDB5", "#C3E8C8", "#E0F3DB", 
                              "#E5F5F9", "#C51B8A", "#DF5D9F", "#FA9FB5", "#FBBFC9", "#FDE0DD", 
                              "#E6550D", "#E9611A", "#EC6E27", "#EF7B35", "#F38742", "#F69450", 
                              "#F9A15D", "#FDAE6B", "#FDB679", "#FDBE87", "#FDC695", "#FDCEA3", 
                              "#FDD6B1", "#FDDEBF", "#FEE6CE", "#8856A7", "#9EBCDA", "#E0ECF4", 
                              "#31A354", "#6FC071", "#ADDD8E", "#D2ECA3", "#F7FCB9", "#DD1C77", 
                              "#E7E1EF", "#2C7FB8", "#63B3BA", "#A3DBB7", "#EDF8B1", "#C1D9ED", 
                              "#C8DCEF", "#CDE0F1", "#D2E3F3", "#D7E6F4", "#DCEAF6", "#E1EDF8", 
                              "#E7F0F9", "#ECF4FB", "#F1F7FD", "#F7FBFF"))
)

new


# orca(new, "Cosmo_N_all/Figures/sunburst.png", width = 800, height = 700)
# Better resolution to just save from plot window

node_data$colour<-NA

node_data$cluster<-gsub("Cosmopolitan ", "", node_data$cluster)

for (i in 1:length(node_data$cluster)){
  node_data$colour[i]<-lineage_info$colour[which(lineage_info$cluster == node_data$cluster[i])]
}

# Plot a nice figure to save
plot_tree<-ggtree(outgroup, colour = "grey50", ladderize = T) %<+% sequence_data +
  geom_tippoint(aes(color=cluster), size=3)  +
  theme(plot.title = element_text(size = 40, face = "bold"))+ 
  scale_color_manual(values=c(lineage_info$colour)) + 
  theme(legend.position = "none")

genotype<-data.frame(lineage = sequence_data$cluster)
rownames(genotype)<-sequence_data$ID

plot_tree<-gheatmap(plot_tree, genotype, offset=-0.01, width=.1, font.size=3, color = NA, 
         colnames_angle=-45, hjust=0) +
  scale_fill_manual(values=c(lineage_info$colour), name="lineage")+ 
  theme(legend.position = "none")

plot_tree

ggsave("Cosmo_N_all/Figures/figure_lineage_tree.png",
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

lineage_table <- data.frame(matrix(ncol = (length(unique(sequence_data$cluster))+1), nrow = length(countries)))
y<-unique(sequence_data$cluster)
x <- c("country", y[order(y)])
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

colnames(lineage_table)<-gsub("Cosmopolitan ", "", colnames(lineage_table))

lineage_table<-lineage_table[,-c(
  which(colnames(lineage_table)%in%lineage_info$cluster[which(lineage_info$country == "Grenada")]))]
lineage_table<-lineage_table[,-c(
  which(colnames(lineage_table)%in%lineage_info$cluster[which(lineage_info$country == "character(0)")]))]

lineage_table<-lineage_table[,-c(which(is.na(colnames(lineage_table))))]

for (i in 1:length(lineage_table$country)) {
  lineage_table$LAT[i]<-world@polygons[which(world$admin == lineage_table$country[i])][[1]]@labpt[2]
  lineage_table$LON[i]<-world@polygons[which(world$admin == lineage_table$country[i])][[1]]@labpt[1]
}


world <- ne_countries(scale = "medium", returnclass = "sf")
plot<-
  ggplot(data = world) +
  geom_sf() 

n <- nrow(lineage_table)
lineage_table$region <- factor(1:n)
lineage_table$radius<-NA
for (i in 1:length(lineage_table$country)) {
  lineage_table$radius[i]<-sum(lineage_table[i,2:138])
}


lineage_table$radius[which(lineage_table$radius %in% 6:10)]<-6
lineage_table$radius[which(lineage_table$radius %in% 11:20)]<-7
lineage_table$radius[which(lineage_table$radius %in% 21:50)]<-8
lineage_table$radius[which(lineage_table$radius %in% 51:100)]<-9
lineage_table$radius[which(lineage_table$radius %in% 101:200)]<-10
lineage_table$radius[which(lineage_table$radius > 200)]<-11

head(lineage_table)

lineage_table$radius<-lineage_table$radius/1.5

colour_table<-data.frame(lineage = colnames(lineage_table[2:138]), colour = NA)
colour_table$lineage<-gsub("Cosmopolitan ", "", colour_table$lineage)

for (i in 1:length(colour_table$lineage)) {
  colour_table$colour[i]<-lineage_info$colour[which(lineage_info$cluster == colour_table$lineage[i])]
}

plot_world<-plot + geom_scatterpie(aes(x=LON, y=LAT, group=region, r=radius),
                                   data=lineage_table, cols=c(colnames(lineage_table)[2:138]), color=NA, alpha=.8)+ 
  theme(legend.position = "none") +scale_fill_manual(values = c(colour_table$colour))+
  coord_sf(xlim = c(-155, 155), ylim = c(-50, 70)) +
  theme(panel.grid.major = element_blank())+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())

plot_zoom<-plot + geom_scatterpie(aes(x=LON, y=LAT, group=region, r=(radius/2)),
                                  data=lineage_table, cols=c(colnames(lineage_table)[2:138]), color=NA, alpha=.8)+ 
  theme(legend.position = "none") + scale_fill_manual(values = c(colour_table$colour))+
  coord_sf(xlim = c(-15,55), ylim = c(-50,62))+
  theme(panel.grid.major = element_blank())+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())

plot_world
plot_zoom
ggsave("Cosmo_N_all/Figures/pie_map.png", plot = plot_world, width = 45, height = 20)
ggsave("Cosmo_N_all/Figures/pie_map_zoom.png", plot = plot_zoom)
write.csv(lineage_table, "Cosmo_N_all/Outputs/country_plot.csv", col.names = F, row.names = F)
write.csv(colour_table, "Cosmo_N_all/Outputs/colours.csv", col.names = F, row.names = F)

