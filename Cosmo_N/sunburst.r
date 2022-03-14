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

lineage_info<-read.csv("Cosmo_N/Outputs/Cosmo_N_updated_lineage_info.csv")
tree<-read.tree("Cosmo_N/Trees/Cosmo_N_aligned.fasta.contree")
outgroup<-read.tree("Cosmo_N/Trees/Cosmo_N_outgroup_aligned.fasta.contree")
outgroup_meta<-read.csv("Cosmo_WGS/Cosmo_WGS_metadata_outgroup.csv")
metadata<-read.csv("Cosmo_N/Cosmo_N_metadata.csv")
sequence_data<-read.csv("Cosmo_N/Outputs/Cosmo_N_sequence_data.csv")

metadata<-metadata[which(metadata$ID %in% sequence_data$ID),]

metadata<-rbind(metadata, outgroup_meta[grep("Asian", outgroup_meta$assignment),])

outgroup<-reroot(outgroup, getMRCA(
  outgroup, tip = which(outgroup$tip.label %in% metadata$ID[(
    grep("Asian", metadata$assignment))])))

metadata<-metadata[-c(grep("Asian", metadata$assignment)),]

outgroup<-tree_subset(outgroup, 2138, levels_back = 0)

lineage_info$colour<-NA

Colours<-c("Reds","Purples","YlOrBr","PuBuGn","YlOrRd","OrRd","PuBu","Pastel1","Greens","Greys",
           "GnBu","BuGn","RdPu","Oranges","BuPu","YlGn","PuRd","YlGnBu")

clades<-c("AF1a","AF1b","AF4","AM1","AM2a","AM3a","AM3b","AM4","CA1","CA2","CA3","CE","Cosmopolitan EE",
          "ME1a","ME2","NEE","WE","Cosmopolitan_")

lineage<-lineage_info$lineage[-c(grep("_", lineage_info$lineage))]
cols<-brewer.pal(9, "Blues")
pal<-colorRampPalette(c(cols))
pal<-rev(pal(length(lineage)))
lineage_info$colour[-c(grep("_", lineage_info$lineage))]<-pal

for (i in 1:length(clades)) {
  lineage<-grep(clades[i], lineage_info$lineage)
  cols<-brewer.pal(3, Colours[i])
  pal<-colorRampPalette(c(cols))
  pal<-rev(pal(length(lineage)))
  lineage_info$colour[(grep(clades[i], lineage_info$lineage))]<-pal
}

lineage_info$lineage<-gsub("Cosmopolitan ", "", lineage_info$lineage)
lineage_info$parent<-gsub("Cosmopolitan ", "", lineage_info$parent)

new<-plot_ly(
  labels = c(lineage_info$lineage),
  parents = c(lineage_info$parent),
  values = c(lineage_info$n_seqs),
  type = "sunburst",
  width = 700,
  height = 700,
  marker = list(colors = lineage_info$colour)
)

new

htmlwidgets::saveWidget(plotly::as_widget(new), "Cosmo_N/Figures/sunburst.html")

# Plot a nice figure to save
plot_tree<-ggtree(outgroup, colour = "grey50", ladderize = T) %<+% sequence_data +
  geom_tippoint(color="grey50", size=4)+
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

ggsave("Cosmo_N/Figures/figure_lineage_tree.png",
       plot = plot_tree)

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
  lineage_table$radius[i]<-sum(lineage_table[i,2:79])
}


lineage_table$radius[which(lineage_table$radius %in% 6:10)]<-6
lineage_table$radius[which(lineage_table$radius %in% 11:20)]<-7
lineage_table$radius[which(lineage_table$radius %in% 21:50)]<-8
lineage_table$radius[which(lineage_table$radius %in% 51:100)]<-9
lineage_table$radius[which(lineage_table$radius %in% 101:200)]<-10
lineage_table$radius[which(lineage_table$radius > 200)]<-11

head(lineage_table)

lineage_table$radius<-lineage_table$radius/1.5

colour_table<-data.frame(lineage = colnames(lineage_table[2:79]), colour = NA)
colour_table$lineage<-gsub("Cosmopolitan ", "", colour_table$lineage)

for (i in 1:length(colour_table$lineage)) {
  colour_table$colour[i]<-lineage_info$colour[which(lineage_info$lineage == colour_table$lineage[i])]
}

plot_world<-plot + geom_scatterpie(aes(x=LON, y=LAT, group=region, r=radius),
                                   data=lineage_table, cols=c(colnames(lineage_table)[2:79]), color=NA, alpha=.8)+ 
  theme(legend.position = "none") +scale_fill_manual(values = c(colour_table$colour))+
  coord_sf(xlim = c(-155, 155), ylim = c(-50, 70)) +
  theme(panel.grid.major = element_blank())+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())

plot_world
ggsave("Cosmo_N/Figures/pie_map.png", plot = plot_world, width = 45, height = 20)

