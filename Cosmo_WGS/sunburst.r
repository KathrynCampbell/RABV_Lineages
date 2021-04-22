rm(list=ls())
library(plotly)
library(treeio)
library(randomcoloR)
library(gridExtra)
if (!require("processx")) install.packages("processx")
library(processx)
library(TreeTools)
library(ggtree)

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

show_col(lineage_info$colour)

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

plot_tree<-plot_tree + theme(legend.position = "left")
ggsave("Cosmo_WGS/Figures/legend_lineage_tree.png",
       plot = plot_tree,
       height = 25, width = 40)
