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

Sys.setenv("PATH" = "/Users/kathryncampbell/miniconda3/bin")

Sys.getenv("PATH")

lineage_info<-read.csv("Tanz_WGS/Outputs/Tanz_WGS_lineage_info.csv")
node_data<-read.csv("Tanz_WGS/Outputs/Tanz_WGS_node_data.csv")
tree<-read.tree("Tanz_WGS/Trees/Tanz_WGS_aligned.fasta.contree")
metadata<-read.csv("Tanz_WGS/Tanz_WGS_metadata.csv")
sequence_data<-read.csv("Tanz_WGS/Outputs/Tanz_WGS_sequence_data.csv")

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

Colours<-c("Reds","Purples","Greens","Greys","Oranges")

clades<-c("AF1b_A","AF1b_B","AF1b_C","AF1b_D","AF1b_E")

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
  marker = list(colors = list("#DE2D26", "#F27058", "#FCAC92", "#FEE0D2", "#756BB1", "#9894C6", 
                              "#BCBDDC", "#D5D5E8", "#EFEDF5", "#31A354", "#A1D99B", "#E5F5E0", 
                              "#F0F0F0", "#FEE6CE", "#F7FBFF"))
)

new



# orca(new, "Tanz_WGS/Figures/sunburst.png", width = 800, height = 700)
# Better resolution to just save from plot window

# Plot a nice figure to save
plot_tree<-ggtree(tree, colour = "grey50", ladderize = T) %<+% sequence_data +
  geom_tippoint(aes(color=cluster), size=3)  +
  ggtitle("Cosmo N Lineage Tree")+
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
plot_tree

ggsave("Tanz_WGS/Figures/figure_lineage_tree.png",
       plot = plot_tree,
       height = 25, width = 40)
