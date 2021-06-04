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

lineage_info<-read.csv("Cosmo_N_all/Outputs/Cosmo_N_all_lineage_info.csv")
node_data<-read.csv("Cosmo_N_all/Outputs/Cosmo_N_all_node_data.csv")
tree<-read.tree("Cosmo_N_all/Trees/Cosmo_N_all_aligned.fasta.contree")
metadata<-read.csv("Cosmo_N_all/Cosmo_N_all_metadata.csv")
sequence_data<-read.csv("Cosmo_N_all/Outputs/Cosmo_N_all_sequence_data.csv")

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
  marker = list(colors = list("#08306B", "#083573", "#083B7B", "#084084", "#08468C", "#084C94", 
                              "#08519C", "#0C57A0", "#115CA5", "#1562A9", "#1967AD", "#1D6CB1", 
                              "#2272B5", "#2878B8", "#2D7DBB", "#3383BE", "#3888C1", "#3E8EC4", 
                              "#DE2D26", "#E0342B", "#E23B30", "#E44236", "#E6493B", "#E85141", 
                              "#EA5846", "#ED5F4C", "#EF6651", "#F16D56", "#F3755C", "#F57C61", 
                              "#F78367", "#F98A6C", "#FC9272", "#FC9778", "#FC9D7F", "#FCA286", 
                              "#FCA88D", "#FCAD94", "#FCB39B", "#FDB9A2", "#FDBEA8", "#FDC4AF", 
                              "#FDC9B6", "#FDCFBD", "#FDD4C4", "#FDDACB", "#FEE0D2", "#756BB1", 
                              "#7D74B6", "#857EBB", "#8E87C0", "#9691C5", "#9E9BCA", "#A7A4CF", 
                              "#AFAED4", "#B7B8D9", "#BFBFDD", "#C5C5E0", "#CBCBE3", "#D1D0E6", 
                              "#D7D6E9", "#DDDCEC", "#E3E1EF", "#E9E7F2", "#EFEDF5", "#FFF7BC", 
                              "#1C9099", "#A6BDDB", "#ECE2F0", "#F03B20", "#F45D2C", "#F87F39", 
                              "#FCA145", "#FEBA58", "#FECB70", "#FEDC88", "#FFEDA0", "#E34A33", 
                              "#E75C40", "#EB6F4E", "#F0825B", "#F49569", "#F8A876", "#FDBB84", 
                              "#FDC28F", "#FDCA9A", "#FDD1A6", "#FDD9B1", "#FDE0BC", "#FEE8C8", 
                              "#2B8CBE", "#A6BDDB", "#ECE7F2", "#CCEBC5", "#FBB4AE", "#31A354", 
                              "#43AB5F", "#56B46B", "#69BE77", "#7BC783", "#8ED08F", "#A1D99B", 
                              "#ACDDA6", "#B7E2B2", "#C3E7BD", "#CEEBC9", "#D9F0D4", "#E5F5E0", 
                              "#636363", "#797979", "#909090", "#A6A6A6", "#BDBDBD", "#C9C9C9", 
                              "#D6D6D6", "#E3E3E3", "#F0F0F0", "#43A2CA", "#6BB9C1", "#93D1B9", 
                              "#B3E1BC", "#C9EACB", "#E0F3DB", "#E5F5F9", "#C51B8A", "#DA4F9B", 
                              "#EF84AC", "#FAACBD", "#FBC6CD", "#FDE0DD", "#E6550D", "#EA661F", 
                              "#EF7832", "#F38A45", "#F89C58", "#FDAE6B", "#FDB97E", "#FDC492", 
                              "#FDCFA6", "#FDDABA", "#FEE6CE", "#8856A7", "#969AC9", "#B4CCE2", 
                              "#E0ECF4", "#31A354", "#62BA6B", "#94D182", "#BBE396", "#D9EFA7", 
                              "#F7FCB9", "#DD1C77", "#E7E1EF", "#2C7FB8", "#63B3BA", "#A3DBB7", 
                              "#EDF8B1", "#4493C7", "#4B98C9", "#529DCC", "#59A2CF", "#60A6D1", 
                              "#67ABD4", "#6FB0D6", "#78B5D8", "#80B9DA", "#89BEDC", "#92C3DE", 
                              "#9AC8E0", "#A2CBE2", "#A9CEE4", "#AFD1E7", "#B6D4E9", "#BDD7EC", 
                              "#C4DAEE", "#C9DDF0", "#CDDFF1", "#D1E2F2", "#D5E5F4", "#D9E7F5", 
                              "#DDEAF6", "#E1EDF8", "#E5F0F9", "#EAF2FA", "#EEF5FC", "#F2F8FD", 
                              "#F7FBFF"))
)

new

previous_colours<-lineage_info$colour[-c(grep("_", lineage_info$cluster))][1]

for (i in 1:length(clades)) {
  previous_colours<-c(previous_colours,
                      lineage_info$colour[(grep(clades[i], lineage_info$cluster))][1])
}

previous$colours<-NA

for (i in 1:length(previous$assignment)) {
  previous$colours[i]<-lineage_info$colour[grep(previous$assignment[i], lineage_info$cluster)][1]
}

previous$colours[which(is.na(previous$colours))]<-c("#FDD835","#F8BBD0","#8D6E63","#FFCCBC","#CDDC39")

show_col(previous$colours)

dput(previous$colours)

previous_sun<-plot_ly(
  labels = c(previous$assignment),
  parents = c(previous$parent),
  values = c(previous$n_seqs),
  type = "sunburst",
  width = 1000,
  height = 700,
  marker = list(colors = list ("#2C7FB8", "#E34A33", "#2B8CBE", "#756BB1", "#31A354", "#1C9099", 
                               "#DD1C77", "#43A2CA", "#FDD835", "#E5F5F9", "#DE2D26", "#636363", 
                               "#F03B20", "#F8BBD0", "#31A354", "#E6550D", "#8D6E63", "#8856A7", 
                               "#FFCCBC", "#FFF7BC", "#CCEBC5", "#C51B8A", "#CDDC39"))
)
previous_sun


orca(new, "Cosmo_N_all/Figures/sunburst.png", width = 800, height = 700)
orca(previous_sun, "Cosmo_N_all/Figures/previous_sunburst.png", width = 800, height = 700)
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

ggsave("Cosmo_N_all/Figures/figure_lineage_tree.png",
       plot = plot_tree,
       height = 25, width = 40)
