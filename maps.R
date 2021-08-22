library(plotly)
library(treeio)
library(randomcoloR)
library(gridExtra)
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
library(ggnewscale)

rm(list=ls())
N_col<-read.csv("Cosmo_N_all/Outputs/colours.csv")
WGS_col<-read.csv("Cosmo_WGS/Outputs/colours.csv")
N_lineage<-read.csv("Cosmo_N_all/Outputs/country_plot.csv")
WGS_lineage<-read.csv("Cosmo_WGS/Outputs/country_plot.csv")


t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

for (i in 1:length(N_col$colour)) {
  N_col$colour[i]<-t_col(N_col$colour[i], percent = 85)
}

world <- ne_countries(scale = "medium", returnclass = "sf")
plot<-
  ggplot(data = world) +
  geom_sf(fill = "lightgrey", colour = "white") +
  theme_bw()

N_lineage$radius[which(N_lineage$radius %in% 1:5)]<-3
N_lineage$radius[which(N_lineage$radius %in% 6:10)]<-5
N_lineage$radius[which(N_lineage$radius %in% 11:25)]<-7
N_lineage$radius[which(N_lineage$radius %in% 26:50)]<-9
N_lineage$radius[which(N_lineage$radius %in% 51:100)]<-11
N_lineage$radius[which(N_lineage$radius %in% 101:200)]<-13
N_lineage$radius[which(N_lineage$radius > 200)]<-15

WGS_lineage$radius[which(WGS_lineage$radius %in% 1:5)]<-3
WGS_lineage$radius[which(WGS_lineage$radius %in% 6:10)]<-5
WGS_lineage$radius[which(WGS_lineage$radius %in% 11:25)]<-7
WGS_lineage$radius[which(WGS_lineage$radius %in% 26:50)]<-9
WGS_lineage$radius[which(WGS_lineage$radius %in% 51:100)]<-11
WGS_lineage$radius[which(WGS_lineage$radius %in% 101:200)]<-13
WGS_lineage$radius[which(WGS_lineage$radius > 200)]<-15


for (i in 2:138) {
  names<-c(names, paste("N", (names(N_lineage[i])), sep = "_"))
}

names<-names[-c(1)]
names(N_lineage)<-c("country", names, "LAT", "LON", "region", "radius")
N_col$lineage<-names(N_lineage[2:138])

new_scale <- function(new_aes) {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}

N_col$white<-"white"

plot_world<-plot + geom_scatterpie(aes(x=LON, y=LAT, group=region, r=radius/3),
                                   data=N_lineage, cols=c(colnames(N_lineage)[2:138]), color=NA)+
  scale_fill_manual(values = c(N_col$white))+
  new_scale_fill()+
  geom_scatterpie(aes(x=LON, y=LAT, group=region, r=radius/3),
                  data=N_lineage, cols=c(colnames(N_lineage)[2:138]), color=NA, alpha=.25)+
  scale_fill_manual(values = c(N_col$colour))+
  new_scale_fill()+
  geom_scatterpie(aes(x=LON, y=LAT, group=region, r=radius/3),
                  data=WGS_lineage, cols=c(colnames(WGS_lineage)[2:74]), color=NA)+
  theme(legend.position = "none") +scale_fill_manual(values = c(WGS_col$colour))+
  coord_sf(xlim = c(-130, 130), ylim = c(-45, 65)) +
  theme(panel.grid.major = element_blank())+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())
  
plot_world
ggsave("Cosmo_WGS/Figures/pie_map.png", plot = plot_world, width = 45, height = 20)



