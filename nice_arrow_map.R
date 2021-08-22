rm(list=ls())

metadata_N<-read.csv("Cosmo_N_all/Cosmo_N_all_metadata.csv")
seq_info_N<-read.csv("Cosmo_N_all/Outputs/Cosmo_N_all_sequence_data.csv")
lineages_N<-read.csv("Cosmo_N_all/Outputs/Cosmo_N_all_lineage_info.csv")

test<-"Cosmopolitan AM3a_A1.3"

lineage_info_N<-lineages_N[which(lineages_N$cluster == test),]
lineage_meta_N<-metadata_N[which(metadata_N$ID %in% seq_info_N$ID[which(seq_info_N$cluster == test)]),]

metadata<-read.csv("Cosmo_WGS/Cosmo_WGS_metadata.csv")
seq_info<-read.csv("Cosmo_WGS/Outputs/Cosmo_WGS_sequence_data.csv")
lineages<-read.csv("Cosmo_WGS/Outputs/Cosmo_WGS_lineage_info.csv")

test<-"D1.1.2"

lineage_info<-lineages[which(lineages$cluster == test),]
lineage_meta<-metadata[which(metadata$ID %in% seq_info$ID[which(seq_info$cluster == test)]),]

lineage_meta<-lineage_meta[-c(which(lineage_meta$ID %in% lineage_meta_N$ID)),]
lineage_meta<-rbind(lineage_meta, lineage_meta_N)
lineage_info<-rbind(lineage_info, lineage_info_N)

rm(lineage_info_N)
rm(lineage_meta_N)

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

world<-ne_countries()

#' **Cleaning the data**
#' Find which country names do not match between data and map file
map_countries <- as.character(world$admin)
countries <- unique(lineage_meta$country)
no_match <- setdiff(countries, map_countries); message(length(no_match), " countries are mis-matched: \n", paste0(no_match, collapse="\n"))

lineage_meta$country[grepl("USA|United States", lineage_meta$country)] <- "United States of America"
lineage_meta$country[grepl("Tanzania", lineage_meta$country)] <- "United Republic of Tanzania"


#' How many entries are there for the remaining mis-matches?
table(lineage_meta$country[which(lineage_meta$country %in% no_match)])

#' ** Tabulating the data**
country_table<-table(lineage_meta$country); country_table
length(country_table)

lineage_meta$LAT<-NA
lineage_meta$LON<-NA

lineage_meta<-lineage_meta[-c(which(lineage_meta$country == "-")),]

for (i in 1:length(lineage_meta$ID)) {
  lineage_meta$LAT[i]<-world@polygons[which(world$admin == lineage_meta$country[i])][[1]]@labpt[2]
  lineage_meta$LON[i]<-world@polygons[which(world$admin == lineage_meta$country[i])][[1]]@labpt[1]
}

countries<-unique(lineage_meta$country)

introductions<-NA

for (i in 1:length(countries)) {
  introductions<-rbind(introductions,     lineage_meta[
    which(lineage_meta$year == min(
      lineage_meta$year[which(lineage_meta$country == countries[i])]))[1]
    ,])
}


introductions<-introductions[-c(1),]

introductions<-introductions[order(introductions$year),]

WGS<-introductions[which(introductions$sequence.gb_length > 10000),]
N<-introductions[-c(which(introductions$sequence.gb_length > 10000)),]

# Change this for different ones
world <- ne_countries(scale = "medium", returnclass = "sf")

plot<-
  ggplot(data = world) +
  geom_sf(fill = "lightgrey", colour = "white") +
  theme_bw() +
  geom_point(data = WGS, aes(x = LON, y = LAT), size = 4, 
           shape = 21, fill = "red", show.legend = T) +
  geom_point(data = N, aes(x = LON, y = LAT), size = 4, 
             shape = 21, fill = "darkred", show.legend = T) +
  coord_sf(xlim = c(-155, -30), ylim = c(-70, 60)) +
  geom_segment(x = introductions$LON[1], y = introductions$LAT[1], 
               xend = introductions$LON[2],yend = introductions$LAT[2],
               arrow = arrow(length = unit(0.03, "npc")), size = 1, colour = "darkred") +
  geom_segment(x = introductions$LON[1], y = introductions$LAT[2], 
               xend = introductions$LON[3],yend = introductions$LAT[3],
               arrow = arrow(length = unit(0.03, "npc")), size = 1, colour = "red")+
  theme(panel.grid.major = element_blank())+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())

for (i in 1:length(introductions$ID)) {
  plot<-plot+
    annotate("text", x=(introductions$LON[i]-3), y=(introductions$LAT[i]-2.5), label= introductions$year[i], size = 4)
}

plot

ggsave("Cosmo_WGS/Figures/figures/case_study_1.png")
