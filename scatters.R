rm(list=ls())

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

#' **Importing the Data**
all_data<-read.csv("Cosmo_WGS/Outputs/Cosmo_WGS_sequence_data.csv", stringsAsFactors=FALSE)
metadata<-read.csv("Cosmo_WGS/Cosmo_WGS_metadata.csv")

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

all_data<-all_data[-c(which(is.na(all_data$cluster))),]

lineage_table <- data.frame(matrix(ncol = (length(unique(all_data$cluster))+1), nrow = length(countries)))
x <- c("country", unique(all_data$cluster))
colnames(lineage_table) <- x

lineage_table$country<-countries

for (i in 1:length(countries)) {
  country_place<-which(lineage_table$country == countries[i])
  lineages<-unique(all_data$cluster[which(all_data$ID %in% metadata$ID[which(metadata$country == countries[i])])])
  
  for (j in 1:length(lineages)) {
    country_lineage<-which(colnames(lineage_table) == lineages[j])
    lineage_table[country_place, country_lineage]<-length(which((all_data$cluster[which(all_data$ID %in% 
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

lineage_table$radius<-lineage_table$radius/2

plot + geom_scatterpie(aes(x=LON, y=LAT, group=region, r=radius),
                       data=lineage_table, cols=c(colnames(lineage_table)[2:74]), color=NA, alpha=.8)+ 
  theme(legend.position = "none")




