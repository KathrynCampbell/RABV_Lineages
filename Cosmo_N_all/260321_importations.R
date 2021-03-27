rm(list=ls())

metadata<-read.csv("Cosmo_N_all/Cosmo_N_all_metadata.csv")
seq_info<-read.csv("Cosmo_N_all/Outputs/Cosmo_N_all_sequence_data.csv")
lineages<-read.csv("Cosmo_N_all/Outputs/Cosmo_N_all_lineage_info.csv")

test<-"NEE_A1.2.1"

lineage_info<-lineages[which(lineages$cluster == test),]
lineage_meta<-metadata[which(metadata$ID %in% seq_info$ID[which(seq_info$cluster == test)]),]

library(rgdal)
library(ggplot2)
library(rworldmap)
library(cleangeo)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(mapproj)

#' Load functions
source("R/shp_to_df.R")
source("R/create_discrete_scale_for_GLUE.R")
#' These scripts were created by Rachel Steenson 31/01/2020

metadata_WE<-lineage_meta


#' Download world map
world_map <- getMap()
plot(world_map)

#' **Cleaning the data**
#' Find which country names do not match between data and map file
map_countries <- as.character(world_map$ADMIN)
countries <- unique(metadata_WE$country)
no_match <- setdiff(countries, map_countries); message(length(no_match), " countries are mis-matched: \n", paste0(no_match, collapse="\n"))

metadata_WE$country[grepl("Ivoire", metadata_WE$country)] <- "Ivory Coast" # The weird characters mean we have to do a string search!
metadata_WE$country[grepl("Tanzania", metadata_WE$country)] <- "United Republic of Tanzania"
metadata_WE$country[grepl("USA|United States", metadata_WE$country)] <- "United States of America"
metadata_WE$country[grepl("Columbia", metadata_WE$country)] <- "Colombia"
metadata_WE$country[grepl("Serbia", metadata_WE$country)] <- "Republic of Serbia"
metadata_WE$country[grepl("Czechia", metadata_WE$country)] <- "Czech Republic"

#' How many entries are there for the remaining mis-matches?
table(metadata_WE$country[which(metadata_WE$country %in% no_match)])

#' ** Tabulating the data**
country_table<-table(metadata_WE$country); country_table
length(country_table)

years1<-metadata_WE %>%
  filter(year %in% 1999) %>%
  group_by(country) %>%
  summarise(n=n()) ; years1

years2<-metadata_WE %>%
  filter(year %in% 2007) %>%
  group_by(country) %>%
  summarise(n=n()) ; years2

years3<-metadata_WE %>%
  filter(year %in% 2011) %>%
  group_by(country) %>%
  summarise(n=n()) ; years3

years4<-metadata_WE %>%
  filter(year %in% 2004:2007) %>%
  group_by(country) %>%
  summarise(n=n()) ; years4

years5<-metadata_WE %>%
  filter(year %in% 2008:2014) %>%
  group_by(country) %>%
  summarise(n=n()) ; years5


#' **PROCESS SHAPEFILE DATA**
#' Clean shapefile - if you run the shp_to_df function without this, you get an error (not normally required!)
world_map <- clgeo_Clean(world_map)

#' Merge summary data with world map
world_map_data_years1 <- merge(world_map, years1, by.x="ADMIN", by.y="country")
world_map_data_years2 <- merge(world_map, years2, by.x="ADMIN", by.y="country")
world_map_data_years3 <- merge(world_map, years3, by.x="ADMIN", by.y="country")
world_map_data_years4 <- merge(world_map, years4, by.x="ADMIN", by.y="country")
world_map_data_years5 <- merge(world_map, years5, by.x="ADMIN", by.y="country")


#' Set number of sequences==NA as 0
world_map_data_years1[is.na(world_map_data_years1$n)] <- 0
table(world_map_data_years1$n, useNA="always")

world_map_data_years2[is.na(world_map_data_years2$n)] <- 0
table(world_map_data_years2$n, useNA="always")

world_map_data_years3[is.na(world_map_data_years3$n)] <- 0
table(world_map_data_years3$n, useNA="always")

world_map_data_years4[is.na(world_map_data_years4$n)] <- 0
table(world_map_data_years4$n, useNA="always")

world_map_data_years5[is.na(world_map_data_years5$n)] <- 0
table(world_map_data_years5$n, useNA="always")

world_map_data_years1_df <- shp_to_df(world_map_data_years1)
world_map_data_years2_df <- shp_to_df(world_map_data_years2)
world_map_data_years3_df <- shp_to_df(world_map_data_years3)
world_map_data_years4_df <- shp_to_df(world_map_data_years4)
world_map_data_years5_df <- shp_to_df(world_map_data_years5)


summary(world_map_data_years1_df$n)
summary(world_map_data_years2_df$n)
summary(world_map_data_years3_df$n)
summary(world_map_data_years4_df$n)
summary(world_map_data_years5_df$n)

# Colours
br <- c(0, 1, 2, 3, 4, 5, 10)
cols = c("white", colorRampPalette(c("#ffe5e5", "#b30000"))(length(br)-1))

# Map breaks to colours
col_pal <- c("0-1"=cols[1], "1-2"=cols[2], "2-3"=cols[3], "3-4"=cols[4], "4-5"=cols[5], "5-10"=cols[6])

#' years 1
world_map_data_years1_df_catscale1 <- create_discrete_scale_for_GLUE(dataframe=world_map_data_years1_df, n_col="n",
                                                                     breaks=br)
years1_graph<-ggplot() +
  geom_polygon(data=world_map_data_years1_df_catscale1, aes(x=long, y=lat, group=group, fill=cat_scale), col="black") +
  scale_fill_manual(name="Number of \nSequences", values=col_pal, drop=FALSE) +
  theme_void() +
  coord_equal() + 
  coord_cartesian(xlim = c(-180, -30), ylim = c(-60, 80))+
  ggtitle("1995") + theme(plot.title = element_text(size=10))
years1_graph

#' years 2
world_map_data_years2_df_catscale1 <- create_discrete_scale_for_GLUE(dataframe=world_map_data_years2_df, n_col="n",
                                                                     breaks=br)
years2_graph<-ggplot() +
  geom_polygon(data=world_map_data_years2_df_catscale1, aes(x=long, y=lat, group=group, fill=cat_scale), col="black") +
  scale_fill_manual(name="Number of \nSequences", values=col_pal, drop=FALSE) +
  theme_void() +
  coord_equal() +
  coord_cartesian(xlim = c(-180, -30), ylim = c(-60, 80))+
  ggtitle("1996-2000") + theme(plot.title = element_text(size=10))
years2_graph

#' years 3
world_map_data_years3_df_catscale1 <- create_discrete_scale_for_GLUE(dataframe=world_map_data_years3_df, n_col="n",
                                                                     breaks=br)
years3_graph<-ggplot() +
  geom_polygon(data=world_map_data_years3_df_catscale1, aes(x=long, y=lat, group=group, fill=cat_scale), col="black") +
  scale_fill_manual(name="Number of \nSequences", values=col_pal, drop=FALSE) +
  theme_void() +
  coord_equal() +
  coord_cartesian(xlim = c(-180, -30), ylim = c(-60, 80))+
  ggtitle("2002") + theme(plot.title = element_text(size=10))
years3_graph

#'years 4
world_map_data_years4_df_catscale1 <- create_discrete_scale_for_GLUE(dataframe=world_map_data_years4_df, n_col="n",
                                                                     breaks=br)
years4_graph<-ggplot() +
  geom_polygon(data=world_map_data_years4_df_catscale1, aes(x=long, y=lat, group=group, fill=cat_scale), col="black") +
  scale_fill_manual(name="Number of \nSequences", values=col_pal, drop=FALSE) +
  theme_void() +
  coord_equal() +
  coord_cartesian(xlim = c(-150, -30), ylim = c(-60, 40))+
  ggtitle("2004-2007") + theme(plot.title = element_text(size=10))
years4_graph

#' years 5
world_map_data_years5_df_catscale1 <- create_discrete_scale_for_GLUE(dataframe=world_map_data_years5_df, n_col="n",
                                                                     breaks=br)
years5_graph<-ggplot() +
  geom_polygon(data=world_map_data_years5_df_catscale1, aes(x=long, y=lat, group=group, fill=cat_scale), col="black") +
  scale_fill_manual(name="Number of \nSequences", values=col_pal, drop=FALSE) +
  theme_void() +
  coord_equal() +
  coord_cartesian(xlim = c(-150, -30), ylim = c(-60, 40))+
  ggtitle("2008-2014") + theme(plot.title = element_text(size=10))
years5_graph


combined<-grid.arrange(years1_graph, years2_graph, years3_graph, 
                       nrow = 1, ncol = 3, top = test)

png(paste("Cosmo_N_all/Figures/Lineage_timemaps/Interesting/", test, ".png", sep = ""), width = 2000, height = 1000)
grid.arrange(years1_graph, years2_graph, years3_graph, 
             nrow = 1, ncol = 3, top = test)
dev.off()

write.csv(lineage_info, file = paste("Cosmo_N_all/Figures/Lineage_timemaps/Interesting/", test, ".csv", sep = ""))

