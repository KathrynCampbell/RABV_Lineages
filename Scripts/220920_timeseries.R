#'---
#'title: Cosmopolitan Lineage Assignment
#'author: Kathryn Campbell
#'date: 22/07/2020
#'---

rm(list=ls())
library(ggplot2)
#############################################
#            IMPORT THE DATA                #
#############################################
sequence_data <- read.csv("Outputs/sequence_data_cosmo.csv")
node_data <- read.csv("Outputs/node_data_cosmo.csv")
# These were created in the lineage assignment script (Scripts/160920_total_automation.R)

attach(sequence_data)
#############################################
#            IMPORT THE DATA                #
#############################################
# Test some different time series resolutions
# 
# All lineages included
# -------------------------------------------

all<- ggplot(sequence_data, aes(x=Year, fill=cluster))+
  geom_histogram(stat="count") ; all

# Major only
# -------------------------------------------
sequence_data$cluster <- gsub("\\..*","",sequence_data$cluster)
major <- ggplot(sequence_data, aes(x=Year, fill=cluster))+
  geom_histogram(stat="count") ; major

help(ggplot)
major

png("figures/LineageTimeSeries.png", width = 1100, height = 600)
dev.off()