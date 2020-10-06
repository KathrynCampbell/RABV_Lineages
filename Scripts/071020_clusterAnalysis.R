#'---------------------------------------------------------
#'title: Tanzanian lineage analysis
#'author: Kathryn Campbell
#'date: 07/10/2020
#'---------------------------------------------------------

rm(list = ls())

#############################################
#            INSTALL PACKAGES               #
#############################################
library(dplyr)

#############################################
#         CLUSTER ANALYSIS TABLE            #
#############################################
sequence_data <- read.csv(file = "Outputs/sequence_data_cosmo.csv")
metadata <- read.csv(file = "Sequences/220720_GLUE_CosmoMeta.csv")
# Import the data about assignment for each sequence
# This is created in the script scripts/160920_total_automation.R

for (i in 1:length(sequence_data$ID)) {
  sequence_data$host[i] <- metadata$sequence.host[which(metadata$sequence.sequenceID == sequence_data$ID[i])]
}
# Pull the information about the host for each sequence from the metadata into the sequence_data

clusters <- sequence_data %>%
  group_by(cluster) %>%
  summarise()
# Create a data frame ready to fill in information about each cluster

for (i in 1:length(clusters$cluster)) {
  clusters$year_first[i] <- sequence_data %>%
    filter(cluster == clusters$cluster[i])%>%
    group_by(Year)%>%
    summarise()%>%
    min()
  clusters$year_last[i] <- sequence_data %>%
    filter(cluster == clusters$cluster[i])%>%
    group_by(Year)%>%
    summarise()%>%
    max()
  clusters$country[i] <- sequence_data %>%
    filter(cluster == clusters$cluster[i])%>%
    group_by(country)%>%
    summarise()
  clusters$host[i] <- sequence_data %>%
    filter(cluster == clusters$cluster[i])%>%
    group_by(host)%>%
    summarise()
}
# For each cluster, find and list the earliest colleciton year, the latest collection year and all the places
# and hosts that cluster has been found in

clusters$n_seqs<-(sequence_data %>%
                    group_by(cluster)%>%
                    summarise(n=n()))$n
# Add another column listing the number of sequences assigned to each cluster