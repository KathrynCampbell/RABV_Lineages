#'---
#'title: Sequence Data Function
#'author: Kathryn Campbell
#'date: 04/09/2020
#'---
#'
#'
#'**SEQUENCE DATA**
#'=================
#'
#'File: sequence_data.r
#'=====================


sequence_data <- function(alignment) {
  # Make a dataframe ready to fill with info about number of gaps and N bases, and length of the alignment and sequence
  seq_data <- data.frame(ID = alignment$nam, N = NA, "-" = NA,
                         Length_before = nchar(alignment$seq[[1]]), Length_after = NA)
  
  for (i in 1:(length(alignment$seq))) {
    seq_data$N[i] <- str_count(alignment$seq[[i]], pattern = 'n')
    seq_data$`-`[i] <- str_count(alignment$seq[[i]], pattern = '-')
    seq_data$Length_after[i] <- (seq_data$Length_before[i] - seq_data$N[i] - seq_data$`-`[i])
  }
  # For each sequence, count the number of n bases and the number of gaps
  # Calculate the length after removing these 
  
  seq_data
}

