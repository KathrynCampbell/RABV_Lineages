create_discrete_scale_for_GLUE <- function(dataframe, n_col, breaks){

  n_breaks = length(breaks)

  # Add a 0 to breaks if not included
  if(breaks[1]!=0){
    breaks <- c(0, breaks)
  }

  # Create dataframe to act as a guide
  temp_df <- data.frame(lower=breaks[1:n_breaks-1],
                        upper=breaks[2:n_breaks])
  temp_df$label <- paste0(temp_df$lower, "-", temp_df$upper)

  # Add discrete categories to dataframe
  for(i in 1:nrow(temp_df)){

    indx <- which(dataframe[[n_col]] >= temp_df$lower[i] & dataframe[[n_col]] < temp_df$upper[i])
    dataframe$cat_scale[indx] <- temp_df$label[i]

  }

  # Set as a factor
  dataframe$cat_scale <- factor(dataframe$cat_scale,
                                        levels=temp_df$label)


  return(dataframe)
}
