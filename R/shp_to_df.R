shp_to_df <- function(shapefile){
  shapefile$id <- rownames(shapefile@data)
  #shapefile_df <- tidy(shapefile, region = "id")
  shapefile_df <- fortify(shapefile, region = "id")
  shapefile_df <- left_join(shapefile_df, shapefile@data, by = "id")
  return(shapefile_df)
}