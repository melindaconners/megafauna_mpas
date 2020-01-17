#' Convert Longitude From 0-360 to -180 +180 Format
#'
#' This function allows you to convert longitudes from the 0-360 format common in climate data
#' to the more standard -180 to +180 format.
#'  lon the input vector of longitudes, from 0 to 360
#' 
#'
Lon360to180 <- function(lon){
  ((lon + 180) %% 360) - 180
}


