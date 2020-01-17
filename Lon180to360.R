#' Convert Longitude From -180 - +180 to 0 0 360 Format
#'
#' This function allows you to convert longitudes ti the 0-360 format common in climate data
#' from the more standard -180 to +180 format.


Lon180to360 <- function(lon){
  lon %% 360
}
