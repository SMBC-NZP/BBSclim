#' GetSppCounts
#'
#' Filter raw BBS counts by species, add 0 counts, and latitude/longitude
#' @param count Data frame containing the raw BBS counts obtained from the function `GetRouteData()`
#' @param AOU Numeric AOU code for the species of interest
#' @param weather Data frame containing the weather data obtained from the function `GetWeather()`
#' @param routes Data frame containing the route data obtained from the function `GetRoutes()`
#' @param years Optional vector containing the years of interest (default in full study period)
#' @param countrynum Optional numeric vector containing the countries of interest (840 = US, 124 = Canada, 484 = Mexico)
#' @param statenum Optional numeric vector containing the states of interest (see BBS website for code values)
#' @return A .csv file containing the following fields:
#' @return   routeID The unique 8 digit route ID for each route
#' @return   Year The year that the count was conducted
#' @return   aou The numeric code for the focal species
#' @return   countN Total individuals of the focal species recorded at stop N (50-stop data) or N-9:N (10-stop data)
#' @return   Latitude The latitude for the route
#' @return   Longitude The longitude for the route
#' @return   Stratum The BBS stratum for the route
#' @return   BCR The Bird Conservation Region for the route
#' @export

GetSppCounts <- function(count = bbs, AOU, run_atrb = weather, route_atrb = routes,
                         years = seq(from = 1997, to = 2014), statenum = NULL, countrynum = NULL){

  spp_counts <- dplyr::filter(count, aou == AOU & Year %in% years)

  if(!is.null(countrynum)){spp_counts <- dplyr::filter(spp_counts, grepl(paste("^", countrynum, sep = ""), routeID))}
  if(!is.null(statenum)){spp_counts <- dplyr::filter(spp_counts, regexpr(as.character(statenum), routeID) == 4)}


  run_atrb <- dplyr::filter(run_atrb, routeID %in% spp_counts$routeID)
  if(!is.null(years)){run_atrb <- dplyr::filter(run_atrb, Year %in% years)}
  run_atrb <- dplyr::select(run_atrb, routeID, Year, RunType)

  ### Add RunType to count data
  spp_counts <- dplyr::left_join(spp_counts, run_atrb)

  ### Run info for 0 counts (run has weather data but no count data)
  count0 <- dplyr::anti_join(run_atrb, spp_counts)

  ### Add 0 counts to data & fill in AOU code
  spp_counts_full <- dplyr::full_join(count0, spp_counts)
  spp_counts_full$aou <- AOU

  ### Include only runs that meet BBS quality standards
  spp_counts_full <- dplyr::filter(spp_counts_full, RunType == 1)
  spp_counts_full <- dplyr::select(spp_counts_full, -RunType)

  ### Fill in 0 counts
  spp_counts_full[is.na(spp_counts_full)] <- 0

  ### Add longitude and latitude
  route_atrb <- dplyr::select(route_atrb, routeID, Latitude, Longitude, Stratum, BCR)
  spp_counts_full <- dplyr::left_join(spp_counts_full, route_atrb)

  write.csv(spp_counts_full,
            paste("inst/output/spp_counts", paste(alpha, "counts.csv", sep = "_"), sep = "/"),
            row.names = FALSE)

  spp_counts_full
}


