#' GetIndices
#'
#' Estimate indices of range dynamics from annual occupancy estimates
#' @param prob_df Data frame containing the annual occupancy probabilities (from GetOcc)
#' @return Data frame containing annual estimates of the following indices:
#' @return    avg.occ = range-wide mean occupancy probability
#' @return    s.lat = southern range limit
#' @return    n.lat = northern range limit
#' @return    avg.lat = occupancy-weighted mean breeding latitude
#' @return    avg.lon = occupancy-weighted mean breeding longitude
#' @export

GetIndices <- function(prob_df){
  prob_grp <- dplyr::group_by(prob_df, Year)
  prob_grp <- dplyr::mutate(prob_grp, cum.occ = cumsum(Prob)/sum(Prob))
  avg.psi <- dplyr::summarise(prob_grp, avg.occ = mean(Prob))

  s.lat <- dplyr::summarise(prob_grp, s.lat = lat[min(which(cum.occ > 0.99))])
  n.lat <- dplyr::summarise(prob_grp, n.lat = lat[min(which(cum.occ > 0.01))])

  avg.lat <- dplyr::summarise(prob_grp, avg.lat = sum(lat * Prob)/sum(Prob))
  avg.lon <- dplyr::summarise(prob_grp, avg.lon = sum(lon * Prob)/sum(Prob))

  indices <- dplyr::left_join(avg.psi, s.lat)
  indices <- dplyr::left_join(indices, n.lat)
  indices <- dplyr::left_join(indices, avg.lat)
  indices <- dplyr::left_join(indices, avg.lon)

  write.csv(indices, file = paste0("inst/output/indices/", alpha, "_ind.csv"), row.name = FALSE)

  indices
}

