#' write_pao
#'
#' Write occupancy and climate data for program Presence
#' @param counts Data frame containing the buffered BBS counts obtained from the function `buffer_BBS()`
#' @param clim Data frame containing the annual bioclim values obtained from the function 'GetBioVars()`
#' @param alpha Four letter alpha code for the species of interest
#' @return A .pao file containing the detection histories, covariates, and summary information to input into Presence
#' @export

write_pao <- function(counts, clim, alpha){

  tot_stops <- length(grep("count", names(counts))) # Stops/route
  n_seasons <- max(counts$Year) - min(counts$Year) + 1
  common <- code_lookup$common[code_lookup$alpha == toupper(alpha)]

  ### Covert count data to long format
  counts <- dplyr::select(counts, routeID, Year, grep("count", names(counts)))
  counts <- tidyr::gather(counts, key = "stop", value = "n", -routeID, -Year)


  ### Add column with presence/absence data
  pres <- dplyr::mutate(counts, occ = ifelse(n > 0, 1, 0))
  pres <- dplyr::select(pres, -n)

  ### Covert back to wide w/ 1 column for each year/stop (i.e., svy)
  pres <- tidyr::unite(pres, svy, Year, stop, sep = "_")
  pres <- pres[!duplicated(pres),]
  pres <- tidyr::spread(pres, key = svy, value = occ)
  pres <- dplyr::arrange(pres, routeID)

  det_hist <- dplyr::select(pres, -routeID)

  spp_clim <- dplyr::arrange(spp_clim, routeID)
  spp_clim2 <- dplyr::rename(spp_clim, Lat = Latitude, Lon = Longitude)
  spp_clim2 <- dplyr::mutate(spp_clim2, Lat = (Lat - mean(Lat))/sd(Lat),
                                        Lon = (Lon - mean(Lon))/sd(Lon))
  spp_clim3 <- dplyr::mutate(spp_clim2, sq_Lat = Lat ^ 2, sq_Lon = Lon ^ 2)

  for(ss in 1:tot_stops) {
    sc <- scale(1:tot_stops)[ss];        #names(sc)  <- paste0('Stop',ss)
    sc2 <- (scale(1:tot_stops)[ss])^2;   #names(sc2) <- paste0('sqStop',ss)
    spp_clim3   <- cbind(spp_clim3, sc, sc2)
    colnames(spp_clim3)[ncol(spp_clim3)-1] <- paste0('Stop',ss)
    colnames(spp_clim3)[ncol(spp_clim3)] <- paste0('sq_Stop',ss)
  }

  sitecovs	<- dplyr::select(spp_clim3, -routeID)

  spp_pao <- create.pao(data = det_hist, nsurveyseason = tot_stops,
                        unitcov = sitecovs, survcov = NULL,
                        title = paste(common, "PRESENCE Analysis", sep = " "),
                        paoname = paste("inst/output/pao", paste(alpha, ".pao", sep = ""), sep = "/"))
  RPresence::write.pao(pao = spp_pao)
  spp_pao
}
