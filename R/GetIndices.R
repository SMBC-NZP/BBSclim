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

GetIndices <- function(prob_df, ...){

  prob_grp <- dplyr::group_by(prob_df, Year)
  prob_grp <- dplyr::mutate(prob_grp, cum.occ = cumsum(Prob)/sum(Prob))
  avg.psi <- dplyr::summarise(prob_grp, avg.occ = mean(Prob))
  avg.psi$se.occ <- delta(index = "avg.psi", est = avg.psi$avg.occ, ...)

  s.lat <- dplyr::summarise(prob_grp, s.lat = range.limit(prob = Prob, coord = lat, limit = "south"))
  s.lat$se.s.lat <- delta(index = "s.lat", est = s.lat$s.lat, ...)

  n.lat <- dplyr::summarise(prob_grp, n.lat = range.limit(prob = Prob, coord = lat, limit = "north"))
  n.lat$se.n.lat <- delta(index = "n.lat", est = n.lat$n.lat, ...)

  avg.lat <- dplyr::summarise(prob_grp, avg.lat = sum(lat * Prob)/sum(Prob))
  avg.lat$se.avg.lat <- delta(index = "avg.lat", est = avg.lat$avg.lat, ...)

  avg.lon <- dplyr::summarise(prob_grp, avg.lon = sum(lon * Prob)/sum(Prob))
  avg.lon$se.avg.lon <- delta(index = "avg.lon", est = avg.lon$avg.lon, ...)

  indices <- dplyr::left_join(avg.psi, s.lat)
  indices <- dplyr::left_join(indices, n.lat)
  indices <- dplyr::left_join(indices, avg.lat)
  indices <- dplyr::left_join(indices, avg.lon)

  write.csv(indices, file = paste0("inst/output/", alpha, "/indices.csv"), row.names = FALSE)

  indices
}

#' delta
#'
#' Estimate se of indices using delta method

delta <- function(index, est, alpha, years, buff, epslon = 0.1e-10) {
    betas <- GetBetas(alpha)

    len.psi <- length(betas$psi.betas)
    len.gam <- length(betas$gam.betas)
    len.eps <- length(betas$eps.betas)
    all.betas <- c(betas$psi.betas, betas$gam.betas, betas$eps.betas)

    grad <- matrix(0, length(est), length(all.betas))

      for (jj in 1:length(all.betas)) {
        all.betas[jj] <- all.betas[jj] + epslon  #  increment one of the beta's
        betas2 <- list(psi.betas = all.betas[1:len.psi],
                       gam.betas = all.betas[(len.psi+1):(len.psi + len.gam)],
                       eps.betas = all.betas[(len.psi + len.gam + 1):length(all.betas)])
        prob_df2 <- GetOccProb(betas2, alpha, years, buffer = buff)
        prob_grp2 <- dplyr::group_by(prob_df2, Year)
        prob_grp2 <- dplyr::mutate(prob_grp2, cum.occ = cumsum(Prob)/sum(Prob))
        if(index == "avg.psi"){
          est2 <- dplyr::summarise(prob_grp2, avg.occ = mean(Prob))$avg.occ
        }

        if(index == "s.lat"){
          est2 <- dplyr::summarise(prob_grp2, s.lat = range.limit(prob = Prob, coord = lat, limit = "south"))$s.lat
        }

        if(index == "n.lat"){
          est2 <- dplyr::summarise(prob_grp2, n.lat = range.limit(prob = Prob, coord = lat, limit = "north"))$n.lat
        }

        if(index == "avg.lat"){
          est2 <- dplyr::summarise(prob_grp2, avg.lat = sum(lat * Prob)/sum(Prob))$avg.lat
        }

        if(index == "avg.lon"){
          est2 <- dplyr::summarise(prob_grp2, avg.lon = sum(lon * Prob)/sum(Prob))$avg.lon
        }

        grad[, jj] <- (est2 - est)/epslon  #  partial derivative = change in real parm / change in betas
        all.betas[[jj]] <- all.betas[[jj]] - epslon  #  increment one of the beta's
      }

    r.vc <- grad %*% betas$vc.mat %*% t(grad) #  real vc matrix (r.vc) = grad * betaVC * grad'
    se <- sqrt(diag(r.vc))
    se
  }

#' range.limit
#'
#' Estimate range limit using cumulative probability method

range.limit <- function(prob, coord, limit){
  xy <- data.frame(x = prob, y = coord)
  x <- arrange(xy, y)$x/sum(xy$x)
  y <- arrange(xy, y)$y

    if(limit == "south"){
      lim <- predict(smooth.spline(cumsum(x), y, spar=0.1), 0.01)$y
    }

    if(limit == "north"){
      lim <- predict(smooth.spline(cumsum(x), y, spar=0.1), 0.99)$y
    }
  lim
  }

