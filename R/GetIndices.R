#' GetIndices
#'
#' Estimate indices of range dynamics from annual occupancy estimates
#' @param prob_df Data frame containing the annual occupancy probabilities (from GetOcc)
#' @return Data frame containing annual estimates of the following indices:
#' @return    avg.psi = proportion of area occupied (i.e., range size)
#' @return    s.lat = southern range limit
#' @return    n.lat = northern range limit
#' @return    avg.lat = occupancy-weighted mean breeding latitude
#' @return    avg.lon = occupancy-weighted mean breeding longitude
#' @export

GetIndices <- function(alpha){
  prob_df <- read.csv(paste0('inst/output/', alpha, '/occ.csv'))
  spp_buff <- read.csv(paste0('inst/output/', alpha, '/count_buff.csv'))

  opts <- read.csv('inst/global_opts.csv')
  years <- seq(from = opts$start_yr, to = opts$end_yr)

  prob_grp <- dplyr::group_by(prob_df, Year)

  avg.psi <- dplyr::summarise(prob_grp, value = mean(Prob))
  avg.psi$sd.err <- delta(index = "avg.psi", est = avg.psi$value, alpha)
  avg.psi$ind <- "avg.psi"

  s.lat <- dplyr::summarise(prob_grp, value = range.limit(cell.probs = Prob, prob = 0.99, coord = lat, limit = "south"))
  s.lat$sd.err <- delta(index = "s.lat", est = s.lat$value, alpha)
  s.lat$ind <- "s.lat"

  s.core <- dplyr::summarise(prob_grp, value = range.limit(cell.probs = Prob, prob = 0.75, coord = lat, limit = "south"))
  s.core$sd.err <- delta(index = "s.core", est = s.core$value, alpha)
  s.core$ind <- "s.core"

  n.lat <- dplyr::summarise(prob_grp, value = range.limit(cell.probs = Prob, prob = 0.99, coord = lat, limit = "north"))
  n.lat$sd.err <- delta(index = "n.lat", est = n.lat$value, alpha)
  n.lat$ind <- "n.lat"

  n.core <- dplyr::summarise(prob_grp, value = range.limit(cell.probs = Prob, prob = 0.75, coord = lat, limit = "north"))
  n.core$sd.err <- delta(index = "n.core", est = n.core$value, alpha)
  n.core$ind <- "n.core"

  avg.lat <- dplyr::summarise(prob_grp, value = sum(lat * Prob)/sum(Prob))
  avg.lat$sd.err <- delta(index = "avg.lat", est = avg.lat$value, alpha)
  avg.lat$ind <- "avg.lat"

  avg.lon <- dplyr::summarise(prob_grp, value = sum(lon * Prob)/sum(Prob))
  avg.lon$sd.err <- delta(index = "avg.lon", est = avg.lon$value, alpha)
  avg.lon$ind <- "avg.lon"

  indices <- dplyr::bind_rows(avg.psi, s.lat)
  indices <- dplyr::bind_rows(indices, s.core)
  indices <- dplyr::bind_rows(indices, n.lat)
  indices <- dplyr::bind_rows(indices, n.core)
  indices <- dplyr::bind_rows(indices, avg.lat)
  indices <- dplyr::bind_rows(indices, avg.lon)

  write.csv(indices, file = paste0("inst/output/", alpha, "/indices.csv"), row.names = FALSE)
}

#' delta
#'
#' Estimate se of indices using delta method

delta <- function(index, est, alpha, epslon = 0.1e-10) {
    betas <- GetBetas(alpha)

    len.psi <- length(betas$psi.betas)
    len.gam <- length(betas$gam.betas)
    len.eps <- length(betas$eps.betas)
    all.betas <- c(betas$psi.betas, betas$gam.betas, betas$eps.betas)

    grad <- matrix(0, length(est), length(all.betas))

      for (jj in 1:length(all.betas)) {
        all.betas[jj] <- all.betas[jj] + epslon  #  increment one of the beta's
        beta2 <- list(psi.betas = all.betas[1:len.psi],
                       gam.betas = all.betas[(len.psi + 1):(len.psi + len.gam)],
                       eps.betas = all.betas[(len.psi + len.gam + 1):length(all.betas)])
        prob_df2 <- GetOccProb(alpha, beta2, Write = FALSE)
        prob_grp2 <- dplyr::group_by(prob_df2, Year)

        if(index == "avg.psi"){
          est2 <- dplyr::summarise(prob_grp2, avg.psi = mean(Prob))$avg.psi
        }

        if(index == "s.lat"){
          est2 <- dplyr::summarise(prob_grp2, s.lat = range.limit(cell.probs = Prob, prob = 0.99, coord = lat, limit = "south"))$s.lat
        }

        if(index == "s.core"){
          est2 <- dplyr::summarise(prob_grp2, s.core = range.limit(cell.probs = Prob, prob = 0.75, coord = lat, limit = "south"))$s.core
        }

        if(index == "n.lat"){
          est2 <- dplyr::summarise(prob_grp2, n.lat = range.limit(cell.probs = Prob, prob = 0.99, coord = lat, limit = "north"))$n.lat
        }

        if(index == "n.core"){
          est2 <- dplyr::summarise(prob_grp2, n.core = range.limit(cell.probs = Prob, prob = 0.75, coord = lat, limit = "north"))$n.core
        }

        if(index == "avg.lat"){
          est2 <- dplyr::summarise(prob_grp2, avg.lat = sum(lat * Prob)/sum(Prob))$avg.lat
        }

        if(index == "avg.lon"){
          est2 <- dplyr::summarise(prob_grp2, avg.lon = sum(lon * Prob)/sum(Prob))$avg.lon
        }

        grad[, jj] <- (est2 - est) / epslon
        all.betas[jj] <- all.betas[jj] - epslon  #  change beta[jj] back to original
      }

    r.vc <- grad %*% betas$vc.mat %*% t(grad)#  real vc matrix (r.vc) = grad * betaVC * grad'
    se <- sqrt(diag(r.vc))
    se
  }

#' range.limit
#'
#' Estimate range limit using cumulative probability method


range.limit <- function(cell.probs, prob, coord, limit){
  xy <- data.frame(x = cell.probs, y = coord)
  x <- dplyr::arrange(xy, y)$x/sum(xy$x)
  y <- dplyr::arrange(xy, y)$y

    if(limit == "south"){
      lim <- predict(smooth.spline(cumsum(x), y, spar=0.1), (1 - prob))$y
    }

    if(limit == "north"){
      lim <- predict(smooth.spline(cumsum(x), y, spar=0.1), prob)$y
    }
  lim
  }

