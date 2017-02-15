#' GetDM
#'
#' Create design matrices for Presence
#' @param pao .pao file for the species of interest
#' @param cov_names Vector containing the names of all covariates
#' @param cov_list List containing the index of covariates to include in each model
#' @param psi_list List containing the index of covariates on psi to include in each model
#' @param index Integer indicating which element in cov_list/psi_list to use
#' @param is.het Logical; Should heterogeniety in detection be modeled? (default = TRUE)
#' @param is.annual Logical; Should detection vary annually? (default = TRUE)
#' @param coord.p Logical; Should detection vary by lat/long (default = TRUE)
#' @param coord.th Logical; Should local availabilty vary by lat/long (default = TRUE)
#' @return A list containing 6 design matrices
#' @return   dm1 Design matrix for psi, theta, theta1
#' @return   dm2 Design matrix for gamma
#' @return   dm3 Design matrix for epsilon
#' @return   dm4 Design matrix for p
#' @return   dm5 Design matrix for theta0
#' @return   dm6 Design matrix for p.pi
#' @export

GetDM	<- function(pao,cov_list,
                  is.het, is.annual, coord.p = TRUE, coord.th = TRUE) {
  opts <- read.csv("inst/global_opts.csv")

  start_yr <- opts$start_yr

  n_surv <- pao$nsurveys  # total num of surveys
  n_seas <- pao$nseasons	# num seasons
  n_count <- pao$nsurveyseason[1]

  years <- seq(from = start_yr, to = start_yr + n_seas - 1)

  num.betas	<- c(1 + length(cov_list$psi.cov),
                 1 + length(cov_list$th0.cov),
                 1 + length(cov_list$th1.cov),
                 1 + length(cov_list$gam.cov),
                 1 + length(cov_list$eps.cov),
                 1 + length(cov_list$p1.cov))


  # 1st design matrix, 1 row for psi, K rows for theta0, K rows for theta1
  dm1 <- matrix('0', n_surv * 2 + 1, sum(num.betas[1:3]))

  # psi
  psi.dm <- c(1, paste(cov_list$psi.cov, as.character(years[1]), sep = "_"))
  if(length(cov_list$psi.cov) == 0) psi.dm <- 1   # if no covs, just need intercept
  dm1[1, 1:num.betas[1]] <- psi.dm

  # th0
  if(!is.null(cov_list$th0.cov)){
    coord.ind.th0 <- grep('Lat|Lon', cov_list$th0.cov)
    coord.mat.th0 <- matrix(rep(cov_list$th0.cov[coord.ind.th0], n_surv), nrow = n_surv, byrow = TRUE)
    if(length(coord.ind.th0) == 0)	coord.mat.th0 <- NULL


    non.coord.th0 <- grep('Lat|Lon', cov_list$th0.cov, invert=T)
    cov.mat.th0 <- matrix(rep(paste(cov_list$th0.cov[non.coord.th0], as.character(years[1]), sep = "_"), n_count), nrow = n_count, byrow = TRUE)
    for (ii in 2:length(years)) {
      cov.mat.th0 <- rbind(cov.mat.th0, matrix(rep(paste(cov_list$th0.cov[non.coord.th0], as.character(years[ii]), sep = "_"), n_count), nrow = n_count, byrow = TRUE))
    }
    if(length(non.coord.th0)==0)	cov.mat.th0 <- NULL

    th0.intercept <- rep(1, n_surv)
    th0.dm <- cbind(th0.intercept, coord.mat.th0, cov.mat.th0)
    if(length(cov_list$th0.cov) == 0) th0.dm  <- rep(1, n_surv)
  }else{
    th0.dm <- rep(1, n_surv)
  }
  dm1[2:(n_surv + 1),(1 + num.betas[1]):sum(num.betas[1:2])] <- th0.dm

  # th1
  if(!is.null(cov_list$th1.cov)){
    coord.ind.th1 <- grep('Lat|Lon', cov_list$th1.cov)
    coord.mat.th1 <- matrix(rep(cov_list$th1.cov[coord.ind.th1], n_surv), nrow = n_surv, byrow = TRUE)
    if(length(coord.ind.th1) == 0)	coord.mat.th1 <- NULL


    non.coord.th1 <- grep('Lat|Lon', cov_list$th1.cov, invert=T)
    cov.mat.th1 <- matrix(rep(paste(cov_list$th1.cov[non.coord.th1], as.character(years[1]), sep = "_"), n_count), nrow = n_count, byrow = TRUE)
    for (ii in 2:length(years)) {
      cov.mat.th1 <- rbind(cov.mat.th1, matrix(rep(paste(cov_list$th1.cov[non.coord.th1], as.character(years[ii]), sep = "_"), n_count), nrow = n_count, byrow = TRUE))
    }
    if(length(non.coord.th1)==0)	cov.mat.th1 <- NULL

    th1.intercept <- rep(1, n_surv)
    th1.dm <- cbind(th1.intercept, coord.mat.th1, cov.mat.th1)
    if(length(cov_list$th1.cov) == 0) th1.dm  <- rep(1, n_surv)
  }else{
    th1.dm  <- rep(1, n_surv)
  }
  dm1[(n_surv + 2):(2 * n_surv + 1), (1 + sum(num.betas[1:2])):sum(num.betas[1:3])] <- th1.dm

  rownames(dm1) <- c('psi', paste0('th0(',1:n_surv,')'), paste0('th1(',1:n_surv,')'))
  colnames(dm1) <- paste0('a',1:sum(num.betas[1:3]))

  # gam
  gam.dm <- c(1, paste(cov_list$gam.cov, as.character(years[2]), sep = "_"))
  for (ii in 3:length(years)) {
    gam.dm <- c(gam.dm, c(1, paste(cov_list$gam.cov, as.character(years[ii]), sep = "_")))
  }
  if(length(cov_list$gam.cov) == 0) gam.dm  <- rep(1, n_seas - 1)
  dm2 <- matrix(gam.dm, n_seas - 1, num.betas[4], byrow=T,
                dimnames = list(paste0('gam', 1:(n_seas-1)), paste0('b', 1:num.betas[4])))

  # eps
  eps.dm <- c(1, paste(cov_list$eps.cov, as.character(years[2]), sep = "_"))
  for (ii in 3:length(years)) {
    eps.dm <- c(eps.dm, c(1, paste(cov_list$eps.cov, as.character(years[ii]), sep = "_")))
  }
  if(length(cov_list$eps.cov) == 0) eps.dm  <- rep(1, n_seas - 1)
  dm3 <- matrix(eps.dm, n_seas - 1, num.betas[5], byrow = T,
                dimnames = list(paste0('eps',1:(n_seas - 1)),paste0('c', 1:num.betas[5])))

  # p
  stop.ind   <- grep('Stop', cov_list$p1.cov)
  scale.stop <- paste0(cov_list$p1.cov[stop.ind], 1)
  for (ii in 2:n_count) {
    scale.stop <- c(scale.stop, paste0(cov_list$p1.cov[stop.ind], ii))
  }
  scale.stop <- matrix(scale.stop, nrow = n_surv, ncol = length(stop.ind), byrow = T)
  if(length(stop.ind) == 0)	scale.stop <- NULL

  time1	<- c(rep(0, n_count),
             rep(c(rep(1, n_count),
                   rep(0, n_surv)), n_seas - 1))
  time2	<- matrix(time1, nrow = n_surv, ncol = n_seas - 1)  # warning is ok
  if(is.annual == FALSE)	time2 <- NULL

  coord.ind.p <- grep('Lat|Lon', cov_list$p1.cov)
  coord.mat.p <- matrix(rep(cov_list$p1.cov[coord.ind.p], n_surv), nrow = n_surv, byrow = TRUE)
  if(length(coord.ind.p) == 0)	coord.mat.p <- NULL

  cov.mat.p <- NULL
  non.stop.p <- grep('Stop|Lat|Lon', cov_list$p1.cov, invert=T)
  for (ii in 1:length(years)) {
    for (jj in 1:(n_surv/n_seas)) {
      cov.mat.p <- rbind(cov.mat.p, paste0(cov_list$p1.cov[non.stop.p], "_", as.character(years[ii])))
    }
  }
  if(length(non.stop.p) == 0)	cov.mat.p <- NULL

  p1.intercept <- rep(1, n_surv)
  p1.dm <- cbind(p1.intercept, time2, scale.stop, coord.mat.p, cov.mat.p)

  zeros <- matrix(0, nrow = dim(p1.dm)[1], ncol = dim(p1.dm)[2])
  dm4 <- rbind(cbind(p1.dm, zeros), cbind(zeros, p1.dm))
  rownames(dm4) <- c(paste0('p1(', 1:n_surv,')'),paste0('p2(', 1:n_surv,')'))     # for het

  if(is.het == FALSE) dm4 <- p1.dm
  if(is.het == FALSE) rownames(dm4) <- c(paste0('p1(',1:n_surv,')'))     	# could get rid of parentheses
  colnames(dm4) <- paste0('d', 1:dim(dm4)[2])

  # theta.pi
  dm5 <- matrix(0, n_seas, 1, dimnames = list(paste0('thta0pi', 1:n_seas), NULL)) 		# NOTE THE ZERO when fixing, note no colname
  # p.pi
  dm6 <- matrix(1, n_seas, 1, dimnames = list(paste0('pi', 1:n_seas),'f1')) 	                # mixture on detection

  dm_list <- list(dm1 = dm1, dm2 = dm2, dm3 = dm3, dm4 = dm4, dm5 = dm5, dm6 = dm6)

  dm_list
}


