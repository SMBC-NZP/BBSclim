#' gof
#'
#' Identify top model, test GOF, save output, and delete output files of all other models
#' @param aic_tab AIC table from RunGamMods
#' @param mods List of models
#' @param cov_data Data frame containing the covariates used to fit the full model
#' @param year_seq Vector containing the years included in the analysis
#' @export
#'

gof <- function(alpha, pao){
  mod_opts <- read.csv("inst/model_opts.csv")
  glob_opts <- read.csv("inst/global_opts.csv")

  annual_aic <- read.csv(paste0("inst/output/", alpha, "/annual_aic.csv"))

  mods1 <- GetGamMods()
  aic_gam <- read.csv(paste0("inst/output/", alpha, "/gam_aic.csv"))
  top <- aic_gam$Model_num[1]
  covs <- mods1[[top]]
  covs.ll <- list(gam_covs = covs$gam.cov, eps_covs = covs$eps.cov)

  mods <- GetPsiMods(covs = covs.ll)

  if(annual_aic$Model[1] == "annual"){
    annual <- TRUE
  }else{
    annual <- FALSE
  }

  year_seq <- seq(from = glob_opts$start_yr, to = glob_opts$end_yr)

  aic_tab <- read.csv(paste0("inst/output/", alpha, "/psi_aic.csv"))

  clim_data <- pao$unitcov

  det_hist <- pao$det.data

  gof.pass <- 0
  model.num <- 0

  while (gof.pass==0) {
    ## Read in .out file for current top model
    modname <- aic_tab$Model[model.num + 1]
    modnum <- aic_tab$Model_num[model.num + 1]
    mod_out <- scan(paste0("inst/output/", alpha, "/pres/", modname,".out"), what='character', sep='\n', quiet=T)

    ## Extract beta coef estimates and se
    jj <- grep('std.error', mod_out)
    jj2 <- grep('Variance-Covariance Matrix of Untransformed', mod_out)

    betas <- mod_out[(jj+1):(jj2-1)]
    coefs <- as.numeric(substr(betas, 41,50))
    std.er <- as.numeric(substr(betas, 54,63))

    ## Covariates included in the current top model
    if(modnum == 31){
      covs_use <- list(psi.cov = c("tmp", "sq_tmp", "dtr", "sq_dtr", "Twet", "sq_Twet",
                                   "Prec", "sq_Prec", "Pwarm", "sq_Pwarm"),
                       th0.cov = mods[[1]]$th0.cov,
                       th1.cov = mods[[1]]$th1.cov,
                       gam.cov = mods[[1]]$gam.cov,
                       eps.cov = mods[[1]]$eps.cov,
                       p1.cov = mods[[1]]$p1.cov)
    }else{
      covs_use <- mods[[modnum]]
    }


    ## Beta coefs
    psi.coefs <- coefs[grep('psi',betas)]
    th0.coefs <- coefs[grep('th0',betas)]
    th1.coefs <- coefs[grep('th1',betas)]
    gam.coefs <- coefs[grep('gam',betas)]
    eps.coefs <- coefs[grep('eps',betas)]
    p1.coefs  <- coefs[grep('\\sp1',betas)]   # require space (\\s) before p1 to avoid tmp1968, etc
    p2.coefs  <- coefs[grep('\\sp2',betas)]
    if(length(p2.coefs) == 0)		p2.coefs <- p1.coefs	# if no het
    pi.coefs <- coefs[grep('pi1',betas)]
    if(length(pi.coefs)==0)		pi.coefs <- 0		# if no het

    ## se
    psi.se <- std.er[grep('psi', betas)]
    gam.se <- std.er[grep('gam', betas)]
    eps.se <- std.er[grep('eps', betas)]

    ## Simulate new detection history from top model
    sim.data	<- sim.bbs.ms(covs = covs_use, cov_data = clim_data,
                            psi.coefs=psi.coefs, th0.coefs=th0.coefs,
                            th1.coefs=th1.coefs, gam.coefs=gam.coefs,
                            eps.coefs=eps.coefs, p1.coefs=p1.coefs,
                            p2.coefs=p2.coefs, pi.coefs=pi.coefs, years=year_seq,
                            is.het = mod_opts$het, is.annual = annual, det_hist, opts = glob_opts)

    ## Create .pao file for simulated data
    write_pao(alpha = alpha, sim = TRUE)

    sim_pao <- RPresence::read.pao(paste0("inst/output/", alpha, "/pres/sim.pao"))

    ## Run Presence on simulated data
      initvals <- c(psi.coefs, th0.coefs, th1.coefs, gam.coefs, eps.coefs, p1.coefs)
      if(length(p2.coefs) == 0)  initvals <- c(initvals, p2.coefs, pi.coefs)

      ## Create design matrices for model
      sim_dm <- GetDM(pao = sim_pao, cov_list = covs_use, is.het = mod_opts$het, is.annual = annual)

      sim_name <- paste0(alpha, "_sim")

      ## Run model
      write_dm_and_run2(pao = sim_pao, cov_list = covs_use, is.het = mod_opts$het, dm_list = sim_dm,
                        modname = sim_name, fixed = TRUE,
                        inits = TRUE, maxfn = '35000 vc lmt=5', alpha = alpha)

      ## Test whether coefs from simulated data are similar to coefs from top model
      gof.pass <- test.presence.gof(modname = sim_name, pao2 = sim_pao, mod = covs_use,
                                    Psi.coefs = psi.coefs, Gam.coefs = gam.coefs, Eps.coefs = eps.coefs,
                                    Psi.se = psi.se, Gam.se = gam.se, Eps.se = eps.se,
                                    is.het = mod_opts$het, is.annual = annual)

      if(gof.pass){
        # If model passes GOF test, change output file name to "top_mod.out"
        file.rename(from = paste0("inst/output/", alpha, "/pres/", modname, ".out"),
                    to = paste0("inst/output/", alpha, "/top_mod.out"))

        # Zip pres folder containing all other .out files
        files2zip <- dir(paste0("inst/output/", alpha, "/pres"), full.names = TRUE)
        zip(zipfile = paste0("inst/output/", alpha, '/presZip'), files = files2zip)

        unlink(paste0("inst/output/", alpha, "/pres"), recursive = TRUE)
      }else{
        # If model does not pass GOF test, delete sim.out file so Presence can run again
        temp.files <- list.files(paste0("inst/output/", alpha, "/pres"))
        sim.file <- temp.files[grep("sim", temp.files)]
        file.remove(paste0("inst/output/", alpha, "/pres/", sim.file))
      }
    }  # end while loop

}

#' sim.bbs.ms
#'
#' Simulate BBS data using covariates and coefficients from top model

sim.bbs.ms <-  function(covs, psi.coefs, th0.coefs, th1.coefs,
                        gam.coefs, eps.coefs, p1.coefs, p2.coefs,
                        pi.coefs, years, cov_data, is.annual, is.het, det_hist, opts) {

  # initial occupancy

  occ		<- matrix(NA, dim(cov_data)[1], length(years))
  occ.prob <- rep(plogis(psi.coefs[1]), dim(cov_data)[1])  # if no covs
  if(length(psi.coefs)>1)   occ.prob	<- plogis(psi.coefs %*% t(cbind(1,cov_data[,paste0(covs$psi.cov, "_", as.character(years[1]))])))   # problem if psi.ind==logical(0)
  occ[,1]		<- rbinom(occ.prob, 1, occ.prob)


  # the scale.stop is the effect of stop number (i.e., time of day)
  if(opts$tenstops){tot.stops <- 5}else{tot.stops <- 50}
  scale.stop <- cbind(scale(1:tot.stops), scale(1:tot.stops)^2)
  stop.ind <- grep('Stop', covs$p1.cov)
  coord.p <- grep('Lat|Lon', covs$p1.cov)  # which covariates are for lat/long
  coord.th0 <- grep('Lat|Lon', covs$th0.cov)  # which covariates are for lat/long
  coord.th1 <- grep('Lat|Lon', covs$th1.cov)  # which covariates are for lat/long
  non.coord.p <- grep('Stop|Lat|Lon', covs$p.cov, invert=T)  # which covariates are for climate
  non.coord.th0 <- grep('Lat|Lon', covs$th0.cov, invert=T)  # which covariates are for climate
  non.coord.th1 <- grep('Lat|Lon', covs$th1.cov, invert=T)  # which covariates are for climate


  # occupancy in following years

  for (yr in 2:length(years)) {
    gam.prob <- rep(plogis(gam.coefs[1]), dim(cov_data)[1]) # if no covs
    if(length(covs$gam.cov)>0) gam.prob	<- plogis(gam.coefs %*% t(cbind(1,cov_data[,paste0(covs$gam.cov,"_",as.character(years[yr]))])))
    eps.prob <- rep(plogis(eps.coefs[1]), dim(cov_data)[1]) # if no covs
    if(length(covs$eps.cov)>0) eps.prob	<- plogis(eps.coefs %*% t(cbind(1,cov_data[,paste0(covs$eps.cov,"_",as.character(years[yr]))])))
    occ[,yr]	<- rbinom(gam.prob, 1, occ[,yr-1]*(1 - eps.prob) + (1-occ[,yr-1])*gam.prob)
  }

  # occupancy and detection at each stop

  stops		<- array(NA, c(dim(cov_data)[1], length(years), tot.stops + 1))
  history		<- array(NA, c(dim(cov_data)[1], length(years), tot.stops))

  for (yr in seq_along(years)) {
    # th0 is prob of presence if NOT PRESENT at previous stop, th1 if PRESENT at previous
    th0 <- plogis(th0.coefs %*% t(cbind(1, cov_data[,covs$th0.cov]))) # if no climate covs
    if(length(non.coord.th0)>0) th0	<- plogis(th0.coefs %*% t(cbind(1, cov_data[,covs$th0.cov[coord.th0]], cov_data[,paste0(covs$th0.cov[non.coord.th0],"_", as.character(years[yr]))])))
    th1 <- plogis(th1.coefs %*% t(cbind(1,cov_data[,covs$th1.cov]))) # if no climate covs
    if(length(non.coord.th0)>0) th1	<- plogis(th1.coefs %*% t(cbind(1,cov_data[,covs$th1.cov], cov_data[,paste0(covs$th1.cov[non.coord.th1],"_",as.character(years[yr]))])))
    # Equilibrium probability = probability that stop 0 is available
    stop.prob	<- th0/(th0+1-th1)
    stops[, yr, 1]	<- rbinom(stop.prob, 1, stop.prob*occ[,yr])

    p.mix		<- rbinom(stop.prob, 1, plogis(pi.coefs)) # 1 means use p1, 0 means use p2

    for (ss in 2:(tot.stops + 1)) {

      stops[,yr,ss]	<- rbinom(stop.prob, 1, stops[,yr,ss-1]*th1*occ[,yr] + (1-stops[,yr,ss-1])*th0*occ[,yr])

      ones <- matrix(0, nrow = nrow(cov_data), ncol = length(years))
      ones[, yr] <- 1
      if(!is.annual) ones <- rep(1, nrow(cov_data))
      stop.mat  <- matrix(scale.stop[ss - 1, stop.ind], nrow=nrow(cov_data), ncol=length(stop.ind), byrow=T)
      if(length(stop.ind) == 0)  stop.mat <- NULL
      if(length(coord.p) > 0)   coord.mat  <- as.matrix(cov_data[ , covs$p1.cov[coord.p]])
      if(length(coord.p) == 0)  coord.mat <- NULL
      if(length(non.coord.p) > 0)   clim.mat  <- as.matrix(cov_data[ , paste0(covs$p1.cov[non.coord.p], "_", as.character(years[yr]))])
      if(length(non.coord.p) == 0)  clim.mat <- NULL
      values    <- cbind(ones, stop.mat, coord.mat, clim.mat)

      p1      <- rep(plogis(p1.coefs[1]), dim(cov_data)[1]) # if no covs
      #if(length(p1.ind)>0)    p1		<- plogis(matrix(p1.coefs,nrow=1) %*% t(cbind(1, cov_data[,paste0(covs$p1.ind],as.character(years[yr]))])))
      if(length(covs$p1.cov)>0)    p1		<- plogis(matrix(p1.coefs,nrow=1) %*% t(values))
      p2      <- rep(plogis(p2.coefs[1]), dim(cov_data)[1]) # if no covs
      #if(length(p2.ind)>0)    p2		<- plogis(matrix(p2.coefs,nrow=1) %*% t(cbind(1, cov_data[,paste0(covs$p1.ind],as.character(years[yr]))])))
      if(length(covs$p1.cov)>0)    p2		<- plogis(matrix(p2.coefs,nrow=1) %*% t(values))
      det.prob	<- p1*p.mix + p2*(1-p.mix)  # this selects p1 or p2
      history[,yr,ss-1]	<- rbinom(stop.prob, 1, det.prob*stops[,yr,ss])
    }
  }


  # format and insert NAs to match original data  #####  DET.HIST NEEDS TO BE RIGHT FORMAT

  history2 <- NULL
  for (ii in seq_along(years)) {
    history[which(is.na(det_hist[,(ii-1)*tot.stops+1])), ii, ] <- NA
    #history[which(is.na(det.hist[,(ii-1)*tot.stops+2])==TRUE),ii,] <- NA
    history2 <- cbind(history2, history[, ii, ])
  }

  write.csv(history2, file = paste0("inst/output/", alpha, "/pres/sim_hist.csv"), row.names = FALSE)
} # end function


#' test.presence.gof
#'
#' Bootstrap GOF test: test if parameters from top model run with simulated data are consistent with original estimates
#' @return 1 if simulated and observed estimates are similar (model not overfit); 0 otherwise

test.presence.gof	<- function(modname, large = 4, pao2, mod, is.annual, is.het,
                              Psi.coefs, Gam.coefs, Eps.coefs, Psi.se, Gam.se, Eps.se) {

  ###	read the output from simulated model
  mod_out2 <- scan(paste0("inst/output/", alpha, "/pres/", modname ,".out"), what='character', sep='\n', quiet=T)


  ## Extract beta coefs and se
  jjx <- grep('std.error', mod_out2)
  jjx2 <- grep('Variance-Covariance Matrix of Untransformed', mod_out2)
  betas2 <- mod_out2[(jjx + 1):(jjx2 - 1)]
  coefs2 <- as.numeric(substr(betas2, 41,50))
  std.er2 <- as.numeric(substr(betas2, 54,63))

  num.betas <- plyr::ldply(mod, function(x) length(x))$V1

  if(is.annual) num.betas[6] <- num.betas[6] + pao2$nseasons

  # check for quadratic terms with significant wrong sign
  test.quad <- 1.96    					# for most models, require all quads to have correct sign
  #if(length(grep("gof", modname))==1) test.quad <- 1.96	# if a gof test, then only reject significantly wrong sign
  if(num.betas[1] > 1) {
    quads <- grep("psi.sq",betas2)
    psi.check 	<- c(T, coefs2[quads]/std.er2[quads]<test.quad)
  }
  if(num.betas[4]>1) {
    quads <- grep("gam1.sq",betas2)
    gam.check 	<- c(T, coefs2[quads]/std.er2[quads]<test.quad)
  }
  if(num.betas[5]>1) {
    quads <- grep("eps1.sq",betas2)
    eps.check 	<- c(T, coefs2[quads]/std.er2[quads]>(-test.quad))
  }


  #### compare bootstrap coefs to original
  psi.coefs2 <- coefs2[grep('psi',betas2)]
  gam.coefs2 <- coefs2[grep('gam',betas2)]
  eps.coefs2 <- coefs2[grep('eps',betas2)]

  psi.se2 <- std.er2[grep('psi', betas2)]
  gam.se2 <- std.er2[grep('gam', betas2)]
  eps.se2 <- std.er2[grep('eps', betas2)]

  z.score <- (c(Psi.coefs, Gam.coefs, Eps.coefs) - c(psi.coefs2, gam.coefs2, eps.coefs2))/sqrt(c(Psi.se, Gam.se, Eps.se)^2 + c(psi.se2, gam.se2, eps.se2)^2)
  #if(mean(abs(z.score)<1.96)>0.8)  gof.pass <- 1
  #if(model.num==nrow(t))	   gof.pass <- 1


  wonky <- min(min(!is.na(std.er2[1:sum(num.betas[1:5])]) > 0), min(abs(!is.na(std.er2))) < large,
               mean(psi.check, na.rm=T), mean(gam.check, na.rm=T), mean(eps.check, na.rm=T),
               mean(abs(z.score)<1.96)>0.8)

  wonky

}	# end test.presence.gof -------------------------------------------------------
