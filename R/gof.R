#' gof
#'
#' Identify top model, test GOF, save output, and delete output files of all other models
#' @param alpha Alpha code for species of interest
#' @export

gof <- function(alpha){
  mod_opts <- read.csv("inst/model_opts.csv")

  pao <- suppressMessages(RPresence::read.pao(paste0("inst/output/", alpha, "/pres/pres_in.pao")))

  annual_aic <- read.csv(paste0("inst/output/", alpha, "/annual_aic.csv"))

  mods1 <- BBSclim::GetPsiMods()
  aic_psi <- read.csv(paste0("inst/output/", alpha, "/psi_aic.csv"))
  top <- aic_psi$Model_num[1]
  covs <- mods1[[top]]

  mods <- BBSclim::GetGamMods(psi_covs = covs$psi.cov)


  if(is.na(annual_aic$LogLik[1])){
    annual <- FALSE
  }else{
    if(annual_aic$Model[1] == "annual"){
      annual <- TRUE
    }else{
      annual <- FALSE
    }
  }

  start_yr <- as.integer(substr(colnames(pao$unitcov)[grep("tmp", colnames(pao$unitcov))][1], 5, 8))

  year_seq <- seq(from = start_yr, to = start_yr + pao$nseasons - 1)

  aic_tab <- read.csv(paste0("inst/output/", alpha, "/gam_aic.csv"))
  aic_tab <- aic_tab[!is.na(aic_tab$AIC),]
  aic_tab$check <- NA

  #if(nrow(aic_tab) == 0) stop("No models passed overfitting test")

  clim_data <- pao$unitcov

  det_hist <- pao$det.data

  if(mod_opts$Parallel){
    cores <- parallel::detectCores()
    if(!is.null(mod_opts$limit.cores)){
      cores <- min(cores, mod_opts$limit.cores, nrow(aic_tab))
    }

    doParallel::registerDoParallel(cores = cores)

    if(nrow(aic_tab) <= cores){
      mod_check <- foreach::foreach(i=1:nrow(aic_tab), .combine = rbind,
                                    .packages = c("dplyr", "BBSclim")) %dopar%{

                                      ## Read in .out file for current top model
                                      modname <- aic_tab$Model[i]
                                      modnum <- aic_tab$Model_num[i]
                                      mod_out <- scan(paste0("inst/output/", alpha, "/pres/", modname,".out"), what='character', sep='\n', quiet=T)
                                      sim_name <- paste0("sim_", modnum)

                                      ## Extract beta coef estimates and se
                                      jj <- grep('std.error', mod_out)
                                      jj2 <- grep('Variance-Covariance Matrix of Untransformed', mod_out)

                                      betas <- mod_out[(jj+1):(jj2-1)]
                                      coefs <- as.numeric(substr(betas, 41,50))
                                      std.er <- as.numeric(substr(betas, 54,63))

                                      ## Covariates included in the current top model
                                      covs_use <- mods[[modnum]]


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
                                      BBSclim::sim.bbs.ms(alpha = alpha, covs = covs_use, cov_data = clim_data,
                                                          psi.coefs = psi.coefs, th0.coefs = th0.coefs,
                                                          th1.coefs = th1.coefs, gam.coefs = gam.coefs,
                                                          eps.coefs = eps.coefs, p1.coefs = p1.coefs,
                                                          p2.coefs = p2.coefs, pi.coefs = pi.coefs, years = year_seq,
                                                          is.het = mod_opts$het, is.annual = annual, det_hist = det_hist,
                                                          opts = mod_opts, name = sim_name)

                                      ## Create .pao file for simulated data
                                      write_pao(alpha = alpha, sim = TRUE, name = sim_name)

                                      sim_pao <- RPresence::read.pao(paste0("inst/output/", alpha, "/pres/", sim_name, ".pao"))

                                      ## Run Presence on simulated data
                                      initvals <- c(psi.coefs, th0.coefs, th1.coefs, gam.coefs, eps.coefs, p1.coefs)
                                      if(length(p2.coefs) == 0)  initvals <- c(initvals, p2.coefs, pi.coefs)

                                      ## Create design matrices for model
                                      sim_dm <- suppressWarnings(BBSclim::GetDM(pao = sim_pao, cov_list = covs_use, is.het = mod_opts$het, is.annual = annual))


                                      fixedpars <- matrix(rep("eq", pao$nseasons), pao$nseasons, 1)
                                      r1 <- dim(sim_dm$dm1)[1] + dim(sim_dm$dm2)[1] + dim(sim_dm$dm3)[1] + dim(sim_dm$dm4)[1]
                                      rownames(fixedpars) <- (r1 + 1):(r1 + pao$nseasons)

                                      ## Run model
                                      RPresence::write_dm_and_run(paoname = pao$paoname, fixed = fixedpars,
                                                                  dms = sim_dm, model = i,
                                                                  noderived = TRUE, limit.real = TRUE,
                                                                  modname = sim_name)

                                      file.rename(from = paste0("pres_", sim_name, ".out"),
                                                  to = paste0("inst/output/", alpha, "/pres/", sim_name, ".out"))

                                      if(file.exists(paste0("inst/output/", alpha, "/pres/", sim_name, ".out"))){
                                        ## Test whether coefs from simulated data are similar to coefs from top model
                                        gof.pass <- test.presence.gof(modname = sim_name, pao2 = sim_pao, mod = covs_use, alpha = alpha,
                                                                      Psi.coefs = psi.coefs, Gam.coefs = gam.coefs, Eps.coefs = eps.coefs,
                                                                      Psi.se = psi.se, Gam.se = gam.se, Eps.se = eps.se,
                                                                      is.het = mod_opts$het, is.annual = annual)
                                      }else{
                                        gof.pass <- 0
                                      }
                                      gof.pass

                                    }
      aic_tab$check <- mod_check
      top_mod <- aic_tab$Model[min(which(aic_tab$check == 1), na.rm = TRUE)]

      if(!is.na(top_mod)){
        aic_tab <- aic_tab[aic_tab$check == 1,]
        write.csv(aic_tab, paste0("inst/output/", alpha, "/gam_aic_check.csv"))

      # If model passes GOF test, change output file name to "top_mod.out"
      file.rename(from = paste0("inst/output/", alpha, "/pres/", top_mod, ".out"),
                  to = paste0("inst/output/", alpha, "/top_mod.out"))

      # Zip pres folder containing all other .out files
      files2zip <- dir(paste0("inst/output/", alpha, "/pres"), full.names = TRUE)
      zip(zipfile = paste0("inst/output/", alpha, '/presZip'), files = files2zip)

      unlink(paste0("inst/output/", alpha, "/pres"), recursive = TRUE)
      }else{
        write.csv(aic_tab, paste0("inst/output/", alpha, "/gam_aic_check.csv"))
        stop("No model passed GOF")
        # file.rename(from = paste0("inst/output/", alpha, "/pres/psi_model_33.out"),
        #             to = paste0("inst/output/", alpha, "/top_mod.out"))
        #
        # files2zip <- dir(paste0("inst/output/", alpha, "/pres"), full.names = TRUE)
        # zip(zipfile = paste0("inst/output/", alpha, '/presZip'), files = files2zip)
        #
        # unlink(paste0("inst/output/", alpha, "/pres"), recursive = TRUE)
      }
    }else{
      pass <- 0
      cycle <- 0
      while(pass == 0){
        aic_tab2 <- aic_tab[(cycle * cores + 1):min(nrow(aic_tab), (cycle * cores + cores)),]

        mod_check <- foreach::foreach(i=1:nrow(aic_tab2), .combine = rbind,
                                      .packages = c("dplyr", "BBSclim")) %dopar%{

                                        ## Read in .out file for current top model
                                        modname <- aic_tab2$Model[i]
                                        modnum <- aic_tab2$Model_num[i]
                                        mod_out <- scan(paste0("inst/output/", alpha, "/pres/", modname,".out"), what='character', sep='\n', quiet=T)
                                        sim_name <- paste0("sim_", modnum)

                                        ## Extract beta coef estimates and se
                                        jj <- grep('std.error', mod_out)
                                        jj2 <- grep('Variance-Covariance Matrix of Untransformed', mod_out)

                                        betas <- mod_out[(jj+1):(jj2-1)]
                                        coefs <- as.numeric(substr(betas, 41,50))
                                        std.er <- as.numeric(substr(betas, 54,63))

                                        ## Covariates included in the current top model
                                        covs_use <- mods[[modnum]]


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
                                        sim.data	<-  BBSclim::sim.bbs.ms(alpha = alpha, covs = covs_use, cov_data = clim_data,
                                                               psi.coefs = psi.coefs, th0.coefs = th0.coefs,
                                                               th1.coefs = th1.coefs, gam.coefs = gam.coefs,
                                                               eps.coefs = eps.coefs, p1.coefs = p1.coefs,
                                                               p2.coefs = p2.coefs, pi.coefs = pi.coefs, years = year_seq,
                                                               is.het = mod_opts$het, is.annual = annual, det_hist = det_hist, opts = mod_opts,
                                                               name = sim_name)

                                        ## Create .pao file for simulated data
                                        write_pao(alpha = alpha, sim = TRUE, name = sim_name)

                                        sim_pao <- RPresence::read.pao(paste0("inst/output/", alpha, "/pres/", sim_name, ".pao"))

                                        ## Run Presence on simulated data
                                        initvals <- c(psi.coefs, th0.coefs, th1.coefs, gam.coefs, eps.coefs, p1.coefs)
                                        if(length(p2.coefs) == 0)  initvals <- c(initvals, p2.coefs, pi.coefs)

                                        ## Create design matrices for model
                                        sim_dm <- suppressWarnings(BBSclim::GetDM(pao = sim_pao, cov_list = covs_use, is.het = mod_opts$het, is.annual = annual))


                                        fixedpars <- matrix(rep("eq", pao$nseasons), pao$nseasons, 1)
                                        r1 <- dim(sim_dm$dm1)[1] + dim(sim_dm$dm2)[1] + dim(sim_dm$dm3)[1] + dim(sim_dm$dm4)[1]
                                        rownames(fixedpars) <- (r1 + 1):(r1 + pao$nseasons)

                                        ## Run model
                                        RPresence::write_dm_and_run(paoname = pao$paoname, fixed = fixedpars,
                                                                    dms = sim_dm, model = i,
                                                                    noderived = TRUE, limit.real = TRUE,
                                                                    modname = sim_name)

                                        suppressMessages(file.rename(from = paste0("pres_", sim_name, ".out"),
                                                    to = paste0("inst/output/", alpha, "/pres/", sim_name, ".out")))

                                        if(file.exists(paste0("inst/output/", alpha, "/pres/", sim_name, ".out"))){
                                          ## Test whether coefs from simulated data are similar to coefs from top model
                                          gof.pass <- test.presence.gof(modname = sim_name, pao2 = sim_pao, mod = covs_use, alpha = alpha,
                                                                        Psi.coefs = psi.coefs, Gam.coefs = gam.coefs, Eps.coefs = eps.coefs,
                                                                        Psi.se = psi.se, Gam.se = gam.se, Eps.se = eps.se,
                                                                        is.het = mod_opts$het, is.annual = annual)
                                        }else{
                                          gof.pass <- 0
                                        }
                                        gof.pass
                                      }
        ifelse(length(mod_check) <- nrow(aic_tab2)) mod_check <- c(mod_check, rep(0, nrow(aic_tab2) - length(mod_check)))
        aic_tab2$check <- mod_check
        top_mod <- aic_tab2$Model[min(which(aic_tab2$check == 1), na.rm = TRUE)]

        if(!is.na(top_mod) & cycle <= ceiling(nrow(aic_tab)/cores)){
          aic_tab2 <- aic_tab2[aic_tab2$check == 1,]
          write.csv(aic_tab2, paste0("inst/output/", alpha, "/gam_aic_check.csv"))

          # If model passes GOF test, change output file name to "top_mod.out"
          file.rename(from = paste0("inst/output/", alpha, "/pres/", top_mod, ".out"),
                      to = paste0("inst/output/", alpha, "/top_mod.out"))

          # Zip pres folder containing all other .out files
          files2zip <- dir(paste0("inst/output/", alpha, "/pres"), full.names = TRUE)
          zip(zipfile = paste0("inst/output/", alpha, '/presZip'), files = files2zip)

          unlink(paste0("inst/output/", alpha, "/pres"), recursive = TRUE)
          pass <- 1
        }else{
          if(cycle < ceiling(nrow(aic_tab)/cores)){
            cycle <- cycle + 1
          }else{
            aic_tab$check <- 0
            write.csv(aic_tab, paste0("inst/output/", alpha, "/gam_aic_check.csv"), row.names = FALSE)

            file.rename(from = paste0("inst/output/", alpha, "/pres/psi_model_33.out"),
                        to = paste0("inst/output/", alpha, "/top_mod.out"))

            files2zip <- dir(paste0("inst/output/", alpha, "/pres"), full.names = TRUE)
            zip(zipfile = paste0("inst/output/", alpha, '/presZip'), files = files2zip)

            unlink(paste0("inst/output/", alpha, "/pres"), recursive = TRUE)
          }
        }

      }

    }

  }else{

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
    covs_use <- list(psi.cov = c("Lat", "sq_Lat", "Lon", "sq_Lon"),
                     th0.cov = NULL,
                     th1.cov = NULL,
                     gam.cov = mods[[1]]$gam.cov,
                     eps.cov = mods[[1]]$eps.cov,
                     p1.cov = mods[[1]]$p1.cov)



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
    sim.data	<-  BBSclim::sim.bbs.ms(alpha = alpha, covs = covs_use, cov_data = clim_data, name = sim_name,
                            psi.coefs = psi.coefs, th0.coefs = th0.coefs,
                            th1.coefs = th1.coefs, gam.coefs = gam.coefs,
                            eps.coefs = eps.coefs, p1.coefs = p1.coefs,
                            p2.coefs = p2.coefs, pi.coefs = pi.coefs, years = year_seq,
                            is.het = mod_opts$het, is.annual = annual, det_hist = det_hist, opts = mod_opts)

    ## Create .pao file for simulated data
    write_pao(alpha = alpha, sim = TRUE, name = sim_name)

    sim_pao <- RPresence::read.pao(paste0("inst/output/", alpha, "/pres/", sim_name, ".pao"))

    ## Run Presence on simulated data
      initvals <- c(psi.coefs, th0.coefs, th1.coefs, gam.coefs, eps.coefs, p1.coefs)
      if(length(p2.coefs) == 0)  initvals <- c(initvals, p2.coefs, pi.coefs)

      ## Create design matrices for model
      sim_dm <- suppressMessages(BBSclim::GetDM(pao = sim_pao, cov_list = covs_use, is.het = mod_opts$het, is.annual = annual))

      sim_name <- paste0(alpha, "_sim")

      fixedpars <- matrix(rep("eq", pao$nseasons), pao$nseasons, 1)
      r1 <- dim(sim_dm$dm1)[1] + dim(sim_dm$dm2)[1] + dim(sim_dm$dm3)[1] + dim(sim_dm$dm4)[1]
      rownames(fixedpars) <- (r1 + 1):(r1 + pao$nseasons)

      ## Run model
      RPresence::write_dm_and_run(paoname = pao$paoname, fixed = fixedpars,
                                  dms = sim_dm, model = i,
                                  noderived = TRUE, limit.real = TRUE,
                                  modname = sim_name)

      file.rename(from = paste0("pres_", sim_name, ".out"),
                  to = paste0("inst/output/", alpha, "/pres/", sim_name, ".out"))

      ## Test whether coefs from simulated data are similar to coefs from top model
      gof.pass <- test.presence.gof(modname = sim_name, pao2 = sim_pao, mod = covs_use, alpha = alpha,
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
}

#' sim.bbs.ms
#'
#' Simulate BBS data using covariates and coefficients from top model
#' @export

sim.bbs.ms <-  function(alpha, covs, psi.coefs, th0.coefs, th1.coefs,
                        gam.coefs, eps.coefs, p1.coefs, p2.coefs,
                        pi.coefs, years, cov_data, is.annual, is.het, det_hist, opts, name) {

  # initial occupancy

  occ	<- matrix(NA, dim(cov_data)[1], length(years))
  occ.prob <- rep(plogis(psi.coefs[1]), dim(cov_data)[1])  # if no covs

  coord.ind.psi <- grep("Lat|Lon", covs$psi.cov)
  psi.covs <- rep(1, dim(cov_data)[1])
  clim.ind.psi <- seq(1:length(covs$psi.cov))[!(seq(1:length(covs$psi.cov)) %in% coord.ind.psi)]
  clim.covs <- NULL
  if(length(clim.ind.psi) > 0){
    clim.covs <- cov_data[,paste0(covs$psi.cov[clim.ind.psi], "_", as.character(years[1]))]
    psi.covs <- cbind(psi.covs, clim.covs)
  }

  coord.ind.psi <- grep('Lat|Lon', covs$psi.cov)
  coord.covs <- NULL
  if(length(coord.ind.psi) > 0){
    coord.covs <- cov_data[, covs$psi.cov[coord.ind.psi]]
    psi.covs <- cbind(psi.covs, coord.covs)
  }


  if(length(psi.coefs) > 1)   occ.prob	<- plogis(psi.coefs %*% t(data.matrix(psi.covs)))   # problem if psi.ind==logical(0)
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
  coord.ind.gam <- grep('Lat|Lon', covs$gam.cov)
  coord.ind.eps <- grep('Lat|Lon', covs$eps.cov)
  clim.ind.gam <- seq(1:length(covs$gam.cov))[!(seq(1:length(covs$gam.cov)) %in% coord.ind.gam)]
  clim.ind.eps <- seq(1:length(covs$eps.cov))[!(seq(1:length(covs$eps.cov)) %in% coord.ind.eps)]

  for (yr in 2:length(years)) {
    gam.covs <- rep(1, dim(cov_data)[1])
    eps.covs <- rep(1, dim(cov_data)[1])

    clim.covs2 <- NULL
    if(length(clim.ind.gam) > 0){
      clim.covs2 <- cov_data[,paste0(covs$gam.cov[clim.ind.gam], "_", as.character(years[yr]))]
      gam.covs <- cbind(gam.covs, clim.covs2)
    }

    coord.cov2 <- NULL
    if(length(coord.ind.gam) > 0){
      coord.cov2 <- cov_data[, covs$gam.cov[coord.ind.gam]]
      gam.covs <- cbind(gam.covs, coord.cov2)
    }

    clim.covs3 <- NULL
    if(length(clim.ind.eps) > 0){
      clim.covs3 <- cov_data[,paste0(covs$eps.cov[clim.ind.eps], "_", as.character(years[yr]))]
      eps.covs <- cbind(eps.covs, clim.covs3)
    }

    coord.cov3 <- NULL
    if(length(coord.ind.eps) > 0){
      coord.cov3 <- cov_data[, covs$eps.cov[coord.ind.eps]]
      eps.covs <- cbind(eps.covs, coord.cov3)
    }


    gam.prob	<- plogis(gam.coefs %*% t(gam.covs))
    eps.prob	<- plogis(eps.coefs %*% t(eps.covs))
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

  write.csv(history2, file = paste0("inst/output/", alpha, "/pres/", name, "_hist.csv"), row.names = FALSE)
} # end function


#' test.presence.gof
#'
#' Bootstrap GOF test: test if parameters from top model run with simulated data are consistent with original estimates
#' @return 1 if simulated and observed estimates are similar (model not overfit); 0 otherwise

test.presence.gof	<- function(modname, large = 4, pao2, mod, is.annual, is.het, alpha,
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
  # test.quad <- 1.96    					# for most models, require all quads to have correct sign
  # #if(length(grep("gof", modname))==1) test.quad <- 1.96	# if a gof test, then only reject significantly wrong sign
  # if(num.betas[1] > 1) {
  #   quads <- grep("psi.sq",betas2)
  #   psi.check 	<- c(T, coefs2[quads]/std.er2[quads]<test.quad)
  # }
  # if(num.betas[4]>1) {
  #   quads <- grep("gam1.sq",betas2)
  #   gam.check 	<- c(T, coefs2[quads]/std.er2[quads]<test.quad)
  # }
  # if(num.betas[5]>1) {
  #   quads <- grep("eps1.sq",betas2)
  #   eps.check 	<- c(T, coefs2[quads]/std.er2[quads]>(-test.quad))
  # }


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
               #min(psi.check, na.rm=T), min(gam.check, na.rm=T), min(eps.check, na.rm=T),
               mean(abs(z.score)<1.96)>0.8)

  wonky

}	# end test.presence.gof -------------------------------------------------------
