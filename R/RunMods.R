#' RunPsiMods
#'
#' Runs full model set for 127 psi models, evalute each, write and return AIC table, delete output file (optional)
#' @param alpha alpha code for species of interest
#' @export

RunPsiMods <- function(alpha){
    opts <- read.csv("inst/model_opts.csv")

    pao <- suppressMessages(RPresence::read.pao(paste0("inst/output/", alpha, "/pres/pres_in_psi.pao")))


    mods <- BBSclim::GetPsiMods()


    if(opts$Parallel){
      cores <- parallel::detectCores()
      if(!is.null(opts$limit.cores)){
        cores <- min(cores, opts$limit.cores)
      }

      doParallel::registerDoParallel(cores = cores)

      if(opts$psi.test){
        mods <- mods[1:cores]
      }

      aic_table <- foreach::foreach(i = 1:length(mods), .combine = rbind,
                                    .packages = c("dplyr", "BBSclim")) %dopar%{
                                      modname <- paste0('psi_model_', i)

                                      ## Create design matrices for model i
                                      spp_dm <- suppressWarnings(BBSclim::GetDM(pao = pao, psi = TRUE, cov_list = mods[[i]], is.annual = FALSE, is.het = FALSE))

                                      fixedpars <- matrix(rep("eq", pao$nseasons), pao$nseasons, 1)
                                      r1 <- dim(spp_dm$dm1)[1] + dim(spp_dm$dm4)[1]
                                      rownames(fixedpars) <- (r1 + 1):(r1 + pao$nseasons)

                                      ## Run model
                                      RPresence::write_dm_and_run(paoname = pao$paoname, fixed = fixedpars,
                                                                  dms = spp_dm, model = i,
                                                                  noderived = TRUE, limit.real = TRUE,
                                                                  modname = modname)

                                      suppressMessages(file.rename(from = paste0("pres_", modname, ".out"),
                                                  to = paste0("inst/output/", alpha, "/pres/", modname, ".out")))

                                      ## Read output file
                                      a <- scan(paste0('inst/output/', alpha, "/pres/", modname, ".out"), what='c',sep='\n',quiet=TRUE)

                                      ## Evaluate model (if model converges, will equal TRUE)
                                      check <- BBSclim::mod_eval(pres_out = a, pao2 = pao, mod = mods[[i]], psi = TRUE, strict.psi = TRUE, strict.gam = FALSE, is.annual = FALSE, is.het = opts$het)

                                      if(check == FALSE){ # If model does not converge, save NA in AIC table
                                        aic_temp <- dplyr::data_frame(Model = modname, Model_num = i, LogLik = NA, nParam = NA,
                                                                      AIC = NA)

                                      }else{ # If model does converge, save results to AIC table
                                        ## Extract log likelihood
                                        j <- grep('-2log', a)
                                        loglike <- as.numeric(unlist(strsplit(a[j],'=',2))[2])

                                        ## Extract AIC
                                        j <- grep('AIC', a)
                                        aic <- as.numeric(unlist(strsplit(a[j],'=',2))[2])

                                        ## Number of parameters
                                        j <- grep('of par', a)
                                        n  <- as.numeric(unlist(strsplit(a[j],'=',2))[2])

                                        aic_temp <- dplyr::data_frame(Model = modname, Model_num = i, LogLik = loglike, nParam = n,
                                                                      AIC = aic)
                                      }
                                      aic_temp
                                    }

      doParallel::stopImplicitCluster()
      aic_table
    }else{
      if(opts$psi.test){
        mods <- mods[length(mods) - 1:length(mods)]
      }

      for(i in 1:length(mods)){

        modname <- paste0('psi_model_', i)

        ## Create design matrices for model i
        spp_dm <- suppressWarnings(BBSclim::GetDM(pao = pao, psi = TRUE, cov_list = mods[[i]], is.annual = annual, is.het = opts$het))

        fixedpars <- matrix(rep("eq", pao$nseasons), pao$nseasons, 1)
        r1 <- dim(spp_dm$dm1)[1] + dim(spp_dm$dm4)[1]
        rownames(fixedpars) <- (r1 + 1):(r1 + pao$nseasons)

        ## Run model
        RPresence::write_dm_and_run(paoname = pao$paoname, fixed = fixedpars,
                                    dms = spp_dm, model = i,
                                    noderived = TRUE, limit.real = TRUE,
                                    modname = modname)

        suppressMessages(file.rename(from = paste0("pres_", modname, ".out"),
                    to = paste0("inst/output/", alpha, "/pres/", modname, ".out")))


        ## Read output file
        a <- scan(paste0('inst/output/', alpha, "/pres/", modname, ".out"), what='c', sep='\n', quiet=TRUE)

      ## Evaluate model (if model converges, will equal TRUE)
      check <- mod_eval(pres_out = a, pao2 = pao, mod = mods[[i]], strict.psi = TRUE, strict.gam = TRUE, is.annual = annual, is.het = opts$het)


        if(check == FALSE){ # If model does not converge, save NA in AIC table
          aic_temp <- dplyr::data_frame(Model = modname, Model_num = i, LogLik = NA, nParam = NA,
                                        AIC = NA)

        }else{ # If model does converge, save results to AIC table
          ## Extract log likelihood
          j <- grep('-2log', a)
          loglike <- as.numeric(unlist(strsplit(a[j],'=',2))[2])

          ## Extract AIC
          j <- grep('AIC', a)
          aic <- as.numeric(unlist(strsplit(a[j],'=',2))[2])

          ## Number of parameters
          j <- grep('of par', a)
          n  <- as.numeric(unlist(strsplit(a[j],'=',2))[2])

          aic_temp <- dplyr::data_frame(Model = modname, Model_num = i, LogLik = loglike, nParam = n,
                                        AIC = aic)
        }

        if(i == 1){
          aic_table <- aic_temp
        }else{
          aic_table <- dplyr::bind_rows(aic_table, aic_temp)
        }

       }
      }


    ## Add delta AIC column and sort by delta AIC
    aic_table <- dplyr::mutate(aic_table, delta_AIC = AIC - min(AIC, na.rm = TRUE))
    aic_table <- dplyr::arrange(aic_table, delta_AIC)

    ## Write AIC table
    write.csv(aic_table, file = paste0("inst/output/", alpha, "/psi_aic.csv"), row.names = FALSE)
}


#' RunGamMods
#'
#' Runs full model set for 961 gamma/epsilon models, evaluate each, write and return AIC table, delete .out files (optional)
#' @param alpha four letter species code
#' @param pao pao object
#' @export


RunGamMods <- function(alpha){
  opts <- read.csv("inst/model_opts.csv")

  annual_aic <- read.csv(paste0("inst/output/", alpha, "/annual_aic.csv"))

  psi_aic <- read.csv(paste0("inst/output/", alpha, "/psi_aic.csv"))
  top_covs <- GetPsiMods()[[psi_aic$Model_num[1]]]

  mods <- GetGamMods(psi_covs = top_covs$psi.cov)

  pao <- suppressMessages(RPresence::read.pao(paste0("inst/output/", alpha, "/pres/pres_in.pao")))

  if(is.na(annual_aic$LogLik[1])){
    annual <- FALSE
  }else{
    if(annual_aic$Model[1] == "annual"){
      annual <- TRUE
    }else{
      annual <- FALSE
    }
  }

  if(opts$Parallel){

    cores <- parallel::detectCores()
    if(!is.null(opts$limit.cores)){
      cores <- min(cores, opts$limit.cores)
    }

    doParallel::registerDoParallel(cores = cores)

    if(opts$gam.test){
      mods <- mods[1:cores]
    }

    aic_table <- foreach::foreach(i=1:length(mods), .combine = rbind,
                                  .packages = c("dplyr", "BBSclim")) %dopar%{

                                    modname <- paste0('gam_model_', i)

                                    ## Create design matrices for model i
                                    spp_dm <- suppressWarnings(BBSclim::GetDM(pao = pao, cov_list = mods[[i]], is.annual = annual, is.het = opts$het))

                                    fixedpars <- matrix(rep("eq", pao$nseasons), pao$nseasons, 1)
                                    r1 <- dim(spp_dm$dm1)[1] + dim(spp_dm$dm2)[1] + dim(spp_dm$dm3)[1] + dim(spp_dm$dm4)[1]
                                    rownames(fixedpars) <- (r1 + 1):(r1 + pao$nseasons)

                                    ## Run model
                                    RPresence::write_dm_and_run(paoname = pao$paoname, fixed = fixedpars,
                                                                dms = spp_dm, model = i,
                                                                noderived = TRUE, limit.real = TRUE,
                                                                modname = modname)

                                    suppressMessages(file.rename(from = paste0("pres_", modname, ".out"),
                                                to = paste0("inst/output/", alpha, "/pres/", modname, ".out")))

                                    ## Read output file
                                    a <- scan(paste0('inst/output/', alpha, "/pres/", modname, ".out"), what='c',sep='\n',quiet=TRUE)

                                    ## Evaluate model (if model converges, will equal TRUE)
                                    check <- BBSclim::mod_eval(pres_out = a, pao2 = pao, mod = mods[[i]], strict.psi = TRUE, strict.gam = TRUE,
                                                               is.het = opts$het, is.annual = annual)

                                    if(check == FALSE){ # If model does not converge, save NA in AIC table
                                      aic_temp <- dplyr::data_frame(Model = modname, Model_num = i, LogLik = NA, nParam = NA,
                                                                    AIC = NA)

                                    }else{ # If model does converge, save results to AIC table
                                      ## Extract log likelihood
                                      j <- grep('-2log', a)
                                      loglike <- as.numeric(unlist(strsplit(a[j],'=',2))[2])

                                      ## Extract AIC
                                      j <- grep('AIC', a)
                                      aic <- as.numeric(unlist(strsplit(a[j],'=',2))[2])

                                      ## Number of parameters
                                      j <- grep('of par', a)
                                      n  <- as.numeric(unlist(strsplit(a[j],'=',2))[2])

                                      aic_temp <- dplyr::data_frame(Model = modname, Model_num = i, LogLik = loglike, nParam = n,
                                                                    AIC = aic)
                                    }
                                    aic_temp
                                  }

    doParallel::stopImplicitCluster()
    aic_table
  }else{
    if(opts$gam.test){
      mods <- mods[1:cores]
    }

    for(i in 1:length(mods)){

      modname <- paste0('gam_model_', i)

      ## Create design matrices for model i
      spp_dm <- suppressWarnings(BBSclim::GetDM(pao = pao, cov_list = mods[[i]], is.het = opts$het, is.annual = annual))

      fixedpars <- matrix(rep("eq", pao$nseasons), pao$nseasons, 1)
      r1 <- dim(spp_dm$dm1)[1] + dim(spp_dm$dm2)[1] + dim(spp_dm$dm3)[1] + dim(spp_dm$dm4)[1]
      rownames(fixedpars) <- (r1 + 1):(r1 + pao$nseasons)

      ## Run model
      RPresence::write_dm_and_run(paoname = pao$paoname, fixed = fixedpars,
                                  dms = spp_dm, model = i,
                                  noderived = TRUE, limit.real = TRUE,
                                  modname = modname)

      suppressMessages(file.rename(from = paste0("pres_", modname, ".out"),
                  to = paste0("inst/output/", alpha, "/pres/", modname, ".out")))

      ## Read output file
      a <- scan(paste0('inst/output/', alpha, "/pres/", modname, ".out"), what='c',sep='\n',quiet=TRUE)

    ## Evaluate model (if model converges, will equal TRUE)
    check <- mod_eval(pres_out = a, pao2 = pao, mod = mods[[i]], strict.psi = TRUE, strict.gam = TRUE, is.het = opts$het, is.annual = annual)

      if(check == FALSE){ # If model does not converge, save NA in AIC table
        aic_temp <- dplyr::data_frame(Model = modname, Model_num = i, LogLik = NA, nParam = NA,
                                      AIC = NA)

      }else{ # If model does converge, save results to AIC table
        ## Extract log likelihood
        j <- grep('-2log', a)
        loglike <- as.numeric(unlist(strsplit(a[j],'=',2))[2])

        ## Extract AIC
        j <- grep('AIC', a)
        aic <- as.numeric(unlist(strsplit(a[j],'=',2))[2])

        ## Number of parameters
        j <- grep('of par', a)
        n  <- as.numeric(unlist(strsplit(a[j],'=',2))[2])

        aic_temp <- dplyr::data_frame(Model = modname, Model_num = i, LogLik = loglike, nParam = n,
                                      AIC = aic)
      }

      if(i == 1){
        aic_table <- aic_temp
      }else{
        aic_table <- dplyr::bind_rows(aic_table, aic_temp)
      }

    }
  }



  ## Add delta AIC column and sort by delta AIC
  aic_table <- dplyr::mutate(aic_table, delta_AIC = AIC - min(AIC, na.rm = TRUE))
  aic_table <- dplyr::arrange(aic_table, delta_AIC)

  ## Write AIC table
  write.csv(aic_table, file = paste0("inst/output/", alpha, "/gam_aic.csv"), row.names = FALSE)

}

