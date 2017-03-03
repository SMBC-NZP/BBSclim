#' RunPsiMods
#'
#' Runs full model set for 31 psi models, evalute each, write and return AIC table, delete output file (optional)
#' @param pao .pao file for species of interest
#' @param alpha alpha code for species of interest
#' @param psi_mods list containing the parameters for each model to evaluate
#' @export

RunPsiMods <- function(alpha, pao){
    opts <- read.csv("inst/model_opts.csv")

    annual_aic <- read.csv(paste0("inst/output/", alpha, "/annual_aic.csv"))

    mods1 <- GetGamMods()
    aic_tab <- read.csv(paste0("inst/output/", alpha, "/gam_aic.csv"))
    top <- aic_tab$Model_num[1]
    covs <- mods1[[top]]
    covs.ll <- list(gam_covs = covs$gam.cov, eps_covs = covs$eps.cov)

    mods <- GetPsiMods(covs = covs.ll)


    if(annual_aic$Model[1] == "annual"){
      annual <- TRUE
    }else{
      annual <- FALSE
    }

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
                                      spp_dm <- BBSclim::GetDM(pao = pao, cov_list = mods[[i]], is.annual = annual, is.het = opts$het)

                                      ## Run model
                                      RPresence::write_dm_and_run(paoname = pao$paoname, fixed = fixedpars,
                                                                  dms = spp_dm, model = i,
                                                                  noderived = TRUE, limit.real = TRUE,
                                                                  modname = modname)

                                      file.rename(from = paste0("pres_", modname, ".out"),
                                                  to = paste0("inst/output/", alpha, "/pres/", modname, ".out"))

                                      ## Read output file
                                      a <- scan(paste0('inst/output/', alpha, "/pres/", modname, ".out"), what='c',sep='\n',quiet=TRUE)

                                      ## Evaluate model (if model converges, will equal TRUE)
                                      check <- BBSclim::mod_eval(pres_out = a, pao2 = pao, mod = mods[[i]], strict = TRUE, is.annual = annual, is.het = opts$het)

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
      aic_table
    }else{
      if(opts$psi.test){
        mods <- mods[length(mods) - 1:length(mods)]
      }

      for(i in 1:length(mods)){

        modname <- paste0('psi_model_', i)

        ## Create design matrices for model i
        spp_dm <- GetDM(pao = pao, cov_list = mods[[i]], is.annual = annual, is.het = opts$het)

        ## Run model
        RPresence::write_dm_and_run(paoname = pao$paoname,
                                    dms = spp_dm, model = i,
                                    noderived = TRUE, limit.real = TRUE,
                                    modname = modname)

        file.rename(from = paste0("pres_", modname, ".out"),
                    to = paste0("inst/output/", alpha, "/pres/", modname, ".out"))


        ## Read output file
        a <- scan(paste0('inst/output/', alpha, "/pres/", modname, ".out"), what='c', sep='\n', quiet=TRUE)

      ## Evaluate model (if model converges, will equal TRUE)
      check <- mod_eval(pres_out = a, pao2 = pao, mod = mods[[i]], strict = TRUE, is.annual = annual, is.het = opts$het)


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


    b <- scan(paste0("inst/output/", alpha, "/pres/psi_model_31.out"), what='c',sep='\n',quiet=TRUE)

    j <- grep('-2log', b)
    loglike <- as.numeric(unlist(strsplit(b[j],'=',2))[2])

    ## Extract AIC
    j <- grep('AIC', b)
    aic <- as.numeric(unlist(strsplit(b[j],'=',2))[2])

    ## Number of parameters
    j <- grep('of par', b)
    n  <- as.numeric(unlist(strsplit(b[j],'=',2))[2])

    aic_last <- dplyr::data_frame(Model = "psi_model_31", Model_num = 31, LogLik = loglike, nParam = n,
                                  AIC = aic)

    aic_table <- dplyr::bind_rows(aic_table, aic_last)


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


RunGamMods <- function(alpha, pao){
  opts <- read.csv("inst/model_opts.csv")

  annual_aic <- read.csv(paste0("inst/output/", alpha, "/annual_aic.csv"))


  mods <- GetGamMods()


  if(annual_aic$Model[1] == "annual"){
    annual <- TRUE
  }else{
    annual <- FALSE
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
                                    spp_dm <- BBSclim::GetDM(pao = pao, cov_list = mods[[i]], is.annual = annual, is.het = opts$het)

                                    ## Run model
                                    RPresence::write_dm_and_run(paoname = pao$paoname,
                                                                dms = spp_dm, model = i,
                                                                modname = modname,
                                                                noderived = TRUE, limit.real = TRUE)

                                    file.rename(from = paste0("pres_", modname, ".out"),
                                                to = paste0("inst/output/", alpha, "/pres/", modname, ".out"))

                                    ## Read output file
                                    a <- scan(paste0('inst/output/', alpha, "/pres/", modname, ".out"), what='c',sep='\n',quiet=TRUE)

                                    ## Evaluate model (if model converges, will equal TRUE)
                                    check <- BBSclim::mod_eval(pres_out = a, pao2 = pao, mod = mods[[i]], strict = FALSE,
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
    aic_table
  }else{
    if(opts$gam.test){
      mods <- mods[1:cores]
    }

    for(i in 1:length(mods)){

      modname <- paste0('gam_model_', i)

      ## Create design matrices for model i
      spp_dm <- GetDM(pao = pao, cov_list = mods[[i]], is.het = opts$het, is.annual = annual)

      ## Run model
      RPresence::write_dm_and_run(paoname = pao$paoname,
                                  dms = spp_dm, model = i,
                                  modname = modname,
                                  noderived = TRUE, limit.real = TRUE)

      file.rename(from = paste0("pres_", modname, ".out"),
                  to = paste0("inst/output/", alpha, "/pres/", modname, ".out"))

      ## Read output file
      a <- scan(paste0('inst/output/', alpha, "/pres/", modname, ".out"), what='c',sep='\n',quiet=TRUE)

    ## Evaluate model (if model converges, will equal TRUE)
    check <- mod_eval(pres_out = a, pao2 = pao, mod = mods[[i]], strict = FALSE, is.het = opts$het, is.annual = annual)

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

  ## Top gamma/epsilon model == last psi model, so rename and save (to avoid running again)
  file.rename(from = paste0("inst/output/", alpha, "/pres/", aic_table$Model[1], ".out"),
              to = paste0("inst/output/", alpha, "/pres/", "psi_model_31.out"))

}

