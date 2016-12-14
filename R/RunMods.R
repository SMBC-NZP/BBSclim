#' RunPsiMods
#'
#' Runs full model set for 31 psi models, evalute each, write and return AIC table, delete output file (optional)
#' @param pao .pao file for species of interest
#' @param alpha alpha code for species of interest
#' @param psi_mods list containing the parameters for each model to evaluate
#' @param del Should .out files be deleted after model output in evaluated?
#' @param time Does p vary annually?
#' @param het Heterogeneous detection?
#' @export

RunPsiMods <- function(pao, alpha, mods = psi_mods, del = TRUE, ...,
                       test = FALSE, Parallel = TRUE){

    if(Parallel){
      cores <- parallel::detectCores()
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)

      if(test){
        mods <- mods[1:cores]
      }

      aic_table <- foreach::foreach(i=1:length(mods), .combine = rbind,
                                    .packages = c("dplyr", "BBSclim")) %dopar%{
                                      modname <- paste0(alpha, '_psi_model_', i)

                                      ## Create design matrices for model i
                                      spp_dm <- BBSclim::GetDM(pao = pao, cov_list = mods[[i]], ...)

                                      ## Run model
                                      BBSclim::write_dm_and_run2(pao = pao, cov_list = mods[[i]], ..., dm_list = spp_dm,
                                                                 modname = modname, fixed = TRUE, out = "temp",
                                                                 inits = TRUE, maxfn = '32000 vc lmt=5', alpha = alpha)

                                      ## Read output file
                                      a <- scan(paste0('inst/output/pres/temp/', modname, ".out"), what='c',sep='\n',quiet=TRUE)

                                      ## Evaluate model (if model converges, will equal TRUE)
                                      check <- BBSclim::mod_eval(pres_out = a, pao2 = pao, mod = mods[[i]], ...)

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
      if(test){
        mods <- mods[1:2]
      }

      for(i in 1:length(mods)){

        modname <- paste0(alpha, '_psi_model_', i)

        ## Create design matrices for model i
        spp_dm <- GetDM(pao = pao, cov_list = mods[[i]], ...)

        ## Run model
        write_dm_and_run2(pao = pao, cov_list = mods[[i]], ..., dm_list = spp_dm,
                          modname = modname, fixed = TRUE, out = "temp",
                          inits = TRUE, maxfn = '32000 vc lmt=5', alpha = alpha)

        ## Read output file
        a <- scan(paste0('inst/output/pres/temp/', modname, ".out"), what='c',sep='\n',quiet=TRUE)

      ## Evaluate model (if model converges, will equal TRUE)
      check <- mod_eval(pres_out = a, pao2 = pao, mod = mods[[i]], ...)


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
          aic_tab <- aic_temp
        }else{
          aic_table <- dplyr::bind_rows(aic_tab, aic_temp)
        }

       }
      }


      ## Add delta AIC column and sort by delta AIC
      aic_table <- dplyr::mutate(aic_table, delta_AIC = AIC - min(AIC))
      aic_table <- dplyr::arrange(aic_table, delta_AIC)

      ## Write AIC table
      write.csv(aic_table, file = paste0("inst/output/aic/psi/", alpha, ".csv"), row.names = FALSE)

      ## Top psi model == last gam model, so rename and save (to avoid running again)
      file.rename(from = paste0("inst/output/pres/temp/", aic_table$Model[1], ".out"),
                  to = paste0("inst/output/pres/temp/", alpha, "_gam_model_961.out"))
      temp.files <- list.files("inst/output/pres/temp")
      temp.psi.files <- temp.files[grep("psi", temp.files)]
      if(del) file.remove(paste0("inst/output/pres/temp/", temp.psi.files))

      ## Return AIC table
      aic_table

}

#' top_covs
#'
#' Get covariates of top model
#' @param aic_tab AIC table from RunPsiMods or RunGamMods
#' @param mods List of models
#' @param psi Return top psi model (TRUE) or gamma/epsilon model (FALSE)
#' @export
#'

top_covs <- function(aic_tab, mods, psi = TRUE){
    if(psi){
      top <- aic_tab$Model_num[1]
      covs <- mods[[top]]$psi.cov
      covs
    }else{
      top <- aic_tab$Model_num[1]
      covs <- mods[[top]]
      covs
    }
}

#' RunGamMods
#'
#' Runs full model set for 961 gamma/epsilon models, evaluate each, write and return AIC table, delete .out files (optional)
#' @param gam_mods list containing the covariates for each model
#' @export


RunGamMods <- function(pao, alpha, mods = gam_mods, del = TRUE, ...,
                       test = FALSE, trim = TRUE, Parallel = TRUE){
  if(Parallel){
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)

    if(test){
      mods <- mods[1:cores]
    }

    aic_table <- foreach::foreach(i=1:length(mods), .combine = rbind,
                                  .packages = c("dplyr", "BBSclim")) %dopar%{

                                    modname <- paste0(alpha, '_gam_model_', i)

                                    ## Create design matrices for model i
                                    spp_dm <- BBSclim::GetDM(pao = pao, cov_list = mods[[i]], ...)

                                    ## Run model
                                    BBSclim::write_dm_and_run2(pao = pao, cov_list = mods[[i]], ..., dm_list = spp_dm,
                                                               modname = modname, fixed = TRUE, out = "temp",
                                                               inits = TRUE, maxfn = '32000 vc lmt=5', alpha = alpha)

                                    ## Read output file
                                    a <- scan(paste0('inst/output/pres/temp/', modname, ".out"), what='c',sep='\n',quiet=TRUE)

                                    ## Evaluate model (if model converges, will equal TRUE)
                                    check <- BBSclim::mod_eval(pres_out = a, pao2 = pao, mod = mods[[i]], ...)

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
    if(test){
      mods <- mods[1:2]
    }

    for(i in 1:length(mods)){

      modname <- paste0(alpha, '_gam_model_', i)

      ## Create design matrices for model i
      spp_dm <- GetDM(pao = pao, cov_list = mods[[i]], ...)

      ## Run model
      write_dm_and_run2(pao = pao, cov_list = mods[[i]], ..., dm_list = spp_dm,
                        modname = modname, fixed = TRUE, out = "temp",
                        inits = TRUE, maxfn = '32000 vc lmt=5', alpha = alpha)

      ## Read output file
      a <- scan(paste0('inst/output/pres/temp/', modname, ".out"), what='c',sep='\n',quiet=TRUE)

    ## Evaluate model (if model converges, will equal TRUE)
    check <- mod_eval(pres_out = a, pao2 = pao, mod = mods[[i]], ...)

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
        aic_tab <- aic_temp
      }else{
        aic_table <- dplyr::bind_rows(aic_tab, aic_temp)
      }

    }
  }

  if(del) file.remove(paste0("inst/output/pres/temp/", modname, ".out"))
  ## Add delta AIC column and sort by delta AIC
  aic_table <- dplyr::mutate(aic_table, delta_AIC = AIC - min(AIC))
  aic_table <- dplyr::arrange(aic_table, delta_AIC)

  ## Write AIC table
  write.csv(aic_table, file = paste0("inst/output/aic/gam/", alpha, ".csv"), row.names = FALSE)

  ## Return AIC table
  if(trim) aic_table <- aic_tab[1:min(nrow(aic_table), 25), ]
  aic_table

}

