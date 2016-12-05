#' RunPsiMods
#'
#' Runs full model set for 31 psi models and writes output to ~/pres/psi
#' @param pao .pao file for species of interest
#' @param alpha alpha code for species of interest
#' @param psi_mods list containing the parameters for each model to evaluate
#' @param del Should .out files be deleted after model output in evaluated?
#' @param time Does p vary annually?
#' @param het Heterogeneous detection?
#' @export

RunPsiMods <- function(pao, alpha, mods = psi_mods, del = TRUE, time = annual, het = het_det,
                       test = FALSE){

    for(i in 1:length(mods)){

      modname <- paste0(alpha, '_psi_model_', i)

      if(!test){
        ## Create design matrices for model i
        spp_dm <- GetDM(pao = pao, cov_list = mods[[i]], time = time, het = het)

        ## Run model
        write_dm_and_run2(pao = pao, cov_list = mods[[i]], het = het_det, dm_list = spp_dm,
                          modname = modname, fixed = TRUE, out = "temp",
                          inits = TRUE, maxfn = '35000 lmt=5', alpha = alpha)
      }

      ## Read output file
      a <- scan(paste0('inst/output/temp/psi/', modname, ".out"), what='c',sep='\n',quiet=TRUE)

      ## Evaluate model (if model converges, will equal TRUE)
      check <- mod_eval(a)

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
        aic_tab <- dplyr::bind_rows(aic_tab, aic_temp)
      }

      if(del) file.remove(paste0("inst/output/pres/temp/", modname, ".out"))
    }

  ## Add delta AIC column and sort by delta AIC
  aic_tab <- dplyr::mutate(aic_tab, delta_AIC = AIC - min(AIC))
  aic_tab <- dplyr::arrange(aic_tab, delta_AIC)

  ## Write AIC table
  write.csv(aic_tab, file = paste0("inst/output/aic/psi/", alpha, ".csv"), row.names = FALSE)

  ## Return AIC table
  aic_tab
}

#' top_psi
#'
#' Get covariates for top psi model
#' @param aic_tab AIC table from RunPsiMods
#' @export
#'

top_psi <- function(aic_tab, mods = psi_mods){
  top <- aic_tab$Model_num[1]
  covs <- mods[[top]]$psi.cov
  covs
}


#' RunGamMods
#'
#' Runs full model set for 961 gamma/epsilon models and writes output to ~/pres/gam_eps
#' @param gam_mods list containing the covariates for each model
#' @export


RunGamMods <- function(pao, alpha, mods = gam_mods, del = TRUE, time = annual, het = het_det,
                       test = FALSE){

  for(i in 1:length(mods)){

    modname <- paste0(alpha, '_gam_model_', i)

    if(!test){
      ## Create design matrices for model i
      spp_dm <- GetDM(pao = pao, cov_list = mods[[i]], time = time, het = het)

      ## Run model
      write_dm_and_run2(pao = pao, cov_list = psi_mods[[i]], het = TRUE, dm_list = spp_dm,
                        modname = modname, fixed = TRUE, out = "temp",
                        inits = TRUE, maxfn = '35000 lmt=5', alpha = alpha)
    }

    ## Read output file
    a <- scan(paste0('inst/output/pres/temp/', modname, ".out"), what='c',sep='\n',quiet=TRUE)

    ## Evaluate model (if model converges, will equal TRUE)
    check <- mod_eval(a)

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
      aic_tab <- dplyr::bind_rows(aic_tab, aic_temp)
    }

    if(del) file.remove(paste0("inst/output/pres/temp/", modname, ".out"))
  }

  ## Add delta AIC column and sort by delta AIC
  aic_tab <- dplyr::mutate(aic_tab, delta_AIC = AIC - min(AIC))
  aic_tab <- dplyr::arrange(aic_tab, delta_AIC)

  ## Write AIC table
  write.csv(aic_tab, file = paste0("inst/output/aic/gam/", alpha, ".csv"), row.names = FALSE)

  ## Return AIC table
  aic_tab
}
