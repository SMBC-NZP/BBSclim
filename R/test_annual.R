#' test_annual
#'
#' Tests whether there is support for annual variation in p
#' @param alpha Alpha code for species of interest
#' @param pao .pao file
#' @return AIC table comparing model with and without annual variation in p
#' @export

test_annual <- function(alpha, pao){
    opts <- read.csv("inst/model_opts.csv")
    mod <- GetGamMods()[[961]]

    if(opts$Parallel){
      cores <- 2
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)

      annual_aic <- foreach::foreach(i=1:2, .combine = rbind,
                                     .packages = c("dplyr", "BBSclim")) %dopar%{

                                      ## Create design matrices for model i
                                      if(i == 1){
                                        modname <- "annual"
                                        spp_dm <- BBSclim::GetDM(pao = pao, cov_list = mod, is.annual = TRUE, is.het = opts$het)
                                        annual <- TRUE
                                      }else{
                                        modname <- "constant"
                                        spp_dm <- BBSclim::GetDM(pao = pao, cov_list = mod, is.annual = FALSE, is.het = opts$het)
                                        annual <- FALSE
                                      }


                                       fixedpars <- matrix(rep("eq", pao$nseasons), pao$nseasons, 1)
                                       r1 <- dim(spp_dm$dm1)[1] + dim(spp_dm$dm2)[1] + dim(spp_dm$dm3)[1] + dim(spp_dm$dm4)[1] 
                                       rownames(fixedpars) <- (r1 + 1):(r1 + pao$nseasons)
                                       
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
                                      check <- BBSclim::mod_eval(pres_out = a, pao2 = pao, mod = mod, is.annual = annual, is.het = opts$het)

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
      annual_aic
    }else{
      for(i in 1:2){

        if(i == 1){
          modname <- "annual"
          spp_dm <- BBSclim::GetDM(pao = pao, cov_list = mod, is.annual = TRUE, is.het = opts$het)
        }else{
          modname <- "constant"
          spp_dm <- BBSclim::GetDM(pao = pao, cov_list = mod, is.annual = FALSE, is.het = opts$het)
        }

        fixedpars <- matrix(rep("eq", pao$nseasons), pao$nseasons, 1)
        r1 <- dim(spp_dm$dm1)[1] + dim(spp_dm$dm2)[1] + dim(spp_dm$dm3)[1] + dim(spp_dm$dm4)[1] 
        rownames(fixedpars) <- (r1 + 1):(r1 + pao$nseasons)
        
        ## Run model
        RPresence::write_dm_and_run(paoname = pao$paoname, fixed = fixedpars,
                                    dms = spp_dm, model = i,
                                    noderived = TRUE, limit.real = TRUE,
                                    modname = modname)

        file.rename(from = paste0("pres_", modname, ".out"),
                    to = paste0("inst/output/", alpha, "/pres/", modname, ".out"))

        ## Read output file
        a <- scan(paste0('inst/output/', alpha, "/pres/", modname, ".out"), what='c', sep='\n', quiet=TRUE)

        ## Evaluate model (if model converges, will equal TRUE)
        check <- BBSclim::mod_eval(pres_out = a, pao2 = pao, mod = mod, is.annual = annual, is.het = opts$het)

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
          annual_aic <- dplyr::bind_rows(aic_table, aic_temp)
        }

      }
    annual_aic
    }


    ## Add delta AIC column and sort by delta AIC
    annual_aic <- dplyr::mutate(annual_aic, delta_AIC = AIC - min(AIC, na.rm = TRUE))
    annual_aic <- dplyr::arrange(annual_aic, delta_AIC)

    ## Write AIC table
    write.csv(annual_aic, file = paste0("inst/output/", alpha, "/annual_aic.csv"), row.names = FALSE)
}
