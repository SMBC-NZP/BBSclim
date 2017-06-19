#' test_annual
#'
#' Tests whether there is support for annual variation in p
#' @return AIC table comparing model with and without annual variation in p
#' @export

test_annual <- function(){
    opts <- read.csv("inst/model_opts.csv")
    spp_list <- read.csv('inst/spp_list_annual.csv')

    spp <- NULL
    ## Check if test_annual has already been run for species
    for(i in 1:length(spp_list$spp)){
      spp_test <- file.exists(paste0("inst/output/", spp_list$spp[i], "/annual_aic.csv"))
      if(!spp_test) spp <- c(spp, as.character(spp_list$spp[i]))
    }

    if(!is.null(spp)){
      ann <- rep(c(TRUE, FALSE), length(spp))
      spp2 <- rep(spp, each = 2)

      ## Fit models
        cores <- parallel::detectCores()
        if(!is.null(opts$limit.cores)){
          cores <- min(cores, opts$limit.cores)
        }
        doParallel::registerDoParallel(cores = cores)

        annual_aic <-  foreach::foreach(i=1:length(spp2), .combine = rbind,
                                        .packages = c("dplyr", "BBSclim")) %dopar%{

            psi_aic <- read.csv(paste0("inst/output/", spp2[i], "/psi_aic.csv"))
            top_covs <- GetPsiMods()[[psi_aic$Model_num[1]]]
            mod <- GetGamMods(psi_covs = top_covs$psi.cov)[[961]]

             spp_pao <- RPresence::read.pao(paste0("inst/output/", spp2[i], "/pres/pres_in.pao"))

             if(ann[i]){modname <- paste0(spp2[i], "_annual")}else{modname <- paste0(spp2[i], "_constant")}
             spp_dm <- suppressWarnings(BBSclim::GetDM(pao = spp_pao, cov_list = mod, is.annual = ann[i], is.het = opts$het))


             fixedpars <- matrix(rep("eq", spp_pao$nseasons), spp_pao$nseasons, 1)
             r1 <- dim(spp_dm$dm1)[1] + dim(spp_dm$dm2)[1] + dim(spp_dm$dm3)[1] + dim(spp_dm$dm4)[1]
             rownames(fixedpars) <- (r1 + 1):(r1 + spp_pao$nseasons)

             ## Run model
             RPresence::write_dm_and_run(paoname = spp_pao$paoname, fixed = fixedpars,
                                         dms = spp_dm, model = i,
                                         noderived = TRUE, limit.real = TRUE,
                                         modname = modname)

             file.rename(from = paste0("pres_", modname, ".out"),
                         to = paste0("inst/output/", spp2[i], "/pres/", sub(".*_", "", modname), ".out"))

             modname <- sub(".*_", "", modname)

             ## Read output file
             a <- scan(paste0('inst/output/', spp2[i], "/pres/", modname, ".out"), what='c',sep='\n',quiet=TRUE)

             ## Evaluate model (if model converges, will equal TRUE)
             check <- BBSclim::mod_eval(pres_out = a, pao2 = spp_pao, mod = mod, is.annual = ann[i], is.het = opts$het,
                                        nuisance = TRUE)

             if(check == FALSE){ # If model does not converge, save NA in AIC table
                 aic_temp <- dplyr::data_frame(Model = modname, LogLik = NA, nParam = NA,
                                               AIC = NA, sp = spp2[i])

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

               aic_temp <- dplyr::data_frame(Model = modname, LogLik = loglike, nParam = n,
                                             AIC = aic, sp = spp2[i])
               }
               aic_temp
             }

        doParallel::stopImplicitCluster()

        for(i in 1:length(spp)){
          aic_temp <- dplyr::filter(annual_aic, sp == spp[i])
          aic_temp <- dplyr::select(aic_temp, -sp)

          ## Add delta AIC column and sort by delta AIC
          aic_temp <- dplyr::mutate(aic_temp, delta_AIC = AIC - min(AIC, na.rm = TRUE))
          aic_temp <- dplyr::arrange(aic_temp, delta_AIC)

          ## Write AIC table
          write.csv(aic_temp, file = paste0("inst/output/", spp[i], "/annual_aic.csv"), row.names = FALSE)
        }
      check <- TRUE
    }else{
      check <- FALSE
    }
     check
}
