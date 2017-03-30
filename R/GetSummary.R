#' GetSummary
#'
#' Generate summary information for report
#' @export

GetSummary <- function(alpha){
  mod_opts <- read.csv("inst/model_opts.csv")

  common <- BBSclim::code_lookup$common[BBSclim::code_lookup$alpha == toupper(alpha)]

  mods1 <- GetGamMods()
  aic_tab <- read.csv(paste0("inst/output/", alpha, "/gam_aic.csv"))
  top <- aic_tab$Model_num[1]
  covs <- mods1[[top]]
  covs.ll <- list(gam_covs = covs$gam.cov, eps_covs = covs$eps.cov)

  mods <- GetPsiMods(covs = covs.ll)

  ## Did p vary annually?
  annual_aic <- read.csv(paste0("inst/output/", alpha, "/annual_aic.csv"))

  if(annual_aic$Model[1] == "annual"){
    annual <- "Yes"
  }else{
    annual <- "No"
  }


  ## Count data
  raw_counts <- read.csv(paste0('inst/output/', alpha, '/raw_counts.csv'))
  used_counts <- read.csv(paste0('inst/output/', alpha, '/no_outlier_counts.csv'))
  buff_counts <- read.csv(paste0('inst/output/', alpha, '/count_buff.csv'))
  outliers <- unique(raw_counts$routeID[!('%in%'(raw_counts$routeID, used_counts$routeID))])

  years <- seq(from = min(raw_counts$Year), to = max(raw_counts$Year))
  nYears <- length(years)

  ### Annual variation in p AIC table
  colnames(annual_aic) <- c("Model", "LogLik", "k", "AIC", "$\\Delta$ AIC")


  ### Psi models AIC table
  psi_aic <- read.csv(paste0('inst/output/', alpha, '/psi_aic.csv'))
  modnum <- psi_aic$Model_num[1]

  ### Gamma/Epsilon AIC table
  gam_aic <- read.csv(paste0('inst/output/', alpha, '/gam_aic.csv'))
  gam_aic <- gam_aic[1:25,]

  ## Extract beta coef estimates and se
  top_mod <- scan(paste0("inst/output/", alpha, "/top_mod.out"), what='character', sep='\n', quiet=T)

  jj <- grep('std.error', top_mod)
  jj2 <- grep('Variance-Covariance Matrix of Untransformed', top_mod)

  betas <- top_mod[(jj+1):(jj2-1)]
  coefs <- round(as.numeric(substr(betas, 41,50)), digits = 2)
  std.er <- round(as.numeric(substr(betas, 54,63)), digits = 2)

  m <- grep('==>name', top_mod)


  ## Covariates included in the top model
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


  ## Beta coefficients for psi, gamma, & epsilon
  psi_beta_tab <- MakeBetatab(coefs = coefs, sd.err = std.er, alpha, nYears = nYears, covs_use = covs_use, years = years)

  ## Beta coefficients for th, th0, p, & omega
  p_beta_tab <- MakeBetatab(coefs = coefs, sd.err = std.er, alpha, nYears = nYears, covs_use = covs_use, years = years, nuisance = TRUE)


  ## Include only psi models that passed GOF test
  psi_aic <- psi_aic[which(psi_aic$Model_num == modnum):nrow(psi_aic),]
  psi_aic <- mutate(psi_aic, delta_AIC = AIC - min(AIC, na.rm = TRUE))
  psi_aic <- psi_aic[1:10,]

  colnames(gam_aic) <- c("Model", "Model #", "LogLik", "k", "AIC", "$\\Delta$ AIC")
  colnames(psi_aic) <- c("Model", "Model #", "LogLik", "k", "AIC", "$\\Delta$ AIC")

  summ <- list(spp_name = common,
       spp_alpha = alpha,
       annual = annual,
       n.routes = length(unique(used_counts$routeID)),
       n.outliers = length(outliers),
       n.buffer = length(unique(buff_counts$routeID)) - length(unique(used_counts$routeID)),
       annual.aic = annual_aic,
       psi.aic = psi_aic,
       gam.aic = gam_aic,
       psi.betas = psi_beta_tab,
       p.betas = p_beta_tab,
       top_mod = modnum)
  summ
}

#' MakeAICtab
#'
#' Make AIC table for including in report
#' @param nuisance If true, returns AIC table of theta, theta', p, & omega; if false, returns AIC table for psi, gamma & epsilon


MakeBetatab <- function(coefs, sd.err, alpha, covs_use, nYears, years, nuisance = FALSE){
  mod_opts <- read.csv("inst/model_opts.csv")

  ## Did p vary annually?
  annual_aic <- read.csv(paste0("inst/output/", alpha, "/annual_aic.csv"))

  if(annual_aic$Model[1] == "annual"){
    is.annual <- TRUE
  }else{
    is.annual <- FALSE
  }

  if(!nuisance){
    ## Covert covs included in psi, gam, & eps models to factor, set levels
    covs_use2 <- factor(unique(c(covs_use$psi.cov, covs_use$gam.cov, covs_use$eps.cov)),
                        levels = c("tmp", "sq_tmp", "Twet", "sq_Twet", "Prec",
                                   "sq_Prec", "Pwarm", "sq_Pwarm", "dtr", "sq_dtr"))


    ## Data frame containing beta coeffecients and se
    beta_est <- matrix(NA, nrow = 11, ncol = 3)
    colnames(beta_est) <- c("$\\psi$", "$\\gamma$", "$\\epsilon$")


    ## Fill in intercept values
    beta_est[1, 1] <- paste(coefs[1], " (", sd.err[1], ")", sep = "")
    beta_est[1, 2] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                     length(covs_use$th1.cov) + 4],
                             " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                            length(covs_use$th1.cov) + 4], ")", sep = "")
    beta_est[1, 3] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                     length(covs_use$th1.cov) + length(covs_use$gam.cov) + 5],
                             " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                            length(covs_use$th1.cov) + length(covs_use$gam.cov) + 5], ")", sep = "")

    ## Fill in coefficients for climate covariates
    for(i in 1:10){
      # Psi
      if(levels(covs_use2)[i] %in% covs_use$psi.cov){
        beta.num <- which(covs_use$psi.cov == levels(covs_use2)[i])
        beta_est[i + 1, 1] <- paste(coefs[1 + beta.num], " (", sd.err[1 + beta.num], ")", sep = "")
      }

      # Gamma
      if(levels(covs_use2)[i] %in% covs_use$gam.cov){
        beta.num <- which(covs_use$gam.cov == levels(covs_use2)[i])
        beta_est[i + 1, 2] <-  paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                              length(covs_use$th1.cov) + 4 + beta.num],
                                      " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                     length(covs_use$th1.cov) + 4 + beta.num], ")", sep = "")
      }

      # Epsilon
      if(levels(covs_use2)[i] %in% covs_use$eps.cov){
        beta.num <- which(covs_use$eps.cov == levels(covs_use2)[i])
        beta_est[i + 1, 3] <-   paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                               length(covs_use$th1.cov) + length(covs_use$gam.cov) + 5 + beta.num],
                                       " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                      length(covs_use$th1.cov) + length(covs_use$gam.cov) + 5 + beta.num], ")", sep = "")
      }
    }


    ## Rename climate covariates
    covs_use2 <- dplyr::recode(covs_use2, tmp = "Temp", sq_tmp = "Temp^2",
                               Twet = "Temp, Wettest Qrt", sq_Twet = "Temp, Wettest Qrt^2",
                               Prec = "Precip", sq_Prec = "Precip^2",
                               Pwarm = "Precip, Warmest Qrt", sq_Pwarm = "Precip, Warmest Qrt^2",
                               dtr = "Diurnal temp range", sq_dtr = "Diurnal temp range^2")

    ## Covert to data frame, add intercept, covert to character, replace NA with "-"
    beta_df <- as.data.frame(beta_est)
    covs <- data.frame(cov = c("Intercept", as.character(levels(covs_use2[order(covs_use2)]))))
    beta_df <- dplyr::bind_cols(covs, beta_df)
    beta_df[, 2:4] <- as.character(unlist(beta_df[, 2:4]))
    beta_df[is.na(beta_df)] <- "-"
    names(beta_df)[1] <- ""

  }else{
    ## Data frame containing beta coeffecients and se for theta, theta', p, omega ----
    covs_use2 <- factor(unique(c(covs_use$th0.cov, covs_use$th1.cov, covs_use$p1.cov)),
                        levels = c("Lat", "sq_Lat", "Lon", "sq_Lon", "Stop", "sq_Stop"))


    if(mod_opts$het){
      if(is.annual){
        beta_est <- matrix(NA, nrow = nYears + length(covs_use2), ncol = 5)
        colnames(beta_est) <- c("$\\theta$", "$\\theta'$", "$p_1$", "$p_2$", "$\\omega$")
      }else{
        beta_est <- matrix(NA, nrow = 1 + length(covs_use2), ncol = 5)
        colnames(beta_est) <- c("$\\theta$", "$\\theta'$", "$p_1$", "$p_2$", "$\\omega$")
      }
    }else{
      if(is.annual){
        beta_est <- matrix(NA, nrow = nYears + length(covs_use2), ncol = 3)
        colnames(beta_est) <- c("$\\theta$", "$\\theta'$", "$p$")
      }else{
        beta_est <- matrix(NA, nrow = 1 + length(covs_use2), ncol = 3)
        colnames(beta_est) <- c("$\\theta$", "$\\theta'$", "$p$")
      }
    }



    # Fill in intercept values (if annual == TRUE, different intercept for each year so don't include)
    beta_est[1, 1] <- paste(coefs[length(covs_use$psi.cov) + 2],
                             " (", sd.err[length(covs_use$psi.cov) + 2], ")", sep = "")
    beta_est[1, 2] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) + 3],
                             " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) + 3], ")", sep = "")
    if(is.annual){
      beta_est[1:nYears, 3] <- paste(coefs[(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                               length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                               length(covs_use$eps.cov) + 6):(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                                                length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                                                                length(covs_use$eps.cov) + 5 + nYears)],
                                      " (", sd.err[(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                      length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                                      length(covs_use$eps.cov) + 6):(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                                                       length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                                                                       length(covs_use$eps.cov) + 5 + nYears)], ")", sep = "")
      if(mod_opts$het){
        beta_est[1:nYears, 4] <- paste(coefs[(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                 length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                                 length(covs_use$eps.cov) + length(covs_use$p1.cov) + nYears + 6):(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                                                                                     length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                                                                                                     length(covs_use$eps.cov) + length(covs_use$p1.cov) + 5 + 2 * nYears)],
                                        " (", sd.err[(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                        length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                                        length(covs_use$eps.cov) + length(covs_use$p1.cov) + nYears + 6):(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                                                                                            length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                                                                                                            length(covs_use$eps.cov) + length(covs_use$p1.cov) + 5 + 2 * nYears)], ")", sep = "")

        beta_est[1, 5] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                         length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                         length(covs_use$eps.cov) + 2 * length(covs_use$p1.cov) +
                                         2 * nYears + 6],
                                 " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                                length(covs_use$eps.cov) + 2 * length(covs_use$p1.cov) +
                                                2 * nYears + 6], ")", sep = "")
      }
    }else{
      beta_est[1, 3] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                       length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                       length(covs_use$eps.cov) + 6],
                               " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                              length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                              length(covs_use$eps.cov) + 6], ")", sep = "")
      if(mod_opts$het){
        beta_est[1, 4] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                         length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                         length(covs_use$eps.cov) + length(covs_use$p1.cov) + 7],
                                 " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                                length(covs_use$eps.cov) + length(covs_use$p1.cov) + 7], ")", sep = "")

        beta_est[1, 5] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                         length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                         length(covs_use$eps.cov) + 2 * length(covs_use$p1.cov) + 8],
                                 " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                                length(covs_use$eps.cov) + 2 * length(covs_use$p1.cov) + 8], ")", sep = "")
      }
    }

    # Fill in coefficients for lat, lon, and stop
    for(i in 1:length(covs_use2)){
      if(is.annual){
        # th0 and th1
        if(levels(covs_use2)[i] %in% covs_use$th0.cov){
          beta.num <- which(covs_use$th0.cov == levels(covs_use2)[i])
          beta_est[i + nYears, 1] <- paste(coefs[length(covs_use$psi.cov) + 2 + beta.num],
                                           " (", sd.err[length(covs_use$psi.cov) + 2 + beta.num], ")", sep = "")
        }

        if(levels(covs_use2)[i] %in% covs_use$th1.cov){
          beta.num <- which(covs_use$th1.cov == levels(covs_use2)[i])
          beta_est[i + nYears, 2] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) + 3 + beta.num],
                                           " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) + 3 + beta.num], ")", sep = "")
        }
        if(mod_opts$het){
          if(levels(covs_use2)[i] %in% covs_use$p1.cov){
            beta.num <- which(covs_use$p1.cov == levels(covs_use2)[i])
            beta_est[i + nYears, 3] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) + length(covs_use$th1.cov) +
                                                   length(covs_use$gam.cov) + length(covs_use$eps.cov) + nYears + 5 + beta.num],
                                             " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) + length(covs_use$th1.cov) +
                                                            length(covs_use$gam.cov) + length(covs_use$eps.cov) + nYears + 5 + beta.num], ")", sep = "")

            beta_est[i + nYears, 4] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) + length(covs_use$th1.cov) +
                                                     length(covs_use$gam.cov) + length(covs_use$eps.cov) + length(covs_use$p1.cov) + 2 * nYears + 5 + beta.num],
                                             " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) + length(covs_use$th1.cov) +
                                                            length(covs_use$gam.cov) + length(covs_use$eps.cov) + length(covs_use$p1.cov) + 2 * nYears + 5 + beta.num], ")", sep = "")
          }
        }else{
            if(levels(covs_use2)[i] %in% covs_use$p1.cov){
              beta.num <- which(covs_use$p1.cov == levels(covs_use2)[i])
              beta_est[i + nYears, 3] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) + length(covs_use$th1.cov) +
                                                       length(covs_use$gam.cov) + length(covs_use$eps.cov) + nYears + 5 + beta.num],
                                               " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) + length(covs_use$th1.cov) +
                                                              length(covs_use$gam.cov) + length(covs_use$eps.cov) + nYears + 5 + beta.num], ")", sep = "")
          }
        }
      }else{
        # th0 and th1
        if(levels(covs_use2)[i] %in% covs_use$th0.cov){
          beta.num <- which(covs_use$th0.cov == levels(covs_use2)[i])
          beta_est[i + 1, 1] <- paste(coefs[length(covs_use$psi.cov) + 2 + beta.num],
                                           " (", sd.err[length(covs_use$psi.cov) + 2 + beta.num], ")", sep = "")
        }

        if(levels(covs_use2)[i] %in% covs_use$th1.cov){
          beta.num <- which(covs_use$th1.cov == levels(covs_use2)[i])
          beta_est[i + 1, 2] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) + 3 + beta.num],
                                           " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) + 3 + beta.num], ")", sep = "")
        }
        if(mod_opts$het){
          if(levels(covs_use2)[i] %in% covs_use$p1.cov){
            beta.num <- which(covs_use$p1.cov == levels(covs_use2)[i])
            beta_est[i + 1, 3] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) + length(covs_use$th1.cov) +
                                                     length(covs_use$gam.cov) + length(covs_use$eps.cov) + 6 + beta.num],
                                             " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) + length(covs_use$th1.cov) +
                                                            length(covs_use$gam.cov) + length(covs_use$eps.cov) + 6 + beta.num], ")", sep = "")

            beta_est[i + 1, 4] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) + length(covs_use$th1.cov) +
                                                     length(covs_use$gam.cov) + length(covs_use$eps.cov) + length(covs_use$p1.cov) + 7 + beta.num],
                                             " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) + length(covs_use$th1.cov) +
                                                            length(covs_use$gam.cov) + length(covs_use$eps.cov) + length(covs_use$p1.cov) + 7 + beta.num], ")", sep = "")
          }
        }else{
          if(levels(covs_use2)[i] %in% covs_use$p1.cov){
            beta.num <- which(covs_use$p1.cov == levels(covs_use2)[i])
            beta_est[i + 1, 3] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) + length(covs_use$th1.cov) +
                                                     length(covs_use$gam.cov) + length(covs_use$eps.cov) + 6 + beta.num],
                                             " (", sd.err[length(covs_use$psi.cov) + length(covs_use$th0.cov) + length(covs_use$th1.cov) +
                                                            length(covs_use$gam.cov) + length(covs_use$eps.cov) + 6 + beta.num], ")", sep = "")
          }
        }
      }

    }

    ## Rename covariates
    covs_use2 <- dplyr::recode(covs_use2, Lat = "Latitude", sq_Lat = "Latitude^2",
                               Lon = "Longitude", sq_Lon = "Longitude^2",
                               sq_Stop = "Stop^2")


    ## Covert to data frame, add intercept, covert to character, replace NA with "-"
    beta_df <- as.data.frame(beta_est)
    if(is.annual){
      covs <- data.frame(cov = c("Intercept", paste0("Year_", seq(from = min(years) + 1, to = max(years))), as.character(levels(covs_use2[order(covs_use2)]))))
    }else{
      covs <- data.frame(cov = c("Intercept", as.character(levels(covs_use2[order(covs_use2)]))))
    }
    beta_df <- dplyr::bind_cols(covs, beta_df)
    beta_df[, 2:6] <- as.character(unlist(beta_df[, 2:6]))
    beta_df[is.na(beta_df)] <- "-"
    names(beta_df)[1] <- ""
  }

  beta_df
}
