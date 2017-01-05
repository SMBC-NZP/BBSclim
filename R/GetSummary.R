#' GetSummary
#'
#' Generate summary information for report
#' @export

GetSummary <- function(alpha, gam_mods){
  model_opts <- read.csv("inst/model_opts.csv")
  global_opts <- read.csv("inst/global_opts.csv")

  nYears <- length(seq(from = global_opts$start_yr, to = global_opts$end_yr))

  common <- code_lookup$common[code_lookup$alpha == toupper(alpha)]

  raw_counts <- read.csv(paste0('inst/output/', alpha, '/raw_counts.csv'))
  used_counts <- read.csv(paste0('inst/output/', alpha, '/no_outlier_counts.csv'))
  buff_counts <- read.csv(paste0('inst/output/', alpha, '/count_buff.csv'))
  outliers <- unique(raw_counts$routeID[!('%in%'(raw_counts$routeID, used_counts$routeID))])
  summ <- data_frame(spp_name = common,
                     n.routes = nrow(used_counts),
                     n.outliers = length(outliers))

  psi_aic <- read.csv(paste0('inst/output/', alpha, '/psi_aic.csv'))
  psi_aic <- psi_aic[1:max(which(!is.na(psi_aic$AIC))),]

  gam_aic <- read.csv(paste0('inst/output/', alpha, '/gam_aic.csv'))
  gam_aic <- gam_aic[1:max(which(!is.na(gam_aic$AIC))),]

  ## Extract beta coef estimates and se
  top_mod <- scan(paste0("inst/output/", alpha, "/top_mod.out"), what='character', sep='\n', quiet=T)

  jj <- grep('std.error', top_mod)
  jj2 <- grep('Variance-Covariance Matrix of Untransformed', top_mod)

  betas <- top_mod[(jj+1):(jj2-1)]
  coefs <- round(as.numeric(substr(betas, 41,50)), digits = 2)
  std.er <- round(as.numeric(substr(betas, 54,63)), digits = 2)

  ## Covariates included in the top model
  covs_use <- gam_mods[[modnum]]

  covs_use2 <- factor(unique(c(covs_use$psi.cov, covs_use$gam.cov, covs_use$eps.cov)),
                      levels = c("tmp", "sq_tmp", "Twet", "sq_Twet", "Prec",
                                 "sq_Prec", "Pwarm", "sq_Pwarm", "dtr", "sq_dtr"))


  ## Data frame containing beta coeffecients and se for psi, gamma, epsilon
  beta_est1 <- matrix(NA, nrow = 11, ncol = 3)
  colnames(beta_est1) <- c("$\\psi$", "$\\gamma$", "$\\epsilon$")


  # Fill in intercept values
  beta_est1[1, 1] <- paste(coefs[1], " (", std.er[1], ")", sep = "")
  beta_est1[1, 2] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                  length(covs_use$th1.cov) + 4],
                          " (", std.er[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                         length(covs_use$th1.cov) + 4], ")", sep = "")
  beta_est1[1, 3] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                  length(covs_use$th1.cov) + length(covs_use$gam.cov) + 5],
                          " (", std.er[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                         length(covs_use$th1.cov) + length(covs_use$gam.cov) + 5], ")", sep = "")

  # Fill in coefficients for climate covariates
  for(i in 1:10){
    # Psi
    if(levels(covs_use2)[i] %in% covs_use$psi.cov){
      beta.num <- which(covs_use$psi.cov == levels(covs_use2)[i])
      beta_est1[i + 1, 1] <- paste(coefs[1 + beta.num], " (", std.er[1 + beta.num], ")", sep = "")
    }

    # Gamma
    if(levels(covs_use2)[i] %in% covs_use$gam.cov){
      beta.num <- which(covs_use$gam.cov == levels(covs_use2)[i])
      beta_est1[i + 1, 2] <-  paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                           length(covs_use$th1.cov) + 4 + beta.num],
                                   " (", std.er[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                  length(covs_use$th1.cov) + 4 + beta.num], ")", sep = "")
    }

    # Epsilon
    if(levels(covs_use2)[i] %in% covs_use$eps.cov){
      beta.num <- which(covs_use$eps.cov == levels(covs_use2)[i])
      beta_est1[i + 1, 3] <-   paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                            length(covs_use$th1.cov) + length(covs_use$gam.cov) + 5 + beta.num],
                                    " (", std.er[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                   length(covs_use$th1.cov) + length(covs_use$gam.cov) + 5 + beta.num], ")", sep = "")
    }
  }


  covs_use2 <- dplyr::recode(covs_use2, tmp = "Temp", sq_tmp = "Temp^2",
                             Twet = "Temp, Wettest Qrt", sq_Twet = "Temp, Wettest Qrt^2",
                             Prec = "Precip", sq_Prec = "Precip^2",
                             Pwarm = "Precip, Warmest Qrt", sq_Pwarm = "Precip, Warmest Qrt^2",
                             dtr = "Diurnal temp range", sq_dtr = "Diurnal temp range^2")

  beta_df1 <- as.data.frame(beta_est1)
  covs <- data.frame(cov = c("Intercept", as.character(levels(covs_use2[order(covs_use2)]))))
  beta_df1 <- dplyr::bind_cols(covs, beta_df1)

  ## Data frame containing beta coeffecients and se for theta, theta', p, omega ----
  covs_use3 <- factor(unique(c(covs_use$th0.cov, covs_use$th1.cov, covs_use$p1.cov)),
                      levels = c("Lat", "sq_Lat", "Lon", "sq_Lon", "Stop", "sq_Stop"))


  if(model_opts$het){
    if(model_opts$annual){
      beta_est2 <- matrix(NA, nrow = nYears + length(covs_use3), ncol = 5)
      colnames(beta_est2) <- c("$\\theta$", "$\\theta'$", "$p_1$", "$p_2$", "$\\omega$")
    }else{
      beta_est2 <- matrix(NA, nrow = 1 + length(covs_use3), ncol = 5)
      colnames(beta_est2) <- c("$\\theta$", "$\\theta'$", "$p_1$", "$p_2$", "$\\omega$")
    }
  }else{
    if(model_opts$annual){
      beta_est2 <- matrix(NA, nrow = nYears + length(covs_use3), ncol = 3)
      colnames(beta_est2) <- c("$\\theta$", "$\\theta'$", "$p$")
    }else{
      beta_est2 <- matrix(NA, nrow = 1 + length(covs_use3), ncol = 3)
      colnames(beta_est2) <- c("$\\theta$", "$\\theta'$", "$p$")
    }
  }



  # Fill in intercept values (if annual == TRUE, different intercept for each year so don't include)
  beta_est2[1, 1] <- paste(coefs[length(covs_use$psi.cov) + 2],
                          " (", std.er[length(covs_use$psi.cov) + 2], ")", sep = "")
  beta_est2[1, 2] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) + 3],
                          " (", std.er[length(covs_use$psi.cov) + length(covs_use$th0.cov) + 3], ")", sep = "")
  if(model_opts$annual){
    beta_est2[1:nYears, 3] <- paste(coefs[(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                     length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                     length(covs_use$eps.cov) + 6):(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                                            length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                                                            length(covs_use$eps.cov) + 5 + nYears)],
                             " (", std.er[(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                             length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                             length(covs_use$eps.cov) + 6):(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                                              length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                                                              length(covs_use$eps.cov) + 5 + nYears)], ")", sep = "")
    if(model_opts$het){
      beta_est2[1:nYears, 4] <- paste(coefs[(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                               length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                               length(covs_use$eps.cov) + length(covs_use$p1.cov) + nYears + 6):(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                                                length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                                                                length(covs_use$eps.cov) + length(covs_use$p1.cov) + 5 + 2 * nYears)],
                                      " (", std.er[(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                      length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                                      length(covs_use$eps.cov) + length(covs_use$p1.cov) + nYears + 6):(length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                                                                                                          length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                                                                                                          length(covs_use$eps.cov) + length(covs_use$p1.cov) + 5 + 2 * nYears)], ")", sep = "")

      beta_est2[1, 5] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                      length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                      length(covs_use$eps.cov) + 2 * length(covs_use$p1.cov) +
                                      2 * nYears + 6],
                              " (", std.er[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                             length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                             length(covs_use$eps.cov) + 2 * length(covs_use$p1.cov) +
                                             2 * nYears + 6], ")", sep = "")
    }
  }else{
    beta_est2[1, 3] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                     length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                     length(covs_use$eps.cov) + 6],
                             " (", std.er[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                            length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                            length(covs_use$eps.cov) + 6], ")", sep = "")
    if(model_opts$het){
      beta_est2[1, 4] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                       length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                       length(covs_use$eps.cov) + length(covs_use$p1.cov) + 7],
                               " (", std.er[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                              length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                              length(covs_use$eps.cov) + length(covs_use$p1.cov) + 7], ")", sep = "")

      beta_est2[1, 5] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                       length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                       length(covs_use$eps.cov) + 2 * length(covs_use$p1.cov) + 8],
                               " (", std.er[length(covs_use$psi.cov) + length(covs_use$th0.cov) +
                                              length(covs_use$th1.cov) + length(covs_use$gam.cov) +
                                              length(covs_use$eps.cov) + 2 * length(covs_use$p1.cov) + 8], ")", sep = "")
    }
  }

  # Fill in coefficients for lat, lon, and stop
  for(i in 1:length(covs_use3)){
    # Psi
   if(model_opts$annual){
     if(levels(covs_use3)[i] %in% covs_use$th0.cov){
       beta.num <- which(covs_use$th0.cov == levels(covs_use3)[i])
       beta_est2[i + nYears, 1] <- paste(coefs[length(covs_use$psi.cov) + 2 + beta.num],
                                    " (", std.er[length(covs_use$psi.cov) + 2 + beta.num], ")", sep = "")
     }

     if(levels(covs_use3)[i] %in% covs_use$th1.cov){
       beta.num <- which(covs_use$th1.cov == levels(covs_use3)[i])
       beta_est2[i + nYears, 2] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) + 3 + beta.num],
                                        " (", std.er[length(covs_use$psi.cov) + length(covs_use$th0.cov) + 3 + beta.num], ")", sep = "")
     }
   }else{
     if(levels(covs_use3)[i] %in% covs_use$th0.cov){
       beta.num <- which(covs_use$th0.cov == levels(covs_use3)[i])
       beta_est2[i + 1, 1] <- paste(coefs[length(covs_use$psi.cov) + 2 + beta.num],
                                         " (", std.er[length(covs_use$psi.cov) + 2 + beta.num], ")", sep = "")
     }

     if(levels(covs_use3)[i] %in% covs_use$th1.cov){
       beta.num <- which(covs_use$th1.cov == levels(covs_use3)[i])
       beta_est2[i + 1, 2] <- paste(coefs[length(covs_use$psi.cov) + length(covs_use$th0.cov) + 3 + beta.num],
                                         " (", std.er[length(covs_use$psi.cov) + length(covs_use$th0.cov) + 3 + beta.num], ")", sep = "")
     }
   }

  }




  list(spp_name = common,
       spp_alpha = alpha,
       n.routes = nrow(used_counts),
       n.outliers = length(outliers),
       n.buffer = nrow(buff_counts) - nrow(used_counts),
       psi.aic = psi_aic,
       gam.aic = gam_aic)
}
