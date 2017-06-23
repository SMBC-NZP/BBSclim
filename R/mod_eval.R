#' mod_eval
#'
#' Checks for problems with convergence
#' @param mod Model output file
#' @param dig Minimum threshold for significant digits
#' @param large Maximum threshold for beta coefficients
#' @return TRUE if model converged; FALSE if not
#' @export

mod_eval <- function(pres_out, mod, pao2, psi = FALSE, dig = 2, large = 8, is.annual, is.het,
                     check.quad = FALSE, strict.psi = FALSE, strict.gam = FALSE,
                     nuisance = FALSE){
  jj <- grep('std.error', pres_out)
  jj2 <- grep('Variance-Covariance Matrix of Untransformed', pres_out)

  betas <- pres_out[(jj + 1):(jj2 - 1)]
  coefs <- as.numeric(substr(betas, 41, 50))
  std.er <- as.numeric(substr(betas, 54, 63))

  num.betas <- plyr::ldply(mod, function(x) length(x))$V1

  if(psi){ # Check that psi parameters are estimable
    psi.check <- NULL

    if(num.betas[1] > 1) {
      ## Check for negative SE
      psi.check 	<- !is.na(std.er[2:(num.betas[1] + 1)])

      ## Check for very coeficients
      psi.check 	<- abs(coefs[2:(num.betas[1] + 1)]) < large

      ## Check for significant wrong sign on quadratic terms
      if(check.quad){
        quads <- grep("psi.sq", betas)
        quads.coord <- quads[(length(quads) - 1):length(quads)]
        quads.clim <- quads[!(quads %in% quads.coord)]
        #psi.check[quads.coord - 1] <- coefs[quads.coord] < 0
        if(strict.psi){
          psi.check[quads.clim - 1] <- coefs[quads.clim] < 0
        }else{
          psi.check[quads.clim - 1] <- coefs[quads.clim]/std.er[quads.clim] < 1.96
        }
      }
    }

    wonky2 <- is.na(sum(std.er[1:(num.betas[1] + 1)]))

    if(mean(psi.check) == 1 & wonky2 == 0){
      ok <- TRUE
    }else{
      ok <- FALSE
    }


  }else{
    if(is.annual) num.betas[6] <- num.betas[6] + pao2$nseasons

    psi.check <- th0.check <- th1.check <- gam.check <- eps.check <- p1.check <- p2.check <- stop.check <- NULL
    ### drop covariates with negative standard errors
    if(nuisance == FALSE){
      if(num.betas[1] > 1) {
        psi.check 	<- !is.na(std.er[2:(num.betas[1] + 1)])
      }
      if(num.betas[4] > 1) {
        gam.check 	<- !is.na(std.er[(5 + sum(num.betas[1:3])):(sum(num.betas[1:4]) + 4)])
      }
      if(num.betas[5] > 1) {
        eps.check 	<- !is.na(std.er[(6 + sum(num.betas[1:4])):(sum(num.betas[1:5]) + 5)])
      }
    }

    if(num.betas[2] > 1) {
      th0.check 	<- !is.na(std.er[(3 + num.betas[1]):(sum(num.betas[1:2]) + 2)])
    }
    if(num.betas[3] > 1) {
      th1.check 	<- !is.na(std.er[(4 + sum(num.betas[1:2])):(sum(num.betas[1:3]) + 3)])
    }

    if(num.betas[6] > 1) {
      p1.check <- p2.check  <- !is.na(std.er[(6 + sum(num.betas[1:5])):(sum(num.betas[1:6]) + 5)])
      if(is.het) p2.check  <- !is.na(std.er[(6 + sum(num.betas[1:6])):(sum(num.betas[1:6]) + 5 + num.betas[6])])
    }



    ### drop covariates with large coefficients (if no neg SE)
    ### change to drop cov with large SE (if no neg SE, i.e., check==1)
    if(mean(c(psi.check, th0.check, th1.check, gam.check, eps.check, p1.check, p2.check, 1))==1) {
      if(nuisance == FALSE){
        if(num.betas[1]>1) {
          psi.check 	<- abs(coefs[2:(num.betas[1] + 1)]) < large
        }
        if(num.betas[4]>1) {
          gam.check 	<- abs(coefs[(5 + sum(num.betas[1:3])):(sum(num.betas[1:4]) + 4)]) < large
        }
        if(num.betas[5]>1) {
          eps.check <- abs(coefs[(6 + sum(num.betas[1:4])):(sum(num.betas[1:5]) + 5)]) < large
        }
      }

      if(num.betas[2]>1) {
        th0.check 	<- abs(coefs[(3 + num.betas[1]):(sum(num.betas[1:2]) + 2)]) < large
      }
      if(num.betas[3]>1) {
        th1.check 	<- abs(coefs[(4 + sum(num.betas[1:2])):(sum(num.betas[1:3]) + 3)]) < large
      }
      if(num.betas[6]>1) {
        p1.check <- p2.check <- abs(coefs[(7 + sum(num.betas[1:5])):(sum(num.betas[1:6]) + 5)]) < large
        if(is.het)	p2.check <- abs(coefs[(7 + sum(num.betas[1:6])):(sum(num.betas[1:6]) + 5 + num.betas[6])])<large
      }
    }	# end the if statement



    ### drop covariates with the wrong sign (if no neg SE)
    if(check.quad){
      if(nuisance == FALSE){
        if(mean(c(psi.check, th0.check, th1.check, gam.check, eps.check, p1.check, p2.check, 1))==1) {
          if(num.betas[1] > 1) {
            quads <- grep("psi.sq", betas)
            quads.coord <- quads[(length(quads) - 1):length(quads)]
            quads.clim <- quads[!(quads %in% quads.coord)]
            #psi.check[quads.coord - 1] <- coefs[quads.coord] < 0
            if(strict.psi){
              psi.check[quads.clim - 1] <- coefs[quads.clim] < 0
            }else{
              psi.check[quads.clim - 1] <- coefs[quads.clim]/std.er[quads.clim] < 1.96
            }
          }
          if(num.betas[4] > 1) {
            quads <- grep("gam1.sq", betas)
            quads.coord <- quads[(length(quads) - 1):length(quads)]
            quads.clim <- quads[!(quads %in% quads.coord)]
            #gam.check[quads.coord - 4 - sum(num.betas[1:3])] <- coefs[quads.coord] < 0
            if(strict.gam){
              gam.check[quads.clim - 4 - sum(num.betas[1:3])]	<- coefs[quads.clim] < 0
            }else{
              gam.check[quads.clim - 4 - sum(num.betas[1:3])]	<- coefs[quads.clim]/std.er[quads.clim] < 1.96
            }
          }
          if(num.betas[5] > 1) {
            quads <- grep("eps1.sq", betas)
            quads.coord <- quads[(length(quads) - 1):length(quads)]
            quads.clim <- quads[!(quads %in% quads.coord)]
            #eps.check[quads.coord - 5 - sum(num.betas[1:4])] <- coefs[quads.coord] > 0
            if(strict.gam){
              eps.check[quads.clim - 5 - sum(num.betas[1:4])]	<- coefs[quads.clim] > 0
            }else{
              eps.check[quads.clim - 5 - sum(num.betas[1:4])]	<- coefs[quads.clim]/std.er[quads.clim] > 1.96
            }
          }
        }
      }
    }

    ### if no problem, finish analysis
    if (mean(c(psi.check, th0.check, th1.check, gam.check, eps.check, p1.check, p2.check, 1))==1){
      coef.ok <- 1
    }else{
      coef.ok <- 0
    }


    ### also check if intercepts have neg var
    if(nuisance == FALSE){
      wonky2 <- is.na(sum(std.er))
    }else{
      wonky2 <- 0
    }


    # this checks for significant digits
    jj <- grep('significant digits', pres_out)
    if(length(jj) == 0){
      wonky3 <- 0
    }else{
      num.conv <- pres_out[jj]
      sig.dig <- as.numeric(unlist(strsplit(num.conv, "\\s+"))[2])
      wonky3 <- sig.dig < dig     # TRUE if sig.dig low
    }



    # if wonky result, move it to discard folder
    if (coef.ok == 1 & wonky2 == 0 & wonky3 == 0) {
      ok <- TRUE
    } else { ok <- FALSE}
  }

  ok
}

