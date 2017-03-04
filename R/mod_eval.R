#' mod_eval
#'
#' Checks for problems with convergence
#' @param mod Model output file
#' @param dig Minimum threshold for significant digits
#' @param large Maximum threshold for beta coefficients
#' @return TRUE if model converged; FALSE if not
#' @export

mod_eval <- function(pres_out, mod, pao2 , dig = 2, large = 8, is.annual, is.het, strict.psi = FALSE, strict.gam = FALSE){
  jj <- grep('std.error', pres_out)
  jj2 <- grep('Variance-Covariance Matrix of Untransformed', pres_out)

  betas <- pres_out[(jj + 1):(jj2 - 1)]
  coefs <- as.numeric(substr(betas, 41, 50))
  std.er <- as.numeric(substr(betas, 54, 63))

  num.betas <- plyr::ldply(mod, function(x) length(x))$V1

  if(is.annual) num.betas[6] <- num.betas[6] + pao2$nseasons

  psi.check <- th0.check <- th1.check <- gam.check <- eps.check <- p1.check <- p2.check <- stop.check <- NULL
  ### drop covariates with negative standard errors
  if(num.betas[1] > 1) {
    psi.check 	<- !is.na(std.er[2:(num.betas[1] + 1)])
  }
  if(num.betas[2] > 1) {
    th0.check 	<- !is.na(std.er[(3 + num.betas[1]):(sum(num.betas[1:2]) + 2)])
  }
  if(num.betas[3] > 1) {
    th1.check 	<- !is.na(std.er[(4 + sum(num.betas[1:2])):(sum(num.betas[1:3]) + 3)])
  }
  if(num.betas[4] > 1) {
    gam.check 	<- !is.na(std.er[(5 + sum(num.betas[1:3])):(sum(num.betas[1:4]) + 4)])
  }
  if(num.betas[5] > 1) {
    eps.check 	<- !is.na(std.er[(6 + sum(num.betas[1:4])):(sum(num.betas[1:5]) + 5)])
  }
  if(num.betas[6] > 1) {
    p1.check <- p2.check  <- !is.na(std.er[(6 + sum(num.betas[1:5])):(sum(num.betas[1:6]) + 5)])
    if(is.het) p2.check  <- !is.na(std.er[(6 + sum(num.betas[1:6])):(sum(num.betas[1:6]) + 5 + num.betas[6])])
  }



  ### drop covariates with large coefficients (if no neg SE)
  ### change to drop cov with large SE (if no neg SE, i.e., check==1)
  if (mean(c(psi.check, th0.check, th1.check, gam.check, eps.check, p1.check, p2.check, 1))==1) {
    if(num.betas[1]>1) {
      psi.check 	<- abs(std.er[2:(num.betas[1] + 1)]) < large
    }
    if(num.betas[2]>1) {
      th0.check 	<- abs(std.er[(3 + num.betas[1]):(sum(num.betas[1:2]) + 2)]) < large
    }
    if(num.betas[3]>1) {
      th1.check 	<- abs(std.er[(4 + sum(num.betas[1:2])):(sum(num.betas[1:3]) + 3)]) < large
    }
    if(num.betas[4]>1) {
      gam.check 	<- abs(std.er[(5 + sum(num.betas[1:3])):(sum(num.betas[1:4]) + 4)]) < large
    }
    if(num.betas[5]>1) {
      eps.check <- abs(std.er[(6 + sum(num.betas[1:4])):(sum(num.betas[1:5]) + 5)]) < large
    }
    if(num.betas[6]>1) {
      p1.check <- p2.check <- abs(std.er[(7 + sum(num.betas[1:5])):(sum(num.betas[1:6]) + 5)]) < large
      if(is.het)	p2.check   	<- abs(std.er[(7 + sum(num.betas[1:6])):(sum(num.betas[1:6]) + 5 + num.betas[6])])<large
    }
  }	# end the if statement



  ### drop covariates with the wrong sign (if no neg SE)
  # before, I only dropped signif coefficients, but now drop all
  if (mean(c(psi.check, th0.check, th1.check, gam.check, eps.check, p1.check, p2.check, 1))==1) {
    if(num.betas[1]>1) {
      quads <- grep("psi.sq", betas)
      if(strict.psi){
        psi.check[quads - 1]	<- coefs[quads] < 0
      }
    }
    if(num.betas[4]>1) {
      quads <- grep("gam1.sq", betas)
      if(strict.gam){
        gam.check[quads - 4 - sum(num.betas[1:3])]	<- coefs[quads] < 0
      }
    }
    if(num.betas[5]>1) {
      quads <- grep("eps1.sq",betas)
      if(strict.gam){
        eps.check[quads - 5 - sum(num.betas[1:4])]	<- coefs[quads] > 0
      }
    }
  }	# end the if statement

  ### if no problem, finish analysis
  if (mean(c(psi.check, th0.check, th1.check, gam.check, eps.check, p1.check, p2.check, 1))==1){
    coef.ok <- 1
  }else{
    coef.ok <- 0
  }


  ### also check if intercepts have neg var
  wonky2 <- is.na(sum(std.er))

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

  ok
}

