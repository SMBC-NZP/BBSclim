#' drop.bad.coef
#'
#' Drops quadratic term if linear term is bad
#' @param mod Model output file
#' @return TRUE if model converged; FALSE if not
#' @export

drop.bad.coef 	<- function(keep, cov) {
  quads <- grep("sq", cov)  	# quadratic terms
  prior <- quads - 1		# linear term prior to a quadratic term (some linear terms don't have a quad)
  need.to.switch <- which(keep[prior] == FALSE)	# which linear terms are bad
  keep[quads[need.to.switch]] = FALSE		# then make the quad term false too
  cov <- cov[keep]
  cov
}

#' mod_eval
#'
#' Checks for problems with convergence
#' @param mod Model output file
#' @param dig Minimum threshold for significant digits
#' @param large Maximum threshold for beta coefficients
#' @return TRUE if model converged; FALSE if not
#' @export

mod_eval <- function(pres_out, mod, pao2 , dig = 3, large = 8, time, het){
  jj <- grep('std.error', pres_out)
  jj2 <- grep('Individual Site estimates of <psi>', pres_out)

  betas <- pres_out[(jj + 1):(jj2 - 1)]
  coefs <- as.numeric(substr(betas, 41, 50))
  std.er <- as.numeric(substr(betas, 54, 63))

  num.betas <- plyr::ldply(mod, function(x) length(x))$V1

  if(time) num.betas[6] <- num.betas[6] + pao2$nseasons

  psi.keep <- th0.keep <- th1.keep <- gam.keep <- eps.keep <- p1.keep <- p2.keep <- stop.keep <- NULL
  ### drop covariates with negative standard errors
  if(num.betas[1] > 1) {
    psi.keep 	<- !is.na(std.er[2:(num.betas[1] + 1)])
    psi.cov	<- drop.bad.coef(psi.keep, mod$psi.cov)
  }
  if(num.betas[2] > 1) {
    th0.keep 	<- !is.na(std.er[(3 + num.betas[1]):(sum(num.betas[1:2]) + 2)])
    th0.cov	<- drop.bad.coef(th0.keep, mod$th0.cov)
  }
  if(num.betas[3] > 1) {
    th1.keep 	<- !is.na(std.er[(4 + sum(num.betas[1:2])):(sum(num.betas[1:3]) + 3)])
    th1.cov	<- drop.bad.coef(th1.keep, mod$th1.cov)
  }
  if(num.betas[4] > 1) {
    gam.keep 	<- !is.na(std.er[(5 + sum(num.betas[1:3])):(sum(num.betas[1:4]) + 4)])
    gam.cov	<- drop.bad.coef(gam.keep, mod$gam.cov)
  }
  if(num.betas[5] > 1) {
    eps.keep 	<- !is.na(std.er[(6 + sum(num.betas[1:4])):(sum(num.betas[1:5]) + 5)])
    eps.cov	<- drop.bad.coef(eps.keep, mod$eps.cov)
  }
  if(num.betas[6] > 1) {
    p1.keep <- p2.keep  <- !is.na(std.er[(7 + sum(num.betas[1:5])):(sum(num.betas[1:6]) + 5)])
    if(het)      p2.keep  <- !is.na(std.er[(7 + sum(num.betas[1:6])):(sum(num.betas[1:6]) + 5 + num.betas[6])])
    if(time){
      p1.cov	<- drop.bad.coef(p1.keep & p2.keep, c(seq(1:pao2$nseasons) * pao2$nsurveyseason[1] - pao2$nsurveyseason[1] + 1, mod$p1.cov))
    }else{
      p1.cov	<- drop.bad.coef(p1.keep & p2.keep, mod$p1.cov)
    }
  }



  ### drop covariates with large coefficients (if no neg SE)
  ### change to drop cov with large SE (if no neg SE, i.e., keep==1)
  if (mean(c(psi.keep, th0.keep, th1.keep, gam.keep, eps.keep, p1.keep, p2.keep, 1))==1) {
    if(num.betas[1]>1) {
      psi.keep 	<- abs(std.er[2:(num.betas[1] + 1)]) < large
      psi.cov	<- drop.bad.coef(psi.keep, mod$psi.cov)
    }
    if(num.betas[2]>1) {
      th0.keep 	<- abs(std.er[(3 + num.betas[1]):(sum(num.betas[1:2]) + 2)]) < large
      th0.cov	<- drop.bad.coef(th0.keep, mod$th0.cov)
    }
    if(num.betas[3]>1) {
      th1.keep 	<- abs(std.er[(4 + sum(num.betas[1:2])):(sum(num.betas[1:3]) + 3)]) < large
      th1.cov	<- drop.bad.coef(th1.keep, mod$th1.cov)
    }
    if(num.betas[4]>1) {
      gam.keep 	<- abs(std.er[(5 + sum(num.betas[1:3])):(sum(num.betas[1:4]) + 4)]) < large
      gam.cov	<- drop.bad.coef(gam.keep, mod$gam.cov)
    }
    if(num.betas[5]>1) {
      eps.keep <- abs(std.er[(6 + sum(num.betas[1:4])):(sum(num.betas[1:5]) + 5)]) < large
      eps.cov	<- drop.bad.coef(eps.keep, mod$eps.cov)
    }
    if(num.betas[6]>1) {
      p1.keep <- p2.keep <- abs(std.er[(7 + sum(num.betas[1:5])):(sum(num.betas[1:6]) + 5)]) < large
      if(het)	p2.keep   	<- abs(std.er[(7 + sum(num.betas[1:6])):(sum(num.betas[1:6]) + 5 + num.betas[6])])<large
      if(time){
        p1.cov	<- drop.bad.coef(p1.keep & p2.keep, c(seq(1:pao2$nseasons) * pao2$nsurveyseason[1] - pao2$nsurveyseason[1] + 1, mod$p1.cov))
      }else{
        p1.cov	<- drop.bad.coef(p1.keep & p2.keep, mod$p1.cov)
      }
    }
  }	# end the if statement



  ### drop covariates with the wrong sign (if no neg SE)
  # before, I only dropped signif coefficients, but now drop all
  if (mean(c(psi.keep, th0.keep, th1.keep, gam.keep, eps.keep, p1.keep, p2.keep, 1))==1) {
    if(num.betas[1]>1) {
      quads <- grep("psi.sq", betas)
      psi.keep[quads - 1] 	<- coefs[quads] < 0	# '-1' is for intercept
      psi.cov	<- drop.bad.coef(psi.keep, mod$psi.cov)
    }
    if(num.betas[4]>1) {
      quads <- grep("gam1.sq", betas)
      gam.keep[quads - 4 - sum(num.betas[1:3])] 	<- coefs[quads] < 0  # quads indicates position in all betas, but gam.keep is just for gam
      gam.cov	<- drop.bad.coef(gam.keep, mod$gam.cov)
    }
    if(num.betas[5]>1) {
      quads <- grep("eps1.sq",betas)
      eps.keep[quads - 5 - sum(num.betas[1:4])] 	<- coefs[quads] > 0
      eps.cov	<- drop.bad.coef(eps.keep, mod$eps.cov)
    }

  }	# end the if statement

  ### if no problem, finish analysis
  if (mean(c(psi.keep, th0.keep, th1.keep, gam.keep, eps.keep, p1.keep, p2.keep, 1))==1){
    coef.ok <- 1
  }else{
    coef.ok <- 0
  }


  ### also check if intercepts have neg var
  wonky2 <- is.na(sum(std.er))

  # this checks for significant digits
  jj <- grep('significant digits', pres_out)
  num.conv <- pres_out[jj]
  sig.dig <- as.numeric(unlist(strsplit(num.conv, "\\s+"))[2])
  wonky3 <- sig.dig < dig     # TRUE if sig.dig low

  # if wonky result, move it to discard folder
  if (coef.ok == 1 & wonky2 == 0 & wonky3 == 0) {
    ok <- TRUE
  } else { ok <- FALSE}

  ok
}

