#' GetGamMods
#'
#' Returns list with full model set for 961 gamma/epsilon models
#' @return A list containing the model covariates for each of the 961 potential gamma/epsilon models
#' @export

GetGamMods <- function(psi_covs){
  p_covs	<- c("Stop", "sq_Stop", "Lat", "sq_Lat", "Lon", "sq_Lon")
  th_covs <- NULL
  all_covs <- c("tmp", "sq_tmp", "dtr", "sq_dtr", "Twet", "sq_Twet", "Prec", "sq_Prec", "Pwarm", "sq_Pwarm")

  n <- length(all_covs)/2

  id <- unlist(
    lapply(1:n,
           function(i)combn((1:n)*2 - 1,i,simplify=F)
    )
    ,recursive=F)

  gam_mods <- sapply(id, function(i) c(all_covs[sort(c(i, i + 1))]))

  gam_mods2 <- list()
  for(i in 1:length(gam_mods)){
    for(j in 1:length(gam_mods)){
      gam_mods2[[length(gam_mods) * (i - 1) + j]] <-  list(psi.cov = psi_covs,
                                                           th0.cov = th_covs, th1.cov = th_covs,
                                                           gam.cov = c(gam_mods[[i]]),
                                                           eps.cov = c(gam_mods[[j]]),
                                                           p1.cov = p_covs)
    }
  }

  gam_mods2[[(length(gam_mods2) + 1)]] <- list(psi.cov = psi_covs,
                                            th0.cov = th_covs, th1.cov = th_covs,
                                            gam.cov = c("Lat", "sq_Lat", "Lon", "sq_Lon"),
                                            eps.cov = c("Lat", "sq_Lat", "Lon", "sq_Lon"),
                                            p1.cov = p_covs)
  gam_mods2
}


#' GetPsiMods
#'
#' Returns list with full model set for 961 gamma/epsilon models
#' @param covs Data.frame containing the climate covariates for gamma & epsilon from the top gamma/epsilon model
#' @return A list containing the model covariates for each of the 31 potential psi models
#' @export


GetPsiMods <- function(){
  all_covs <- c("tmp", "sq_tmp", "dtr", "sq_dtr", "Twet", "sq_Twet", "Prec", "sq_Prec", "Pwarm", "sq_Pwarm")
  p_covs	<- NULL
  th_covs <- NULL

  n <- length(all_covs)/2

  id <- unlist(
    lapply(1:n,
           function(i)combn((1:n)*2 - 1,i,simplify=F)
    )
    ,recursive=F)

  psi_mods <- sapply(id, function(i) c(all_covs[sort(c(i, i + 1))]))

  psi_mods2 <- list()

  for(i in 1:length(psi_mods)){
      psi_mods2[[i]] <-  list(psi.cov = c(psi_mods[[i]]),
                             th0.cov = th_covs, th1.cov = th_covs,
                             p1.cov = p_covs)
  }

  psi_mods2[[length(psi_mods) + 1]] <- list(psi.cov = c("Lat", "sq_Lat", "Lon", "sq_Lon"),
                                            th0.cov = th_covs, th1.cov = th_covs,
                                            p1.cov = p_covs)

  psi_mods2
}
