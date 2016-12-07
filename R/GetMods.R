#' GetPsiMods
#'
#' Returns list with full model set for 31 psi models
#' @return A list containing the model covariates for each of the 31 potential psi models
#' @export

GetPsiMods <- function(){
  p_covs	<- c("Stop", "sq_Stop", "Lat", "sq_Lat", "Lon", "sq_Lon")
  th_covs <- c("Lat", "sq_Lat", "Lon", "sq_Lon")
  psi_covs <- c("tmp", "sq_tmp", "dtr", "sq_dtr", "Twet", "sq_Twet", "Prec", "sq_Prec", "Pwarm", "sq_Pwarm")

  n <- length(psi_covs)/2

  id <- unlist(
    lapply(1:n,
           function(i)combn((1:n)*2 - 1,i,simplify=F)
    )
    ,recursive=F)

  psi_mods <- sapply(id, function(i) c(psi_covs[sort(c(i, i + 1))]))

  psi_mods2 <- lapply(psi_mods, function(x) list(psi.cov = x,
                                                 th0.cov = th_covs, th1.cov = th_covs,
                                                 gam.cov = psi_covs, eps.cov = psi_covs,
                                                 p1.cov = p_covs))
  psi_mods2
}


#' GetGamMods
#'
#' Returns list with full model set for 961 gamma/epsilon models
#' @param psi_cov Vector containing the climate covariates from the top psi model
#' @return A list containing the model covariates for each of the 961 potential psi models
#' @export


GetGamMods <- function(psi_covs){
  gam_covs <- c("tmp", "sq_tmp", "dtr", "sq_dtr", "Twet", "sq_Twet", "Prec", "sq_Prec", "Pwarm", "sq_Pwarm")
  p_covs	<- c("Stop", "sq_Stop", "Lat", "sq_Lat", "Lon", "sq_Lon")
  th_covs <- c("Lat", "sq_Lat", "Lon", "sq_Lon")

  n <- length(gam_covs)/2

  id <- unlist(
    lapply(1:n,
           function(i)combn((1:n)*2 - 1,i,simplify=F)
    )
    ,recursive=F)

  gam_mods <- sapply(id, function(i) c(gam_covs[sort(c(i, i + 1))]))

  gam_mods2 <- list()
  for(i in 1:length(gam_mods)){
    for(j in 1:length(gam_mods)){
      gam_mods2[[length(gam_mods) * (i - 1) + j]] <-  list(psi.cov = psi_covs,
                                                           th0.cov = th_covs, th1.cov = th_covs,
                                                           gam.cov = gam_mods[[i]], eps.cov = gam_mods[[j]],
                                                           p1.cov = p_covs)
    }
  }

  gam_mods2
}
