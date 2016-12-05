#' RunPsiMods
#'
#' Runs full model set for 31 psi models and writes output to ~/pres/psi
#' @export

RunPsiMods <- function(pao, alpha){
    for(i in 1:length(psi_mods)){
      spp_dm <- GetDM(pao = pao, cov_list = psi_mods[[i]])

      modname <- paste0(alpha, '_psi_model_', i)

      write_dm_and_run2(pao = pao, cov_list = psi_mods[[i]], het = TRUE, dm_list = spp_dm,
                        modname = modname, fixed = TRUE, out = "psi",
                        inits = TRUE, maxfn = '35000 lmt=5', alpha = alpha)
    }
}


#' RunGamMods
#'
#' Runs full model set for 961 gamma/epsilon models and writes output to ~/pres/gam_eps
#' @param gam_mods list containing the covariates for each model
#' @export


RunGamMods <- function(gam_mods, pao, alpha){
  for(i in 1:length(psi_mods)){
    spp_dm <- GetDM(pao = pao, cov_list = gam_mods[[i]])

    modname <- paste0(alpha, '_psi_model_', i)

    write_dm_and_run2(pao = pao, cov_list = gam_mods[[i]], het = TRUE, dm_list = spp_dm,
                      modname = modname, fixed = TRUE, out = "gam_eps",
                      inits = TRUE, maxfn = '35000 lmt=5', alpha = alpha)
  }
}
