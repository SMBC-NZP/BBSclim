# library(httr)
# set_config(config(ssl_verifypeer = 0L))
# devtools::install_github('crushing05/rBBS')
library(rBBS)
library(BBSclim)

tenstops <- TRUE

### Import raw 10-stop BBS data
bbs <- GetRouteData(TenStops = tenstops)

### Import weather data
weather <- GetWeather()

### Import route information
routes <- GetRoutes()

het_det <- TRUE
annual <- TRUE

### Get count data for species
alpha <- "lowa"

spp_AOU <- GetAOU(alpha)

spp_counts <- GetSppCounts(AOU = spp_AOU, Write = TRUE, path = 'inst/output/spp_counts')

### Remove outliers
# 1) Compute knn matrix (for each point, mean distance to knn)
# 2) Estimate mean and sd of knn matrix
# 3) remove routes with knn > mu.knn + 3sd

### Read raw BBS counts
spp_counts2 <- dplyr::filter(spp_counts, Longitude > -110)
spp_buff <- buffer_BBS(spp_count = spp_counts2)


### Load rasters containing bioclim estimates
data("NA_biovars")


### Extract annual bioclim estimates for Wood Thrush BBS routes
spp_clim <- GetBioVars(counts = spp_buff)


### Save occupancy and climate data as .pao file for input into Presence

write_pao(counts = spp_buff, covs = spp_clim, alpha = alpha, TenStops = tenstops)


### Run psi models
Test <- TRUE

psi_mods <- GetPsiMods()

spp_pao <- RPresence::read.pao(paste0("inst/output/pao/", alpha, ".pao"))


spp_psi_aic <- RunPsiMods(pao = spp_pao, mods = psi_mods, alpha = alpha, test = Test,
                          time = annual, het = het_det, del = FALSE)


### Run gam/eps models with top covariates from psi models
spp_psi_covs <- top_covs(aic_tab = spp_psi_aic, mods = psi_mods)

gam_mods <- GetGamMods(psi_covs = spp_psi_covs)


spp_gam_aic <- RunPsiMods(pao = spp_pao, mods = gam_mods, alpha = alpha, test = Test,
                          time = annual, het = het_det, del = FALSE)




### Move top model output to 'top/' and delete others
# spp_top_mod <- top_covs(aic_tab = spp_gam_aic, mods = gam_mods, psi = FALSE)
#
#
#
# ### Goodness-of-fit of top model
# aic_temp <- dplyr::data_frame(Model = "mod1", Model_num = 1, LogLik = 259231, nParam = 2,
#                               AIC = 2165464)
#
# aic_temp2 <- dplyr::data_frame(Model = "mod2", Model_num = 2, LogLik = 252346, nParam = 3,
#                               AIC = 2161654)
# aic_tab <- dplyr::bind_rows(aic_temp, aic_temp2)
# aic_tab <- dplyr::mutate(aic_tab, delta_AIC = AIC - min(AIC))
# aic_tab <- dplyr::arrange(aic_tab, delta_AIC)
#
# spp_top_mod <- top_mod(aic_tab = aic_tab, mods = gam_mods, psi = FALSE)

spp_gam_aic <- RunGamMods(pao = spp_pao, mods = gam_mods, alpha = alpha, test = Test,
                          time = annual, het = het_det, del = FALSE)




### Move top model output to 'top/' and delete others
# spp_top_mod <- top_covs(aic_tab = spp_gam_aic, mods = gam_mods, psi = FALSE)
#
#
#
# ### Goodness-of-fit of top model
aic_temp <- dplyr::data_frame(Model = "lowa_psi_model_3", Model_num = 1, LogLik = 259231, nParam = 2,
                              AIC = 2165464)

aic_temp2 <- dplyr::data_frame(Model = "lowa_psi_model_4", Model_num = 2, LogLik = 252346, nParam = 3,
                              AIC = 2161654)
aic_tab <- dplyr::bind_rows(aic_temp, aic_temp2)
aic_tab <- dplyr::mutate(aic_tab, delta_AIC = AIC - min(AIC))
aic_tab <- dplyr::arrange(aic_tab, delta_AIC)

spp_top_mod <- top_mod(aic_tab = aic_tab, mods = gam_mods, psi = FALSE)





