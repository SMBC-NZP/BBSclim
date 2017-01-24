# library(httr)
# set_config(config(ssl_verifypeer = 0L))
# devtools::install_github('crushing05/rBBS')
# devtools::install_github('crushing05/BBS.tenstop')
# devtools::install_github('crushing05/BBS.fiftystop')
library(rBBS)
library(BBSclim)

# Set global options
global_opts(tenstops = FALSE)

model_opts(psi.test = TRUE, gam.test = TRUE)

### Import raw 10-stop BBS data
bbs <- GetBBS(is.tenstops = FALSE)


### Get count data for species

CreateSpp(alpha)

spp_AOU <- GetAOU(alpha)

GetSppCounts(AOU = spp_AOU, Write = TRUE, path = paste0('inst/output/', alpha))

### Remove outliers

RemoveOutliers(path = paste0('inst/output/', alpha))

# need to save output in a way that # outliers can be added to report

### Add buffer around routes with counts > 0
buffer_BBS(bbs = bbs, alpha = alpha)


### Make map of routes
make_route_plot(alpha)


### Load rasters containing bioclim estimates
data("NA_biovars")


### Extract annual bioclim estimates for Wood Thrush BBS routes
spp_clim <- GetBioVars(alpha = alpha)


### Save occupancy and climate data as .pao file for input into Presence

write_pao(alpha = alpha)

### Test for annual variation in p
gam_mods <- GetGamMods()

spp_pao <- RPresence::read.pao(paste0("inst/output/", alpha, "/pres/pres_in.pao"))

spp_annual <- test_annual(alpha = alpha, pao = spp_pao)

### Run gamma/epsilon models

RunGamMods(alpha, pao = spp_pao, mods = gam_mods)


### Run psi models with top covariates from gamma/epsilon models
spp_gam_covs <- top_covs(alpha)

spp_psi_mods <- GetPsiMods(covs = spp_gam_covs)

RunPsiMods(alpha, pao = spp_pao)


### Goodness-of-fit of top model
spp_top_mod <- top_covs(aic_tab = spp_gam_aic, mods = gam_mods, psi = FALSE)


spp_gof <- gof(aic_tab = spp_gam_aic, mods = gam_mods, covs = spp_pao$unitcov,
               year_seq = years, Tenstops = tenstops, alpha = alpha,
               time = annual, het = het_det, det_hist = spp_pao$det.data)


### Estimate annual occupancy probability for all raster cells in range
spp_betas <- GetBetas(alpha)
GetOccProb(alpha = alpha, betas = spp_betas)


### Estimate annual indices of range dynamics
GetIndices(alpha)


lowa <- GetSummary(alpha)
lowa
