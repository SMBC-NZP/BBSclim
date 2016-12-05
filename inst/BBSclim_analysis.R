
library(rBBS)
library(RPresence)
library(BBSclim)

### Import raw 10-stop BBS data
bbs <- GetRouteData()

### Import weather data
weather <- GetWeather()

### Import route information
routes <- GetRoutes()

het_det <- TRUE
annual <- TRUE

### Get count data for WOTH
alpha <- "lowa"

spp_AOU <- GetAOU(alpha)

spp_counts <- GetSppCounts(AOU = spp_AOU)


### Read raw BBS counts
spp_counts2 <- dplyr::filter(spp_counts, Longitude > -110)
spp_buff <- buffer_BBS(spp_count = spp_counts2)


### Load rasters containing bioclim estimates
data("NA_biovars")


### Extract annual bioclim estimates for Wood Thrush BBS routes
spp_clim <- GetBioVars(counts = spp_buff)


### Save occupancy and climate data as .pao file for input into Presence

write_pao(counts = spp_buff, clim = spp_clim, alpha = alpha)


### Run psi models

psi_mods <- GetPsiMods()
psi_mods <- psi_mods[1:2]

spp_pao <- read.pao(paste0("inst/output/pao/", alpha, ".pao"))

spp_psi_aic <- RunPsiMods(pao = spp_pao, mods = psi_mods, alpha = alpha, test = FALSE,
                      time = annual, het = het_det)


### Run gam/eps models with top covariates from psi models
spp_psi_covs <- top_psi(aic_tab = spp_psi_aic, mods = psi_mods)

gam_mods <- GetGamMods(psi_covs = spp_psi_covs)
gam_mods <- gam_mods[1:2]

spp_gam_aic <- RunPsiMods(pao = spp_pao, mods = gam_mods, alpha = alpha, test = TRUE,
                          time = annual, het = het_det)
