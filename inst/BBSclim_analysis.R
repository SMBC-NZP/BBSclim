# library(httr)
# set_config(config(ssl_verifypeer = 0L))
# devtools::install_github('crushing05/rBBS')
# devtools::install_github('crushing05/BBS.tenstop')
# devtools::install_github('crushing05/BBS.fiftystop')
library(rBBS)
library(BBSclim)

# Set global options
source("R/glob_opts.R")

### Import raw 10-stop BBS data
bbs <- GetBBS(is.tenstops = tenstops)


### Get count data for species
alpha <- "kewa"

CreateSpp(alpha)

spp_AOU <- GetAOU(alpha)

spp_counts <- GetSppCounts(AOU = spp_AOU, Write = TRUE, path = paste0('inst/output/', alpha))

### Remove outliers

spp_counts2 <- RemoveOutliers(counts = spp_counts)

# need to save output in a way that # outliers can be added to report

### Add buffer around routes with counts > 0
spp_buff <- buffer_BBS(spp_count = spp_counts2, bbs = bbs)


### Make map of routes
make_route_plot(spp_buff, spp_counts2, spp_counts, alpha)


### Load rasters containing bioclim estimates
data("NA_biovars")


### Extract annual bioclim estimates for Wood Thrush BBS routes
spp_clim <- GetBioVars(counts = spp_buff)


### Save occupancy and climate data as .pao file for input into Presence

write_pao(counts = spp_buff, covs = spp_clim, alpha = alpha, TenStops = tenstops)


### Run psi models
psi_mods <- GetPsiMods()

spp_pao <- RPresence::read.pao(paste0("inst/output/", alpha, "/spp.pao"))


spp_psi_aic <- RunPsiMods(pao = spp_pao, mods = psi_mods, alpha = alpha,
                          time = annual, het = het_det, del = FALSE)

# spp_psi_aic <- read.csv("inst/output/aic/psi/lowa.csv")

### Run gam/eps models with top covariates from psi models
spp_psi_covs <- top_covs(aic_tab = spp_psi_aic, mods = psi_mods)

gam_mods <- GetGamMods(psi_covs = spp_psi_covs)

spp_gam_aic <- RunGamMods(pao = spp_pao, mods = gam_mods, alpha = alpha, test = Test,
                          time = annual, het = het_det, del = FALSE)

# spp_gam_aic <- read.csv("inst/output/aic/gam/lowa.csv")

# ### Goodness-of-fit of top model
spp_top_mod <- top_covs(aic_tab = spp_gam_aic, mods = gam_mods, psi = FALSE)


spp_gof <- gof(aic_tab = spp_gam_aic, mods = gam_mods, covs = spp_pao$unitcov,
               year_seq = years, Tenstops = tenstops, alpha = alpha,
               time = annual, het = het_det, det_hist = spp_pao$det.data)


### Estimate annual occupancy probability for all raster cells in range
spp_betas <- GetBetas(alpha)
spp_occ <- GetOccProb(spp_betas, alpha, years, buffer = spp_buff)

ggplot(spp_occ, aes(x = lon, y = lat, fill = Prob)) + geom_raster() + facet_wrap(~Year)

### Estimate annual indices of range dynamics
spp_ind <- GetIndices(prob_df = spp_occ, alpha, years, spp_buff)
