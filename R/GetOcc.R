#' GetOccProb
#'
#' Estimate annual probability of occupancy using parameter estimates from top model
#' @param alpha Four letter alpha code for species of interest
#' @param buff_method Optional: Which buffer method to use? Options include: "rec", "bcr"
#' @param buffer Dataframe containing the buffered spp count data
#' @export

GetOccProb <- function(alpha, betas = NULL, buff_method = "rec", Write = TRUE){
  gam_aic <- read.csv(paste0('inst/output/', alpha, '/gam_aic_check.csv'))
  pass <- gam_aic$check[1] == 1

  if(!is.na(pass)){
    if(is.null(betas)) betas <- BBSclim::GetBetas(alpha)

    buffer <- read.csv(paste0('inst/output/', alpha, '/count_buff.csv'))
    years <- seq(from = min(buffer$Year), to = max(buffer$Year))
    climate <- raster.to.array(alpha, years)

    r.psi <- matrix(0, dim(climate)[1], length(years))

    # Initial occupancy prob
    psi.mat <- matrix(betas$psi.betas, nrow = dim(climate)[1], ncol = length(betas$psi.betas), byrow = TRUE)
    psi.cov <- names(betas$psi.betas)
    r.psi[, 1] <- plogis(rowSums(psi.mat * climate[, names(betas$psi.betas), 1]))

    # Extinction/colonization prob
    gam.mat <- matrix(betas$gam.betas, nrow = dim(climate)[1], length(betas$gam.betas), byrow=T)
    gam.cov <- names(betas$gam.betas)
    eps.mat <- matrix(betas$eps.betas, nrow = dim(climate)[1], length(betas$eps.betas), byrow=T)
    eps.cov <- names(betas$eps.betas)
    for (yy in 2:length(years)) {
      r.gam <- matrix(plogis(rowSums(gam.mat * climate[, names(betas$gam.betas), yy])), dim(climate)[1], 1)  #  real colonization for each site
      r.eps <- matrix(plogis(rowSums(eps.mat * climate[, names(betas$eps.betas), yy])), dim(climate)[1], 1)  #  real extinction for each site
      #   compute psi for years 2 ... years
      r.psi[, yy]<- r.psi[, yy - 1] * (1 - r.eps) + (1 - r.psi[, yy - 1]) * r.gam
    }

    all.values <- raster::getValues(NA_biovars$biovars1997[[1]])   # extract climate values
    coord.all <- sp::coordinates(NA_biovars$biovars1997)
    coord     <- coord.all[!is.na(all.values),]


    psi.df <- as.data.frame(r.psi)
    names(psi.df) <- as.character(years)
    psi.df$lat <- coord[, "y"]
    psi.df$lon <- coord[, "x"]

    psi.df2 <- tidyr::gather(psi.df, -lat, -lon, key = "Year", value = "Prob")

    if(Write){
      if(is.null(buff_method)){
        write.csv(psi.df2, file = paste0("inst/output/", alpha, "/occ.csv"), row.names = FALSE)
      }else{
        psi.df3 <- dplyr::filter(psi.df2, lat < max(buffer$Latitude) & lat > min(buffer$Latitude) &
                                   lon > min(buffer$Longitude) & lon < max(buffer$Longitude))
        write.csv(psi.df3, file = paste0("inst/output/", alpha, "/occ.csv"), row.names = FALSE)
      }
    }else{
      if(is.null(buff_method)){
        psi.df2
      }else{
        psi.df3 <- dplyr::filter(psi.df2, lat < max(buffer$Latitude) & lat > min(buffer$Latitude) &
                                   lon > min(buffer$Longitude) & lon < max(buffer$Longitude))
        psi.df3
      }
    }
  }else{
    psi.df <- data.frame(occ = NA)
    write.csv(psi.df, file = paste0("inst/output/", alpha, "/occ.csv"), row.names = FALSE)
  }
}

#' raster.to.array
#'
#' Convert bioclim rasters in array w/ dim n.cell x n.vars x n.years

raster.to.array <- function(alpha, years) {
    index <- c(1,2,8,12,18)
    scale.values <- read.csv(paste0("inst/output/", alpha, "/clim_scale.csv"))

    xy.values <- read.csv(paste0("inst/output/", alpha, "/route_clim.csv"))
    mu.lat <- mean(xy.values$Latitude)
    sd.lat <- sd(xy.values$Latitude)

    mu.lon <- mean(xy.values$Longitude)
    sd.lon <- sd(xy.values$Longitude)

    for (ii in seq_along(years)) {

      # get climate data within masked area
      for (jj in seq_along(index)) {
        assign(paste0("bio.", index[jj]), NA_biovars[[paste0("biovars",as.character(years[ii]))]][[index[jj]]])    # mask the climate data
        assign(paste0("all.values.",index[jj]), raster::getValues(get(paste0("bio.",index[jj]))))                                                 # extract climate values
        assign(paste0("values.",index[jj]), get(paste0("all.values.",index[jj]))[!is.na(get(paste0("all.values.",index[jj])))])                   # remove NAs
        assign(paste0("center.values.",index[jj]), (get(paste0("values.",index[jj])) - scale.values[jj,"mean"])/scale.values[jj,"sd"])     # center and scale
        if (ii+jj==2)  climate <- array(1, dim=c(length(get(paste0("values.",index[jj]))), 15, length(years)))  # first dim is num of sites, 11 is for intercept + 10 covariates
        # alternate linear and quadratic terms, i.e., tmp, sqtmp, dtr, sqdtr, etc
        climate[,jj*2,ii] <- get(paste0("center.values.",index[jj]))     # must be same order as params, need intercept!
        climate[,jj*2+1,ii] <- get(paste0("center.values.",index[jj]))^2     # must be same order as params, need intercept!
      }  # end jj loop

      ## Scale lat/long based on mean/sd used to fit model
      xy <- as.data.frame(raster::rasterToPoints(bio.1))
      sy <- (xy$y - mu.lat)/sd.lat
      sx <- (xy$x - mu.lon)/sd.lon
      climate[, length(index)*2 + 2, ii] <- sx
      climate[, length(index)*2 + 3, ii] <- sx^2
      climate[, length(index)*2 + 4, ii] <- sy
      climate[, length(index)*2 + 5, ii] <- sy^2

    } # end ii loop
    dimnames(climate) <- list(1:(dim(climate)[1]),
                              c("int","tmp","sq_tmp","dtr","sq_dtr","Twet","sq_Twet","Prec","sq_Prec","Pwarm","sq_Pwarm", "Lat", "sq_Lat", "Lon", "sq_Lon"),
                              c(years))
   climate

  }   # end function


#' get.betas
#'
#' Retrieve beta values and var-cov matrix from top model
#' @export

GetBetas <- function(alpha) {

  # read top model and pull out betas and vc matrix
  top.model.out <- scan(paste0("inst/output/", alpha, "/top_mod.out"), what = 'character', sep = '\n', quiet = T)

  jj <- grep('std.error', top.model.out)
  jj.end <- grep('Variance-Covariance Matrix of Untransformed', top.model.out)
  betas <- top.model.out[(jj + 1):(jj.end - 1)]

  psi.betas <- betas[grep('psi', betas)]
  loc.per <- regexpr("psi", psi.betas)
  psi.names <- substr(psi.betas, loc.per + 4, loc.per + 11)
  psi.names <- gsub("[0-9]", "", psi.names)[-1]
  psi.names <- gsub("_$", "", psi.names)
  psi.betas <- as.numeric(substr(psi.betas, 41, 50))
  names(psi.betas) <- trimws(c("int", psi.names))

  gam.betas <- betas[grep('gam', betas)]
  loc.per <- regexpr("gam", gam.betas)
  gam.names <- substr(gam.betas, loc.per + 5, loc.per + 12)
  gam.names <- gsub("[0-9]","",gam.names)[-1]
  gam.names <- gsub("_$", "", gam.names)
  gam.betas <- as.numeric(substr(gam.betas, 41, 50))
  names(gam.betas) <- trimws(c("int", gam.names))

  eps.betas <- betas[grep('eps', betas)]
  loc.per <- regexpr("eps", eps.betas)
  eps.names <- substr(eps.betas, loc.per + 5, loc.per + 12)
  eps.names <- gsub("[0-9]","",eps.names)[-1]
  eps.names <- gsub("_$", "", eps.names)
  eps.betas <- as.numeric(substr(eps.betas, 41, 50))
  names(eps.betas) <- trimws(c("int", eps.names))

  jj <- grep('Variance-Covariance Matrix of Untransformed', top.model.out)
  jj.end <- grep('Individual Site estimates of <psi>', top.model.out)
  raw.vc <- top.model.out[(jj + 2):(jj.end - 2)]
  raw.vc2 <- (strsplit(raw.vc, " +"))
  raw.vc3 <- do.call("rbind", raw.vc2)
  first.dim <- dim(raw.vc3)[1]
  vc.mat <- matrix(as.numeric(raw.vc3[, -(1:2)]), nrow=first.dim)

  rownames(vc.mat) <- raw.vc3[,2]
  col.nms <- top.model.out[(jj + 1)]
  col.nms <- unlist(strsplit(col.nms, " +"))
  colnames(vc.mat) <- col.nms[-1]

  # the betas.vc$vc includes all params. need just psi, gam, eps, and not theta, p
  gam.eps <- grep("B|C", colnames(vc.mat))
  vc <- vc.mat[c(1:length(psi.betas), gam.eps), c(1:length(psi.betas), gam.eps)]

  return(list(psi.betas = psi.betas, gam.betas = gam.betas, eps.betas = eps.betas,
              vc.mat = vc))
}
