library(raster)

# Read in data and create annual summaries
## Precipitation
  prec_files <- list.files("data-raw/prec/")[!grepl("gz", list.files("data-raw/prec/"))]

  for (ii in 1:length(prec_files)) {
    if(ii==1) {
      prec <- brick(paste("data-raw/prec/", prec_files[ii], sep = ""))
    } else {
      prec2 <- brick(paste("data-raw/prec/",  prec_files[ii], sep = ""))
      prec <- stack(prec, prec2)
    }
  }


## Minimum temperature
  tmn_files <- list.files("data-raw/tmn/")[!grepl("gz", list.files("data-raw/tmn/"))]

  for (ii in 1:length(tmn_files)) {
    if(ii==1) {
      tmn <- brick(paste("data-raw/tmn/", tmn_files[ii], sep = ""))
    } else {
      tmn2 <- brick(paste("data-raw/tmn/",  tmn_files[ii], sep = ""))
      tmn <- stack(tmn, tmn2)
    }
  }



## Maximum temperature
  tmx_files <- list.files("data-raw/tmx/")[!grepl("gz", list.files("data-raw/tmx/"))]

  for (ii in 1:length(tmx_files)) {
    if(ii==1) {
      tmx <- brick(paste("data-raw/tmx/", tmx_files[ii], sep = ""))
    } else {
      tmx2 <- brick(paste("data-raw/tmx/",  tmx_files[ii], sep = ""))
      tmx <- stack(tmx, tmx2)
    }
  }


# Set lat/long boundaries
  n.amer		<- extent(-180, -50, 20, 85) # only analyze this area

# crop to North America only
  prec <- crop(prec, n.amer)
  tmn <- crop(tmn, n.amer)
  tmx <- crop(tmx, n.amer)	# cut down data to just N Amer

  start.yr	<- 1991
  start.mo	<- 6		  # 6 is June, etc. Start of 12 mo. year. If 6, then year is June-May.
  end.yr		<- 2014

### Subset annual (June - May) rasters
# find the layer for the start year and month
  first.layer	<- grep(paste0(start.yr, ".", formatC(start.mo, width = 2, format = "d", flag = "0")),   names(tmn))

# assemble monthly tmn data into years
  for (ii in 1:(end.yr-start.yr)) {
    assign(paste("prec",ii+start.yr, sep=""), subset(prec, (first.layer-12+ii*12):(first.layer-1+ii*12)))
    assign(paste("tmin", ii + start.yr, sep = ""), subset(tmn, (first.layer - 12 + ii * 12):(first.layer - 1 + ii * 12)))
    assign(paste("tmax",ii+start.yr, sep=""), subset(tmx, (first.layer-12+ii*12):(first.layer-1+ii*12)))
  }


# Run biovars function to create summary climate layers ----
library(dismo)
  start.bbs <- 1996 # Should be 1 year before first BBS year

  for (ii in 1:(end.yr - start.bbs)) {
    assign(paste("biovars", ii + start.bbs, sep=""),
           biovars(get(paste("prec", ii + start.bbs, sep="")),
                   get(paste("tmin", ii + start.bbs, sep="")),
                   get(paste("tmax", ii + start.bbs, sep=""))))
  }

  biovar_files <- paste("biovars", seq(from = start.bbs + 1, to = end.yr), sep = "")

  vars <- names(get(biovar_files[1]))

  NA_biovars <- lapply(biovar_files, get)
  names(NA_biovars) <- biovar_files

# BCR shapefile
  if (!require(gpclib)) install.packages("gpclib", type="source")
  maptools::gpclibPermit()
  new.crs <- maptools::CRS("+proj=longlat +datum=WGS84")
  bcr.shapefile <- rgdal::readShapePoly("data-raw/bcr/BCR", verbose=TRUE, proj4string=new.crs)
  bcr.shapefile@data$id <- rownames(bcr.shapefile@data)
  bcr.points <- ggplot2::fortify(bcr.shapefile, region = "id")
  bcr <- merge(bcr.points, bcr.shapefile@data, by = "id")code_lookup <- read.csv("data-raw/code_lookup.csv")

# Alpha code look-up table
  code_lookup <- read.csv("data-raw/code_lookup.csv")
  code_lookup <- dplyr::select(code_lookup, -group)


  devtools::use_data(NA_biovars, bcr, code_lookup, overwrite = TRUE)

