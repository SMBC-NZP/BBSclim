#' PlotRoutes
#'
#' Plot map with BBS route locations
#' @param alpha 4-letter alpha code for species of interest
#' @param psi Should occupancy probability (in year 1) be plotted? (*doesn't work yet*)
#' @export

PlotRoutes <- function(alpha){
  buff_routes <- read.csv(here::here(paste0('inst/output/', alpha, '/count_buff.csv')))
  raw_routes <- read.csv(here::here(paste0('inst/output/', alpha, '/raw_counts.csv')))
  used_routes <- read.csv(here::here(paste0('inst/output/', alpha, '/no_outlier_counts.csv')))

  usa <- ggplot2::map_data("state")
  canada <- ggplot2::map_data("worldHires", "Canada")
  mexico <- ggplot2::map_data("worldHires", "Mexico")

  spp_buff2 <- dplyr::distinct(buff_routes, routeID, .keep_all = TRUE)
  spp_rts <- dplyr::distinct(used_routes, routeID, .keep_all = TRUE)
  spp_out <- dplyr::anti_join(raw_routes, used_routes)
  spp_out <- dplyr::distinct(spp_out, routeID, .keep_all = TRUE)

  xmin <- min(min(spp_buff2$Longitude), min(spp_out$Longitude)) - 2
  xmax <- max(max(spp_buff2$Longitude), max(spp_out$Longitude)) + 2

  ymin <- min(min(spp_buff2$Latitude), min(spp_out$Latitude)) - 2
  ymax <- max(max(spp_buff2$Latitude), max(spp_out$Latitude)) + 2


  p <- ggplot() + theme(panel.background = element_rect(fill = "#CDD2D4", color = NA), axis.title = element_blank())
  p <- p + coord_map(projection = "lambert", lat0 = 25, lat1 = 50, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  p <- p + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "#F0F0F1",  color="#909090") +
    geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "#F0F0F1", color="#909090") +
    geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = "#F0F0F1", color="#909090")
  p <- p + geom_point(data = spp_buff2, aes(x = Longitude, y = Latitude), size = 0.75, color = "#777777", alpha = 0.9)
  p <- p + geom_point(data = spp_rts, aes(x = Longitude, y = Latitude), size = 1.4, color = "#dc322f")
  p <- p + geom_point(data = spp_out, aes(x = Longitude, y = Latitude), size = 5, shape = 21, color = "#778899")
  p <- p + geom_point(data = spp_out, aes(x = Longitude, y = Latitude), size = 1.4, color = "#dc322f")
  p <- p + scale_x_continuous(breaks = seq(from = 10*(xmin%/%10 + as.logical(xmin%%10)),
                                           to = 10*(xmax%/%10 + as.logical(xmax%%10)), by = 10))
  p <- p + scale_y_continuous(breaks = seq(from = 5*(ymin%/%5 + as.logical(ymin%%5)),
                                           to = 5*(ymax%/%5 + as.logical(ymax%%5)), by = 10))
  p
}


#' PlotLat
#'
#' Plot graph of changes in latitude indices
#' @param alpha 4-letter alpha code for species of interest
#' @export

PlotLat <- function(alpha, ci = FALSE){
  indices <- read.csv(here::here(paste0('inst/output/', alpha, '/indices.csv')))
  lat.indices <- dplyr::filter(indices, ind == "n.lat" | ind == "s.lat" |
                               ind == "n.core" | ind == "s.core" | ind == "avg.lat")


  p <- ggplot(lat.indices, aes(x = Year, y = value, group = ind, color = ind))
  p <- p + geom_line(aes(linetype = ind))
  p <- p + scale_y_continuous("Latitude")
  p <- p + scale_linetype_manual(values = c("solid", "dashed", "longdash", "dashed", "longdash"))
  p <- p + scale_color_manual(values = c("black", "grey35", "grey50", "grey35", "grey50"))
  p <- p + scale_x_continuous(breaks = seq(from = min(lat.indices$Year),
                                           to = max(lat.indices$Year),
                                           by = 2))
  p <- p + theme(legend.position = "none")

  if(ci){
    p <- p + geom_ribbon(aes(ymin = (value - 1.96 * sd.err), ymax = (value + 1.96 * sd.err)), alpha = 0.2)
    p
  }else{
    p
  }
}

#' PlotLon
#'
#' Plot graph of changes in longitude indices
#' @param alpha 4-letter alpha code for species of interest
#' @export

PlotLon <- function(alpha, ci = TRUE){
  indices <- read.csv(here::here(paste0('inst/output/', alpha, '/indices.csv')))
  lon.indices <- dplyr::filter(indices, ind == "avg.lon")

  p <- ggplot(lon.indices, aes(x = Year, y = value, group = ind, color = ind))
  p <- p + geom_line(aes(linetype = ind))
  p <- p + scale_y_continuous("(West)        Longitude        (East)")
  p <- p + scale_linetype_manual(values = c("solid", "dashed", "longdash", "dashed", "longdash"))
  p <- p + scale_color_manual(values = c("black", "grey35", "grey50", "grey35", "grey50"))
  p <- p + scale_x_continuous(breaks = seq(from = min(lon.indices$Year),
                                           to = max(lon.indices$Year),
                                           by = 2))
  p <- p + theme(legend.position = "none")

  if(ci){
    p <- p + geom_ribbon(aes(ymin = (value - 1.96 * sd.err), ymax = (value + 1.96 * sd.err)), alpha = 0.2)
    p
  }else{
    p
  }
}

#' PlotPsi
#'
#' Plot graph of changes in range-wide occupancy
#' @param alpha 4-letter alpha code for species of interest
#' @export

PlotPsi <- function(alpha){
  indices <- read.csv(here::here(paste0('inst/output/', alpha, '/indices.csv')))
  psi.indices <- dplyr::filter(indices, ind == "avg.psi")
  y.min <- 0
  y.max <- plyr::round_any(max(psi.indices$value + 1.96 * psi.indices$sd.err), 0.1, f = ceiling)

  p <- ggplot(psi.indices, aes(x = Year, y = value, group = ind))
  p <- p + geom_point()
  p <- p + geom_line()
  p <- p + geom_point(color = "white", size = 6)
  p <- p + geom_point(size = 4)
  p <- p + geom_ribbon(aes(ymin = (value - 1.96 * sd.err), ymax = (value + 1.96 * sd.err)), alpha = 0.2)
  p <- p + scale_y_continuous("Occupancy", limits = c(y.min, y.max))
  p <- p + scale_x_continuous(breaks = seq(from = min(psi.indices$Year),
                                           to = max(psi.indices$Year),
                                           by = 2))
  p
}

#' MapPsi
#'
#' Annual occupancy probability maps
#' @param alpha 4-letter alpha code for species of interest
#' @param proj Should maps be projected? (Probably but much slower)
#' @export


MapPsi <- function(alpha, proj = TRUE, show.points = FALSE){
  psi <- read.csv(here::here(paste0('inst/output/', alpha, '/occ.csv')))
  indices <- read.csv(here::here(paste0('inst/output/', alpha, '/indices.csv')))
  limits <- dplyr::filter(indices, ind %in% c("s.lat", "n.lat"))
  core <- dplyr::filter(indices, ind %in% c("s.core", "n.core"))
  center <- dplyr::filter(indices, ind %in% c("avg.lat"))

  xmin <- min(psi$lon) - 1
  xmax <- max(psi$lon) + 1

  ymin <- min(psi$lat) - 1
  ymax <- max(psi$lat) + 1

  if(show.points) routes <- read.csv(here::here(paste0('inst/output/', alpha, '/count_buff.csv')))

  if(proj){
    p <- ggplot() + facet_wrap(~Year)
    p <- p + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "#F0F0F1",  color="#909090") +
      geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "#F0F0F1", color="#909090") +
      geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = "#F0F0F1", color="#909090")
    p <- p + geom_tile(data = psi, aes(x = lon, y = lat, fill = Prob))
    p <- p + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA,  color="#909090") +
      geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = NA, color="#909090") +
      geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = NA, color="#909090")
    p <- p + coord_map(projection = "lambert", lat0 = 25, lat1 = 50, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
    p <- p + geom_hline(data = limits, aes(yintercept = value, group = ind), linetype = "dashed", color = "grey70")
    p <- p + geom_hline(data = core, aes(yintercept = value, group = ind), linetype = "longdash", color = "grey70")
    p <- p + geom_hline(data = center, aes(yintercept = value), color = "grey70")
    if(show.points){
     p <- p + geom_point(data = routes, aes(x = Longitude, y = Latitude, size = speciestotal), color = "black", alpha = 0.5)
     p <- p + scale_size_continuous(guide = FALSE)
    }
    p <- p + scale_fill_continuous(low = "#F0F0F1", high = "#ef6e0b")
    p <- p + scale_x_continuous("Longitude",breaks = seq(from = 10*(xmin%/%10 + as.logical(xmin%%10)),
                                                         to = 10*(xmax%/%10 + as.logical(xmax%%10)), by = 10))
    p <- p + scale_y_continuous("Latitude", breaks = seq(from = 5*(ymin%/%5 + as.logical(ymin%%5)),
                                                         to = 5*(ymax%/%5 + as.logical(ymax%%5)), by = 10))
    p <- p + theme(panel.background = element_rect(fill = "#CDD2D4", color = NA), axis.title = element_blank())
    p
  }else{
    p <- ggplot()
    p <- p + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "#F0F0F1",  color="#909090") +
      geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "#F0F0F1", color="#909090") +
      geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = "#F0F0F1", color="#909090")
    p <- p + geom_raster(data = psi, aes(x = lon, y = lat, fill = Prob)) + facet_wrap(~Year)
    p <- p + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA,  color="#909090") +
      geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = NA, color="#909090") +
      geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = NA, color="#909090")
    p <- p + geom_hline(data = limits, aes(yintercept = value, group = ind), linetype = "dashed", color = "grey70")
    p <- p + geom_hline(data = core, aes(yintercept = value, group = ind), linetype = "longdash", color = "grey70")
    p <- p + geom_hline(data = center, aes(yintercept = value), color = "grey70")
    if(show.points){
      p <- p + geom_point(data = routes, aes(x = Longitude, y = Latitude, size = speciestotal), color = "black", alpha = 0.5)
      p <- p + scale_size_continuous(guide = FALSE)
    }
    p <- p + scale_fill_continuous(low = "#F0F0F1", high = "#ef6e0b")
    p <- p + theme(panel.background = element_rect(fill = "#CDD2D4", color = NA), axis.title = element_blank())
    p <- p + scale_x_continuous("Longitude", breaks = seq(from = 10*(xmin%/%10 + as.logical(xmin%%10)),
                                                          to = 10*(xmax%/%10 + as.logical(xmax%%10)), by = 10),
                                limits = c(xmin, xmax))
    p <- p + scale_y_continuous("Latitude", breaks = seq(from = 5*(ymin%/%5 + as.logical(ymin%%5)),
                                                         to = 5*(ymax%/%5 + as.logical(ymax%%5)), by = 10),
                                limits = c(ymin, ymax))
    p
  }


}

#' md2html
#'
#' Function to covert .md file to html (ht https://github.com/richfitz/modeladequacy)
#' @export

md2html <- function(filename, dest = NULL) {
  if(is.null(dest)){ dest <- paste0(tools::file_path_sans_ext(filename), ".html")}
  opts <- setdiff(markdownHTMLOptions(TRUE), 'base64_images')
  markdownToHTML(filename, dest, options = opts)
}

