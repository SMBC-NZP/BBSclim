#' psi_2_gif
#'
#' Saves annual psi maps as .png files for conversion to gif
#' @param alpha 4-letter alpha code
#' @export

psi_2_gif <- function(alpha){

  dir.create(paste0('inst/output/', alpha, '/occ_gif'))
  psi <- read.csv(paste0('inst/output/', alpha, '/occ.csv'))

  years <- unique(psi$Year)

  usa <- ggplot2::map_data("state")
  canada <- ggplot2::map_data("worldHires", "Canada")
  mexico <- ggplot2::map_data("worldHires", "Mexico")

  xmin <- min(psi$lon)
  xmax <- max(psi$lon)

  ymin <- min(psi$lat)
  ymax <- max(psi$lat)

  for(i in 1:length(years)){
    psi2 <- dplyr::filter(psi, Year == years[i])

    p <- ggplot()
    p <- p + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "#F0F0F1",  color="#909090") +
      geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "#F0F0F1", color="#909090") +
      geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = "#F0F0F1", color="#909090")
    p <- p + geom_tile(data = psi2, aes(x = lon, y = lat, fill = Prob))
    p <- p + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA,  color="#909090") +
      geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = NA, color="#909090") +
      geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = NA, color="#909090")
    p <- p + coord_map(projection = "lambert", lat0 = 25, lat1 = 50, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
    p <- p + scale_fill_continuous(low = "#F0F0F1", high = "#ef6e0b", guide = FALSE)
    p <- p + scale_x_continuous("Longitude",breaks = seq(from = 10*(xmin%/%10 + as.logical(xmin%%10)),
                                                         to = 10*(xmax%/%10 + as.logical(xmax%%10)), by = 10))
    p <- p + scale_y_continuous("Latitude", breaks = seq(from = 5*(ymin%/%5 + as.logical(ymin%%5)),
                                                         to = 5*(ymax%/%5 + as.logical(ymax%%5)), by = 10))
    p <- p + theme(panel.background = element_rect(fill = "#CDD2D4", color = NA), axis.title = element_blank())
    p <- p + labs(title = years[i])
    png(file = paste0('inst/output/', alpha, '/occ_gif/occ_', years[i], '.png'))
    print(p)
    dev.off()
  }
}
