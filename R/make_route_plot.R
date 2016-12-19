#' make_route_plot
#' 
#' plot all routes, including counts, buffer, and outliers
#' @export

make_route_plot <- function(buff_routes, used_routes, raw_routes, alpha){
  usa <- map_data("state")
  canada <- map_data("worldHires", "Canada")
  mexico <- map_data("worldHires", "Mexico")

  spp_buff2 <- dplyr::distinct(buff_routes, routeID, .keep_all = TRUE)
  spp_rts <- dplyr::distinct(used_routes, routeID, .keep_all = TRUE)
  spp_out <- dplyr::anti_join(raw_routes, used_routes)
  spp_out <- dplyr::distinct(spp_out, routeID, .keep_all = TRUE)
  
  xmin <- min(min(spp_buff2$Longitude), min(spp_out$Longitude)) - 4
  xmax <- max(max(spp_buff2$Longitude), max(spp_out$Longitude)) + 4
  
  ymin <- min(min(spp_buff2$Latitude), min(spp_out$Latitude)) - 4
  ymax <- max(max(spp_buff2$Latitude), max(spp_out$Latitude)) + 4
  
  if(nrow(spp_out) == 1){
    summ <- paste(nrow(spp_rts), " routes (", nrow(spp_out), " outlier)", sep = "")
  }else{
    summ <- paste(nrow(spp_rts), " routes (", nrow(spp_out), " outliers)", sep = "")
  }

  p <- ggplot() + theme_minimal() + theme(panel.background = element_rect(fill = "#CDD2D4", color = NA), axis.title = element_blank())
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
  # p <- p + annotate("text", x = -71, 
                    # y = 5*(ymin%/%5 + as.logical(ymin%%5)) + 5, label = summ, size = 6, color = "#494949")
                    # 
  to.pdf(print(p), filename = paste0("inst/output/", alpha, "/route_map.pdf"), width = 10, height = 6)
}

#' to.pdf
#' 
#' function to save figures as pdf objects 

to.pdf <- function(expr, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}