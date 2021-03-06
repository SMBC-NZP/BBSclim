#' buffer_BBS
#'
#' Add all routes within rectangular buffer around the maximum and minimum observed latitude and longitude
#' @param spp_count Data frame containing the count data for the focal species obtained using the `GetSppCounts()` function
#' @param routes Data frame containing the route data obtained from the function `GetRoutes()`
#' @param buffer Numeric value indicating the number of degrees to buffer the count data (default = 2)
#' @param method Which buffer method to use? Options include: "rec", "bcr", "bcr + rec" (default = "rec")
#' @return A .csv file containing the following fields:
#' @return   routeID The unique 8 digit route ID for each route
#' @return   Year The year that the count was conducted
#' @return   aou The numeric code for the focal species
#' @return   countN Total individuals of the focal species recorded at stop N (50-stop data) or N-9:N (10-stop data)
#' @return   Latitude The latitude for the route
#' @return   Longitude The longitude for the route
#' @export


buffer_BBS	<- function(alpha, bbs, buffer = 2, method = "rec", raw = FALSE) {

  if(raw){
    spp_count <- read.csv(paste0("inst/output/", alpha, "/raw_counts.csv"))
  }else{
    spp_count <- read.csv(paste0("inst/output/", alpha, "/no_outlier_counts.csv"))
  }

  route_atrb <- bbs$routes

  if(method == "rec"){
    start_year <- min(spp_count$Year)
    end_year <- max(spp_count$Year)
    n_yr <- end_year - start_year + 1

    lat_occ	<- range(spp_count$Latitude)	# latitude of occupied routes
    long_occ	<- range(spp_count$Longitude)		# longitude of occupied routes
    long_occ[long_occ==-Inf] <- 0

    max_lat	<- max(lat_occ) + buffer
    min_lat	<- min(lat_occ) -  buffer
    min_long	<- min(long_occ) - buffer
    max_long	<- max(long_occ) + buffer

    ## Identify routes that buffer observed counts
    buff_routes	<- dplyr::filter(route_atrb, Latitude < max_lat & Latitude > min_lat & Longitude < max_long & Longitude > min_long)
    buff_routes <- dplyr::anti_join(buff_routes, spp_count, by = "routeID")
    buff_routes <- dplyr::select(buff_routes, routeID, Latitude, Longitude, Stratum, BCR)
    
    ## Create data frame containing years that buffered routes were run (0 count)
    buff_run <- dplyr::filter(bbs$weather,  routeID %in% buff_routes$routeID & Year %in% seq(from = start_year, to = end_year) & RunType == 1)
    buff_run <- dplyr::distinct(buff_run, routeID, Year, .keep_all = FALSE)
    
    ## Add lat, long, stratum, & BCR
    buff_run <- dplyr::left_join(buff_run, buff_routes)
    
    ### Create data frame with 0 counts for buffered routes
    col_counts <- grep("count|stop", names(spp_count), value = TRUE)
    count_buff <- dplyr::as_data_frame(matrix(0, nrow = nrow(buff_run),
                                              ncol = length(col_counts)))
    names(count_buff) <- col_counts
    count_buff <- dplyr::bind_cols(buff_run, count_buff)
    count_buff$aou <- unique(spp_count$aou)
    count_buff$speciestotal <- 0
    
    spp_count_buff <- dplyr::bind_rows(spp_count, count_buff)

    write.csv(spp_count_buff,
              paste0("inst/output/", alpha, "/count_buff.csv"),
              row.names = FALSE)
    }

  if(method == "bcr"){
    start_year <- min(spp_count$Year)
    end_year <- max(spp_count$Year)
    n_yr <- end_year - start_year + 1

    bcr <- unique(spp_count$BCR)
    bcr		<- bcr[bcr < 38]				# cuts out Mexico

    buff_routes	<- dplyr::filter(route_atrb, BCR %in% bcr)
    buff_routes <- dplyr::anti_join(buff_routes, spp_count, by = "routeID")
    buff_routes <- dplyr::select(buff_routes, routeID, Latitude, Longitude, Stratum, BCR)
    buff_routes <- buff_routes[rep(seq_len(nrow(buff_routes)), each = n_yr),]

    ### Create data frame with 0 counts for buffered routes
    col_counts <- grep("count", names(spp_count), value = TRUE)
    count_buff <- dplyr::as_data_frame(matrix(0, nrow = nrow(buff_routes), ncol = length(col_counts)))
    names(count_buff) <- col_counts
    count_buff$aou <- unique(spp_count$aou)
    count_buff$Year <- rep(seq(start_year, end_year), length(unique(buff_routes$routeID)))
    count_buff <- dplyr::bind_cols(buff_routes, count_buff)

    spp_count_buff <- dplyr::bind_rows(spp_count, count_buff)

    write.csv(spp_count_buff,
              paste0("inst/output/", alpha, "/count_buff.csv"),
              row.names = FALSE)
    }

  if(method == "bcr + rec"){
    start_year <- min(spp_count$Year)
    end_year <- max(spp_count$Year)
    n_yr <- end_year - start_year + 1

    bcr <- unique(spp_count$BCR)
    bcr		<- bcr[bcr < 38]				# cuts out Mexico

    lat_occ	<- range(spp_count$Latitude)	# latitude of occupied routes
    long_occ	<- range(spp_count$Longitude)		# longitude of occupied routes
    long_occ[long_occ==-Inf] <- 0

    max_lat	<- max(lat_occ) + buffer
    min_lat	<- min(lat_occ) -  buffer
    min_long	<- min(long_occ) - buffer
    max_long	<- max(long_occ) + buffer

    buff_routes	<- dplyr::filter(route_atrb, BCR %in% bcr & Latitude < max_lat & Latitude > min_lat & Longitude < max_long & Longitude > min_long)
    buff_routes <- dplyr::anti_join(buff_routes, spp_count, by = "routeID")
    buff_routes <- dplyr::select(buff_routes, routeID, Latitude, Longitude, Stratum, BCR)
    buff_routes <- buff_routes[rep(seq_len(nrow(buff_routes)), each = n_yr),]

    ### Create data frame with 0 counts for buffered routes
    col_counts <- grep("count", names(spp_count), value = TRUE)
    count_buff <- dplyr::as_data_frame(matrix(0, nrow = nrow(buff_routes), ncol = length(col_counts)))
    names(count_buff) <- col_counts
    count_buff$aou <- unique(spp_count$aou)
    count_buff$Year <- rep(seq(start_year, end_year), length(unique(buff_routes$routeID)))
    count_buff <- dplyr::bind_cols(buff_routes, count_buff)

    spp_count_buff <- dplyr::bind_rows(spp_count, count_buff)

    write.csv(spp_count_buff,
              paste0("inst/output/", alpha, "/count_buff.csv"),
              row.names = FALSE)
  }
}
