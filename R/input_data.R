#' input_data
#'
#' Formats BBS data, climate data, and presence input files for all species in spp_list
#' @export

input_data <- function(start_yr = 1972, end_yr = 2014){
    opts <- read.csv("inst/model_opts.csv")
    spp_list <- read.csv('inst/spp_list.csv')

    year_range <- seq(from = start_yr, to = end_yr)

    bbs_data <- BBSclim::GetBBS()

    spp <- NULL
    ## Check if test_annual has already been run for species
    for(i in 1:length(spp_list$spp)){
      spp_test <- file.exists(paste0("inst/output/", spp_list$spp[i], "/pres/pres_in.pao"))
      if(!spp_test) spp <- c(spp, as.character(spp_list$spp[i]))
    }

    if(!is.null(spp)){
      if(length(spp) > 1){
        cores <- parallel::detectCores()
        if(!is.null(opts$limit.cores)){
          cores <- min(cores, opts$limit.cores, length(spp))
        }
        doParallel::registerDoParallel(cores = cores)

        spp_data <- foreach::foreach(i=1:length(spp), .packages = c("dplyr", "BBSclim")) %dopar%{

          ### If new species in list, create directory
          BBSclim::CreateSpp(spp[i])

          ### Get AOU numeric code and raw BBS count data
          spp_AOU <- rBBS::GetAOU(spp[i])
          suppressMessages(rBBS::GetSppCounts(bbs_raw = bbs_data, AOU = spp_AOU,
                                              Write = TRUE, years = year_range,
                                              path = paste0('inst/output/', spp[i])))

          ### Remove outliers
          suppressMessages(rBBS::RemoveOutliers(path = paste0('inst/output/', spp[i])))

          ### Add buffer around routes with counts > 0
          suppressMessages(BBSclim::buffer_BBS(bbs = bbs_data, alpha = spp[i]))

          ### Extract climate variables for route locations
          BBSclim::GetBioVars(alpha = spp[i])

          ### Write .pao input file for Presence
          BBSclim::write_pao(alpha = spp[i])

        }
        doParallel::stopImplicitCluster()
        check <- TRUE
      }else{
        ### If new species in list, create directory
        BBSclim::CreateSpp(spp)

        ### Get AOU numeric code and raw BBS count data
        spp_AOU <- rBBS::GetAOU(spp)
        suppressMessages(rBBS::GetSppCounts(bbs_raw = bbs_data, AOU = spp_AOU,
                                            Write = TRUE, years = year_range,
                                            path = paste0('inst/output/', spp)))

        ### Remove outliers
        suppressMessages(rBBS::RemoveOutliers(path = paste0('inst/output/', spp)))

        ### Add buffer around routes with counts > 0
        suppressMessages(BBSclim::buffer_BBS(bbs = bbs_data, alpha = spp))

        ### Extract climate variables for route locations
        BBSclim::GetBioVars(alpha = spp)

        ### Write .pao input file for Presence
        BBSclim::write_pao(alpha = spp)
      }
    }else{
      check <- TRUE
    }
     check
}
