#' global_opts
#'
#' Set global options for data preparation
#' @export

global_opts <- function(tenstops = TRUE, start_yr = 1997, end_yr = 2014){
  opts <- data.frame(tenstops = tenstops, # Use 10-stop or 50-stop BBS data?\
                     start_yr = start_yr, # Start year for BBS data
                     end_yr = end_yr)   # End year for BBS data

  write.csv(opts, "inst/global_opts.csv", row.names = FALSE)
}


#' model_opts
#'
#' Set model options for running presence
#' @export

model_opts <- function(psi.test = FALSE, gam.test = FALSE, het = TRUE,
                       Parallel = TRUE, limit.cores = 50){
    opts <- data.frame(het = het,       # Heterogeneous detection?
                       psi.test = psi.test,  # Run subset of psi models?
                       gam.test =gam.test,  # Run subset of gam models?
                       Parallel = Parallel,
                       limit.cores = limit.cores)      # Delete .out files once top model is selected?

  write.csv(opts, "inst/model_opts.csv", row.names = FALSE)
}

