#' model_opts
#'
#' Set model options for running presence
#' @export

model_opts <- function(tenstops = TRUE, psi.test = FALSE, gam.test = FALSE, het = TRUE,
                       Parallel = TRUE, limit.cores = 50){
    opts <- data.frame(tenstops = tenstops,
                       het = het,
                       psi.test = psi.test,  # Run subset of psi models?
                       gam.test =gam.test,  # Run subset of gam models?
                       Parallel = Parallel,
                       limit.cores = limit.cores)      # Delete .out files once top model is selected?

  write.csv(opts, "inst/model_opts.csv", row.names = FALSE)
}

