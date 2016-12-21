## Global options for data preparation

global_opts <- data.frame(tenstops = TRUE, # Use 10-stop or 50-stop BBS data?\
                          start_yr = 1997, # Start year for BBS data
                          end_yr = 2014)   # End year for BBS data

write.csv(global_opts, "inst/global_opts.csv", row.names = FALSE)


## Model options for running presence
model_opts <- data.frame(het = TRUE,       # Heterogeneous detection?
                         annual = TRUE,    # Annual variation in p?
                         psi.test = TRUE,  # Run subset of psi models?
                         gam.test = TRUE,  # Run subset of gam models?
                         Parallel = TRUE)      # Delete .out files once top model is selected?

write.csv(model_opts, "inst/model_opts.csv", row.names = FALSE)