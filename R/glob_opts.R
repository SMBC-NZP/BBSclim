## Global options for models
tenstops <- TRUE
start_yr <- 1997
end_yr <- 2014
years <- seq(from = start_yr, to = end_yr)
het <- TRUE
annual <- TRUE
test <- TRUE
alpha <- "kewa"

save(tenstops, start_yr, end_yr, years, het, annual, test, alpha, file = "inst/glob_opts.RData")
