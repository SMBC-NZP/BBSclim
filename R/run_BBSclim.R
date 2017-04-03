#' run_BBSclim
#'
#' Run full analysis
#' @param ... Arguments to be passed to model_opts()
#' @return If you're lucky, .html files and model results for each species in spp_list.csv
#' @export


run_BBSclim <- function(...){
  httr::set_config(httr::config(ssl_verifypeer = 0L))


  ### Install packages from
  suppressMessages(devtools::install_github('crushing05/BBS.tenstop'))
  suppressMessages(devtools::install_github('crushing05/BBS.fiftystop'))
  suppressMessages(devtools::install_github('richfitz/datastorr'))
  suppressMessages(devtools::install_github('crushing05/rBBS'))
  suppressMessages(devtools::install_github('richfitz/remake'))
  suppressMessages(remake::install_missing_packages("remake.yml"))
  suppressMessages(devtools::install_github('krlmlr/here'))


  BBSclim::model_opts(...)

  spp <- read.csv("inst/spp_list.csv")
  spp_names <- as.character(spp$spp)
  vals <- list(spp_names = whisker::iteratelist(spp_names, value="spp_name"))

  str <- whisker::whisker.render(readLines("config/remake.yml.whisker"), vals)
  writeLines(str, "remake.yml")

  remake::make()
}
