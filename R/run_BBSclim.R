#' run_BBSclim
#'
#' Run full analysis
#' @param ... Arguments to be passed to model_opts()
#' @return If you're lucky, .html files and model results for each species in spp_list.csv
#' @export


run_BBSclim <- function(...){
  httr::set_config(httr::config(ssl_verifypeer = 0L))

  ### Install packages from
  devtools::install_github('crushing05/BBS.tenstop')
  devtools::install_github('crushing05/BBS.fiftystop')
  devtools::install_github('richfitz/datastorr')
  devtools::install_github('crushing05/rBBS')
  devtools::install_github('richfitz/remake')
  remake::install_missing_packages("remake.yml")
  devtools::install_github('krlmlr/here')


  BBSclim::model_opts(...)

  spp <- read.csv("inst/spp_list.csv")
  spp_names <- as.character(spp$spp)
  vals <- list(spp_names = whisker::iteratelist(spp_names, value="spp_name"))

  str <- whisker::whisker.render(readLines("config/remake.yml.whisker"), vals)
  writeLines(str, "remake.yml")

  remake::make()
}
