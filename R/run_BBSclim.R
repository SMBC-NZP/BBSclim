#' run_BBSclim
#'
#' Run full analysis
#' @param ... Arguments to be passed to model_opts()
#' @param annual TRUE/FALSE Should annual models be fit for all species in spp_list_annual?
#' @param runmods TRUE/FALSE Should psi & gamma/epsilon models be fit for all species in spp_list_runmods?
#' @export


run_BBSclim <- function(annual, runmods, install.pkgs = FALSE, ...){
  if(install.pkgs){
    httr::set_config(httr::config(ssl_verifypeer = 0L))

    ### Install packages from
    suppressMessages(devtools::install_github('crushing05/BBS.tenstop'))
    suppressMessages(devtools::install_github('crushing05/BBS.fiftystop'))
    suppressMessages(devtools::install_github('richfitz/datastorr'))
    suppressMessages(devtools::install_github('crushing05/rBBS'))
    suppressMessages(devtools::install_github('richfitz/remake'))
    suppressMessages(remake::install_missing_packages("remake.yml"))
    suppressMessages(devtools::install_github('krlmlr/here'))
  }



  BBSclim::model_opts(...)

  if(annual){
    spp <- read.csv("inst/spp_list_annual.csv")
    spp_names <- as.character(spp$spp)
    spp_names <- spp_names[!duplicated(spp_names)]
    vals <- list(spp_names = whisker::iteratelist(spp_names, value="spp_name"))

    str <- whisker::whisker.render(readLines("config/remake_annual.yml.whisker"), vals)
    writeLines(str, "remake_annual.yml")

    remake::make(remake_file = "remake_annual.yml")
  }


  if(runmods){
    spp <- read.csv("inst/spp_list_runmods.csv")
    spp_names <- as.character(spp$spp)
    spp_names <- spp_names[!duplicated(spp_names)]
    vals <- list(spp_names = whisker::iteratelist(spp_names, value="spp_name"))

    str <- whisker::whisker.render(readLines("config/remake_runmods.yml.whisker"), vals)
    writeLines(str, "remake_runmods.yml")

    remake::make(remake_file = "remake_runmods.yml")
  }
}
