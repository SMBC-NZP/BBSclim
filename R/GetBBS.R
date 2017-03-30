#' GetBBS
#'
#' Retrieve BBS data files
#' @export

GetBBS <- function(is.tenstops = NULL){
  if(is.null(is.tenstops)){
    opts <- read.csv("inst/model_opts.csv")
    if(opts$tenstops){
      bbs <- BBS.tenstop::get_BBS10()
    }else{
      bbs <- BBS.fiftystop::get_BBS50()
    }
  }else{
    if(is.tenstops){
      bbs <- BBS.tenstop::get_BBS10()
    }else{
      bbs <- BBS.fiftystop::get_BBS50()
    }
  }
  bbs
}
