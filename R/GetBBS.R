#' GetBBS
#' 
#' Retrieve BBS data files
#' @export

GetBBS <- function(Tenstops){
  if(Tenstops){
    bbs <- BBS.tenstop::get_BBS10()
  }else{
    bbs <- BBS.fiftystop::get_BBS50()
  }
  bbs
}