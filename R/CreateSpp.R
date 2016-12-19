#' CreateSpp
#' 
#' Create folder to store results
#' @export

CreateSpp <- function(alpha){
  dir.create(paste0("inst/output/", alpha))
  
  
  file.copy(from = "inst/output/.gitignore", to = paste0("inst/output/", alpha, "/.gitignore"))
}

CreateSpp(alpha = "scta")
