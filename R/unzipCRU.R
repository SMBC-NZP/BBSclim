#' Unzip CRU climate data.
#'
#' \code{unzipCRU} unzips the compressed CRU climate data and
#' adds the decompressed files to the "data" folder for each variable.
#'
#' @param varb The climate variable of interest
#' @param in_dir The directory containing the compressed files
#' @param out_dir The director where the uncompressed files will be stored
#' @export


unzipCRU <- function(varb, in_dir = "data-raw/", out_dir = "data-raw/"){
  var_in_dir <- paste(in_dir, varb, "/", sep = "")
  var_files <- list.files(var_in_dir)

  var_out_dir <- paste(out_dir, varb, "/", sep = "")

  for(i in 1:length(var_files)){
    R.utils::gunzip(filename = paste(var_in_dir, var_files[i], sep = ""),
           destname = paste(var_out_dir, gsub(".gz", "", var_files[i]), sep = ""),
           remove = FALSE, overwrite = TRUE)
  }
}
