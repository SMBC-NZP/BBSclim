#' GetAOU
#'
#' Get numeric AOU species code using common name or four letter alpha code
#' @param alpha Data frame containing the (buffered) count data for the focal species
#' @param common Integer vector containing the bioclim variables of interest
#' @return The four digit numeric code for the species of interest
#' @export

GetAOU <- function(alpha = NULL, common = NULL){
  if(!is.null(alpha)){
    alpha = toupper(alpha)
    AOU = code_lookup$AOU[code_lookup$alpha == alpha]}
  if(!is.null(common)){
    common = toupper(common)
    AOU = code_lookup$AOU[code_lookup$alpha == common]
  }
  AOU
}

