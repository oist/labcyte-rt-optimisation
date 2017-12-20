#' makePlate384
#' 
#' Creates a data frame representing a 384-well plate
#' 
#' At the moment the row names are not homogeneous in length
#' ("A1" instead of "A01").  This may change later.
#' 
#' @note Class of columns may be changed later without notice !
#' 
#' @examples 
#' head(makePlate384)
#' summary(makePlate384())
#' 
#' @export makePlate384

makePlate384 <- function() {
  
  plate <- data.frame(
      Row = rep(LETTERS[1:16], 24)
    , Col = unlist(lapply(1:24, rep, 16))
  )
  
  rownames(plate) <- paste0( plate$Row
                           , plate$Col)
  plate
}