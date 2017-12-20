#' readResults
#' 
#' Read a qPCR result file in SDS 2.4 AQ Results 1.0 format
#' 
#' @return Returns a data frame with one well per row, sorted in
#' the same order as the plate output by the function \code{makePlate384}.
#' 
#' @param file Path to the file
#' 
#' @examples 
#' \dontrun{
#' readResults("H3903W3V.txt")
#' }
#' 
#' @export readResults

readResults <- function(file) {
  
  df <- read.table( file
                  , skip = 10
                  , head = T
                  , nrows = 384
                  , sep="\t"
                  , row.names = "Sample.Name"
                  , stringsAsFactors = FALSE)

  df$Ct[df$Ct == "Undetermined"] <- NA
 
  df$Ct <- as.numeric(df$Ct)
  
  # Sort the results
  df[rownames(makePlate384()),]
}