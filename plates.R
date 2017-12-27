#' Well coordinates
#' 
#' Objects of this class represent a position in a microwell plate.  This class
#' exists so that the validity of the coordinates can be tested (for instance
#' \code{M16} would be rejected for a 96-well plate).
#' 
#' @slot well A position in the microwell plate.
#' @slot plateFormat A plate format.
#' 
#' @examples 
#' well <- Well(well = "A01", plateFormat = "384")
#' well
#' 
#' @export

Well <- setClass("Well",
  representation = representation(
    well        = "character",
    plateFormat = "character"),
  prototype = list(
    well        = "",
    plateFormat = "undefined"))

setMethod("as.character", "Well", function(x)
  paste0(x@well, " (", x@plateFormat, "-well format)"))

setMethod("show", "Well", function(object) cat(as.character(object)))

setGeneric("Row", function(object) standardGeneric("Row"))

setMethod("Row", "Well", function(object)
  gsub("[[:digit:]]+", "", object@well))

setGeneric("Column", function(object) standardGeneric("Column"))

setMethod("Column", "Well", function(object)
  gsub("[[:alpha:]]+", "", object@well))

validPlateFormats <- c("undefined", "96", "384", "1536")

setValidity("Well", function(object) {
  if (! object@plateFormat %in% validPlateFormats)
    return(paste0( "Supported plate formats: "
                 , paste(sapply(validPlateFormats, dQuote), collapse = ", ")
                 , "."))
  if (object@well == "")
    return("Missing well coordinates.")
})


#' Move to the next well
#' 
#' From left to right, and up to down, give the position of the next well.
#' 
#' @param well A \code{\link{Well}} object.
#' 
#' @examples 
#' 
#' nextWell(Well(well="A12", plateFormat="96"))
#' nextWell(Well(well="A12", plateFormat="384"))
#' 
#' @return A \code{\link{Well}} object.
#' 
#' @export

setGeneric("nextWell", function(well) standardGeneric("nextWell"))
setMethod ("nextWell", "Well", function(well) {
  if (well@plateFormat == "undefined") stop("Can not use undefined format.")
  wells <- num_to_well(1:well@plateFormat, plate = well@plateFormat)
  position <- which(wells == well@well)
  if (position >= well@plateFormat) stop("Next well out of range.")
  Well(well=wells[position + 1], plateFormat = well@plateFormat)
})


#' A multiwell plate
#' 
#' Objects of this class represent a multiwell plate.
#' 
#' @slot plate A table (\code{data.frame}, \code{DataFrame}, \code{data.table}
#'       or \code{tibble}).
#' @slot deadVolume The dead volume of wells in this plate.
#' @slot maxVolume The maxiumum colume of wells in this plate.
#'       
#' @examples 
#' 
#' plate <- Plate( plate      = tibble::tibble(well = num_to_well(1:384, plate = "384"))
#'               , deadVolume = 10
#'               , maxVolume  = 100)
#' plate
#' 
#' @importFrom platetools num_to_well
#' @export

Plate <- setClass( "Plate",
  slots     =    c( plate      = "ANY"
                  , deadVolume = "numeric"
                  , maxVolume  = "numeric"),
  prototype = list( plate = list()
                  , deadVolume = 10
                  , maxVolume  = 100)
)

setMethod("as.character", "Plate", function(x)
  paste0( "A Plate with data about ", nrow(x@plate), " wells "
        , "(dead volume: ", x@deadVolume
        , "; max volume: ", x@maxVolume, ")."))

setMethod("show", "Plate", function(object) cat(as.character(object)))

setMethod("colnames", "Plate", function(x) colnames(x@plate))


#' Get reagent name
#' 
#' In a source plate, get the name of the reagent contained in a given well.
#' 
#' @param plate A \code{\link{Plate}} object.
#' @param well A \code{\link{Well}} object.
#' 
#' @examples 
#' 
#' sourcePlate <- Plate(plate = tibble::tibble(well = platetools::num_to_well(1:384, plate = "384")))
#' sourcePlate %<>%
#'   setWell(Well(well = "A01", plateFormat = "384"), "dNTP", 100) %>%
#'   setWell(Well(well = "A02", plateFormat = "384"), "dNTP", 100) %>%
#'   setWell(Well(well = "A03", plateFormat = "384"), "buffer", 100)
#' 
#' sourcePlate %>% sourceReagent(Well(well = "A01"))
#' sourcePlate %>% sourceReagent(Well(well = "A03"))
#' sourcePlate %>% sourceReagent()
#' 
#' @family Plate method
#' 
#' @export

setGeneric("sourceReagent", function(plate, well) standardGeneric("sourceReagent"))

setMethod("sourceReagent", c("Plate", "Well"), function(plate, well) {
  wellName <- well@well
  plateTable <- plate@plate
  plateRow <- plateTable[plateTable$well == wellName, ]
  plateRow <- plateRow[, -1] # Removing the "well" column.
  index <- which(!is.na(plateRow))
  if (length(index) == 0)
    stop("No reagent in well ", wellName, ".")
  if (length(index) > 1)
    stop( "More than one reagent in well", wellName
        , ". This should not happen in source plates")
  colnames(plateRow)[index]
})

setMethod("sourceReagent", c("Plate", "missing"), function(plate, well) {
  colnames(plate@plate)[-1]
})

#' Get reagent volume
#' 
#' In a given plate, get the volume of the reagent contained in a given well.
#' 
#' @param plate A \code{\link{Plate}} object.
#' @param well A \code{\link{Well}} object.
#' 
#' @examples 
#' 
#' destPlate <- Plate(plate = tibble::tibble(well = platetools::num_to_well(1:384, plate = "384")))
#' destPlate %<>%
#'   setWell(Well(well = "A01", plateFormat = "384"), "dNTP", 50) %>%
#'   setWell(Well(well = "A02", plateFormat = "384"), "dNTP", 100) %>%
#'   setWell(Well(well = "A01", plateFormat = "384"), "buffer", 50)
#' 
#' destPlate %>% plateWellVolume(Well(well = "A01"))
#' destPlate %>% plateWellVolume(Well(well = "A01"), "buffer")
#' 
#' @export

setGeneric("plateWellVolume", function(plate, well, what) standardGeneric("plateWellVolume"))

setMethod("plateWellVolume", c("Plate", "Well", "missing"), function(plate, well, what) {
  if (ncol(plate@plate) == 1) return(0)
  wellName   <- well@well
  plateTable <- plate@plate
  plateRow   <- plateTable[plateTable$well == wellName, ]
  plateRow   <- plateRow[, -1] # Removing the "well" column.
  sum(plateRow, na.rm = TRUE)
})

setMethod("plateWellVolume", c("Plate", "Well", "character"), function(plate, well, what) {
  if (ncol(plate@plate) == 1) return(0)
  wellName   <- well@well
  plateTable <- plate@plate
  plateRow   <- plateTable[plateTable$well == wellName, ]
  plateRow   <- plateRow[, -1] # Removing the "well" column.
  plateRow[[what]]
})


#' Find a well that can provide enough reagent
#' 
#' Given a \code{\link{Plate}} object, check if the reagent is available, and
#' if yes, in which well.
#' 
#' @param plate A \code{\link{Plate}} object.
#' @param reagent A reagent name.
#' @param start A \code{\link{Well}} object (to avoid backtracking).
#' 
#' @return a \code{\link{Well}} object.
#' 
#' @examples 
#' sourcePlate <- Plate(plate = tibble::tibble(well = platetools::num_to_well(1:384, plate = "384")))
#' sourcePlate %<>%
#'   setWell(Well(well = "A01", plateFormat = "384"), "dNTP", 100) %>%
#'   setWell(Well(well = "A02", plateFormat = "384"), "dNTP", 100) %>%
#'   setWell(Well(well = "A03", plateFormat = "384"), "buffer", 100)
#' 
#' seekReagent(sourcePlate, "buffer")
#' 
#' @export

setGeneric("seekReagent", function(object, reagent, start) standardGeneric("seekReagent"))

setMethod ("seekReagent", c("Plate", "character", "Well"), function(object, reagent, start) {
  plateTable <- object@plate
  # Truncate the table to the start position
  plateTable <- plateTable[which(plateTable$well == start@well):nrow(plateTable),]
  # Remove empty wells
  plateTable <- plateTable[!as.vector(is.na(plateTable[,reagent])),]
  well <- plateTable$well[1]
  if (is.na(well)) stop("Reagent not found")
  Well(well=well)
})

setMethod ("seekReagent", c("Plate", "character", "missing"), function(object, reagent, start) {
  seekReagent(object, reagent, Well(well="A01", plateFormat = "undefined"))
})


#' Set the contents of a well
#' 
#' @param plate A \code{\link{Plate}} object.
#' @param well A \code{\link{Well}} object.
#' @param what A reagent name.
#' @param volume A volume (in nanoliters).
#' 
#' @return Returns a \code{\link{Plate}} object.
#' 
#' @examples 
#' 
#' sourcePlate <- Plate(plate = tibble::tibble(well = platetools::num_to_well(1:384, plate = "384")))
#' sourcePlate@plate
#' sourcePlate <- setWell(sourcePlate, Well(well = "A01", plateFormat = "384"), "dNTP", 100)
#' sourcePlate <- setWell(sourcePlate, Well(well = "A02", plateFormat = "384"), "dNTP", 100)
#' sourcePlate <- setWell(sourcePlate, Well(well = "A03", plateFormat = "384"), "buffer", 100)
#' sourcePlate@plate
#' 
#' @export

setGeneric("setWell", function(plate, well, what, volume) standardGeneric("setWell"))

setMethod( "setWell"
         , c("Plate", "Well", "character", "numeric")
         , function(plate, well, what, volume) {
  tbl <- plate@plate
  tbl[tbl$well == well@well, what] <- volume
  plate@plate <- tbl
  if (validObject(plate)) plate
})