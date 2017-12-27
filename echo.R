#' @include plates.R

#' Transducer state
#' 
#' Objects of this class represent the state of the transduction system in the
#' Echo machine: the position of the tranducer in the coordinate system of the
#' source plate, and the position of the target plate on top of the transducer.
#' 
#' @slot source A \code{\link{Well}} in the source plate.
#' @slot destination A \code{\link{Well}} in the destination plate.
#' 
#' @examples 
#' 
#' Transducer(
#'   source      = Well(well = "A01", plateFormat = "384"),
#'   destination = Well(well = "A01", plateFormat = "96"))
#' 
#' @export

Transducer <- setClass("Transducer",
  slots = c(
    source      = "Well",
    destination = "Well"))

setMethod("as.character", "Transducer", function(x)
  paste0("Transducer at the following coordinates:\n",
      "  Source:      ", as.character(x@source), "\n",
      "  Destination: ", as.character(x@destination)))

setMethod("show", "Transducer", function(object) cat(as.character(object)))


#' Echo state
#' 
#' Objects of this class represent the state of an Echo machine loaded with
#' a source and a destination plate, and with its transducer at a given
#' position.
#' 
#' The model of the echo machine determines the resolution of the liquid
#' transfers: 2.5 nL increments for the \dQuote{555} and 25 nL increments for
#' the \dQuote{525}.
#' 
#' @slot source A source \code{\link\{Plate}}.
#' @slot destination A destination \code{\link\{Plate}}.
#' @slot transducer A \code{\link{Transducer}} state.
#' @slot model The model of the Echo machine.
#' @slot log A log of the transfers
#' 
#' @examples 
#' 
#' plate <- Plate(plate = tibble::tibble(well = platetools::num_to_well(1:384, plate = "384")))
#' transducer <- Transducer( source      = Well(well = "A01", plateFormat = "384")
#'                         , destination = Well(well = "A01", plateFormat = "96"))
#' echo <- Echo( source      = plate
#'             , destination = plate
#'             , transducer  = transducer
#'             , model       = "525")
#' echo
#' 
#' @export

Echo <- setClass("Echo",
  slots = c(
    source      = "Plate",
    destination = "Plate",
    transducer  = "Transducer",
    model       = "character",
    log         = "list"
  ))

setMethod("show", "Echo", function(object) {
  cat("A ", object@model, " Echo machine:\n",
      "  Source plate:      ", as.character(object@source), "\n",
      "  Destination plate: ", as.character(object@destination), "\n",
      as.character(object@transducer), "\n",
      length(object@log), " elements in the log.", sep = "")
})

#' Change a well position
#' 
#' Makes sure that the updated well position is in the same plate format as the
#' original one.
#' 
#' @param from A \code{\link{Well}} object.
#' @param to A \code{\link{Well}} object.
#' 
#' @return Returns a \code{\link{Well}} object.
#' 
#' @examples 
#' 
#' changeWell( Well(well = "A01", plateFormat = "384")
#'           , Well(well = "B03", plateFormat = "undefined"))
#' 
#' @export

setGeneric("changeWell", function(from, to) standardGeneric("changeWell"))

setMethod("changeWell", c("Well", "Well"), function (from, to) {
  plateFormat <- from@plateFormat
  if (    plateFormat != "undefined" &
       to@plateFormat != "undefined" &
       to@plateFormat != plateFormat)
    stop( "Plate formats not compatible (from: ", sQuote(plateFormat)
        , ", to: ", sQuote(to@plateFormat), ").")
  to@plateFormat <- plateFormat
  if (validObject(to)) to
})

#' Change the state of the transducer
#' 
#' Moves the tranducer to a new source position, or moves the destination plate,
#' or both
#' 
#' @param transducer A \code{\link{Transducer}} object.
#' @param source A \code{\link{Well}} object.
#' @param destination A \code{\link{Well}} object.
#' 
#' @return Returns a transducer object in the new state.
#' 
#' @examples 
#' 
#' transducer <- Transducer( source      = Well(well = "A01", plateFormat = "384")
#'                         , destination = Well(well = "A01", plateFormat = "96"))
#' transducer
#' 
#' moveTransducer(transducer, Well(well = "A02", plateFormat = "384"))
#' moveTransducer(transducer, destination = Well(well = "A02", plateFormat = "96"))
#' moveTransducer( transducer
#'               , source      = Well(well = "A02", plateFormat = "384")
#'               , destination = Well(well = "A02", plateFormat = "96"))
#' 
#' plate <- Plate(plate = tibble::tibble(well = platetools::num_to_well(1:384, plate = "384")))
#'               
#' echo <- Echo( source      = plate
#'             , destination = plate
#'             , transducer  = transducer
#'             , model       = "525")
#' 
#' echo
#' 
#' moveTransducer( echo, source      = Well(well = "A02", plateFormat = "384"))
#' 
#' moveTransducer( echo, destination = Well(well = "A02", plateFormat = "96"))
#'
#' moveTransducer( echo
#'               , source      = Well(well = "A02", plateFormat = "384")
#'               , destination = Well(well = "A02", plateFormat = "96"))

setGeneric( "moveTransducer"
          , function(object, source, destination) standardGeneric("moveTransducer"))

setMethod( "moveTransducer"
        , c("Transducer", "Well", "missing")
        , function(object, source, destination) {
  object@source <- changeWell(object@source, source)
  if (validObject(object)) object
})

setMethod( "moveTransducer"
        , c("Transducer", "missing", "Well")
        , function(object, source, destination) {
  object@destination <- changeWell(object@destination, destination)
  if (validObject(object)) object
})

setMethod( "moveTransducer"
        , c("Transducer", "Well", "Well")
        , function(object, source, destination) {
  object <- moveTransducer(object, source = source)
  object <- moveTransducer(object, destination = destination)
  if (validObject(object)) object
})

setMethod( "moveTransducer"
        , c("Echo", "Well", "missing")
        , function(object, source, destination) {
  object@transducer <- moveTransducer(object@transducer, source)
  if (validObject(object)) object
})

setMethod( "moveTransducer"
        , c("Echo", "missing", "Well")
        , function(object, source, destination) {
  object@transducer <- moveTransducer(object@transducer, destination = destination)
  if (validObject(object)) object
})

setMethod( "moveTransducer"
        , c("Echo", "Well", "Well")
        , function(object, source, destination) {
  object@transducer <- moveTransducer(object@transducer, source, destination)
  if (validObject(object)) object
})


#' Show transfer logs
#' 
#' Extract the logs from an \code{\link{Echo}} object and formats them as a
#' table.
#' 
#' @param echo \code{\link{Echo}}
#' 
#' @examples 
#' 
#' showLogs(echo)
#' 
#' @export

setGeneric("showLogs", function(echo) standardGeneric("showLogs"))

setMethod("showLogs", "Echo", function(echo)
  do.call(rbind, echo@log)
)


#' Transfer volume from a source to a destination plate
#' 
#' Reduce the volume of a reagent in the source plate at the coordinates of
#' the transducer, add that volume to the destination plate,
#' and register the transfer in the log.
#' 
#' @param echo An \code{\link{Echo}} object.
#' @param volume Volume to transfer.
#' 
#' The \code{what} parameter must correspond to a column name in the source
#' plate.
#' 
#' @return Returns An \code{\link{Echo}} object with the \code{source},
#' \code{destination} and \code{log} slots updated.
#' 
#' @examples 
#' 
#' sourcePlate <- Plate(plate = tibble::tibble(well = platetools::num_to_well(1:384, plate = "384")))
#' sourcePlate %<>%
#'   setWell(Well(well = "A01", plateFormat = "384"), "dNTP", 100) %>%
#'   setWell(Well(well = "A02", plateFormat = "384"), "dNTP", 100) %>%
#'   setWell(Well(well = "A03", plateFormat = "384"), "buffer", 100)
#' sourcePlate@plate
#' 
#' destinationPlate <- Plate(plate = tibble::tibble(well = platetools::num_to_well(1:384, plate = "384")))
#' 
#' transducer <- Transducer( source      = Well(well = "A01", plateFormat = "384")
#'                         , destination = Well(well = "A01", plateFormat = "96"))
#'               
#' echo <- Echo( source      = sourcePlate
#'             , destination = destinationPlate
#'             , transducer  = transducer
#'             , model       = "525")
#' 
#' echo %>% showLogs
#' 
#' echo %<>%
#'   transfer(10) %>%
#'   transfer(10) %>%
#'   moveTransducer(destination = Well(well = "A02", plateFormat = "96")) %>%
#'   transfer(5)
#'   
#' echo %>% showLogs
#' 
#' @importFrom magritter %<>%
#' @importFrom magrittr add
#' @export

setGeneric("transfer", function(echo, volume) setGeneric("transfer"))

setMethod( "transfer"
         , c("Echo", "numeric")
         , function (echo, volume) {
  what            <-   sourceReagent(echo@source, echo@transducer@source)
  availableVolume <- plateWellVolume(echo@source, echo@transducer@source)
  if (!enoughVolume(echo, volume))
    stop( "Not enough ", what, " remaining (want: ", volume
        , ", has: ", availableVolume
        , ", min. well volume: ", echo@source@deadVolume, ".")
  destinationVolume <- plateWellVolume(echo@destination, echo@transducer@destination)
  destinationReagentVolume <- plateWellVolume(echo@destination, echo@transducer@destination, what)
  if (destinationVolume + volume > echo@destination@maxVolume)
    stop( "Overflow in destination well", echo@transducer@destination@well
        , " (transfer: ", volume
        , ", current total volume: ", destinationVolume
        , ", max. well volume: ", echo@destination@maxVolume, ".")
  # Take the volume from the source
  echo@source %<>% setWell(echo@transducer@source, what, availableVolume - volume)
  # Shoot it to the destination
  echo@destination %<>% setWell(echo@transducer@destination, what, destinationReagentVolume + volume)
  # Log the transfer
  echo@log[[length(echo@log) + 1]] <-
    list( from = echo@transducer@source@well
        , to   = echo@transducer@destination@well
        , vol  = volume
        , what = what)
  if (validObject(echo)) echo
})

setGeneric("enoughVolume", function(echo, volume) standardGeneric("enoughVolume"))
setMethod ("enoughVolume", c("Echo", "numeric"), function(echo, volume) {
    availableVolume <- plateWellVolume(echo@source, echo@transducer@source)
    availableVolume - volume > echo@source@deadVolume
})

#' Plan Echo transfers
#' 
#' Given a source plate, a target plate and the representation of a filled
#' target plate, calculate the liquid tranfers necessary perform the filling.
#' 
#' Each reagent of the plan plate is transferred one after the other.  Once
#' a reagent is selected, a source well is sought for the transfer and the
#' transducer is moved to it.
#' 
#' @param source A source \code{\link\{Plate}}.
#' @param destination A destination \code{\link\{Plate}}.
#' @param plan A \code{\link\{Plate}} representing the additions to the
#'       destination plate.
#' @param model Which Echo model to use.
#' 
#' @return A transfer log
#' 
#' @seealso \code{\link{showLogs}}
#' 
#' @examples 
#' sourcePlate <- Plate(plate = tibble::tibble(well = platetools::num_to_well(1:384, plate = "384")))
#' sourcePlate %<>%
#'   setWell(Well(well = "A01", plateFormat = "384"), "dNTP", 100) %>%
#'   setWell(Well(well = "A02", plateFormat = "384"), "dNTP", 100) %>%
#'   setWell(Well(well = "A03", plateFormat = "384"), "buffer", 100)
#' sourcePlate@plate
#' 
#' planPlate <- Plate(plate = tibble::tibble(well = platetools::num_to_well(1:384, plate = "384")))
#' planPlate %<>%
#'   setWell(Well(well = "A01", plateFormat = "384"), "dNTP", 50) %>%
#'   setWell(Well(well = "A02", plateFormat = "384"), "dNTP", 100) %>%
#'   setWell(Well(well = "A01", plateFormat = "384"), "buffer", 50)
#' planPlate@plate
#' 
#' destinationPlate <- Plate(plate = tibble::tibble(well = platetools::num_to_well(1:384, plate = "384")))
#' 
#' 
#' @export

planTransfers <- function(source, destination, plan, model = "525") {
  echo <- Echo( source      = sourcePlate
              , destination = destinationPlate
              , transducer  = Transducer( source = Well(well = "A01",      plateFormat = "384")
                                        , destination = Well(well = "A01", plateFormat = "384"))
              , model       = model)
  
# To do: better support plate formats.  Plate objects should have a plateFormat
  # slot; in most cases it should not be retreived from the transducer.
# To do: reform the sourceReagent function (the name is misleading, since now
# it can also be run on destination plates)
  
  for (reagent in sourceReagent(plan)) {
    # Move transducer in upper-left corner
    echo %<>% moveTransducer(source = Well(well = "A01", plateFormat = echo@transducer@source@plateFormat))
    # Find reagent in source plate
    echo %<>% seekReagent(reagent)
  }
}

setMethod("seekReagent", c("Echo", "character", "missing"), function(object, reagent, start) {
  echo %>% moveTransducer(source = seekReagent(echo@source, reagent, echo@transducer@source))
})
