#' simpleBrush
#' 
#' Interactive plot for brushing pairwise comparisons.
#' 
#' @param df The data frame containing the data to plot.  Its first two columns
#'        contain the data and the next columns, if any, may contain annotation.
#' @param plotFunction The function to plot the data.  Its two first arguments must
#'        be numerical vectors of horizontal and vertical coordinates.
#' @param ... Other parameters to be passed to plotFunction.
#' 
#' @examples 
#' \dontrun{
#' simpleBrush(cars)
#' }
#' 
#' @importFrom magrittr '%>%'
#' 
#' @export simpleBrush

simpleBrush <- function(df, plotFunction = plot, ...) {
  
  if( !require("shiny") )
    stop("Package 'shiny' not available, please install it.")
  
  if( !require("magrittr") )
    stop("Package 'magrittr' not available, please install it.")
  
  ui <-
    mainPanel( plotOutput( "plot", brush = "brush")
             , tableOutput("table")) %>%
    fluidPage %>%
    shinyUI

  server <- function(input, output) {
  
    output$plot <- renderPlot({
      plotFunction(df[,1:2], ...)
    })
  
    output$table <- renderTable({
      res <- brushedPoints( df
                          , input$brush
                          , colnames(df)[1]
                          , colnames(df)[2])
      if (nrow(res) == 0)
        return()
      res
    })
  }
  
  runApp( list( ui=ui, server=server)
        , host="0.0.0.0")
}