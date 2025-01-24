## Only run this example in interactive R sessions
  ui <- fluidPage(
    sidebarPanel(
      checkboxInput("XXX","XXX"),
      selectInput("plotType", "Plot Type",
                  c(Scatter = "scatter", Histogram = "hist")
      ),
      # Only show this panel if the plot type is a histogram
      conditionalPanel(
        condition = "input.XXX == 1",
        selectInput(
          "breaks", "Breaks",
          c("Sturges", "Scott", "Freedman-Diaconis", "[Custom]" = "custom")
        ),
        # Only show this panel if Custom is selected
        conditionalPanel(
          condition = "input.breaks == 'custom'",
          sliderInput("breakCount", "Break Count", min = 1, max = 50, value = 10)
        )
      )
    ),
    mainPanel(
      plotOutput("plot")
    )
  )
  
  server <- function(input, output) {
    x <- rnorm(100)
    y <- rnorm(100)
    
    output$plot <- renderPlot({
      if (input$plotType == "scatter") {
        plot(x, y)
      } else {
        breaks <- input$breaks
        if (breaks == "custom") {
          breaks <- input$breakCount
        }
        
        hist(x, breaks = breaks)
      }
    })
  }
  
  
  shinyApp(ui, server)
  