#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(mcmapper)
library(bayescpm)
library(knitr)


source("include.R")



global <- list()



# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("Bayesain power, assurance, and VoI calculator (CONFIDENTIAL - please do not share the link)"),
    tabsetPanel(id="input_tabset",
      tabPanel("Introduction",HTML("
               <H3><B>Welcome to the power, assurance, and VoI calculator for Bayesian and decision-theoretic approach for sample size considerations in risk prediction models</B></H3><BR/>
               <P>This Web app helps you to determine the power, sample size, or expected value of learning for the external validation of your risk prediction model based on uncertainty around its net benefit (NB).</P>
               <P><B>To use this app, you need the following categories of information: </B><BR/>
               
                  1. Outcome prevalence in the target population. <BR/>
                  
                  2. The performance of the model in terms of discrimination and calibration.<BR/>

                  3. Metrics of interest and precision / assurance targets.<BR/>

                  4. General setup for calculations (range of target sample size, number of Monte Carlo simulations).<BR/></P
                  
                  </P>

                  <HR/>
                  <DIV style='display: flex;	 align-items: flex-end;  align-text: center'>App version 2025.01.28. For questions and bug reports contact msafavi@mail.ubc.ca</DIV>
      ")),
      tabPanel("Evidence on model performance",
        HTML("Here we solicit your assessment of the model performance and your uncertainty around this asessment."),
        hr(),
        
        sidebarPanel("Evidence on outcome prevalence",
           selectInput("evidence_prev_dist_type", "Distribution type", choices=c("logitnorm", "beta", "probitnorm"), selected="beta"),
           selectInput("evidence_prev_input_type", "How do you want to parameterize", choices=unname(choices$evidence_prev_input_type), selected="Mean and SD"),
           
           fluidRow(
              column(4,
                numericInput("evidence_prev_input_parm1","Parm 1",value=0.427967,min=0,max=1,width="75px")),
              column(4,
                numericInput("evidence_prev_input_parm2","Parm 2",value=0.0295,min=0,max=1,width="75px")))
           ),
        
        sidebarPanel("Evidence on c-statistic",
          selectInput("evidence_cstat_dist_type", "Distribution type", choices=c("logitnorm", "beta", "probitnorm"), selected="beta"),
          selectInput("evidence_cstat_input_type", "How do you want to parameterize", choices=unname(choices$evidence_cstat_input_type), selected=unname(choices$evidence_cstat_input_type)[2] ),
          fluidRow(
            column(4,
              numericInput("evidence_cstat_input_parm1","Parm 1",value=0.7607,min=0,max=1,width="75px")),
            column(4,
              numericInput("evidence_cstat_input_parm2","Parm 2",value=0.0062,min=0,max=1,width="75px")))
        ),
        
        sidebarPanel("Evidence on calibration slope",
          selectInput("evidence_cal_slp_dist_type", "Distribution type", choices=c("norm")),
          selectInput("evidence_cal_slp_input_type", "How do you want to parameterize", choices=unname(choices$evidence_cal_slp_input_type), selected=unname(choices$evidence_cal_slp_input_type)[2] 
                         ),
          fluidRow(
            column(4,
                   numericInput("evidence_cal_slp_input_parm1","Parm 1",value=0.9950,min=0.5,max=2,width="75px")),
            column(4,
                   numericInput("evidence_cal_slp_input_parm2","Parm 2",value=0.0237,min=0,max=1,width="75px"))),
          
          selectInput("evidence_cal_other_type", "And one of:", choices=unname(choices$evidence_cal_other_type), selected="Mean calibration"),
          conditionalPanel("input.evidence_cal_other_type!=''",
            selectInput("evidence_cal_other_dist_type", "Distribution type", choices=c("norm","lognorm")),
            radioButtons("evidence_cal_other_input_type", "How do you want to parameterize", choices=unname(choices$evidence_cal_other_input_type), selected=unname(choices$evidence_cal_other_input_type)[2]),
            fluidRow(
              column(4,
                     numericInput("evidence_cal_other_input_parm1","Parm 1",value=-0.00934,min=0,max=1,width="75px")),
              column(4,
                     numericInput("evidence_cal_other_input_parm2","Parm 2",value=0.1245,min=0,max=1,width="75px")))
          ),
           
        ), sidebarPanel(numericInput("xxx", "xxx", value=1)),
        actionButton("btn_test_evidence", label="Test evidence"), uiOutput("test_evidence_results")
      ),
      tabPanel("Targets",
        HTML("Please choose outcome types, metrics, and target values."),
        hr(),
        fluidRow(
             column(12,
               sidebarPanel(
                 checkboxInput("b_stats","I want to specify precision targets for statistical metrics (discrimination and calibration)", value=TRUE),
                 hr(),
                 conditionalPanel(
                    HTML("Please specify the type of outcomes"),
                    condition = "input.b_stats==1",
                    checkboxInput("b_fciw","Conventional (frequentist) CI width", value=TRUE),
                    checkboxInput("b_eciw","Expected CI width", value=TRUE),
                    checkboxInput("b_qciw","Quantile of CI width"),
                    conditionalPanel(
                     condition = "input.b_qciw==1",
                     sliderInput("qciw", "Quantile (assurance) value", 0.0, 1, value=0.9),
                    ),hr(),
                     HTML("Please specify metrics of interest"),
                     conditionalPanel(condition="input.b_fciw+input.b_eciw+input.b_qciw==0", HTML("<font color='red'>At least one item needs to be selected.</font>")),
                     checkboxInput("b_target_cstat","c-statistc", value=TRUE),
                      conditionalPanel(
                        condition = "input.b_target_cstat == 1",
                        sliderInput("ciw_cstat", "95%CI width", 0.0, 0.3, value=0.1)
                      ),
                    checkboxInput("b_target_cal_slp","Calibration slope"), 
                    conditionalPanel(
                      condition = "input.b_target_cal_slp == 1",
                      sliderInput("ciw_cal_slp", "95%CI width", 0.0, 0.5, value=0.2)
                    ),
                    checkboxInput("b_target_cal_oe","Calibration O/E", value=TRUE), 
                    conditionalPanel(
                      condition = "input.b_target_cal_oe == 1",
                      sliderInput("ciw_cal_oe", "95%CI width", 0.0, 0.5, value=0.2)
                    ),
                    checkboxInput("b_target_cal_int","Calibration intercept"), 
                    conditionalPanel(
                      condition = "input.b_target_cal_int == 1",
                      sliderInput("ciw_cal_int", "95%CI width", 0.0, 0.5, value=0.2)
                    ),
                    checkboxInput("b_target_cal_mean","Calibration in the large (mean calibration)"), 
                    conditionalPanel(
                      condition = "input.b_target_cal_mean == 1",
                      sliderInput("ciw_cal_mean", "95%CI width", 0.0, 0.5, value=0.2)
                    )
                 )
             ),
             sidebarPanel(
               checkboxInput("b_nb","I want to specify decision-theoretic targets for net benefit (NB)", value=TRUE), 
               conditionalPanel(
                 condition = "input.b_nb == 1",
                 sliderInput("threshold", "Risk threshold", min=0, max=1, value=0.2),
                 checkboxInput("b_nb_voi","VoI (EVPI and EVSI)", value=TRUE),
                 checkboxInput("b_nb_assurance","Assurance", value = 1),
                 conditionalPanel(
                   condition = "input.b_nb_assurance == 1",
                   sliderInput("nb_assurance", "Value", 0.0, 1, value=0.9)
                 )
               )
             )
          ),
        )
      ),
      tabPanel("Setup & run",
        sidebarPanel(
          textAreaInput("N", "Please enter sample sizes of interest (comma separated)","250,500,1000,2000,4000,8000"),
          selectInput("dist_type", "Distribution of calibrated risks", choices=c("logitnorm","beta","probitnorm")),
          numericInput("n_sim", "Monte Carlo sample size", min=100, max=10^6, value=1000),
          numericInput("seed", "Random number seed (optional)", value = NULL, width="50px"),
          checkboxInput("b_impute_cor","Impute correlation", value=1),
          selectInput("method","Computation method", choices = unname(choices$method_type))
        ),
        actionButton("btn_run", "Run"), 
        actionButton("btn_clear_console", "Clear output"),
        actionButton("btn_show_args", "Show args"),
        pre(uiOutput("console")), #getwd(), paste(dir(getwd()),collapse=","),
        p("For any non-trivial computations, please use the R package or the local version of this app on your computer")
      ),
      tabPanel("Report",
               actionButton("btn_gen_report", "Generate report"),
               actionButton("btn_show_all", "Show results"),
               uiOutput("report1")
      )
    )
)





















# Define server logic required to draw a histogram
server <- function(input, output)
{
  require(dplyr)
  
  output$test_evidence_results<- reactive({
    e1 <- collect_evidence(input)
    e2 <- process_evidence(e1)
    paste(deparse(list(L1=e1,L2=e2)),collapse="")
  }) %>%  bindEvent(input$btn_test_evidence)
  
  observeEvent(input$btn_run,
   {
     args <- gen_args(input)
     
     global$args <<- args
     
     progress <- shiny::Progress$new()
     # Make sure it closes when we exit this reactive, even if there's an error
     on.exit(progress$close())
     progress$set(message = "Computing", value = 0)
     args$ex_args$f_progress <- function(txt){progress$inc(message=txt)}
     
     if(is.numeric(input$seed)) {set.seed(input$seed)}
     
     console <- capture.output({res <- do.call(BayesCPM, args=args)}, type="message")
     
     global$results <<- res
     
     output$console <- renderUI(HTML(paste("Run complete. Console text:", paste(console, collapse = "\n"))))
     
     .GlobalEnv$output <- global #For debugging
     
   })
  
  
#   output$console <- reactive({
#     
#     
#     
#   }) %>%  bindEvent(input$btn_run)
#   
#   
  
  observeEvent(input$btn_clear_console,
              {
                 output$console <- renderText("")
              })
  
  observeEvent(input$btn_show_args,
               {
                 output$console <- renderText(paste(deparse(isolate(gen_args(input))), collapse="\n"))
               })
  
  observeEvent(input$btn_gen_report, {
    #rmarkdown::render(input="Report.rmd", params=list(data=global$results), clean=T, runtime="shiny")
    #output$report1 <- renderUI(includeHTML("Report.html"))
    fl <- paste0(tempdir(),"/report.html")
    rmarkdown::render(input="Report.Rmd", output_file=fl, params=list(data=global), clean=T)
    output$report1 <- renderUI(tags$iframe(src=base64enc::dataURI(file=fl, mime="text/html"), style='width:80vw;height:80vh;'))
  })
  
  observeEvent(input$btn_show_all, {
    output$report1 <- renderUI(pre(paste(deparse(global),collapse = "\n")))
  })
  
}




shinyApp(ui = ui, server = server)
