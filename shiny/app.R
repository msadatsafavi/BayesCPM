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
library(evsiexval)
library(mcmapper)
library(bayescpm)
library(knitr)


source("include.R")



global <- list()



# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("Bayesain power calculator (CONFIDENTIAL - please do not share the link)"),
    tabsetPanel(id="input_tabset",
      tabPanel("Introduction",HTML("
               <H3><B>Welcome to the Bayes CPM package for Bayesian and decision-theoretic approach for sample size considerations in risk prediction models</B></H3><BR/>
               <P>This Web app helps you determine the optimal sample size for the external validation of your risk prediction model based on uncertainty around its net benefit (NB).</P>
               <P><B>To use this app, you need the following categories of information: </B><BR/>
                  1. The performance of the model in terms of its sensitivity and specificity at the risk threshold of interest (alongside outcome prevalence).<BR/>

                  2. Precision targets

                  4. General setup for calculations (range of target sample size, number of Monte Carlo simulations).<BR/></P>

                  
                  </P>

                  <HR/>
                  <DIV style='display: flex;	 align-items: flex-end;  align-text: center'>App version 2023.12.20. For questions and bug reports contact msafavi@mail.ubc.ca</DIV>
      ")),
      tabPanel("Evidence on model performance",
        HTML("Here we solicit your assessment of the model performance and your uncertainty around this asessment. Your knowledge of model performance can be specified in different ways."),
        hr(),
        
        sidebarPanel("Prevalence",
           selectInput("evidence_prev_dist_type", "Distribution type", choices=c("logitnorm", "beta", "probitnorm")),
           radioButtons("evidence_prev_input_type", "Input parameters", choices=unname(choices$evidence_prev_input_type)),
           
           fluidRow(
              column(4,
                numericInput("evidence_prev_input_parm1","Parm 1",value=0.47,min=0,max=1,width="75px")),
              column(4,
                numericInput("evidence_prev_input_parm2","Parm 2",value=0.50,min=0,max=1,width="75px")))
           ),
        
        sidebarPanel("c-statistic",
          selectInput("evidence_cstat_dist_type", "Distribution type", choices=c("logitnorm", "beta", "probitnorm"), selected="beta"),
          radioButtons("evidence_cstat_input_type", "Input parameters", choices=unname(choices$evidence_cstat_input_type), selected=unname(choices$evidence_cstat_input_type)[2] ),
          fluidRow(
            column(4,
              numericInput("evidence_cstat_input_parm1","Parm 1",value=0.75,min=0,max=1,width="75px")),
            column(4,
              numericInput("evidence_cstat_input_parm2","Parm 2",value=0.05,min=0,max=1,width="75px")))
        ),
        
        sidebarPanel("Calibration",
          selectInput("evidence_cal_slp_dist_type", "Distribution type", choices=c("norm")),
          radioButtons("evidence_cal_slp_input_type", "Input parameters", choices=unname(choices$evidence_cal_slp_input_type), selected=unname(choices$evidence_cal_slp_input_type)[2] 
                         ),
          fluidRow(
            column(4,
                   numericInput("evidence_cal_slp_input_parm1","Parm 1",value=1,min=0,max=1,width="75px")),
            column(4,
                   numericInput("evidence_cal_slp_input_parm2","Parm 2",value=0.01,min=0,max=1,width="75px"))),
          
          selectInput("evidence_cal_other_type", "And one of:", choices=unname(choices$evidence_cal_other_type)),
          conditionalPanel("input.evidence_cal_other_type!=''",
            selectInput("evidence_cal_other_dist_type", "Distribution type", choices=c("norm","lognorm")),
            radioButtons("evidence_cal_other_input_type", "Input parameters", choices=unname(choices$evidence_cal_other_input_type), selected=unname(choices$evidence_cal_other_input_type)[2]),
            fluidRow(
              column(4,
                     numericInput("evidence_cal_other_input_parm1","Parm 1",value=0,min=0,max=1,width="75px")),
              column(4,
                     numericInput("evidence_cal_other_input_parm2","Parm 2",value=0.1,min=0,max=1,width="75px")))
          ),
           
        ), 
        actionButton("btn_test_evidence", label="Test evidence"), uiOutput("test_evidence_results")
      ),
      tabPanel("Targets",
        fluidRow(
             column(12,
               sidebarPanel(
                 checkboxInput("b_stats","Statistical metrics (discrimination and calibration)", value=TRUE),
                 conditionalPanel(
                   condition = "input.b_stats==1",
                   checkboxInput("b_fciw","Conventional (frequentist) CI width", value=TRUE),
                   checkboxInput("b_eciw","Expected CI width", value=TRUE),
                   checkboxInput("b_qciw","Quantile of CI width"),
                   conditionalPanel(
                     condition = "input.b_qciw==1",
                     sliderInput("qciw", "Assurance value", 0.0, 1, value=0.5)
                   ),
                   conditionalPanel(condition="input.b_fciw+input.b_eciw+input.b_qciw==0", HTML("<font color='red'>At least one item needs to be selected.</font>")),
                    br(),
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
               checkboxInput("b_nb","Net benefit", value=TRUE), 
               conditionalPanel(
                 condition = "input.b_nb == 1",
                 numericInput("threshold", "Risk threshold", min=0, max=1, value=0.5),
                 checkboxInput("b_nb_voi","VoI (EVPI and EVSI)", value=TRUE),
                 checkboxInput("b_nb_assurance","Assurance"),
                 conditionalPanel(
                   condition = "input.b_nb_assurance == 1",
                   sliderInput("nb_assurance", "Value", 0.0, 1, value=0.5)
                 )
               )
             )
          ),
        )
      ),
      tabPanel("Analysis setup",
        textAreaInput("N", "Please enter sample sizes of interest (comma separated)","250,500,1000,2000,4000,8000"),
        selectInput("dist_type", "Distribution of calibrated risks", choices=c("logitnorm","beta","probitnorm")),
        numericInput("n_sim", "Monte Carlo sample size", min=100, max=10^6, value=1000),
        numericInput("seed", "Random number seed (optional)", value = NULL, width="50px"),
        checkboxInput("b_impute_cor","Impute correlation"),
        selectInput("method","Computation method", choices = c("sample","2s")),
        p("For any non-trivial computations, please use the R package or the local version of this app on your computer")
      ),
      tabPanel("Run",
        actionButton("btn_gen_code", "Generate code"), actionButton("clear_results", "Clear results"),
        uiOutput("results1"),
        pre(id="console","")
      ),
      tabPanel("Report",
               actionButton("btn_gen_report", "Generate report"),
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
  
  output$results1 <- reactive({
    
    args <- gen_args(input)
    
    global$args <<- args
    
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Computing", value = 0)
    args$ex_args$f_progress <- function(txt){progress$inc(message=txt)}
    
    #paste(deparse(args),collapse = "")
    
    console <- capture.output({res <- do.call(BayesCPM, args=args)}, type="message")
    
    global$results <<- res
    
    paste("Run complete. Console text:", paste(console, collapse = ""))
    
  }) %>%  bindEvent(input$btn_gen_code)
  
  observeEvent(input$btn_gen_report, {
    #rmarkdown::render(input="Report.rmd", params=list(data=global$results), clean=T, runtime="shiny")
    #output$report1 <- renderUI(includeHTML("Report.html"))
    rmarkdown::render(input="Report.rmd", params=list(data=global$results), clean=T)
    output$report1 <- renderUI(tags$iframe(src=base64enc::dataURI(file="Report.html", mime="text/html"), style='width:80vw;height:80vh;'))
  })
}




shinyApp(ui = ui, server = server)
