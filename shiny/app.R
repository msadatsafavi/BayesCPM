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

source("include.R")


make95CrIFromBeta <- function(a,b)
{
  paste0("95%CI ",
         round(qbeta(0.025,a,b)*100,1),
         "%-",
         round(qbeta(0.975,a,b)*100,1),"%"
  )
}


withConsoleRedirect <- function(containerId, expr) {
  # Change type="output" to type="message" to catch stderr
  # (messages, warnings, and errors) instead of stdout.
  txt <- capture.output(results <- expr, type = "output")
  if (length(txt) > 0) {
    insertUI(paste0("#", containerId), where = "beforeEnd",
             ui = paste0(txt, "\n", collapse = "")
    )
  }
  results
}



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
           selectInput("evidence_prev_dist_type", "Distribution type", choices=c("logitnormal", "beta", "probitnormal")),
           radioButtons("evidence_prev_input_type", "Input parameters", choices=unname(choices$evidence_prev_input_type)),
           
           fluidRow(
              column(4,
                numericInput("evidence_prev_input_parm1","Parm 1",value=NA,min=0,max=1,width="75px")),
              column(4,
                numericInput("evidence_prev_input_parm2","Parm 2",value=NA,min=0,max=1,width="75px")))
           ),
        
        sidebarPanel("c-statistic",
          selectInput("evidence_cstat_dist_type", "Distribution type", choices=c("logit-normal", "beta", "probit-normal")),
          radioButtons("evidence_cstat_input_type", "Input parameters", choices=unname(choices$evidence_cstat_input_type)),
          fluidRow(
            column(4,
                   numericInput("evidence_cstat_input_parm1","Parm 1",value=NA,min=0,max=1,width="75px")),
            column(4,
                   numericInput("evidence_cstat_input_parm2","Parm 2",value=NA,min=0,max=1,width="75px")))
        ),
        
        sidebarPanel("Calibration",
          selectInput("evidence_cal_slp_dist_type", "Distribution type", choices=c("normal")),
          radioButtons("evidence_cal_slp_input_type", "Input parameters", choices=unname(choices$evidence_cal_slp_input_type)),
          fluidRow(
            column(4,
                   numericInput("evidence_cal_slp_input_parm1","Parm 1",value=NA,min=0,max=1,width="75px")),
            column(4,
                   numericInput("evidence_cal_slp_input_parm2","Parm 2",value=NA,min=0,max=1,width="75px"))),
          
          selectInput("evidence_cal_other_type", "And one of:", choices=unname(choices$evidence_cal_other_type)),
          conditionalPanel("input.evidence_cal_other_type!=''",
            radioButtons("evidence_cal_other_input_type", "Input parameters", choices=c("Point estimate and upper 95%CI bound", 
                                                                               "Mean and SD")),
            fluidRow(
              column(4,
                     numericInput("evidence_cal_other_input_parm1","Parm 1",value=NA,min=0,max=1,width="75px")),
              column(4,
                     numericInput("evidence_cal_other_input_parm2","Parm 2",value=NA,min=0,max=1,width="75px")))
          )
        )
      ),
      tabPanel("Criteria",
        fluidRow(
          
             column(12,
                    
             checkboxInput("b_stats","Statistical metrics (discrimination and calibration)"),
             conditionalPanel(
               condition = "input.b_stats == 1",
               
               checkboxInput("e_ciw","Expected CI width", value=TRUE),
               checkboxInput("q_ciw","Quantile of CI width"),
               conditionalPanel(
                 condition = "input.q_ciw == 1",
                 sliderInput("ass_stat", "Assurance value", 0.0, 1, value=0.5)
               ),
               
               conditionalPanel(condition="input.e_ciw+input.q_ciw==0", HTML("<font color='red'>At least one item needs to be selected.</font>"))
             ),
             
             hr(),
             
             checkboxInput("b_nb","Net benefit"), 
             conditionalPanel(
               condition = "input.b_nb == 1",
               checkboxInput("b_voi_nb","VoI (EVPI and EVSI"),
               checkboxInput("b_a_nb","Assurance"),
               conditionalPanel(
                 condition = "input.b_a_nb == 1",
                 sliderInput("ass_nb", "Value", 0.0, 1, value=0.5)
               )
             )
          ),
        )
      ),
      
      tabPanel("Target precision",
        
         conditionalPanel(condition="input.b_stats",
                          
           checkboxInput("b_cstat","c-statistc"),
           conditionalPanel(
             condition = "input.b_cstat == 1",
             sliderInput("ciw_cstat", "95%CI width", 0.0, 0.3, value=0.1)
           ),
           
           checkboxInput("b_cal_slp","Calibration slope"), 
           conditionalPanel(
             condition = "input.b_cal_slp == 1",
             sliderInput("ciw_cal_slp", "95%CI width", 0.0, 0.5, value=0.2)
           ),
           
           checkboxInput("b_cal_oe","Calibration O/E"), 
           conditionalPanel(
             condition = "input.b_cal_oe == 1",
             sliderInput("ciw_cal_oe", "95%CI width", 0.0, 0.5, value=0.2)
           ),
           
           checkboxInput("b_cal_int","Calibration intercept"), 
           conditionalPanel(
             condition = "input.b_cal_int == 1",
             sliderInput("ciw_cal_int", "95%CI width", 0.0, 0.5, value=0.2)
           ),
           
           checkboxInput("b_cal_mean","Calibration in the large (mean calibration)"), 
           conditionalPanel(
             condition = "input.b_cal_mean == 1",
             sliderInput("ciw_cal_mean", "95%CI width", 0.0, 0.5, value=0.2)
           )
         ),
         
         conditionalPanel(condition="input.b_nb==1",
            numericInput("NB threshold", "Threshold", min=0, max=1, value=0.5),
         )
      ),               
      tabPanel("Analysis setup",
        textAreaInput("N", "Please enter sample sizes of interest (comma separated)","0,500,1000,2000,4000,8000,16000"),
        fluidRow(
          column(4,
            numericInput("n_sim", "Monte Carlo sample size", min=100, max=10^6, value=10^5)),
          br(),
          column(4,
            numericInput("seed", "Random number seed (optional)", value = NULL))
        ),
        br(),
        p("For any non-trivial computations, please use the R package or the local version of this app on your computer")
      ),
      tabPanel("Results",
        actionButton("btn_gen_code", "Generate the code"), actionButton("clear_results", "Clear results"),
        uiOutput("results1"),
        pre(id="console","")
      )
    ),
)





















# Define server logic required to draw a histogram
server <- function(input, output)
{
  
  observeEvent(input$btn_gen_code, {generate_code(input)})
#   #browser()
#   global_vars <<- list(
#     params=list(
#       evidence_type="", #ind_beta, prev_cs_A_B, sample
#       evidence=list()
#     ),
#     default_values=list(
#       exchange_rates=list(z=0.02, lambda=100),
#       evidence=list(
#         ind_beta=list(prev = c(43L, 457L), se = c(41L, 2L), sp = c(147, 310)),
#         prev_cs_A_B=list(prev=c(point=0.086000, upper_ci=0.110659), cs=c(point=0.8474887, upper_ci=0.9147262), A=c(point=0.5379061, sd=0.3824605), B=c(point=1.1793535, sd=0.1692672))
#       )
#     ),
#     result_level=0,
#     sample=data.frame(),
#     str_val="", #Will contain text generated by process_input()
#     str_desc ="" #General textual description of the parameters
#   )
# 
#   observeEvent(input$z, {
#     z <- input$z/100
#     output$z_desc <- renderText(paste0("The risk threshold is the threshold on predicted risk at which point the care provider / patient is indifferent between treatment or no treatment.\n
#           A risk threshold of ",input$z, "% indicates that the benefit of a true positive case is equal (in the opposite direction) to the harm of ",(1-z)/z," false positive cases."))
#   })
#   
#   observeEvent(input$b_cstat, {
#     
#   })
# 
#   observeEvent(input$lambda, {
#     output$lambda_desc <- renderText(paste("The Number Needed to Study (NNS) represents the trade-off between sampling efforts and clinical utility. For example, a NNS of 100 means the investigator believes the efforts required to procure 100 more samples is equal to the benefit of one true positive diagnosis.\n
#           To choose this value, think of the context: How important the clinical event is? For example, for a catastrophic event like a stroke, one might justify procuring a sample in hundreds, or thousands, to prevent one more event\n
#           This is also related to the ease at which samples can be obtained. Will this external validation study be based on primary or secondary data collection? The former might require significantly more efforts. \n
#           In thinking about sampling efforts, do not consider one-time costs and efforts required for setting up a validation study. Those one-time costs are not affected by sample size. Instead, think of  'incremental' effort of procuring samples.
#           "))
#   })
# 
#   observeEvent(input$infer_lambda, {
#     z <- input$z/100
#     updateSliderInput(inputId="lambda", value=max(1, (1-z)/(z)))
#     output$infer_lambda_desc <- renderText(paste("The inferred minimum NNS is based on the notion that in most contexts, it is justifiable to obtain one more observation insofar as that observation is expected to prevent one incorrect medical decision.\n
#                                                  Incorrect decisions can be in terms of not detecting an individual who will experience the event (missing a true positive), or unnecessarily treating an individuals who will not experience the event (causing a false positive).\n"
#                                                  ,ifelse(z<0.5,
#                                                     paste0("Because at risk threshold of ", input$z, "% each true positive detection is equal to", (1-z)/z, " false positives, and because the EVSI is in true positive units, the desired NNS should be at least", (1-z)/z,"."),
#                                                     ""
#                                                     ),
#                                                   "In general, this reasoning results in NNS of max(1,(1-z)/z) where z is the clinical threshold."))
#   })
# 
#   observeEvent(input$clear_results, {
#     global_vars$result_level<<-0
#     output$results1=renderUI(HTML(""))
#   })
# 
#   observeEvent(input$evidence_type, {
#     choice <- substring(input$evidence_type,1,1)
#     output$n_sim_outer_note <- renderText("")
#     output$evidence_inputs <- renderUI("")
# 
#     if(choice=="0")
#     {
#       global_vars$params$evidence_type <<- ""
#     }
#     if(choice=="1")
#     {
#       global_vars$params$evidence_type <<- "ind_beta"
#       output$evidence_inputs <- renderUI(list(
#         renderText("This method is based on trivariate specification (outcome prevalence, sensitivity, and specificity), with Beta distribution for modeling uncertainty for each component."),
#         hr(),
#         fluidRow(
#           column(4,sliderInput("prev",label="Expected outcome risk (outcome prevalence) (%)", min=0, max=100, step=0.1, value=global_vars$default_values$evidence$ind_beta$prev[1]/sum(global_vars$default_values$evidence$ind_beta$prev)*100, width="100%")),
#           column(4,numericInput("prev_n",label="Your estimate of the average risk of the clinical outcome is based on how many observations?", value=sum(global_vars$default_values$evidence$ind_beta$prev))),
#           column(4, textOutput("prev_dist"))
#         ),hr(),
#         fluidRow(
#           column(4,sliderInput("se",label="Sensitivity of the model at the chosen risk threshold (%)", min=0, max=100, step=0.1, value=global_vars$default_values$evidence$ind_beta$se[1]/sum(global_vars$default_values$evidence$ind_beta$se)*100, width="100%")),
#           column(4,numericInput("se_n",label="Your estimate of sensitivity is based on how many observations?", value=sum(global_vars$default_values$evidence$ind_beta$se))),
#           column(4, textOutput("se_dist"))
#         ),hr(),
#         fluidRow(
#           column(4,sliderInput("sp",label="Specificity the model at the chosen risk threshold (%)", min=0, max=100, step=0.1, value=global_vars$default_values$evidence$ind_beta$sp[1]/sum(global_vars$default_values$evidence$ind_beta$sp)*100, width="100%")),
#           column(4,numericInput("sp_n",label="Your estimate of sensitivity is based on how many observations?", value=sum(global_vars$default_values$evidence$ind_beta$sp))),
#           column(4, textOutput("sp_dist"))
#         )
#       ))
#       observeEvent(input$prev, {
#         a <- input$prev/100*input$prev_n
#         b <- input$prev_n-input$prev/100*input$prev_n
#         output$prev_dist <- renderText(paste0("prev~Beta(", a,",", b ,") | ", make95CrIFromBeta(a,b)))
#       })
#       observeEvent(input$se, {
#         a <- input$se/100*input$se_n
#         b <- input$se_n-input$se/100*input$se_n
#         output$se_dist <- renderText(paste0("se~Beta(", a,",", b ,") | ", make95CrIFromBeta(a,b)))
#       })
#       observeEvent(input$sp, {
#         a <- input$sp/100*input$sp_n
#         b <- input$sp_n-input$sp/100*input$prev_n
#         output$sp_dist <- renderText(paste0("sp~Beta(", a,",", b ,") | ", make95CrIFromBeta(a,b)))
#       })
#       shinyjs::disable("n_sim_outer")
#       output$n_sim_outer_note <- renderText("Outer simulation is not applicable to the selected type of evidence specification on model performance because of conjugate probability distirbutions.")
#     }
# 
#     if(choice=="2")
#     {
#       global_vars$params$evidence_type <<- "prev_cs_A_B"
#       output$evidence_inputs <- renderUI(list(
#         renderText("This method is based on the specification of your informaiton on typical performance metrics of a model (c-statistic, and calibration intercept and slope)."),
#         hr(),
#         fluidRow(
#           column(4, sliderInput("prev_m_uci",label="Expected outcome risk (outcome prevalence) (%)", min=0, max=100, step=0.1, value=global_vars$default_values$evidence$prev_cs_A_B$prev*100, width="100%")),
#           column(4, textOutput("prev_dist"))
#         ),hr(),
#         fluidRow(
#           column(4, sliderInput("cstat_m_uci",label="Expected value and upper 95%CI bound of the c-statistic", min=0.51, max=0.99, step=0.001, value=global_vars$default_values$evidence$prev_cs_A_B$cs, width="100%")),
#           column(4, textOutput("cstat_dist"))
#         ),hr(),
#         fluidRow(
#           column(4, sliderInput("cal_intercept",label="Expected calibration intercept", min=-1, max=1, step=0.01, value=global_vars$default_values$evidence$prev_cs_A_B$A[1], width="100%")),
#           column(4, numericInput("cal_intercept_sd",label="SE of the intercept", value=global_vars$default_values$evidence$prev_cs_A_B$A[2])),
#           column(4, textOutput("cal_intercept_dist"))
#         ),
#         fluidRow(
#           column(4, sliderInput("cal_slope",label="Expected calibration slope", min=-2, max=2, step=0.01, value=global_vars$default_values$evidence$prev_cs_A_B$B[1], width="100%")),
#           column(4, numericInput("cal_slope_sd",label="SE of the slope", value=global_vars$default_values$evidence$prev_cs_A_B$B[2])),
#           column(4, textOutput("cal_slope_dist"))
#         )
#       ))
#       observeEvent(input$cstat_m_uci, {
#         res <- solve_beta_given_mean_quantile(input$cstat_m_uci[1],input$cstat_m_uci[2],0.975)
#         output$cstat_dist <- renderText(paste0("c~Beta(", round(res[1],2),",",round(res[2],2),") | ", make95CrIFromBeta(res[1],res[2])))
#       })
#       observeEvent(input$prev_m_uci, {
#         res <- solve_beta_given_mean_quantile(input$prev_m_uci[1]/100,input$prev_m_uci[2]/100,0.975)
#         output$prev_dist <- renderText(paste0("prev~Beta(", round(res[1],2),",",round(res[2],2),") | ", make95CrIFromBeta(res[1],res[2])))
#       })
#       shinyjs::enable("n_sim_outer")
#       output$n_sim_outer_note <- renderText("This value cannot be set to >100 for the web app. For larger values use the R package directly")
#     }
# 
#     if(choice=="3")
#     {
#       global_vars$params$evidence_type <<- "sample"
#       output$evidence_inputs <- renderUI(list(fileInput("sample_upload", NULL, buttonLabel = "Upload...", multiple = FALSE),
#                                               tableOutput("post_sample_upload")
#                                               )
#                                        )
#       observeEvent(input$sample_upload,{
#         ext <- tools::file_ext(input$sample_upload$name)
#         dt <- switch(ext,
#                  csv = vroom::vroom(input$sample_upload$datapath, delim = ","),
#                  tsv = vroom::vroom(input$sample_upload$datapath, delim = "\t"),
#                  validate("Invalid file; Please upload a .csv or .tsv file")
#         )
# 
#         global_vars$params$evidence <<- dt
#         if(is.na(sum(match(colnames(dt),c("se","sp","prev")))))
#         {
#           output$post_sample_upload <- renderText("Error: the data should have the following three columns: prev, se, sp")
#         }
#         else
#         {
#           output$post_sample_upload <- renderUI(list(renderTable(head(dt)), actionButton("sample_process","Process"), uiOutput("sample_process_results")))
#         }
#       })
# 
#       observeEvent(input$sample_process,{
#         output$sample_process_results <- renderUI(paste(dim(global_vars$params$evidence),collapse=","))
#       })
#     }
#   })
# 
# 
#   observeEvent(input$run1, {
#     res <- process_input(input)
# 
#     stuff_to_render <- list(renderUI(pre(global_vars$str_desc)), renderUI(pre(paste(deparse(global_vars$params),collapse=""))))
# 
#     if(res==T)
#     {
#       if(global_vars$params$evidence_type=="ind_beta")
#       {
#         VoI <- EVSI_ag(global_vars$params$evidence, global_vars$params$z, n_sim=global_vars$params$n_sim_inner, future_sample_sizes=c())
#       }
# 
#       if(global_vars$params$evidence_type=="prev_cs_A_B")
#       {
#         global_vars$params$sample <<- gen_triplets(global_vars$params$n_sim_outer, z=global_vars$params$z, prev=global_vars$params$evidence$prev, cs=global_vars$params$evidence$cs, A=global_vars$params$evidence$A, B=global_vars$params$evidence$B)
#         stuff_to_render <- c(stuff_to_render, renderText(paste("Results are based on a sample of", nrow(global_vars$params$sample), "draws.")))
#         VoI <- evsiexval::EVSI_g(global_vars$params$sample[,c('prev','se','sp')], global_vars$params$z, n_sim=1, future_sample_sizes=c())
#       }
# 
#       if(global_vars$params$evidence_type=="sample")
#       {
#         stuff_to_render <- c(stuff_to_render, renderText(paste("Results are based on a sample of", nrow(global_vars$params$sample), "draws.")))
#         VoI <- evsiexval::EVSI_g(global_vars$params$evidence[,c('prev','se','sp')], global_vars$params$z, n_sim=global_vars$params$n_sim_inner, future_sample_sizes=c(1))
#       }
# 
#       EVPI <- VoI$EVPI
#       summary <- VoI$summary
# 
#       require("knitr")
#       stuff_to_render <- c(stuff_to_render,list(
#         renderTable(summary, rownames=T, digits=6),
#         HTML(ifelse(max(summary['P_best',])>0.99, "<B style='color:red; font-weight:bold;'>In more than 99% of simulations the same stratgy had the highest NB. This indicates there is not much uncertainty around this decision. VoI analysis might be degenrate and non-informative.</B>","")),
#         hr(),
#         renderText(paste("EVPI=",format(EVPI, nsmall=5))),
#         renderText(paste("Population EVPI=",format(EVPI*global_vars$params$N, nsmall=2))),
#         hr(),
#         actionButton("evsi_run","Run EVSI analysis"),
#         uiOutput("results2")
#         ))
#       global_vars$result_level <<- 1
#     }
#     else
#     {
#       global_vars$result_level <<- 0
#       stuff_to_render <- c(stuff_to_render, renderText(paste("Error:", global_vars$str_val)))
#     }
#     output$results1 <- renderUI(stuff_to_render)
#     output$results2 <- renderUI("")
#   })
# 
# 
#   observeEvent(input$evsi_run, {
# 
#     # input_list <- reactiveValuesToList(input)
#     # toggle_inputs(input_list,F,F)
#     showModal(modalDialog("Computation in progress... Please do not change any input values before the computation is completed."))
#     res <- ENBS(global_vars$params)
#     removeModal()
# 
# 
#     require("knitr")
# 
#     output$results2 <- renderUI(list(
#       renderTable(res),
#       renderPlot({
#         plot(res$n_star, res$EVSI, type='l', xlab="Sample size of the future study", ylab="EVSI")
#         title("Expected Value of Sample Information (in true positive units)")
#         y2 <- pretty(c(0,global_vars$params$N*res$EVPI[1]))
#         axis(4, at=y2/global_vars$params$N, labels=y2)
#         mtext("Population EVSI", side = 4)
#         lines(c(0,max(res$n_star)), rep(res$EVPI[1],2), col="gray")
#       }),
# 
#       renderPlot({
#         plot(res$n_star, res$ENBS, type='l', xlab="Sample size of the future study", ylab="ENBS")
#         title("Expected Net benefit of Sampling (in true positive units)")
#         winner <- which.max(res$ENBS)
#         if(winner!=1 & winner!=length(res$ENBS))
#         {
#           lines(c(res$n_star[winner], res$n_star[winner]), c(0, max(res$ENBS)), col='red')
#           text(res$n_star[winner]*1.1,max(res$ENBS)/2,paste0("Optimal sample size:",res$n_stars[winner]), col='red')
#         }
#         else
#         {
#           text(0,max(res$ENBS)/2,"Edge case!", col='red')
#         }
#       })
#     ))
#     global_vars$result_level <<- 2
#   })
# }
# 
# 
# 
# toggle_inputs <- function(input_list,enable_inputs=T,only_buttons=FALSE)
# {
#   # Subset if only_buttons is TRUE.
# 
#   if(only_buttons){
#     buttons <- which(sapply(input_list,function(x) {any(grepl('Button',attr(x,"class")))}))
#     input_list = input_list[buttons]
#   }
# 
#   # Toggle elements
#   for(x in names(input_list))
#     if(enable_inputs){
#       shinyjs::enable(x)} else {
#         shinyjs::disable(x) }
}



shinyApp(ui = ui, server = server)
