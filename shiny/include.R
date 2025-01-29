choices <- list(evidence_prev_input_type=c(  pe_cih="Point estimate and upper 95%CI bound",
                                             pe_ss="Point estimate and sample size",
                                             m_sd="Mean and SD",
                                             parms="Distribution parameters"),
                evidence_cstat_input_type=c( pe_cih="Point estimate and upper 95%CI bound", 
                                             m_sd="Mean and SD",
                                             parms="Distribution parameters"),
                evidence_cal_slp_input_type=c(pe_cih="Point estimate and upper 95%CI bound", 
                                              m_sd="Mean and SD"),
                evidence_cal_other_type=c("",
                                          cal_oe="O/E ratio",
                                          cal_mean="Mean calibration",
                                          cal_int="Calibration intercept"),
                evidence_cal_other_input_type=c(pe_cih="Point estimate and upper 95%CI bound", 
                                                m_sd="Mean and SD",
                                                parms="Distribution parameters"),
                method_type=c(sample="Simulate validation data", 
                              "2s"="Two-level SE equations")
            )



collect_evidence <- function(input)
{
  evidence <- list()
  
  #Prevalence
  evidence$prev$type <- input$evidence_prev_dist_type
  evidence_prev_input_type <- names(choices$evidence_prev_input_type)[match(input$evidence_prev_input_type,choices$evidence_prev_input_type)]
  p1 <- input$evidence_prev_input_parm1
  p2 <- input$evidence_prev_input_parm2
  if(evidence_prev_input_type=="pe_cih")
  {
    evidence$prev <- c(evidence$prev, mean=p1, cih=p2)
  }
  if(evidence_prev_input_type=="pe_ss") #This one is unconventional and does not have a counterpart in the package
  {
    evidence$prev <- c(evidence$prev, mean=p1, var=p1*(1-p1)/p2)
  }
  if(evidence_prev_input_type=="m_sd")
  {
    evidence$prev <- c(evidence$prev, mean=p1, sd=p2)
  }
  if(evidence_prev_input_type=="parms")
  {
    evidence$prev <- c(evidence$prev, parm1=p1, parm2=p2)
  }
  
  #cstat
  evidence$cstat$type <- input$evidence_cstat_dist_type
  evidence_cstat_input_type <- names(choices$evidence_cstat_input_type)[match(input$evidence_cstat_input_type,choices$evidence_cstat_input_type)]
  p1 <- input$evidence_cstat_input_parm1
  p2 <- input$evidence_cstat_input_parm2
  if(evidence_cstat_input_type=="pe_cih")
  {
    evidence$cstat <- c(evidence$cstat, mean=p1, cih=p2)
  }
  if(evidence_cstat_input_type=="m_sd")
  {
    evidence$cstat <- c(evidence$cstat, mean=p1, sd=p2)
  }
  if(evidence_cstat_input_type=="parms")
  {
    evidence$cstat <- c(evidence$cstat, parm1=p1, parm2=p2)
  }
  
  #cal_slp
  evidence$cal_slp$type <- input$evidence_cal_slp_dist_type
  evidence_cal_slp_input_type <- names(choices$evidence_cal_slp_input_type)[match(input$evidence_cal_slp_input_type,choices$evidence_cal_slp_input_type)]
  p1 <- input$evidence_cal_slp_input_parm1
  p2 <- input$evidence_cal_slp_input_parm2
  if(evidence_cal_slp_input_type=="pe_cih")
  {
    evidence$cal_slp<-c(evidence$cal_slp, mean=p1, cih=p2)
  }
  if(evidence_cal_slp_input_type=="m_sd")
  {
    evidence$cal_slp<-c(evidence$cal_slp, mean=p1, sd=p2)
  }
  if(evidence_cal_slp_input_type=="parms")
  {
    evidence$cal_slp<-c(evidence$cal_slp, parm1=p1, parm2=p2)
  }
  
  
  #cal_other
  cal_other_name <- names(choices$evidence_cal_other_type)[match(input$evidence_cal_other_type,choices$evidence_cal_other_type)]
  evidence[[cal_other_name]] <- list()
  evidence[[cal_other_name]]$type <- input$evidence_cal_other_dist_type
  evidence_cal_other_input_type <- names(choices$evidence_cal_other_input_type)[match(input$evidence_cal_other_input_type,choices$evidence_cal_other_input_type)]
  p1 <- input$evidence_cal_other_input_parm1
  p2 <- input$evidence_cal_other_input_parm2
  if(evidence_cal_other_input_type=="pe_cih")
  {
    evidence[[cal_other_name]]<-c(evidence[[cal_other_name]], mean=p1, cih=p2)
  }
  if(evidence_cal_other_input_type=="m_sd")
  {
    evidence[[cal_other_name]]<-c(evidence[[cal_other_name]], mean=p1, sd=p2)
  }
  if(evidence_cal_other_input_type=="parms")
  {
    evidence[[cal_other_name]]<-c(evidence[[cal_other_name]],parm1=p1, parm2=p2)
  }
  
  evidence
}



collect_targets <- function(input)
{
  targets <- list()
  if(input$b_target_cstat)
  {
    targets$cstat <- input$ciw_cstat
  }
  if(input$b_target_cal_oe)
  {
    targets$cal_oe <- input$ciw_cal_oe
  }
  if(input$b_target_cal_mean)
  {
    targets$cal_mean <- input$ciw_cal_mean
  }
  if(input$b_target_cal_slp)
  {
    targets$cal_slp <- input$ciw_cal_slp
  }
  if(input$b_target_cal_int)
  {
    targets$cal_int <- input$ciw_cal_int
  }
  targets
}





collect_rules <- function(input)
{
  rules <- list()
  if(input$b_fciw)
  {
    rules$fciw <- TRUE
  }
  if(input$b_eciw)
  {
    rules$eciw <- TRUE
  }
  if(input$b_qciw)
  {
    rules$qciw<-input$qciw
  }
  if(input$b_nb_voi)
  {
    rules$nb_voi <- TRUE
  }
  if(input$b_nb_assurance)
  {
    rules$nb_assurance <- input$nb_assurance
  }
  rules
}



gen_args <- function(input)
{
  out <- list()
  
  out$N <- eval(parse(text = paste("c(",input$N,")")))
  
  out$dist_type <- input$dist_type
  
  out$impute_cor <- input$b_impute_cor
  
  out$n_sim <- input$n_sim
  
  out$method=names(choices$method_type)[match(input$method,choices$method_type)]
  
  out$evidence <- collect_evidence(input)
  
  out$target_ciws <- collect_targets(input)
  
  out$rules <- collect_rules(input)
  
  if(!is.null(input$threshold))
  {
    out$threshold <- input$threshold
  }
  
  out
}