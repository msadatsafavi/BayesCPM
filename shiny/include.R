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
                                          cal_int="Calibration intercept")
)



generate_code <- function(input)
{
  browser()
  evidence=list()
  evidence$prev$type <- input$evidence_prev_dist_type
  evidence_prev_input_type <- names(choices$evidence_prev_input_type)[match(input$evidence_prev_input_type,choices$evidence_prev_input_type)]
  p1 <- input$evidence_prev_input_parm2
  p2 <- input$evidence_prev_input_parm2
  if(evidence_prev_input_type=="pc_cih")
  {
    
  }
}
