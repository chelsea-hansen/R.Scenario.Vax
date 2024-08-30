#' Scenario Projections
#'
#' Run scenario projections for different immunization coverage
#'
#' @param fitted_parms A list of fitted parameters estimated using the fit_model function
#' @param parmset A list of fixed parameters
#' @param yinit A matrix of initial compartment values
#' @param yinit.vector A vector of initial compartment values
#' @param data_start The start data of the RSV time series data
#' @param projection_start The user defined start date of the projection period
#' @param projection_end The user defined end date of the projection period
#' @param adult_start The start date of senior vaccination for the projection period
#' @param adult_end The end date of senior vaccination for the projection period
#' @param adult75_doses The number of vaccine doses administered to adults 75+ during the period of vaccine administration
#' @param adult75_doses_last_year The number of vaccine doses administered to adults 75+ in the previous season
#' @param adult65_74_doses The number of vaccine doses administered to adults 65-74 during the period of vaccine administration
#' @param adult65_74_doses_last_year The number of vaccine doses administered to adults 65-74 in the previous season
#' @param maternal_start The start date of the maternal vaccination for the projection period
#' @param maternal_end The end date of maternal vaccination for the projection period
#' @param maternal_doses The number of doses administered to pregnant women during the administration period
#' @param monoclonal_catchup_start The start date of monoclonal antibody administration to infants <8 months
#' @param monoclonal_catchup_end The end date of monoclonal antibody administration to infants <8 months
#' @param monoclonal_catchup_doses The number of monoclonal antibody doses administered to infants <8 months (excluding birth doses)
#' @param monoclonal_birth_start The start date of monoclonal antibody birth doses
#' @param monoclonal_birth_end The end date of monoclonal antibody birth doses
#' @param monoclonal_birth_doses The number of monoclonal antibody doses administered as birth doses
#' @param scenario_name A user defined scenario name. Suggest to name Scenarios A,B,C, D, etc.
#' @param confidence_intervals Whether to calculate confidence intervals.Default = TRUE
#'
#' @return A data frame of weekly RSV hospitalizations by age group for the user defined projection period.
#' @export
#' @import deSolve
#' @import magrittr
#' @import dplyr
#' @import MMWRweek
#' @import ggplot2
#' @import cowplot
#' @importFrom pracma sigmoid
#' @importFrom stats median quantile
#' @importFrom tidyr pivot_longer
#'
#
#'
#' @examples
#' dat = get_data(state_or_county="state",state_abbr="CA",county_name=NULL)
#' parmset=dat[[1]]
#' yinit=dat[[2]]
#' yinit.vector=dat[[3]]
#' scenario_a = scenario_projection(fitted_parms=fitLL,
#'                             parmset = parmset,
#'                             yinit=yinit,
#'                             yinit.vector=yinit.vector,
#'                             data_start = '2017-07-08',
#'                             projection_start = '2024-10-01',
#'                             projection_end = '2025-06-01',
#'                             adult_start = '2024-08-01',
#'                             adult_end = '2025-05-01',
#'                             adult75_doses = 100000,
#'                             adult75_doses_last_year=250000,
#'                             adult65_74_doses = 100000,
#'                             adult65_74_doses_last_year=250000,
#'                             maternal_start = '2024-09-01',
#'                             maternal_end = '2025-04-01',
#'                             maternal_doses = 50000,
#'                             monoclonal_catchup_start = '2024-10-01',
#'                             monoclonal_catchup_end = '2025-04-01',
#'                             monoclonal_catchup_doses = 50000,
#'                             monoclonal_birth_start = '2024-10-01',
#'                             monoclonal_birth_end = '2025-04-01',
#'                             monoclonal_birth_doses = 50000,
#'                             scenario_name="Scenario A",
#'                             confidence_intervals=TRUE)
scenario_projection = function(fitted_parms,
                                parmset,
                                yinit,
                                yinit.vector,
                                 data_start,
                                 projection_start,
                                 projection_end,
                                 adult_start,
                                 adult_end,
                                 adult75_doses,
                                 adult75_doses_last_year,
                                 adult65_74_doses,
                                 adult65_74_doses_last_year,
                                 maternal_start,
                                 maternal_end,
                                 maternal_doses,
                                 monoclonal_catchup_start,
                                 monoclonal_catchup_end,
                                 monoclonal_catchup_doses,
                                 monoclonal_birth_start,
                                 monoclonal_birth_end,
                                 monoclonal_birth_doses,
                                 scenario_name,
                                 confidence_intervals=TRUE){


  dates = data.frame(dates=seq(from=as.Date(data_start), to=as.Date(projection_end), by="weeks")) %>%
    mutate(MMWRweek(dates),
           new_date=MMWRweek2Date(MMWRweek=.data$MMWRweek,MMWRyear=.data$MMWRyear)) %>%
    select("new_date")

  fit_times = seq(1, nrow(dates)+104, by=1)

  sen75.1 = seq(from=as.Date(data_start), to=as.Date(adult_start), by="weeks")
  sen75.2 = seq(from=as.Date(adult_start), to=as.Date(adult_end), by="weeks")
  sen75.3 = seq(from=as.Date(adult_end), to=as.Date(projection_end), by="weeks")
  sen75.4 = c(sigmoid(1:length(sen75.2), a = .26679, b = 12.50297),1)
  sen75.5  =c(0,lag(sen75.4,1)[c(-1)])
  sen75.6 = sen75.4-sen75.5
  #sen6 = sen6[c(-1)]
  adult75_vax = c(rep(0, length(sen75.1)+104), adult75_doses_last_year,sen75.6*adult75_doses, rep(0, length(sen75.3)))[1:length(fit_times)]
  #senior_vax
  #plot(senior_vax)
  adult75_cum = cumsum(adult75_vax)
 # plot(adult75_cum)
  #senior_cum

  sen65.1 = seq(from=as.Date(data_start), to=as.Date(adult_start), by="weeks")
  sen65.2 = seq(from=as.Date(adult_start), to=as.Date(adult_end), by="weeks")
  sen65.3 = seq(from=as.Date(adult_end), to=as.Date(projection_end), by="weeks")
  sen65.4 = c(sigmoid(1:length(sen65.2), a = .26679, b = 12.50297),1)
  sen65.5  =c(0,lag(sen65.4,1)[c(-1)])
  sen65.6 = sen65.4-sen65.5
  #sen6 = sen6[c(-1)]
  adult65_vax = c(rep(0, length(sen65.1)+104), adult65_74_doses_last_year,sen65.6*adult65_74_doses, rep(0, length(sen75.3)))[1:length(fit_times)]
  #senior_vax
  #plot(senior_vax)
  adult65_cum = cumsum(adult65_vax)
  #plot(adult65_cum)
  #senior_cum



  mat1 = seq(from=as.Date(data_start), to=as.Date(maternal_start), by="weeks")
  mat2 = seq(from=as.Date(maternal_start), to=as.Date(maternal_end), by="weeks")
  mat3 = seq(from=as.Date(maternal_end), to=as.Date(projection_end), by="weeks")
  mat4 = c(sigmoid(1:length(mat2), a = .26679, b = 12.50297),1)
  mat5  =c(0,lag(mat4,1)[c(-1)])
  mat6 = mat4-mat5
  #mat6 = mat6[c(-1)]
  maternal_vax = c(rep(0, length(mat1)+110),mat6*maternal_doses, rep(0, length(mat3)))[1:length(fit_times)]
  mat_cum = cumsum(maternal_vax)
 # plot(mat_cum)
  #mat_cum
  #plot(tail(mat_cum,50))


  mon1 = seq(from=as.Date(data_start), to=as.Date(monoclonal_catchup_start), by="weeks")
  mon2 = seq(from=as.Date(monoclonal_catchup_start), to=as.Date(monoclonal_catchup_end), by="weeks")
  mon3 = seq(from=as.Date(monoclonal_catchup_end), to=as.Date(projection_end), by="weeks")
  mon4 = c(sigmoid(1:length(mon2), a = .26679, b = 12.50297),1)
  mon5  =c(0,lag(mon4,1)[c(-1)])
  mon6 = mon4-mon5
  #mon6 = mon6[c(-1)]
  monoclonal_catchup = c(rep(0, length(mon1)+104),mon6*monoclonal_catchup_doses, rep(0, length(mon3)))[1:length(fit_times)]
  mon_cum = cumsum(monoclonal_catchup)
 # plot(mon_cum)
  #mon_cum
  #plot(tail(mon_cum,50))

  bir1 = seq(from=as.Date(data_start), to=as.Date(monoclonal_birth_start), by="weeks")
  bir2 = seq(from=as.Date(monoclonal_birth_start), to=as.Date(monoclonal_birth_end), by="weeks")
  bir3 = seq(from=as.Date(monoclonal_birth_end), to=as.Date(projection_end), by="weeks")
  bir4 = c(sigmoid(1:length(bir2)+1, a = .26679, b = 12.50297),1)
  bir5  =c(0,lag(bir4,1)[c(-1)])
  bir6 = bir4-bir5
  #bir6 = bir6[c(-1)]
  monoclonal_birth = c(rep(0, length(bir1)+104),bir6*monoclonal_birth_doses, rep(0, length(bir3)))[1:length(fit_times)]
  bir_cum = cumsum(monoclonal_birth)
  #plot(bir_cum)
  #bir_cum
  #plot(tail(bir_cum,50))

  #check = cbind(senior_vax, maternal_vax, monoclonal_catchup, monoclonal_birth)[53:465,1:4]
  #check2 = cbind(dates, check)


  #parmset$monoclonal_birth = monoclonal_birth
  #parmset$monoclonal_catchup = monoclonal_catchup
  parmset$monoclonal_01 = monoclonal_birth
  parmset$monoclonal_23 = monoclonal_catchup*.34
  parmset$monoclonal_45 = monoclonal_catchup*.33
  parmset$monoclonal_67 = monoclonal_catchup*.33
  parmset$maternal_vax = maternal_vax
  parmset$senior_vax_75 = adult75_vax
  parmset$senior_vax_65_74 = adult65_vax

  if(confidence_intervals==TRUE){
    newH=data.frame(date=NA, sample=NA,Age=NA,value=NA)
    for(l in 1:100){
      parmset$baseline.txn.rate=fitted_parms[[9]][l,"beta"]
      parmset$b1=fitted_parms[[9]][l,"b1"]
      parmset$phi=fitted_parms[[9]][l,"phi"]
      report_infants = fitted_parms[[4]]
      report_children = fitted_parms[[5]]
      report_adults = fitted_parms[[6]]
      report_seniors65 = fitted_parms[[7]]
      report_seniors75 = fitted_parms[[8]]

      output <- ode(y=yinit.vector, times=fit_times,method = "ode45",
                    func=MSIRS_immunization_dynamics,
                    parms=parmset)

      t0=nrow(dates)
      al <- nrow(yinit)
      output <- tail(output,t0)
      St <- output[,-1]
      I1 <- St[,grep('I1', colnames(St))]
      I2 <- St[,grep('I2', colnames(St))]
      I3 <- St[,grep('I3', colnames(St))]
      I4 <- St[,grep('I4', colnames(St))]
      S1 <- St[,grep('S1', colnames(St))]
      S2 <- St[,grep('S2', colnames(St))]
      S3 <- St[,grep('S3', colnames(St))]
      S0 <- St[,grep('S0', colnames(St))]
      R1 <- St[,grep('R1', colnames(St))]
      R2 <- St[,grep('R2', colnames(St))]
      R3 <- St[,grep('R3', colnames(St))]
      R4 <- St[,grep('R4', colnames(St))]
      M0<- St[,grep('M0', colnames(St))]
      Si<- St[,grep('Si', colnames(St))]
      Mn<- St[,grep('Mn', colnames(St))]
      Mv<- St[,grep('Mv', colnames(St))]
      N<- St[,grep('N', colnames(St))]
      Vs1<- St[,grep('Vs1', colnames(St))]
      Vs2<- St[,grep('Vs2', colnames(St))]

      beta <-  parmset$baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2


      lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
      for (t in 1:t0){
        lambda1[t,] <- as.vector((1+parmset$b1*cos(2*pi*(t-parmset$phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+parmset$seed)%*%beta)/sum(St[t,]))}


      hosp1 <- c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2, report_infants*0.15, report_infants*0.15, rep(report_children, 2), rep(.001, 6))
      hosp2 <- hosp1 * 0.4
      hosp3 <- c(rep(0.001, 8), rep(report_adults, 4), report_seniors65,report_seniors75)


      H1 <- matrix(0, nrow = t0, ncol = al)
      for (i in 1:al) {
        H1[, i] <-
          hosp1[i] * parmset$RRHm * parmset$sigma3 * M0[, i] * lambda1[, i] +
          hosp1[i] * parmset$RRHm * parmset$sigma3 * Mn[, i] * lambda1[, i] +
          hosp1[i] * parmset$RRHm * parmset$sigma3 * Mv[, i] * lambda1[, i] +
          hosp1[i] * N[, i] * lambda1[, i] +
          hosp1[i] * S0[, i] * lambda1[, i] +
          hosp1[i] * Si[, i] * lambda1[, i] +
          hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
          hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
          hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
          hosp3[i] * parmset$sigma3 * Vs1[, i] * lambda1[, i] +
          hosp3[i] * parmset$sigma3 * Vs2[, i] * lambda1[, i]
      }

      H2 = cbind(rowSums(H1[,1:3]),
                 rowSums(H1[,4:6]),
                 rowSums(H1[,7:8]),
                 rowSums(H1[,9:12]),
                 H1[,13],H1[,14])



      H = data.frame(H2)
      names(H)=c("<6m","6-11m","1-4yrs","5-64yrs","65-74yrs","75+yrs")
      H$All = rowSums(H2)
      H$date = dates$new_date

      H = H %>% filter(.data$date>=projection_start) %>%
        mutate(sample=l) %>%
        pivot_longer(cols=c("<6m":"All"),names_to="Age",values_to="value")


      newH=rbind(newH,H)


    }

    newH=newH %>%  filter(!is.na(.data$date)) %>% mutate(date=as.Date(.data$date))

    results = newH %>%
      group_by(.data$date, .data$Age) %>%
      summarize(median = median(.data$value),
                lower = quantile(.data$value, probs=0.025),
                upper = quantile(.data$value, probs=0.975)) %>%
      mutate(scenario = scenario_name)


    totals = results %>%
      group_by(.data$Age) %>%
      summarize(median = sum(.data$median),
                lower = sum(.data$lower),
                upper = sum (.data$upper)) %>%
      mutate(scenario=scenario_name,
             Age=factor(.data$Age,levels=c("<6m","6-11m","1-4yrs","5-64yrs","65-74yrs","75+yrs","All")))



    plot1 = ggplot(data=results %>% filter(.data$Age=="All"))+
      theme_bw()+
      geom_ribbon(aes(x=date, ymin=.data$lower, ymax=.data$upper),color="steelblue",alpha=0.3)+
      geom_line(aes(x=date, y=.data$median),color="navy",linewidth=1.5)+
      labs(x=NULL, y="Weekly RSV Hospitalizations/ED Visits")

    plot2 = ggplot(data=totals)+
      theme_bw()+
      geom_bar(aes(x=.data$Age, y=.data$median),stat="identity",fill="navy")+
      geom_errorbar(aes(x=.data$Age,ymin=.data$lower,ymax=.data$upper),width=0.2)+
      geom_text(aes(x=.data$Age, y=.data$upper+200,label=round(.data$median)))+
      labs(x=NULL,y="Total RSV Hospitalizations/ED Visits")
    plot3 = plot_grid(plot1,plot2, ncol=2)
    print(plot3)

  }else{
  parmset$baseline.txn.rate = fitted_parms[[1]]
  parmset$b1 = fitted_parms[[2]]
  parmset$phi = fitted_parms[[3]]
  report_infants = fitted_parms[[4]]
  report_children = fitted_parms[[5]]
  report_adults = fitted_parms[[6]]
  report_seniors65 = fitted_parms[[7]]
  report_seniors75 = fitted_parms[[8]]

  output <- ode(y=yinit.vector, times=fit_times,method = "ode45",
                func=MSIRS_immunization_dynamics,
                parms=parmset)

  t0=nrow(dates)
  al <- nrow(yinit)
  output <- tail(output,t0)
  St <- output[,-1]
  I1 <- St[,grep('I1', colnames(St))]
  I2 <- St[,grep('I2', colnames(St))]
  I3 <- St[,grep('I3', colnames(St))]
  I4 <- St[,grep('I4', colnames(St))]
  S1 <- St[,grep('S1', colnames(St))]
  S2 <- St[,grep('S2', colnames(St))]
  S3 <- St[,grep('S3', colnames(St))]
  S0 <- St[,grep('S0', colnames(St))]
  R1 <- St[,grep('R1', colnames(St))]
  R2 <- St[,grep('R2', colnames(St))]
  R3 <- St[,grep('R3', colnames(St))]
  R4 <- St[,grep('R4', colnames(St))]
  M0<- St[,grep('M0', colnames(St))]
  Si<- St[,grep('Si', colnames(St))]
  Mn<- St[,grep('Mn', colnames(St))]
  Mv<- St[,grep('Mv', colnames(St))]
  N<- St[,grep('N', colnames(St))]
  Vs1<- St[,grep('Vs1', colnames(St))]
  Vs2<- St[,grep('Vs2', colnames(St))]

  beta <-  parmset$baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2


  lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
  for (t in 1:t0){
    lambda1[t,] <- as.vector((1+parmset$b1*cos(2*pi*(t-parmset$phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+parmset$seed)%*%beta)/sum(St[t,]))}


  hosp1 <- c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2, report_infants*0.15, report_infants*0.15, rep(report_children, 2), rep(.001, 6))
  hosp2 <- hosp1 * 0.4
  hosp3 <- c(rep(0.001, 8), rep(report_adults, 4), report_seniors65,report_seniors75)


  H1 <- matrix(0, nrow = t0, ncol = al)
  for (i in 1:al) {
    H1[, i] <-
      hosp1[i] * parmset$RRHm * parmset$sigma3 * M0[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$sigma3 * Mn[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$sigma3 * Mv[, i] * lambda1[, i] +
      hosp1[i] * N[, i] * lambda1[, i] +
      hosp1[i] * S0[, i] * lambda1[, i] +
      hosp1[i] * Si[, i] * lambda1[, i] +
      hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * Vs1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * Vs2[, i] * lambda1[, i]
  }

  H2 = cbind(rowSums(H1[,1:3]),
             rowSums(H1[,4:6]),
             rowSums(H1[,7:8]),
             rowSums(H1[,9:12]),
             H1[,13],H1[,14])



  H = data.frame(H2)
  names(H)=c("<6m","6-11m","1-4yrs","5-64yrs","65-74yrs","75+yrs")
  H$All = rowSums(H2)
  H$date = dates$new_date

results = H %>% filter(.data$date>=projection_start)


 totals = results %>%
   pivot_longer(cols=c("<6m":"All"),names_to="Age",values_to="value") %>%
   group_by(.data$Age) %>%
   summarize(total = sum(.data$value)) %>%
   mutate(scenario = scenario_name,
          Age=factor(.data$Age,levels=c("<6m","6-11m","1-4yrs","5-64yrs","65-74yrs","75+yrs","All")))

  plot1 = ggplot(data=results)+
    theme_bw()+
    geom_line(aes(x=.data$date, y=.data$All),color="navy",linewidth=1.5)+
    labs(x=NULL, y="Weekly RSV Hospitalizations/ED Visits")

  plot2 = ggplot(data=totals)+
    theme_bw()+
    geom_bar(aes(x=.data$Age, y=.data$total),stat="identity",fill="navy")+
    geom_text(aes(x=.data$Age, y=.data$total+200,label=round(.data$total)))+
    labs(x=NULL,y="Total RSV Hospitalizations/ED Visits")

  plot3 = plot_grid(plot1,plot2, ncol=2)
  print(plot3)
}
  return(results)

}

