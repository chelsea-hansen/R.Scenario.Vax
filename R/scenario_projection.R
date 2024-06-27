#' Scenario Projection - No uncertainty intervals
#'
#' @param fitted_parms A list of fitted parameters estimated using the fit_model function
#' @param parmset A list of fixed parameters
#' @param yinit A matrix of initial compartment values
#' @param yinit.vector A vector of initial compartment values
#' @param data_start The start data of the RSV time series data
#' @param projection_start The user defined start date of the projection period
#' @param projection_end The user defined end date of the projection period
#' @param senior_start The start date of senior vaccination for the projection period
#' @param senior_end The end date of senior vaccination for the projection period
#' @param senior_doses The number of vaccine doses administered to seniors during the period of vaccine administration
#' @param senior_doses_last_year The number of vaccine doses administered to seniors in the previous season
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
#'
#' @return A data frame of weekly RSV hospitalizations by age group for the user defined projection period.
#' @export
#' @import deSolve
#' @import magrittr
#' @import dplyr
#' @import ggplot2
#' @import cdcfluview
#' @import cowplot
#' @importFrom pracma sigmoid
#
#'
#' @examples
#' dat = get_data(state_or_county="state",state_abbr="CA",county_name=NULL)
#' parmset=dat[[1]]
#' yinit=dat[[2]]
#' yinit.vector=dat[[3]]
#' scenario_a = scenario_projection1(fitted_parms=fitLL,
#'                             parmset = parmset,
#'                             yinit=yinit,
#'                             yinit.vector=yinit.vector,
#'                             data_start = '2017-07-08',
#'                             projection_start = '2024-10-01',
#'                             projection_end = '2025-06-01',
#'                             senior_start = '2024-08-01',
#'                             senior_end = '2025-05-01',
#'                             senior_doses = .15,
#'                             senior_doses_last_year=.15,
#'                             maternal_start = '2024-09-01',
#'                             maternal_end = '2025-04-01',
#'                             maternal_doses = 50000,
#'                             monoclonal_catchup_start = '2024-10-01',
#'                             monoclonal_catchup_end = '2025-04-01',
#'                             monoclonal_catchup_doses = 50000,
#'                             monoclonal_birth_start = '2024-10-01',
#'                             monoclonal_birth_end = '2025-04-01',
#'                             monoclonal_birth_doses = 50000,
#'                             scenario_name="Scenario A")
scenario_projection = function(fitted_parms,
                                parmset,
                                yinit,
                                yinit.vector,
                                 data_start,
                                 projection_start,
                                 projection_end,
                                 senior_start,
                                 senior_end,
                                 senior_doses,
                                 senior_doses_last_year,
                                 maternal_start,
                                 maternal_end,
                                 maternal_doses,
                                 monoclonal_catchup_start,
                                 monoclonal_catchup_end,
                                 monoclonal_catchup_doses,
                                 monoclonal_birth_start,
                                 monoclonal_birth_end,
                                 monoclonal_birth_doses,
                                 scenario_name){


  dates = data.frame(dates=seq(from=as.Date(data_start), to=as.Date(projection_end), by="weeks")) %>%
    mutate(mmwr_week(dates),
           new_date=mmwr_week_to_date(week=.data$mmwr_week,year=.data$mmwr_year)) %>%
    select("new_date")

  fit_times = seq(1, nrow(dates)+52, by=1)

  sen1 = seq(from=as.Date(data_start), to=as.Date(senior_start), by="weeks")
  sen2 = seq(from=as.Date(senior_start), to=as.Date(senior_end), by="weeks")
  sen3 = seq(from=as.Date(senior_end), to=as.Date(projection_end), by="weeks")
  sen4 = c(sigmoid(1:length(sen2), a = .26679, b = 12.50297),1)
  sen5  =c(0,lag(sen4,1)[c(-1)])
  sen6 = sen4-sen5
  #sen6 = sen6[c(-1)]
  senior_vax = c(rep(0, length(sen1)+52), senior_doses_last_year,sen6*senior_doses, rep(0, length(sen3)))[1:length(fit_times)]
  #senior_vax
  #plot(senior_vax)
  senior_cum = cumsum(senior_vax)
  #senior_cum

  mat1 = seq(from=as.Date(data_start), to=as.Date(maternal_start), by="weeks")
  mat2 = seq(from=as.Date(maternal_start), to=as.Date(maternal_end), by="weeks")
  mat3 = seq(from=as.Date(maternal_end), to=as.Date(projection_end), by="weeks")
  mat4 = c(sigmoid(1:length(mat2), a = .26679, b = 12.50297),1)
  mat5  =c(0,lag(mat4,1)[c(-1)])
  mat6 = mat4-mat5
  #mat6 = mat6[c(-1)]
  maternal_vax = c(rep(0, length(mat1)+58),mat6*maternal_doses, rep(0, length(mat3)))[1:length(fit_times)]
  mat_cum = cumsum(maternal_vax)
  #mat_cum
  #plot(tail(mat_cum,50))


  mon1 = seq(from=as.Date(data_start), to=as.Date(monoclonal_catchup_start), by="weeks")
  mon2 = seq(from=as.Date(monoclonal_catchup_start), to=as.Date(monoclonal_catchup_end), by="weeks")
  mon3 = seq(from=as.Date(monoclonal_catchup_end), to=as.Date(projection_end), by="weeks")
  mon4 = c(sigmoid(1:length(mon2), a = .26679, b = 12.50297),1)
  mon5  =c(0,lag(mon4,1)[c(-1)])
  mon6 = mon4-mon5
  #mon6 = mon6[c(-1)]
  monoclonal_catchup = c(rep(0, length(mon1)+52),mon6*monoclonal_catchup_doses, rep(0, length(mon3)))[1:length(fit_times)]
  mon_cum = cumsum(monoclonal_catchup)
  #mon_cum
  #plot(tail(mon_cum,50))

  bir1 = seq(from=as.Date(data_start), to=as.Date(monoclonal_birth_start), by="weeks")
  bir2 = seq(from=as.Date(monoclonal_birth_start), to=as.Date(monoclonal_birth_end), by="weeks")
  bir3 = seq(from=as.Date(monoclonal_birth_end), to=as.Date(projection_end), by="weeks")
  bir4 = c(sigmoid(1:length(bir2)+1, a = .26679, b = 12.50297),1)
  bir5  =c(0,lag(bir4,1)[c(-1)])
  bir6 = bir4-bir5
  #bir6 = bir6[c(-1)]
  monoclonal_birth = c(rep(0, length(bir1)+52),bir6*monoclonal_birth_doses, rep(0, length(bir3)))[1:length(fit_times)]
  bir_cum = cumsum(monoclonal_birth)
  #bir_cum
  #plot(tail(bir_cum,50))

  check = cbind(senior_vax, maternal_vax, monoclonal_catchup, monoclonal_birth)[53:465,1:4]
  check2 = cbind(dates, check)


  parmset$monoclonal_birth = monoclonal_birth
  parmset$monoclonal_catchup = monoclonal_catchup
  parmset$maternal_vax = maternal_vax
  parmset$senior_vax = senior_vax

  baseline.txn.rate = fitted_parms[[1]]
  b1 = fitted_parms[[2]]
  phi = fitted_parms[[3]]
  hosp_prop = fitted_parms[[5]]

  output <- ode(y=yinit.vector, times=fit_times,method = "ode45",
                func=MSIRS_immunization_dynamics,
                parms=c(parmset,
                        baseline.txn.rate=baseline.txn.rate,
                        b1=b1,
                        phi=phi))

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

  beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2


  lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
  for (t in 1:t0){
    lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+parmset$seed)%*%beta)/sum(St[t,]))}


  hosp = c(rep(hosp_prop[1],3),rep(hosp_prop[2],3),rep(hosp_prop[3],2),rep(hosp_prop[4],4),hosp_prop[5])
  H1=matrix(0,nrow=t0,ncol=al)#Number of hospitalizations by age
  for (i in 1:al){
    H1[,i]=
      hosp[i]*parmset$RRHm*parmset$sigma3*M0[,i]*lambda1[,i]+
      hosp[i]*parmset$RRHn*parmset$RRIn*parmset$sigma3*Mn[,i]*lambda1[,i]+
      hosp[i]*parmset$RRHm*parmset$RRHv*parmset$RRIv*parmset$sigma3*Mv[,i]*lambda1[,i]+
      hosp[i]*parmset$RRHn*parmset$RRIn*N[,i]*lambda1[,i]+
      hosp[i]*S0[,i]*lambda1[,i]+
      hosp[i]*Si[,i]*lambda1[,i]+
      hosp[i]*parmset$sigma1*S1[,i]*lambda1[,i]+
      hosp[i]*parmset$sigma2*S2[,i]*lambda1[,i]+
      hosp[i]*parmset$sigma3*S3[,i]*lambda1[,i]+
      hosp[i]*parmset$RRIs*parmset$RRHs*parmset$sigma3*Vs1[,i]*lambda1[,i]+
      hosp[i]*parmset$RRIs*parmset$RRHs*parmset$sigma3*Vs2[,i]*lambda1[,i]}


  H2 = cbind(rowSums(H1[,1:3]),
             rowSums(H1[,4:6]),
             rowSums(H1[,7:8]),
             rowSums(H1[,9:12]),
             H1[,13])

  # H3 = t(t(H2) * hosp_prop)


  H = data.frame(H2)
  H$All = rowSums(H2)
  age_list = c("<6m","5-11m","1-4yrs","5-59yrs","60+yrs","All")
  names(H)=age_list
  H$date = dates$new_date
  H = H %>% filter(date>=projection_start) %>%
    mutate(scenario=scenario_name)

  totals = H %>% select(-"date", -"scenario")
  totals=data.frame(total = colSums(totals)) %>%
    mutate(age = age_list,
           age=factor(.data$age,levels=c("<6m","6-11m","1-4yrs","5-59yrs","60+yrs","All")))

  plot1 = ggplot(data=H)+
    theme_bw()+
    geom_line(aes(x=.data$date, y=.data$All),color="navy",linewidth=1.5)+
    labs(x=NULL, y="Weekly RSV Hospitalizations/ED Visits")

  plot2 = ggplot(data=totals)+
    theme_bw()+
    geom_bar(aes(x=.data$age, y=.data$total),stat="identity",fill="navy")+
    geom_text(aes(x=.data$age, y=.data$total+200,label=round(.data$total)))+
    labs(x=NULL,y="Total RSV Hospitalizations/ED Visits")

  plot3 = plot_grid(plot1,plot2, ncol=2)
  print(plot3)

  return(H)

}

