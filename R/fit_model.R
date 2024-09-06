#' Fit model parameters (function)
#'
#' This function uses maximum likelihood estimation to fit model parameters. The function uses the optim function from the stats package in R.
#'
#' @param time_series An RSV hospitalization time series. Make sure this is rounded to whole numbers.
#' @param age_dist The proportion of RSV hospitalizations in each age group (age groups include: <6 months,6-11 months,1-4 years,5-64 years, 65-74 years, 75+ years)
#' @param parmset A list of fixed parameters (retrieved using the get_data() function)
#' @param yinit A matrix of initial values for each model compartment(retrieved from get_data() function)
#' @param yinit.vector A vector of initial values for each model compartment (retrieved from get_data() function)
#'
#' @return A list of fitted parameter values. Prints a plot showing the fitted model.
#' @export
#' @import deSolve
#' @import magrittr
#' @import dplyr
#' @import ggplot2
#' @import lhs
#' @importFrom utils tail head
#' @importFrom stats dpois optim dmultinom
#
#' @examples
#' dat = get_data(state_or_county="state",state_abbr="CA",county_name=NULL)
#' parmset=dat[[1]]
#' yinit=dat[[2]]
#' yinit.vector=dat[[3]]
#'
#' weekly_rsv = round(timeseries[which(timeseries$state=="California"),"value"])
#'
#' fitLL = fit_model(time_series = weekly_rsv$value, age_dist = c(.19,.08,.23,.17,.10,.23),
#' parmset=parmset,yinit=yinit,yinit.vector=yinit.vector)

fit_model = function(time_series, age_dist, parmset, yinit,yinit.vector){

  fit_times = seq(1,length(time_series)+104,by=1)
  #parmset$monoclonal_01 = rep(0,length(time_series)+104)
  #parmset$monoclonal_23 = rep(0,length(time_series)+104)
  #parmset$monoclonal_45 = rep(0,length(time_series)+104)
  #parmset$monoclonal_67 = rep(0,length(time_series)+104)
  parmset$monoclonal_birth = rep(0,length(time_series)+104)
  parmset$monoclonal_catchup = rep(0,length(time_series)+104)
  parmset$maternal_vax = rep(0,length(time_series)+104)
  parmset$senior_vax_65_74 = rep(0,length(time_series)+104)
  parmset$senior_vax_75 = rep(0,length(time_series)+104)

  fitmodel <-  function(parameters,dat) {
    protrans <- parameters[1] # parameter for baseline transmission rate
    baseline.txn.rate = 5+(10*(exp(protrans))) / (1+exp(protrans))
    amp <- parameters[2] # parameter for seasonal amplitude
    b1 <-  .08+(.35*(exp(amp))) / (1+exp(amp))
    trans <- parameters[3] # parameter for seasonal phase
    phi <-  (2*pi*(exp(trans))) / (1+exp(trans)) # transform to between 0 and 2pi
    report_infants <- 1 / (1 + exp(-parameters[4]))
    report_children <- 1 / (1 + exp(-parameters[5]))
    report_adults <- 1 / (1 + exp(-parameters[6]))
    report_seniors65 <- 1 / (1 + exp(-parameters[7]))
    report_seniors75 <- 1 / (1 + exp(-parameters[8]))


    # Simulate the model with initial conditions and timesteps defined above, and parameter values from function call
    results <- ode(y=yinit.vector, method = "ode45", times=fit_times,
                   func=MSIRS_immunization_dynamics,
                   parms=c(parmset,baseline.txn.rate=baseline.txn.rate,b1=b1,phi=phi))

    t0 <- length(time_series)
    al <- nrow(yinit)
    output <- head(results,t0)
    St <- output[,-1]
    I1 <- St[,grep('I1', colnames(St))]
    I2 <- St[,grep('I2', colnames(St))]
    I3 <- St[,grep('I3', colnames(St))]
    I4 <- St[,grep('I4', colnames(St))]
    R1 <- St[,grep('R1', colnames(St))]
    R2 <- St[,grep('R2', colnames(St))]
    R3 <- St[,grep('R3', colnames(St))]
    R4 <- St[,grep('R4', colnames(St))]
    S1 <- St[,grep('S1', colnames(St))]
    S2 <- St[,grep('S2', colnames(St))]
    S3 <- St[,grep('S3', colnames(St))]
    S0 <- St[,grep('S0', colnames(St))]
    M0 <- St[,grep('M0', colnames(St))]
    Si<- St[,grep('Si', colnames(St))]
    Mn1<- St[,grep('Mn1', colnames(St))]
    Mn2<- St[,grep('Mn2', colnames(St))]
    Mv1<- St[,grep('Mv1', colnames(St))]
    Mv2<- St[,grep('Mv2', colnames(St))]
    N1<- St[,grep('N1', colnames(St))]
    N2<- St[,grep('N2', colnames(St))]
    Vs1<- St[,grep('Vs1', colnames(St))]
    Vs2<- St[,grep('Vs2', colnames(St))]

    beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2# q depends on transmission type (whether depends on population density or not); c2 is the contact matrix

    lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
    for (t in 1:t0){
      lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+parmset$seed)%*%beta)/sum(St[t,]))}

    hosp1 <- c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2, report_infants*0.15, report_infants*0.15, rep(report_children, 2), rep(.001, 6))
    hosp2 <- hosp1 * 0.4
    hosp3 <- c(rep(0.001, 8), rep(report_adults, 4), report_seniors65, report_seniors75)

    H1 <- matrix(0, nrow = t0, ncol = al)
    for (i in 1:al) {
      H1[, i] <-
        hosp1[i] * parmset$RRHm * parmset$sigma3 * M0[, i] * lambda1[, i] +
        hosp1[i] * parmset$RRHm * parmset$RRHn1 * parmset$sigma3 * Mn1[, i] * lambda1[, i] +
        hosp1[i] * parmset$RRHm * parmset$RRHn2 *  Mn2[, i] * lambda1[, i] +
        hosp1[i] * parmset$RRHm * parmset$RRHv1 * parmset$sigma3 * Mv1[, i] * lambda1[, i] +
        hosp1[i] * parmset$RRHm * parmset$RRHv2 * Mv1[, i] * lambda1[, i] +
        hosp1[i] * parmset$RRHn1 * N1[, i] * lambda1[, i] +
        hosp1[i] * parmset$RRHn2 * N2[, i] * lambda1[, i] +
        hosp1[i] * S0[, i] * lambda1[, i] +
        hosp1[i] * Si[, i] * lambda1[, i] +
        hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
        hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
        hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
        hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs1[, i] * lambda1[, i] +
        hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs2[, i] * lambda1[, i]
    }

    H <- rowSums(H1)

    H2 <- cbind(rowSums(H1[, 1:3]), rowSums(H1[, 4:6]), rowSums(H1[, 7:8]), rowSums(H1[, 9:12]), H1[, 13],H1[,14])
    age_dist2 <- colSums(H2)


    LLall <- sum(dpois(x = time_series, lambda = H, log = TRUE))
    LLmulti <- dmultinom(x = age_dist2, prob = age_dist, log = TRUE)

    LL <- LLall + LLmulti



   #  H1=matrix(0,nrow=t0,ncol=al)#Number of hospitalizations by age
   # for (i in 1:al){
   #   H1[,i]=
   #     parmset$RRHm*parmset$sigma3*M0[,i]*lambda1[,i]+
   #     parmset$RRIn*parmset$sigma3*Mn[,i]*lambda1[,i]+
   #     parmset$RRIv*parmset$sigma3*Mv[,i]*lambda1[,i]+
   ##     parmset$RRIn*N[,i]*lambda1[,i]+
    #    S0[,i]*lambda1[,i]+
    #    Si[,i]*lambda1[,i]+
    #    parmset$sigma1*S1[,i]*lambda1[,i]+
    #    parmset$sigma2*S2[,i]*lambda1[,i]+
   #     parmset$sigma3*S3[,i]*lambda1[,i]+
    #    parmset$RRIs*parmset$sigma3*Vs1[,i]*lambda1[,i]+
    #    parmset$RRIs*parmset$sigma3*Vs2[,i]*lambda1[,i]}

#
  #  H = rowSums(H1)*reporting_rate #combine into single time series
  #  H2 <- cbind(H[, 1:3], H1[,4:6],rowSums(H[, 7:8]), rowSums(H[, 9:12]), H1[, 13])
  #  age_dist2 <- colSums(H2[1:144,])

   # LL <- sum(dpois(x = dat, lambda =H, log = TRUE)) # fit to timeseries
   # LLmulti <- dmultinom(x = age_dist2, prob = age_dist, log = TRUE)

   # LL <- LLall + LLmulti

    return(LL)
  }


  # Run optimization function  --------------------------------------

  fitLL <- optim(par = c(-0.5,-2,2,-2,-4,-8,-5,-5),
                # lower = c(2,-2.4,1.5,-10,-10,-10,-10),
                # upper = c(3,-1.2,2.5,2,-2,-2,-2),
                 fn = fitmodel,        # the distance function to optimize
                 dat = time_series,  # the dataset to fit to (dpois function)
                control = list(fnscale=-1))# the log likelihood is negative; here we maximize the log likelihood)
  #save your parameters
  #fitLL = readRDS("C:/Users/hansencl/OneDrive - National Institutes of Health/Desktop/RSV Final/Model/2. Calibration/parameters_6Mar24.rds")
 # baseline.txn.rate=exp(fitLL$par[1])
 # b1=exp(fitLL$par[2])
  baseline.txn.rate = 5+(10*(exp(fitLL$par[1]))) / (1+exp(fitLL$par[1]))
  b1 <-  .08+(.35*(exp(fitLL$par[2]))) / (1+exp(fitLL$par[2]))
  phi=(2*pi*(exp(fitLL$par[3]))) / (1+exp(fitLL$par[3]))
  report_infants <- 1 / (1 + exp(-fitLL$par[4]))
  report_children <- 1 / (1 + exp(-fitLL$par[5]))
  report_adults <- 1 / (1 + exp(-fitLL$par[6]))
  report_seniors65 <- 1 / (1 + exp(-fitLL$par[7]))
  report_seniors75 <- 1 / (1 + exp(-fitLL$par[8]))


  set.seed(123)
  h=100
  lhs<-maximinLHS(h,3)

  beta_lower = baseline.txn.rate*.95
  beta_upper = baseline.txn.rate*1.05
  b1_lower = b1*.95
  b1_upper = b1*1.05
  phi_lower = phi*.985
  phi_upper = phi*1.015

  lhs_parms <- cbind(
    beta = lhs[,1]*(beta_upper-beta_lower)+beta_lower,
    b1 = lhs[,2]*(b1_upper-b1_lower)+b1_lower,
    phi = lhs[,3]*(phi_upper-phi_lower)+phi_lower)


  # Plot results with fit parameters   ---------------------------------------------

  output <- ode(y=yinit.vector, times=fit_times,method = "ode45",
                func=MSIRS_immunization_dynamics,
                parms=c(parmset,
                        baseline.txn.rate=baseline.txn.rate,
                        b1=b1,
                        phi=phi))

  t0=length(time_series)
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
  Mn1<- St[,grep('Mn1', colnames(St))]
  Mn2<- St[,grep('Mn2', colnames(St))]
  Mv1<- St[,grep('Mv1', colnames(St))]
  Mv2<- St[,grep('Mv2', colnames(St))]
  N1<- St[,grep('N1', colnames(St))]
  N2<- St[,grep('N2', colnames(St))]
  Vs1<- St[,grep('Vs1', colnames(St))]
  Vs2<- St[,grep('Vs2', colnames(St))]

  beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2


  lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
  for (t in 1:t0){
    lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+parmset$seed)%*%beta)/sum(St[t,]))}


  hosp1 <- c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2, report_infants*0.15, report_infants*0.15, rep(report_children, 2), rep(.001, 6))
  hosp2 <- hosp1 * 0.4
  hosp3 <- c(rep(0.001, 8), rep(report_adults, 4), report_seniors65, report_seniors75)

  H1 <- matrix(0, nrow = t0, ncol = al)
  for (i in 1:al) {
    H1[, i] <-
      hosp1[i] * parmset$RRHm * parmset$sigma3 * M0[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHn1 * parmset$sigma3 * Mn1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHn2 *  Mn2[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHv1 * parmset$sigma3 * Mv1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHv2 * Mv1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHn1 * N1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHn2 * N2[, i] * lambda1[, i] +
      hosp1[i] * S0[, i] * lambda1[, i] +
      hosp1[i] * Si[, i] * lambda1[, i] +
      hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 *  parmset$RRHs * Vs1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 *  parmset$RRHs * Vs2[, i] * lambda1[, i]
  }

  H2 <- rowSums(H1)
  prop2 <- cbind(rowSums(H1[, 1:3]), rowSums(H1[, 4:6]), rowSums(H1[, 7:8]), rowSums(H1[, 9:12]), H1[, 13],H1[, 14])
  age_dist2 <- colSums(prop2)/sum(prop2)


  plot1 = ggplot2::ggplot()+
    theme_bw()+
    geom_area(aes(x=1:length(time_series), y=time_series,fill="Data"))+
    geom_line(aes(x=1:length(time_series), y=H2,color="Model"),linewidth=1.5)+
    scale_fill_manual(name=NULL,values=c("steelblue3"))+
    scale_color_manual(name=NULL,values=c("maroon4"))+
    labs(x=NULL, y="RSV Hospitalizations")

  ages = c("<6m","6-11m","1-4yrs","5-64yrs","65-74yrs","75+yrs")
  ages_order = factor(ages, levels=ages,ordered =TRUE)
  plot2 = ggplot2::ggplot()+
    theme_bw()+
    geom_bar(aes(x=ages_order,y=age_dist,fill="Data"),stat="identity")+
    geom_point(aes(x=ages_order,y=age_dist2,color="Model"),size=5)+
    scale_fill_manual(name=NULL,values=c("steelblue3"))+
    scale_color_manual(name=NULL,values=c("maroon4"))+
    labs(x=NULL, y="Age Distribution")

  plot = plot_grid(plot1,plot2, nrow=2)
  print(plot)

 # fitted_parms = list(baseline.txn.rate,b1,phi,reporting_rate, hosp_props, lhs_parms)
  fitted_parms = list(baseline.txn.rate,b1,phi,report_infants, report_children, report_adults, report_seniors65,report_seniors75,lhs_parms)

  return(fitted_parms)

}
