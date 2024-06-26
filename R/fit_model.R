#' Fit model parameters
#'
#' @param time_series An RSV hospitalization/ED visit time series
#' @param age_dist The proportion of RSV hospitalizations in each age group (age groups include: <6 months,6-11 months,1-4 years,5-59 years, 60+ years)
#' @param parmset A list of fixed parameters (retrieved in the get_data function)
#' @param yinit A matrix of compartment initial values (retrieved from get_data function)
#' @param yinit.vector A vector of compartment initial values (retrieved from get_data function)
#'
#' @return A list of fitted parameter values. Prints a plot showing the fitted model.
#' @export
#' @import deSolve
#' @import magrittr
#' @import dplyr
#' @import ggplot2
#' @import lhs
#' @importFrom utils tail head
#' @importFrom stats dpois optim
#
#' @examples
#' dat = get_data(state_or_county="state",state_abbr="CA",county_name=NULL)
#' parmset=dat[[1]]
#' yinit=dat[[2]]
#' yinit.vector=dat[[3]]
#' fitLL = fit_model(time_series = weekly_rsv, age_dist = c(.19,.11,.15,.08,.17,.30),
#' parmset=parmset,yinit=yinit,yinit.vector=yinit.vector)
fit_model = function(time_series, age_dist, parmset, yinit,yinit.vector){

  fit_times = seq(1,length(time_series)+52,by=1)
  parmset$monoclonal_birth = rep(0,length(time_series)+52)
  parmset$monoclonal_catchup = rep(0,length(time_series)+52)
  parmset$maternal_vax = rep(0,length(time_series)+52)
  parmset$senior_vax = rep(0,length(time_series)+52)

  fitmodel <-  function(parameters,dat) {
    protrans <- parameters[1] # parameter for baseline transmission rate
    baseline.txn.rate = 6+(3*(exp(protrans))) / (1+exp(protrans)) #transform between 6 and 9
    amp <- parameters[2] # parameter for seasonal amplitude
    b1 <-   .05+(.5*(exp(amp))) / (1+exp(amp)) #transform between.05 and .5
    trans <- parameters[3] # parameter for seasonal phase
    phi <-  (2*pi*(exp(trans))) / (1+exp(trans)) # transform to between 0 and 2pi
    reporting_rate = 1/(1+exp(-parameters[4])) #reporting rate = how many infections are reported as hospitalizations in the data

    # Simulate the model with initial conditions and timesteps defined above, and parameter values from function call
    results <- ode(y=yinit.vector, method = "ode45", t=fit_times,
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
    Mn<- St[,grep('Mn', colnames(St))]
    Mv<- St[,grep('Mv', colnames(St))]
    N<- St[,grep('N', colnames(St))]
    Vs1<- St[,grep('Vs1', colnames(St))]
    Vs2<- St[,grep('Vs2', colnames(St))]

    beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2# q depends on transmission type (whether depends on population density or not); c2 is the contact matrix

    lambda1=matrix(0,nrow=t0,ncol=13)#Force of infection
    for (t in 1:t0){
      lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+parmset$seed)%*%beta)/sum(St[t,]))}

    H1=matrix(0,nrow=t0,ncol=al)#Number of hospitalizations by age
    for (i in 1:al){
      H1[,i]=
        parmset$RRHm*parmset$sigma3*M0[,i]*lambda1[,i]+
        parmset$RRIn*parmset$sigma3*Mn[,i]*lambda1[,i]+
        parmset$RRIv*parmset$sigma3*Mv[,i]*lambda1[,i]+
        parmset$RRIn*N[,i]*lambda1[,i]+
        S0[,i]*lambda1[,i]+
        Si[,i]*lambda1[,i]+
        parmset$sigma1*S1[,i]*lambda1[,i]+
        parmset$sigma2*S2[,i]*lambda1[,i]+
        parmset$sigma3*S3[,i]*lambda1[,i]+
        parmset$RRIs*parmset$sigma3*Vs1[,i]*lambda1[,i]+
        parmset$RRIs*parmset$sigma3*Vs2[,i]*lambda1[,i]}


    H = rowSums(H1)*reporting_rate #combine into single time series

    LL <- sum(dpois(x = dat, lambda =H, log = TRUE)) # fit to timeseries

    return(LL)
  }


  # Run optimization function  --------------------------------------

  fitLL <- optim(par = c(0,-2,0,0),
                 fn = fitmodel,        # the distance function to optimize
                 dat = time_series,  # the dataset to fit to (dpois function)
                 control = list(fnscale=-1))# the log likelihood is negative; here we maximize the log likelihood)
  #save your parameters
  #fitLL = readRDS("C:/Users/hansencl/OneDrive - National Institutes of Health/Desktop/RSV Final/Model/2. Calibration/parameters_6Mar24.rds")
  baseline.txn.rate=6+(3*(exp(fitLL$par[1]))) / (1+exp(fitLL$par[1]))
  b1=.05+(.5*(exp(fitLL$par[2]))) / (1+exp(fitLL$par[2]))
  phi=(2*pi*(exp(fitLL$par[3]))) / (1+exp(fitLL$par[3]))
  reporting_rate = 1/(1+exp(-fitLL$par[4]))


  set.seed(123)
  h=100
  lhs<-maximinLHS(h,3)

  #  beta=fitLL$par[1]
  #  beta_lower = fitLL$par[1] - abs(fitLL$par[1]*.5)
  #  beta_upper = fitLL$par[1] + abs(fitLL$par[1]*.5)
  #  b1_lower = fitLL$par[2] - abs(fitLL$par[2]*.5)
  #  b1_upper = fitLL$par[2] + abs(fitLL$par[2]*.5)
  # phi_lower = fitLL$par[3] - abs(fitLL$par[3]*.5)
  #  phi_upper = fitLL$par[3] + abs(fitLL$par[3]*.5)

  beta_lower = baseline.txn.rate*.985
  beta_upper = baseline.txn.rate*1.015
  b1_lower = b1*.985
  b1_upper = b1*1.015
  phi_lower = phi*.985
  phi_upper = phi*1.015

  lhs_parms <- cbind(
    beta = lhs[,1]*(beta_upper-beta_lower)+beta_lower,
    b1 = lhs[,2]*(b1_upper-b1_lower)+b1_lower,
    phi = lhs[,3]*(phi_upper-phi_lower)+phi_lower)


  #transform amplitude and phase
  # lhs_parms[,"beta"]= 6+(3*(exp(lhs_parms[,"beta"]))) / (1+exp(lhs_parms[,"beta"]))
  # lhs_parms[,"b1"] = .05+(.5*(exp(lhs_parms[,"b1"]))) / (1+exp(lhs_parms[,"b1"]))
  # lhs_parms[,"phi"] = (2*pi*(exp(lhs_parms[,"phi"]))) / (1+exp(lhs_parms[,"phi"]))



  # Plot results with fit parameters   ---------------------------------------------

  output <- ode(y=yinit.vector, t=fit_times,method = "ode45",
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
  Mn<- St[,grep('Mn', colnames(St))]
  Mv<- St[,grep('Mv', colnames(St))]
  N<- St[,grep('N', colnames(St))]
  Vs1<- St[,grep('Vs1', colnames(St))]
  Vs2<- St[,grep('Vs2', colnames(St))]

  beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2


  lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
  for (t in 1:t0){
    lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+parmset$seed)%*%beta)/sum(St[t,]))}

  H1=matrix(0,nrow=t0,ncol=al)#Number of hospitalizations by age
  for (i in 1:al){
    H1[,i]=
      parmset$RRHm*parmset$sigma3*M0[,i]*lambda1[,i]+
      parmset$RRIn*parmset$sigma3*Mn[,i]*lambda1[,i]+
      parmset$RRIv*parmset$sigma3*Mv[,i]*lambda1[,i]+
      parmset$RRIn*N[,i]*lambda1[,i]+
      S0[,i]*lambda1[,i]+
      Si[,i]*lambda1[,i]+
      parmset$sigma1*S1[,i]*lambda1[,i]+
      parmset$sigma2*S2[,i]*lambda1[,i]+
      parmset$sigma3*S3[,i]*lambda1[,i]+
      parmset$RRIs*parmset$sigma3*Vs1[,i]*lambda1[,i]+
      parmset$RRIs*parmset$sigma3*Vs2[,i]*lambda1[,i]}

  H2 = rowSums(H1)*reporting_rate
  H2.1= H2*age_dist[1]
  H2.2= H2*age_dist[2]
  H2.3= H2*age_dist[3]
  H2.4= H2*age_dist[4]
  H2.5= H2*age_dist[5]

  H3 = cbind(rowSums(H1[,1:3]),
             rowSums(H1[,4:6]),
             rowSums(H1[,7:8]),
             rowSums(H1[,9:12]),
             H1[,13])

  H3.1 = sum(H2.1)/sum(H3[,1])
  H3.2 = sum(H2.2)/sum(H3[,2])
  H3.3 = sum(H2.3)/sum(H3[,3])
  H3.4 = sum(H2.4)/sum(H3[,4])
  H3.5 = sum(H2.5)/sum(H3[,5])

  hosp_props = c(H3.1,H3.2,H3.3,H3.4,H3.5)


  plot1 = ggplot()+
    theme_bw()+
    geom_area(aes(x=1:201, y=time_series),fill="seashell3")+
    geom_line(aes(x=1:201, y=H2),color="navy",linewidth=1.5)+
    labs(x=NULL, y="RSV Hospitalizations")
  print(plot1)

  fitted_parms = list(baseline.txn.rate,b1,phi,reporting_rate, hosp_props, lhs_parms)

  return(fitted_parms)

}
