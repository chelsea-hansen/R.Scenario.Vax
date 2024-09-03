#' MSIRS dynamic transmission model for RSV (function)
#'
#' This function fits the equations which dictate the flow of individuals between compartments.
#' This function is used in the fit_model() and scenario_projection() functions, but does not need to be used directly by the package user.
#'
#' @param times A numeric vector of times.
#' @param y A vector of starting values for each model compartment.
#' @param parms A list of fixed parameter values.
#'
#' @return A matrix of values.
#' @export
#' @import deSolve
#' @examples
#'

#' dat = get_data(state_or_county="state",state_abbr="CA",county_name=NULL)
#' parmset=dat[[1]] #fixed parameters
#' yinit.vector=dat[[3]]
#'
#' #time steps
#' fit_times = seq(1,100,by=1)
#'
#' parms=c(parmset,baseline.txn.rate=7,b1=.12,phi=3.2)#add the fitted parameters
#'
#' results <- MSIRS_immunization_dynamics(times=fit_times, y=yinit.vector,parms=parms)

MSIRS_immunization_dynamics <- function(times,y,parms){

  States<-array(y, dim=dim(parms$yinit.matrix))
  dimnames(States) <- dimnames(parms$yinit.matrix)

  length.step <- switch(parms$time.step,
                        'month' = 30.44,
                        'week' = 7,
                        'daily' = 1)

  period <- switch(parms$time.step,
                   'month' = 12,
                   'week' = 52.1775,
                   'daily' = 365.25)


  um=parms$um
  omega = 1/(parms$DurationMatImmunityDays/length.step)
  waning1 = 1/(parms$recover1/length.step)
  waning2 = 1/(parms$recover2/length.step)
  waning3 = 1/(parms$recover3/length.step)
  waning4 = 1/(parms$recover3/length.step)

  #parameters for monoclonals
  mono_01 = parms$monoclonal_01[times]
  mono_23 = parms$monoclonal_23[times]
  mono_45 = parms$monoclonal_45[times]
  mono_67 = parms$monoclonal_67[times]
  waningN=1/(parms$waningN/length.step) #duration of nirsevimab protection
  RRIn=parms$RRIn #relative risk of infection for infants receiving nirsevimab (default is to set to 1)

  #parameters for maternal vaccines
  birth_V = parms$maternal_vax[times]#infants born to vaccinated mothers
  waningV=1/(parms$waningV/length.step)#duration of maternal vaccine
  RRIv=parms$RRIv

  #parameters for seniors
  V_s65=parms$senior_vax_65_74[times]
  V_s75=parms$senior_vax_75[times]
  waningS=1/(parms$waningS/2/length.step) #duration of protection from vaccination (divided by 2 because 2 compartments)
  RRIs=parms$RRIs #relative risk of infection for vaccinated seniors

  mu <- 1 / (parms$WidthAgeClassMonth * ifelse(parms$time.step == 'week', 4.345, ifelse(parms$time.step == 'daily', 30.44, 1)))



  gamma1= 1/(parms$dur.days1/length.step)  #converts 1/days to 1/length.step
  gamma2= 1/(parms$dur.days2/length.step)
  gamma3= 1/(parms$dur.days3/length.step)
  gamma4= gamma3



  #Pull out the states  for the model as vectors
  M0 <-  States[,'M0']
  Mv <-  States[,'Mv']
  Mn <-  States[,'Mn'] # newborns who receive monoclonals at birth
  N <-  States[,'N'] #infants who receive nirsevimab ahead of the season but not at birth

  Si <-  States[,'Si']
  S0 <-  States[,'S0']
  I1 <-  States[,'I1']
  R1 <-  States[,'R1']


  S1 <-  States[,'S1']
  I2 <-  States[,'I2']
  R2 <-  States[,'R2']

  S2 <-  States[,'S2']
  I3 <-  States[,'I3']
  R3 <-  States[,'R3']

  S3 <-  States[,'S3']
  I4 <-  States[,'I4']
  R4 <-  States[,'R4']


  Vs1 <-  States[,'Vs1'] # vaccination of seniors
  Vs2 <-  States[,'Vs2']

  N.ages <- length(M0)

  ####################
  seasonal.txn <- (1+parms$b1*cos(2*pi*(times-parms$phi*period)/period))# seasonality waves

  transmission_unittime <-  parms$baseline.txn.rate/(parms$dur.days1/length.step)

  beta=transmission_unittime*parms$c2

  beta_a_i <- seasonal.txn * beta/(sum(States)^parms$q)

  infectiousN <- I1 + parms$rho1*I2 + parms$rho2*I3 + parms$rho2*I4 + parms$seed

  lambda <- infectiousN %*% beta_a_i

  lambda <- as.vector(lambda)
  ##########transmission dynamics##########################

  dy <- matrix(NA, nrow=N.ages, ncol=ncol(States))
  colnames(dy) <- colnames(States)

  period.birth.rate <-parms$PerCapitaBirthsYear #weekly number of births
  Aging.Prop <- c(0,mu[1:(N.ages-1)])


  dy[,'M0'] <- period.birth.rate - c(birth_V,rep(0,13)) -
    c(mono_01,mono_23,rep(0,12))-
    parms$sigma3*lambda*M0 -
    omega*M0-
    (mu+um)*M0 +
    Aging.Prop*c(0,M0[1:(N.ages-1)])

  # newborns who receive monoclonals
  dy[,'Mn'] <- c(mono_01,mono_23,rep(0,12))-
    RRIn*parms$sigma3*lambda*Mn -
    waningN*Mn-
    (mu+um)*Mn +
    Aging.Prop*c(0,Mn[1:(N.ages-1)]) #aging in

  dy[,'Mv'] <- c(birth_V,rep(0,13)) -
    RRIn*parms$sigma3*lambda*Mv -
    waningV*Mv-
    (mu+um)*Mv +
    Aging.Prop*c(0,Mv[1:(N.ages-1)]) #aging in

  # infants <8 months who did not receive nirsevimab at birth but receive a catch-up dose ahead of the season
  # these infants can be in the M0 or S0 compartments
  dy[,'N'] <- c(0,0,mono_45,mono_67,rep(0,10)) - #from both M0 and S0 compartments
    RRIn*lambda*N -
    waningN*N-
    (mu + um)*N +
    Aging.Prop*c(0,N[1:(N.ages-1)])


  dy[,'Si'] <-  waningN*N + waningN*Mn + waningV*Mv-
    lambda*Si -
    (mu + um)*Si +
    Aging.Prop*c(0,Si[1:(N.ages-1)])

  dy[,'S0'] <- omega*M0 -
    c(0,0,mono_45,mono_67,rep(0,10)) -
    lambda*S0 -
    (mu + um)*S0 +
    Aging.Prop*c(0,S0[1:(N.ages-1)])

  dy[,'I1'] <-  lambda*S0+ lambda*Si+ RRIn*lambda*N +
    parms$sigma3*lambda*M0 +  RRIn*parms$sigma3*lambda*Mn + RRIv*parms$sigma3*lambda*Mv-
    gamma1*I1-
    (mu + um)*I1 +
    Aging.Prop*c(0,I1[1:(N.ages-1)])

  dy[,'R1'] <-  gamma1*I1 -
    waning1*R1 -
    (mu + um)*R1 +
    Aging.Prop*c(0,R1[1:(N.ages-1)])

  dy[,'S1'] <-waning1*R1 -
    parms$sigma1*lambda*S1 -
    (mu+um)*S1 +
    Aging.Prop*c(0,S1[1:(N.ages-1)])

  dy[,'I2'] <- parms$sigma1*lambda*S1 -
    gamma2*I2 -
    (mu + um)*I2 +
    Aging.Prop*c(0,I2[1:(N.ages-1)])

  dy[,'R2'] <-  gamma2*I2 -
    waning2*R2 -
    (mu + um)*R2 +
    Aging.Prop*c(0,R2[1:(N.ages-1)])

  dy[,'S2'] <- waning2*R2  -
    parms$sigma2*lambda*S2 -
    (mu+um)*S2 +
    Aging.Prop*c(0,S2[1:(N.ages-1)])

  dy[,'I3'] <- parms$sigma2*lambda*S2 -
    gamma3*I3 -
    (mu+um)*I3 +
    Aging.Prop*c(0,I3[1:(N.ages-1)])

  dy[,'R3'] <-  gamma3*I3 -
    waning3*R3-
    (mu + um)*R3 +
    Aging.Prop*c(0,R3[1:(N.ages-1)])

  dy[,'S3'] <- waning3*R3 + waning4*R4  + waningS*Vs2-
    c(rep(0,12),V_s65,V_s75)-
    parms$sigma3*lambda*S3 -
    (mu + um)*S3 +
    Aging.Prop*c(0,S3[1:(N.ages-1)])


  #make a vaccination compartment for the >65
  dy[,'Vs1'] <- c(rep(0,12),V_s65,V_s75) -
    RRIs*parms$sigma3*lambda*Vs1 -
    waningS*Vs1  -
    (mu + um)*Vs1 +
    Aging.Prop*c(0,Vs1[1:(N.ages-1)])

  #add another vaccination compartment
  dy[,'Vs2'] <- waningS*Vs1 -
    RRIs*parms$sigma3*lambda*Vs2 -
    waningS*Vs2 -
    (mu + um)*Vs2 +
    Aging.Prop*c(0,Vs2[1:(N.ages-1)])

  dy[,'I4'] <- parms$sigma3*lambda*S3 +
    RRIs*parms$sigma3*lambda*Vs1+
    RRIs*parms$sigma3*lambda*Vs2-
    gamma4*I4 -
    (mu + um)*I4 +
    Aging.Prop*c(0,I4[1:(N.ages-1)])


  dy[,'R4'] <-  gamma4*I4 -
    waning4*R4 -
    (mu + um)*R4 +
    Aging.Prop*c(0,R4[1:(N.ages-1)])

  derivs <- as.vector(dy)

  res <- list(derivs)

  return(res)
}



