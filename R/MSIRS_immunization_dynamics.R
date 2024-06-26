#' MSIRS dynamic transmission model for RSV
#'
#' @param times A numeric vector of times
#' @param y A vector of starting values for model compartments
#' @param parms A list of fixed parameter values
#'
#' @return A matrix of values.
#' @export
#'
#' @examples
#'

#' dat = get_data(state_or_county="state",state_abbr="CA",county_name=NULL)
#' parmset=dat[[1]]
#' yinit.vector=dat[[3]]
#' fit_times = seq(1,100,by=1)
#' parms=c(parmset,baseline.txn.rate=7,b1=.12,phi=3.2)
#' results <- MSIRS_immunization_dynamics(times=fit_times, y=yinit.vector,parms=parms)


MSIRS_immunization_dynamics <- function(times,y,parms){

  States<-array(y, dim=dim(parms$yinit.matrix))
  dimnames(States) <- dimnames(parms$yinit.matrix)
  um=parms$um

  if(parms$time.step=='month'){
    period=12
    length.step=30.44 #days
  }else if(parms$time.step=='week'){
    period=52.1775
    length.step=7 #days
  }

  omega = 1/(parms$DurationMatImmunityDays/length.step)
  waning1 = 1/(parms$recover1/length.step)
  waning2 = 1/(parms$recover2/length.step)
  waning3 = 1/(parms$recover3/length.step)
  waning4 = 1/(parms$recover3/length.step)



  #parameters for monoclonals
  birth_N = parms$monoclonal_birth[times]#infants who receive birth dose of nirsevimab
  cup_N = parms$monoclonal_catchup[times] #infants who receive nirsevimab catch-up dose
  waningN=1/(parms$waningN/length.step) #duration of nirsevimab protection
  RRIn=parms$RRIn #relative risk of infection for infants receiving nirsevimab (default is to set to 1)

  #parameters for maternal vaccines
  birth_V = parms$maternal_vax[times]#infants born to vaccinated mothers
  waningV=1/(parms$waningV/length.step)#duration of maternal vaccine
  RRIv=parms$RRIv


  #parameters for seniors
  V_s=parms$senior_vax[times] # indicator function for when the vaccine will be administered * vaccine coverage
  waningS=1/(parms$waningS/2/length.step) #duration of protection from vaccination (divided by 2 because 2 compartments)
  RRIs=parms$RRIs #relative risk of infection for vaccinated seniors


  mu= 1/parms$WidthAgeClassMonth
  if(parms$time.step=='week'){
    mu= 1/(parms$WidthAgeClassMonth*4.345)
  }

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

  # baseline.txn.rate is the probability of transmission given contact per capita
  # (parms$dur.days1/length.step) is the duration of infectiousness of primary infection
  # q depends on transmission type (whether depends on population density or not)
  # density (q=0) vs frequency-dependent (q=1) transmission
  # c2 is the contact matrix
  # beta is transmisibility per unit time
  transmission_unittime <-  parms$baseline.txn.rate/(parms$dur.days1/length.step)
  beta=transmission_unittime*parms$c2

  beta_a_i <- seasonal.txn * beta/(sum(States)^parms$q)

  infectiousN <- I1 + parms$rho1*I2 + parms$rho2*I3 + parms$rho2*I4 + parms$seed
  lambda <- infectiousN %*% beta_a_i
  lambda <- as.vector(lambda)
  ##########transmission dynamics##########################

  dy <- matrix(NA, nrow=N.ages, ncol=ncol(States))
  colnames(dy) <- colnames(States)

  period.birth.rate <- log(parms$PerCapitaBirthsYear+1)/period #adjust annual birthrates to weekly scale
  #birth_N = log(birth_N+1)
  #birth_V = log(birth_V+1)
  #cup_N = log(cup_N+1)
  #mu represents aging to the next class
  #um is death rate/emigration; adjust it to reproduce population growth

  Aging.Prop <- c(0,mu[1:(N.ages-1)])
  cup_prop = cup_N/sum(States[1:4,1:2])

  dy[,'M0'] <- period.birth.rate*sum(States) - c(birth_N,rep(0,12)) - c(birth_V,rep(0,12)) -
    #(1-(birth_N+birth_V))*period.birth.rate*sum(States) -
    parms$sigma3*lambda*M0 -
    c(cup_prop*sum(States[1,"M0"]),cup_prop*sum(States[2,"M0"]),cup_prop*sum(States[3,"M0"]),cup_prop*sum(States[4,"M0"]),rep(0,9))-
    (omega+(mu+um))*M0 +
    Aging.Prop*c(0,M0[1:(N.ages-1)])

  # newborns who receive monoclonals
  dy[,'Mn'] <- c(birth_N,rep(0,12)) + c(cup_prop*sum(States[1,"M0"]),cup_prop*sum(States[2,"M0"]),cup_prop*sum(States[3,"M0"]),cup_prop*sum(States[4,"M0"]),rep(0,9)) -
    #birth_N*period.birth.rate*sum(States) + cup_N*sum(States[1:4,"M0"]) -
    RRIn*parms$sigma3*lambda*Mn -
    waningN*Mn-
    (mu+um)*Mn +
    Aging.Prop*c(0,Mn[1:(N.ages-1)]) #aging in

  dy[,'Mv'] <- c(birth_V,rep(0,12)) -
    #birth_V*period.birth.rate*sum(States) -
    RRIn*parms$sigma3*lambda*Mv -
    waningV*Mv-
    (mu+um)*Mv +
    Aging.Prop*c(0,Mv[1:(N.ages-1)]) #aging in

  # infants <8 months who did not receive nirsevimab at birth but receive a catch-up dose ahead of the season
  # these infants can be in the M0 or S0 compartments
  dy[,'N'] <- c(cup_prop*sum(States[1,"S0"]),cup_prop*sum(States[2,"S0"]),cup_prop*sum(States[3,"S0"]),cup_prop*sum(States[4,"S0"]),rep(0,9)) -
    #cup_N*sum(States[1:4,"S0"])  -
    waningN*N-
    RRIn*lambda*N -
    (mu + um)*N +
    Aging.Prop*c(0,N[1:(N.ages-1)])


  dy[,'Si'] <-  waningN*N + waningN*Mn + waningV*Mv-
    lambda*Si -
    (mu + um)*Si +
    Aging.Prop*c(0,Si[1:(N.ages-1)])

  dy[,'S0'] <- omega*M0 -
    c(cup_prop*sum(States[1,"S0"]),cup_prop*sum(States[2,"S0"]),cup_prop*sum(States[3,"S0"]),cup_prop*sum(States[4,"S0"]),rep(0,9))-
    #cup_N -
    #cup_N*sum(States[1:4,"S0"]) -
    lambda*S0 -
    (mu + um)*S0 +
    Aging.Prop*c(0,S0[1:(N.ages-1)])

  dy[,'I1'] <-  lambda*S0+ lambda*Si+ RRIn*lambda*N +
    parms$sigma3*lambda*M0 +  RRIn*parms$sigma3*lambda*Mn + RRIv*parms$sigma3*lambda*Mv-
    (gamma1 + mu + um)*I1 +
    Aging.Prop*c(0,I1[1:(N.ages-1)])

  dy[,'R1'] <-  gamma1*I1 -
    (waning1 + mu + um)*R1 +
    Aging.Prop*c(0,R1[1:(N.ages-1)])

  dy[,'S1'] <-waning1*R1 -
    parms$sigma1*lambda*S1 -
    (mu+um)*S1 +
    Aging.Prop*c(0,S1[1:(N.ages-1)])

  dy[,'I2'] <- parms$sigma1*lambda*S1 -
    (gamma2 + mu + um)*I2 +
    Aging.Prop*c(0,I2[1:(N.ages-1)])

  dy[,'R2'] <-  gamma2*I2 -
    (waning2 + mu + um)*R2 +
    Aging.Prop*c(0,R2[1:(N.ages-1)])

  dy[,'S2'] <- waning2*R2  -
    parms$sigma2*lambda*S2 -
    (mu+um)*S2 +
    Aging.Prop*c(0,S2[1:(N.ages-1)])

  dy[,'I3'] <- parms$sigma2*lambda*S2 -
    (gamma3 + mu+um)*I3 +
    Aging.Prop*c(0,I3[1:(N.ages-1)])

  dy[,'R3'] <-  gamma3*I3 -
    (waning3 + mu + um)*R3 +
    Aging.Prop*c(0,R3[1:(N.ages-1)])

  dy[,'S3'] <- waning3*R3 + waning4*R4  + waningS*Vs2-
    parms$sigma3*lambda*S3 -
    c(rep(0,12),V_s*sum(States[13,"S3"])) - #seniors who get vaccinated
    (mu + um)*S3 +
    Aging.Prop*c(0,S3[1:(N.ages-1)])


  #make a vaccination compartment for the >65
  dy[,'Vs1'] <- c(rep(0,12),V_s*sum(States[13,"S3"])) +
    c(rep(0,12),V_s*sum(States[13,"R4"])) -
    RRIs*parms$sigma3*lambda*Vs1 -
    (mu + um)*Vs1 -
    waningS*Vs1 +
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
    c(rep(0,12),V_s*sum(States[13,"R4"]))-
    (waning4 + mu + um)*R4 +
    Aging.Prop*c(0,R4[1:(N.ages-1)])

  derivs <- as.vector(dy)

  res <- list(derivs)

  return(res)
}


