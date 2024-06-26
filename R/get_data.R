#' Download and format data for RSV dynamic transmission model
#'
#' @param state_or_county describe if you are modeling for a "state" or a "county"
#' @param state_abbr The two letter state abbreviation (you need this even if you are modeling a county)
#' @param county_name The name of the county (if modeling a state list NULL here)
#'
#' @return a list of values
#' @export
#' @import tidycensus
#' @import magrittr
#' @import dplyr
#' @importFrom utils tail
#'
#'
#' @examples
#' data_california = get_data(state_or_county="state",state_abbr="CA",county_name=NULL)
#' data_kingcounty = get_data(state_or_county="county",state_abbr="WA",county_name="King")
get_data = function(state_or_county, state_abbr, county_name){

  contact = matrix(c(0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.3840000, 1.1520000, 0.4066667, 0.3710833, 1.5046667, 0.6354167, 0.3431250,
                     0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.3840000, 1.1520000, 0.4066667, 0.3710833, 1.5046667, 0.6354167, 0.3431250,
                     0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.3840000, 1.1520000, 0.4066667, 0.3710833, 1.5046667, 0.6354167, 0.3431250,
                     0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.3840000, 1.1520000, 0.4066667, 0.3710833, 1.5046667, 0.6354167, 0.3431250,
                     0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.3840000, 1.1520000, 0.4066667, 0.3710833, 1.5046667, 0.6354167, 0.3431250,
                     0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.0640000, 0.3840000, 1.1520000, 0.4066667, 0.3710833, 1.5046667, 0.6354167, 0.3431250,
                     0.3840000, 0.3840000, 0.3840000, 0.3840000, 0.3840000, 0.3840000, 0.3840000, 1.1520000, 0.4066667, 0.3710833, 1.5046667, 0.6354167, 0.3431250,
                     1.1520000, 1.1520000, 1.1520000, 1.1520000, 1.1520000, 1.1520000, 1.1520000, 1.1520000, 0.4066667, 0.3710833, 1.5046667, 0.6354167, 0.3431250,
                     0.4066667, 0.4066667, 0.4066667, 0.4066667, 0.4066667, 0.4066667, 0.4066667, 0.4066667, 6.6400000, 1.7350000, 3.3550000, 1.8550000, 1.0200000,
                     0.3710833, 0.3710833, 0.3710833, 0.3710833, 0.3710833, 0.3710833, 0.3710833, 0.3710833, 1.7350000, 8.0550000, 2.7275000, 2.5675000, 1.1625000,
                     1.5046667, 1.5046667, 1.5046667, 1.5046667, 1.5046667, 1.5046667, 1.5046667, 1.5046667, 3.3550000, 2.7275000, 4.9075000, 3.1775000, 1.5050000,
                     0.6354167, 0.6354167, 0.6354167, 0.6354167, 0.6354167, 0.6354167, 0.6354167, 0.6354167, 1.8550000, 2.5675000, 3.1775000, 3.8225000, 2.3000000,
                     0.3431250, 0.3431250, 0.3431250, 0.3431250, 0.3431250, 0.3431250, 0.3431250, 0.3431250, 1.0200000, 1.1625000, 1.5050000, 2.3000000, 3.3000000),
                   nrow=13, ncol=13)

  burn_in=matrix(c(6.463157e-01, 4.148358e-01, 2.631371e-01, 1.647599e-01, 1.022259e-01, 6.328545e-02, 1.376335e-02, 1.183210e-03, 5.932823e-05, 1.393462e-06, 1.160965e-08, 0.000000e+00, 0.000000e+00,
                   3.384971e-01, 5.240746e-01, 5.996743e-01, 6.080935e-01, 5.843443e-01, 5.494106e-01, 3.540142e-01, 1.142966e-01, 1.217465e-02, 4.856697e-04, 8.324120e-06, 1.458630e-07, 0.000000e+00,
                   1.445836e-04, 1.877570e-04, 1.966434e-04, 1.906305e-04, 1.771783e-04, 1.637321e-04, 1.064949e-04, 3.803028e-05, 1.026271e-05, 4.805043e-07, 1.160965e-08, 0.000000e+00, 0.000000e+00,
                   6.208153e-03, 2.351284e-02, 5.062077e-02, 8.074394e-02, 1.069442e-01, 1.255334e-01, 1.164279e-01, 4.868788e-02, 1.236930e-02, 1.108115e-03, 3.823058e-05, 1.122023e-06, 3.253169e-08,
                   8.548321e-03, 3.598081e-02, 8.251572e-02, 1.381997e-01, 1.918796e-01, 2.376542e-01, 3.455607e-01, 2.609483e-01, 5.524949e-02, 4.534159e-03, 1.537466e-04, 5.217406e-06, 2.277218e-07,
                   1.490553e-06, 5.960541e-06, 1.191778e-05, 2.085021e-05, 2.828896e-05, 3.423490e-05, 5.275214e-05, 4.492995e-05, 2.463050e-05, 2.378496e-06, 9.287721e-08, 0.000000e+00, 0.000000e+00,
                   1.266970e-04, 6.094653e-04, 1.617839e-03, 3.228804e-03, 5.546126e-03, 8.770089e-03, 3.851550e-02, 5.161391e-02, 2.955543e-02, 4.358630e-03, 2.292558e-04, 9.727939e-06, 3.686924e-07,
                   1.550175e-04, 7.778506e-04, 2.173505e-03, 4.613854e-03, 8.495623e-03, 1.438015e-02, 1.009849e-01, 2.282458e-01, 1.018308e-01, 1.509968e-02, 8.266536e-04, 4.244613e-05, 2.656754e-06,
                   0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.488474e-06, 8.668192e-06, 2.209535e-05, 2.560790e-05, 4.444665e-06, 2.902413e-07, 1.122023e-08, 0.000000e+00,
                   1.490553e-06, 1.043095e-05, 3.426362e-05, 9.531526e-05, 2.233339e-04, 4.703578e-04, 1.316401e-02, 6.317381e-02, 8.497254e-02, 2.320240e-02, 2.141899e-03, 1.525727e-04, 9.553472e-06,
                   1.490553e-06, 4.470406e-06, 1.787667e-05, 5.361483e-05, 1.340004e-04, 2.917409e-04, 1.596334e-02, 1.959572e-01, 4.589380e-01, 5.539364e-01, 5.443667e-01, 5.667089e-01, 6.380972e-01,
                   0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 9.906505e-07, 1.264938e-05, 7.711692e-05, 1.083297e-04, 1.246644e-04, 1.118769e-04, 8.515711e-05,
                   0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.488893e-06, 4.465422e-06, 1.437186e-03, 3.577557e-02, 2.447129e-01, 3.971579e-01, 4.521101e-01, 4.329679e-01, 3.618048e-01,
                   0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                   0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                   0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                   0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                   0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                   0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00),
                 nrow=13,ncol=19)

  colnames(burn_in)= c("M0","S0","I1","R1","S1","I2","R2","S2","I3","R3","S3","I4","R4","Mn","Mv","N","Si","Vs1","Vs2")

  pop = tidycensus::get_estimates(
    geography = state_or_county,
    product = "characteristics",
    breakdown = "AGEGROUP",
    breakdown_labels = TRUE,
    state = state_abbr,
    county = county_name)%>%
    select("AGEGROUP","value")

  pop = tail(pop,18)

  new_agegrps = c("0to4","5to9",rep("10to19",2),rep("20-39",4),rep("40to59",4),rep("60+",6))

  pop = pop %>%
    dplyr::mutate(new_age = factor(new_agegrps,levels=unique(new_agegrps))) %>%
    dplyr::group_by(.data$new_age) %>%
    dplyr::summarize(population = sum(.data$value))

  adults=c(pop$population[2:6])
  infants = rep(pop$population[1]*.034,6)
  children = c(pop$population[1]*.2, rep(pop$population[2]*.6))
  all_ages = c(infants,children, adults)
  yinit = burn_in*all_ages


  #format the data for feeding into the model
  N_ages <- nrow(yinit)
  agenames <- c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y","5-9Y","10-19Y","20-39Y","40-59Y","60Y+")
  al <- N_ages
  rownames(yinit) <- agenames
  yinit.vector <- as.vector(unlist(yinit))
  name.array <- array(NA, dim=dim(yinit))
  for(i in 1:dim(name.array)[1]){
    for(j in 1:dim(name.array)[2]){
      name.array[i,j] <- paste(dimnames(yinit)[[1]][i],dimnames(yinit)[[2]][j]  )
    }
  }
  name.vector <- as.vector(name.array)
  names(yinit.vector) <- name.vector

  yinit.vector=yinit.vector


  #additional parameters
  seed = sum(pop$population)/100000 #seeding 1 infection per 100K population

  vital_stats = tidycensus::get_estimates(
    geography = state_or_county,
    product = "components",
    breakdown_labels = TRUE,
    state = state_abbr,
    county = county_name)

  birth_rate = vital_stats %>%  dplyr::filter(.data$variable=="RBIRTH") %>%  dplyr::select("value")/1000
  birth=c(birth_rate$value, rep(0,12))

  mig_rate = vital_stats %>%  dplyr::filter(.data$variable=="RNETMIG") %>%  dplyr::select("value")/1000
  um = mig_rate$value*-1/52.1775


  parmset<-list(PerCapitaBirthsYear=birth,
                WidthAgeClassMonth=c(rep(2,times=6), 12,12*3, 60, 120, 240, 240, 240),#time spend in each age class (months)
                DurationMatImmunityDays = 112,#duration of maternal immunity (days)
                RRHm = 0.7,#relative risk of hospitalization given infection for those with maternal immunity
                recover1 = 365.25*.5, #days spent in R1 compartment
                recover2 = 365.25*.5, #days spent in R2 compartment
                recover3 = 365.25,#days spent in R3 compartment
                recover4 = 365.25,#days spent in R4 compartment
                um=um, #net death and migration rate
                rho1=0.75,#relative infectiousness following first infection
                rho2=0.51,#relative infectiousness following 2+ infections
                dur.days1=10, #duration of first infection (days)
                dur.days2=7, #duration of second infection (days)
                dur.days3=5, #duration of third+ infection (days)
                yinit.matrix=yinit, #initial compartments
                q=1,
                c2=contact,
                sigma1=0.76,#reduced susceptibility following first infection
                sigma2=0.6,#reduced susceptibility following second infection
                sigma3=0.4,#reduced susceptibility following third infection (and with maternal immunity)
                length.step = 7,
                time.step='week',
                seed=seed,
                waningN = 180,
                waningV=180,
                waningS = 730.5,
                RRIn=1,
                RRIv=1,
                RRIs=1,
                RRHn=.2,
                RRHv=.3,
                RRHs=.1)

  return(list(parmset,yinit, yinit.vector))
}
