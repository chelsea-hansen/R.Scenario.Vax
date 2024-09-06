#' Download and format data for RSV dynamic transmission model (function)
#'
#' Download and format data for RSV dynamic transmission model (function).
#' The function returns a list of 3 items. These will be used in the fit_model() function
#'
#' @param state_or_county describe if you are modeling a "state" or a "county"
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

contact = matrix(c(0.064,	0.064,	0.064,	0.064,	0.064,	0.064,	0.384,	1.152,	0.406666667,	0.371083333,	1.504666667,	0.752333333,	0.157583333,	0.068625,
                      0.064, 0.064,	0.064,	0.064,	0.064,	0.064,	0.384,	1.152,	0.406666667,	0.371083333,	1.504666667,	0.752333333,	0.157583333,	0.068625,
                      0.064,	0.064,	0.064,	0.064,	0.064,	0.064,	0.384,	1.152,	0.406666667,	0.371083333,	1.504666667,	0.752333333,	0.157583333,	0.068625,
                      0.064,	0.064,	0.064,	0.064,	0.064,	0.064,	0.384,	1.152,	0.406666667,	0.371083333,	1.504666667,	0.752333333,	0.157583333,	0.068625,
                      0.064,	0.064,	0.064,	0.064,	0.064,	0.064,	0.384,	1.152,	0.406666667,	0.371083333,	1.504666667,	0.752333333,	0.157583333,	0.068625,
                      0.064,	0.064,	0.064,	0.064,	0.064,	0.064,	0.384,	1.152,	0.406666667,	0.371083333,	1.504666667,	0.752333333,	0.157583333,	0.068625,
                      0.384,	0.384,	0.384,	0.384,	0.384,	0.384,	0.384,	1.152,	0.406666667,	0.371083333,	1.504666667,	0.752333333,	0.157583333,	0.068625,
                      1.152,	1.152,	1.152,	1.152,	1.152,	1.152,	1.152,	1.152,	0.406666667,	0.371083333,	1.504666667,	0.752333333,	0.157583333,	0.068625,
                      0.406666667,	0.406666667,	0.406666667,	0.406666667,	0.406666667,	0.406666667,	0.406666667,	0.406666667,	6.64,	1.735,	3.355,	2.23,	0.47,	0.175,
                      0.371083333,	0.371083333,	0.371083333,	0.371083333,	0.371083333,	0.371083333,	0.371083333,	0.371083333,	1.735,	8.055,	2.7275,	2.7925,	0.62,	0.3175,
                      1.504666667,	1.504666667,	1.504666667,	1.504666667,	1.504666667,	1.504666667,	1.504666667,	1.504666667,	3.355,	2.7275,	4.9075,	3.715,	0.705,	0.2625,
                      0.752333333,	0.752333333,	0.752333333,	0.752333333,	0.752333333,	0.752333333,	0.752333333,	0.752333333,	2.23,	2.7925,	3.715,	4.228,	1.068,	0.581,
                      0.157583333,	0.157583333,	0.157583333,	0.157583333,	0.157583333,	0.157583333,	0.157583333,	0.157583333,	0.47,	0.62,	0.705,	1.068,	1.755,	1.07,
                      0.068625,	0.068625,	0.068625,	0.068625,	0.068625,	0.068625,	0.068625,	0.068625,	0.175,	0.3175,	0.2625,	0.581,	1.07,	1.47),nrow=14,ncol=14)

  burn_in=matrix(c(0.64428,	0.32564,	0.00028,	0.01625,	0.01288,	0.00000,	0.00040,	0.00027,	0.00000,	0.00001,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,0,0,0,
                   0.41100,	0.49064,	0.00039,	0.05069,	0.04445,	0.00002,	0.00161,	0.00115,	0.00000,	0.00003,	0.00001,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,0,0,0,
                   0.25972,	0.55316,	0.00043,	0.09263,	0.08730,	0.00003,	0.00375,	0.00286,	0.00000,	0.00009,	0.00003,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,0,0,0,
                   0.16343,	0.56219,	0.00043,	0.13038,	0.13077,	0.00005,	0.00681,	0.00563,	0.00000,	0.00023,	0.00009,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,  0.00000,0,0,0,
                   0.10303,	0.54786,	0.00041,	0.15800,	0.16902,	0.00006,	0.01097,	0.00990,	0.00000,	0.00052,	0.00021,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,0,0,0,
                   0.06531,	0.52368,	0.00039,	0.17476,	0.20152,	0.00008,	0.01658,	0.01620,	0.00000,	0.00104,	0.00043,  0.00000,	0.00001,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,0,0,0,
                   0.01484,	0.34654,	0.00026,	0.15149,	0.29169,	0.00012,	0.05890,	0.09625,	0.00002,	0.01974,	0.01879,	0.00000,	0.00135,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,0,0,0,
                   0.00132,	0.11379,	0.00010,	0.06344,	0.21684,	0.00010,	0.07157,	0.19973,	0.00005,	0.08108,	0.22240,	0.00002,	0.02958,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,0,0,0,
                   0.00007,	0.01064,	0.00002,	0.01646,	0.04014,	0.00005,	0.03815,	0.07573,	0.00005,	0.09686,	0.51239,	0.00011,	0.20933,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,  0.00000,0,0,0,
                   0.00000,	0.00042,	0.00000,	0.00141,	0.00317,	0.00000,	0.00519,	0.01052,	0.00001,	0.02397,	0.62474,	0.00014,	0.33042,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,0,0,0,
                   0.00000,	0.00001,	0.00000,	0.00005,	0.00010,	0.00000,	0.00026,	0.00053,	0.00000,	0.00208,	0.60429,	0.00017,	0.39251,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,0,0,0,
                   0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00001,	0.00002,	0.00000,	0.00011,	0.62848,	0.00016,	0.37121,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,0,0,0,
                   0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00001,	0.00000,	0.00001,	0.85130,	0.00005,	0.14863,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,0,0,0,
                   0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.91779,	0.00003,	0.08218,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,	0.00000,0,0,0),
                 nrow=14,ncol=22, byrow=TRUE)

  colnames(burn_in)= c("M0","S0","I1","R1","S1","I2","R2","S2","I3","R3","S3","I4","R4","Mn1","Mn2","Mv1","Mv2","N1","N2","Si","Vs1","Vs2")



pop = tidycensus::get_estimates(
    geography = state_or_county,
    product = "characteristics",
    breakdown = c("AGEGROUP"),
    breakdown_labels = TRUE,
    vintage=2022,
    state = state_abbr,
    county = county_name)%>%
    select("AGEGROUP","value")

  pop = tail(pop,18)

  new_agegrps = c("0-4","5-9",rep("10-19",2),rep("20-39",4),rep("40-64",5),rep("65-74",2),rep("75+",3))

  pop = pop %>%
    dplyr::mutate(new_age = factor(new_agegrps,levels=unique(new_agegrps))) %>%
    dplyr::group_by(.data$new_age) %>%
    dplyr::summarize(population = sum(.data$value))

  adults=c(pop$population[2:7])
  infants = rep(pop$population[1]*.034,6)
  children = c(pop$population[1]*.2, rep(pop$population[2]*.6))
  all_ages = c(infants,children, adults)
  yinit = burn_in*all_ages


  #format the data for feeding into the model
  N_ages <- nrow(yinit)
  agenames <- c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y","5-9Y","10-19Y","20-39Y","40-64Y","65-74Y","75+Y")
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
    vintage=2022,
    state = state_abbr,
    county = county_name)

  birth_rate = vital_stats %>%  dplyr::filter(.data$variable=="BIRTHS") %>%  dplyr::select("value")/52.1775
  birth=c(birth_rate$value, rep(0,13))

  mig_rate = vital_stats %>%  dplyr::filter(.data$variable=="RNETMIG") %>%  dplyr::select("value")/1000
  um = mig_rate$value*-1/52.1775


  parmset<-list(PerCapitaBirthsYear=birth,
                WidthAgeClassMonth=c(rep(2,times=6), 12,12*3, 60, 120, 240, 300,120,96),#time spend in each age class (months)
                DurationMatImmunityDays = 90,#duration of maternal immunity (days)
                RRHm = 1,#relative risk of hospitalization given infection for those with maternal immunity
                recover1 = 365.25*.5, #days spent in R1 compartment
                recover2 = 365.25*.5, #days spent in R2 compartment
                recover3 = 358.9,#days spent in R3 compartment
                recover4 = 358.9,#days spent in R4 compartment
                um=um, #net death and migration rate
                rho1=0.75,#relative infectiousness following first infection
                rho2=0.51,#relative infectiousness following 2+ infections
                dur.days1=10, #duration of first infection (days)
                dur.days2=7, #duration of second infection (days)
                dur.days3=5, #duration of third+ infection (days)
                yinit.matrix=yinit, #initial compartments
                q=1,
                c2=contact,
                sigma1=0.89,#reduced susceptibility following first infection
                sigma2=0.89*0.81,#reduced susceptibility following second infection
                sigma3=0.89*0.81*0.33,#reduced susceptibility following third infection (and with maternal immunity)
                length.step = 7,
                time.step='week',
                seed=seed,
                waningN = 90,
                waningN2 = 90,
                waningV=90,
                waningV2=90,
                waningS = 730.5,
                RRIn=1,
                RRIv=1,
                RRIs=1,
                RRHn1=.2,
                RRHn2=.2,
                RRHv1=.45,
                RRHv2=.45,
                RRHs=.2)

  return(list(parmset,yinit, yinit.vector))
}
