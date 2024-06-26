## code to prepare `DATASET` dataset goes here

library(lubridate)
library(cdcfluview)
library(imputeTS)
library(totalcensus)
library(tidyverse)
library(lubridate)

`%notin%` = Negate(`%in%`)


# Pull data from RSV SMH github -------------------------------------------


dat1 = read_csv("https://raw.githubusercontent.com/midas-network/rsv-scenario-modeling-hub/main/target-data/2024-06-21_rsvnet_hospitalization.csv") %>%
  mutate(state = convert_fips_to_names(location),
         state = ifelse(is.na(state),"US",state)) %>%
  filter(age_group %in% c("0-0.49","0.5-0.99","1-1.99","2-4","5-17","18-49","50-64","65-130"),
         target=="inc hosp") %>%
  arrange(date) %>%
  mutate(month = month(date),
         year = year(date),
         week = epiweek(date),
         season = ifelse(month<=6,paste0(year-1,"-",year),paste0(year,"-",year+1)), #July to June seasons
         age_group = factor(age_group, levels=c("0-0.49","0.5-0.99","1-1.99","2-4","5-17","18-49","50-64","65-130"),
                            labels=c("<6m","6-11m","1-2y","2-4y","5-17y","18-49y","50-64y","65+y")),flag=ifelse(age_group %in% c("<6m","6-11m","1-2y","2-4y","5-17y")&
                                                                                                                  season %in% c("2014-2015","2015-2016","2016-2017","2017-2018"),1,0)) %>%
  filter(flag==0) %>%
  select(-flag)


age_distribution = dat1 %>%
  filter(season %notin% c("2014-2015","2015-2016","2016-2017","2017-2018")) %>%
  mutate(group2 = factor(age_group, levels=c("<6m","6-11m","1-2y","2-4y","5-17y","18-49y","50-64y","65+y"),
                         labels=c("<6m","6-11m","1-4y","1-4y","5-59y","5-59y","5-59y","60+y"))) %>%
    group_by(group2,state) %>%
    summarize(total=sum(value,na.rm=TRUE)) %>%
  ungroup() %>%
  group_by(state) %>%
  mutate(prop = total/sum(total)) %>%
  arrange(state, group2) %>%
  select(state,"age_group"=group2,"proportion"=prop)

# Prepare time series data ------------------------------------------------
#state have different lengths of data

state_data = data.frame(date=NA,state=NA, hosp=NA)

seasons_6 = c("GA","MD","MN","NY","OR","TN")

for(s in 1:length(seasons_6)){
  dat2 = dat1 %>%
    filter(state==seasons_6[s])

  #fill in 2014-15 - 2017-18 seasons with available data
  fill1 = dat2 %>%
    filter(season %in% c("2018-2019","2019-2020"),
           age_group %in% c("<6m","6-11m","1-2y","2-4y","5-17y")) %>%
    mutate(week = case_when(week>=40 ~week+1,
                            week<40 ~week),
           year = case_when(year==2018~2014,
                            year==2019~2015,
                            year==2020~2016),
           season = case_when(season=="2018-2019"~"2014-2015",
                              season=="2019-2020"~"2015-2016"),
           date = mmwr_week_to_date(year=year, week=week,day=7),
           dup = case_when(date=='2016-01-09' & week==53 ~1,
                           date!='2016-01-09'| week!=53 ~0)) %>%
    filter(dup==0) %>%
    select(-dup)


  fill2 = dat2 %>%
    filter(season %in% c("2018-2019","2019-2020"),
           age_group %in% c("<6m","6-11m","1-2y","2-4y","5-17y")) %>%
    mutate(year = case_when(year==2018~2016,
                            year==2019~2017,
                            year==2020~2018),
           season = case_when(season=="2018-2019"~"2016-2017",
                              season=="2019-2020"~"2017-2018"),
           date = mmwr_week_to_date(year=year, week=week,day=7))

  dat3 = rbind(dat2, fill1, fill2)


  ts = dat3 %>%
    filter(date>='2014-10-11') %>%
    group_by(date) %>%
    summarize(hosp = sum(value)) %>%
    ungroup() %>%
    mutate(state=seasons_6[s],
           hosp = na_interpolation(hosp,option="linear"))
state_data=rbind(state_data,ts) %>% filter(!is.na(date)) %>% mutate(date=as.Date(date))
}

#states that only have data starting in 2018 for all ages
seasons_2 = c("CO","CT","MI","UT")

for(s in 1:length(seasons_2)){
  dat2 = dat1 %>%
    filter(state==seasons_2[s])

  ts = dat2 %>%
    filter(date>='2018-10-06') %>%
    group_by(date) %>%
    summarize(hosp = sum(value)) %>%
    ungroup() %>%
    mutate(state=seasons_2[s],
           hosp = na_interpolation(hosp,option="linear"))

  state_data=rbind(state_data,ts) %>% filter(!is.na(date)) %>% mutate(date=as.Date(date))
}

#California data starts in 2015
dat2 = dat1 %>%
  filter(state=="CA")

#fill in 2014-15 - 2017-18 seasons with available data
fill1 = dat2 %>%
  filter(season %in% c("2019-2020"),
         age_group %in% c("<6m","6-11m","1-2y","2-4y","5-17y")) %>%
  mutate(week = case_when(week>=40 ~week+1,
                          week<40 ~week),
         year = case_when(year==2019~2015,
                          year==2020~2016),
         season = case_when(season=="2019-2020"~"2015-2016"),
         date = mmwr_week_to_date(year=year, week=week,day=7),
         dup = case_when(date=='2016-01-09' & week==53 ~1,
                         date!='2016-01-09'| week!=53 ~0)) %>%
  filter(dup==0) %>%
  select(-dup)


fill2 = dat2 %>%
  filter(season %in% c("2018-2019","2019-2020"),
         age_group %in% c("<6m","6-11m","1-2y","2-4y","5-17y")) %>%
  mutate(year = case_when(year==2018~2016,
                          year==2019~2017,
                          year==2020~2018),
         season = case_when(season=="2018-2019"~"2016-2017",
                            season=="2019-2020"~"2017-2018"),
         date = mmwr_week_to_date(year=year, week=week,day=7))

dat3 = rbind(dat2, fill1, fill2)


ts = dat3 %>%
  filter(date>='2015-10-03') %>%
  group_by(date) %>%
  summarize(hosp = sum(value)) %>%
  ungroup() %>%
  mutate(state="CA",
         hosp = na_interpolation(hosp,option="linear"))
state_data=rbind(state_data,ts) %>% filter(!is.na(date)) %>% mutate(date=as.Date(date))

#new mexico data starts in 2017

dat2 = dat1 %>%
  filter(state=="NM")

fill2 = dat2 %>%
  filter(season %in% c("2019-2020"),
         age_group %in% c("<6m","6-11m","1-2y","2-4y","5-17y")) %>%
  mutate(year = case_when(year==2019~2017,
                          year==2020~2018),
         season = case_when(season=="2019-2020"~"2017-2018"),
         date = mmwr_week_to_date(year=year, week=week,day=7))

dat3 = rbind(dat2,fill2)


ts = dat3 %>%
  filter(date>='2017-10-07') %>%
  group_by(date) %>%
  summarize(hosp = sum(value)) %>%
  ungroup() %>%
  mutate(state="NM",
         hosp = na_interpolation(hosp,option="linear"))
state_data=rbind(state_data,ts) %>% filter(!is.na(date)) %>% mutate(date=as.Date(date))

ggplot(data=state_data)+
  geom_line(aes(x=date,y=hosp))+
  facet_wrap(~state,scales='free')

timeseries=state_data %>% filter(date<'2020-05-01'|date>='2021-05-01') %>%
  mutate(date = ifelse(date <'2020-05-01', date %m+% years(1), date),
         date=as.Date(date)) %>%
  rename("value"=hosp)

ggplot(data=timeseries)+
  geom_line(aes(x=date,y=value))+
  facet_wrap(~state,scales='free')

usethis::use_data(timeseries, overwrite = TRUE)
usethis::use_data(age_distribution, overwrite = TRUE)
