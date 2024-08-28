## code to prepare `DATASET` dataset goes here
#rm(list=ls())
library(lubridate)
library(cdcfluview)
library(imputeTS)
library(tidyverse)
library(lubridate)
library(zoo)

`%notin%` = Negate(`%in%`)


# Data from RSV-NET -------------------------------------------
dat = read.csv("data-raw/Weekly_Rates_of_Laboratory-Confirmed_RSV_Hospitalizations_from_the_RSV-NET_Surveillance_System_20240826.csv") %>%
  filter(Race=="All",Sex=="All",State!="RSV-NET",
         Age.Category %in% c("0-<6 months","6mo-<12 months","1-4 years","5-17 years","18-49 years","50-64 years","65+ years")) %>%
  mutate(Age.Category = ifelse(Age.Category=="0-<6 months","<6m",
                               ifelse(Age.Category=="6mo-<12 months","6-11m",Age.Category)),
         mmwr_week(Week.ending.date),
         date= as.Date(Week.ending.date)) %>%
  select(-Sex,-Race,-Cumulative.Rate,-Type) %>%
  filter(!is.na(Rate))


#check starting season for each state and age group

start = dat %>%
  group_by(State) %>%
  summarize(start = min(Season))


#fill in child data for missing years
age_groups_of_interest <- c('<6m', '6-11m', '1-4 years', '5-17 years')

#start with states that have begin in 2016

group1 = dat %>% filter(State %in% c("California","Georgia","Maryland","Minnesota","New York","Oregon"))

grid1 = expand.grid(unique(group1$date),unique(group1$State),unique(group1$Age.Category))

average_rates1 <- group1 %>%
  filter(Age.Category %in% age_groups_of_interest, Season %in% c('2018-19', '2019-20')) %>%
  group_by(State, Age.Category, mmwr_week) %>%
  summarize(Average.Rate = mean(Rate, na.rm = TRUE)) %>%
  ungroup()

join1  = group1 %>%
  full_join(grid1, by=c("date"="Var1","State"="Var2","Age.Category"="Var3")) %>%
  select(State,date,Age.Category,Rate) %>%
  mutate(mmwr_week(date)) %>%
  full_join(average_rates1, by=c("State","Age.Category","mmwr_week")) %>%
  mutate(adj_rate = ifelse(is.na(Rate),Average.Rate,Rate)) %>%
  filter(date>='2016-10-08') %>%
  select(date,State,Age.Category,adj_rate)

ggplot(join1)+
  geom_line(aes(x=date,y=adj_rate))+
  facet_grid(cols=vars(State),rows=vars(Age.Category),scales="free")

dates1 = seq(from=as.Date('2016-10-08'),to=as.Date('2024-08-17'),by="week")
grid1.2 = expand.grid(dates1,unique(join1$State),unique(join1$Age.Category))

complete1 = join1 %>%
  full_join(grid1.2, by=c("date"="Var1","State"="Var2","Age.Category"="Var3")) %>%
  mutate(adj_rate=na_interpolation(adj_rate))

ggplot(complete1)+
  geom_line(aes(x=date,y=adj_rate))+
  facet_grid(cols=vars(State),rows=vars(Age.Category),scales="free")

#States that start in 2017

group2 = dat %>% filter(State %in% c("New Mexico"))

grid2 = expand.grid(unique(group2$date),unique(group2$State),unique(group2$Age.Category))

average_rates2 <- group2 %>%
  filter(Age.Category %in% age_groups_of_interest, Season %in% c('2018-19', '2019-20')) %>%
  group_by(State, Age.Category, mmwr_week) %>%
  summarize(Average.Rate = mean(Rate, na.rm = TRUE)) %>%
  ungroup()

join2  = group2 %>%
  full_join(grid2, by=c("date"="Var1","State"="Var2","Age.Category"="Var3")) %>%
  select(State,date,Age.Category,Rate) %>%
  mutate(mmwr_week(date)) %>%
  full_join(average_rates2, by=c("State","Age.Category","mmwr_week")) %>%
  mutate(adj_rate = ifelse(is.na(Rate),Average.Rate,Rate)) %>%
  select(date,State,Age.Category,adj_rate)


ggplot(join2)+
  geom_line(aes(x=date,y=adj_rate))+
  facet_grid(cols=vars(State),rows=vars(Age.Category),scales="free")

dates2 = seq(from=as.Date('2017-10-07'),to=as.Date('2024-08-17'),by="week")
grid2.2 = expand.grid(dates2,unique(join2$State),unique(join2$Age.Category))

complete2 = join2 %>%
  full_join(grid2.2, by=c("date"="Var1","State"="Var2","Age.Category"="Var3")) %>%
  mutate(adj_rate=na_interpolation(adj_rate))

ggplot(complete2)+
  geom_line(aes(x=date,y=adj_rate))+
  facet_grid(cols=vars(State),rows=vars(Age.Category),scales="free")


#States that only start in 2018 - don't need to fill in child data for these
group3 = dat %>% filter(State %in% c("Colorado","Connecticut","Michigan","Utah")) %>%
  select(date,State,Age.Category,"adj_rate"=Rate)

dates3 = seq(from=as.Date('2018-10-06'),to=as.Date('2024-08-17'),by="week")
grid3.2 = expand.grid(dates3,unique(group3$State),unique(group3$Age.Category))

complete3 = group3 %>%
  full_join(grid3.2, by=c("date"="Var1","State"="Var2","Age.Category"="Var3")) %>%
  mutate(adj_rate=na_interpolation(adj_rate))

ggplot(complete3)+
  geom_line(aes(x=date,y=adj_rate))+
  facet_grid(cols=vars(State),rows=vars(Age.Category),scales="free")



complete_data = rbind(complete1, complete2,complete3)


pop1 = read.delim("data-raw/Single-Race Population Estimates 2010-2020 by State and Single-Year Age.txt")
pop2 = read.delim("data-raw/Single-Race Population Estimates 2020-2022 by State and Single-Year Age.txt")

pop=rbind(pop1,pop2) %>%
  filter(!is.na(States.Code)) %>%
  mutate(age_group = case_when(Five.Year.Age.Groups.Code=="1"~"1",
                               Five.Year.Age.Groups.Code=="1-4"~"1-4 years",
                               Five.Year.Age.Groups.Code %in% c("5-9","10-14","15-19")~"5-17 years",
                               Five.Year.Age.Groups.Code %in% c("20-24","25-29","30-34","35-39","40-44","45-49")~"18-49 years",
                               Five.Year.Age.Groups.Code %in% c("50-54","55-59","60-64")~"50-64 years",
                               Five.Year.Age.Groups.Code %in% c("65-69","70-74","75-79","80-84","85+")~"65+ years")) %>%
  group_by(States,age_group,Yearly.July.1st.Estimates) %>%
  summarize(population = sum(Population))


child1 = pop %>% filter(age_group=="1") %>%
  mutate(age_group="<6m",
         population = population*0.5)

child2 = pop %>% filter(age_group=="1") %>%
  mutate(age_group="6-11m",
         population = population*0.5)

pop_new = pop %>% filter(age_group!="1") %>%
  bind_rows(child1,child2)


dat2 = complete_data %>%
  mutate(year = year(date)) %>%
  left_join(pop_new, by=c("Age.Category"="age_group","State"="States","year"="Yearly.July.1st.Estimates")) %>%
  arrange(Age.Category,State,date) %>%
  mutate(population = na_locf(population),
         count = adj_rate/100000*population,
         new_age = factor(Age.Category, levels=c("<6m","6-11m","1-4 years","5-17 years","18-49 years","50-64 years","65+ years"),
                          labels=c("<6m","6-11m","1-4 years","5-64 years","5-64 years","5-64 years","65+ years"))) %>%
  group_by(date,State,new_age) %>%
  summarize(count=sum(count),
            population=sum(population))


ggplot(dat2)+
  geom_line(aes(x=date,y=count))+
  facet_grid(rows=vars(new_age),cols=vars(State),scales="free")


dist_overall = dat2 %>%
  group_by(State,new_age) %>%
  summarize(total = sum(count)) %>%
  ungroup() %>%
  group_by(State) %>%
  mutate(prop = total/sum(total))


dist_by_time = dat2 %>%
  mutate(mmwr_week(date),
         season = ifelse(mmwr_week>=40, paste0(mmwr_year,"-",mmwr_year+1),paste0(mmwr_year-1,"-",mmwr_year))) %>%
  filter(season %notin% c("2016-2017","2017-2018")) %>%
  mutate(period = ifelse(season %in% c("2018-2019","2019-2020"),"pre-pandemic",
                         ifelse(season %in% c("2021-2022","2022-2023"),"rebound",
                                             ifelse(season %in% c("2023-2024"),"current","pandemic")))) %>%
  group_by(State,period,new_age) %>%
  summarize(total = sum(count)) %>%
  ungroup() %>%
  group_by(State,period) %>%
  mutate(prop = total/sum(total)) %>%
  ungroup()


age_distribution = dist_by_time %>% filter(period=="pre-pandemic") %>%
  select("state"=State,"age_group"=new_age,"proportion"=prop)


scaling = dist_by_time %>%
  filter(period %in% c("pre-pandemic","current")) %>%
  mutate(adjust = ifelse(period=="pre-pandemic",total/2,total)) %>%
  pivot_wider(id_cols=c(State, new_age),names_from=period,values_from=adjust) %>%
  mutate(diff = current/`pre-pandemic`)



timeseries = dat2 %>%
  group_by(date, State) %>%
  summarize(count=sum(count),
            population=sum(population)) %>%
  ungroup() %>%
  group_by(State) %>%
  mutate(smooth = round(rollmean(count,k=3,align="center",fill="extend"))) %>%
  select("state"=State,date,"value"=smooth)


ggplot(dat3)+
  geom_line(aes(x=date,y=value))+
  facet_grid(rows=vars(State),scales="free")

usethis::use_data(timeseries, overwrite = TRUE)
usethis::use_data(age_distribution, overwrite = TRUE)




#Age distributions



#Scaling pre-2020 data?
