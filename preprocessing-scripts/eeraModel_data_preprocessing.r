# ----------------------------------
# libraries
# ----------------------------------

library(dplyr)
library(readr)
library(tidyr)
library(readxl)
library(lubridate)



# ----------------------------------
# load functions
# ----------------------------------

source("eeraModel_data_preprocessing_functions.R")


#======================
# Settings
#======================
do.write = F
upto = "2020-06-09"


# ----------------------------------
# create data use for model
# ----------------------------------

#create the case and death data

# cases pre change in recording practices
WATTY62DATA <- "https://raw.githubusercontent.com/watty62/Scot_covid19/master/data/processed/regional_cases.csv"
WATTY62POP <- "https://raw.githubusercontent.com/watty62/Scot_covid19/master/data/processed/HB_Populations.csv"
WATTY62DEATHS <- "https://raw.githubusercontent.com/watty62/Scot_covid19/master/data/processed/regional_deaths.csv"

pop_df <- readr::read_csv(file =WATTY62POP , col_types = cols())
case_df = readr::read_csv(file = WATTY62DATA, col_types = cols())
death_df <- readr::read_csv(file = WATTY62DEATHS, col_types = cols()) 

# cases post change in recording practices
url_scot_country <- "https://www.gov.scot/binaries/content/documents/govscot/publications/statistics/2020/04/coronavirus-covid-19-trends-in-daily-data/documents/trends-in-number-of-people-in-hospital-with-confirmed-or-suspected-covid-19/trends-in-number-of-people-in-hospital-with-confirmed-or-suspected-covid-19/govscot%3Adocument/Trends%2Bin%2Bdaily%2BCOVID-19%2Bdata%2B-%2B110520.xlsx?forceDownload=true"
httr::GET(url_scot_country, httr::write_disk(scot_data_gov <- tempfile(fileext = "govcoviddata.xlsx")))

df_cases_tot <- readxl::read_excel(scot_data_gov, sheet = "Table 5 - Testing", skip = 3)%>%
  rename("Date" = "...1", 
         "confirmed_cases" = "Positive")%>%
  mutate(new_cases = `Daily...5`)%>%
  select(Date, confirmed_cases, new_cases)

df_deaths_tot <- readxl::read_excel(scot_data_gov, sheet = "Table 8 - Deaths", skip =2) %>%
  select(Date, `Number of COVID-19 confirmed deaths registered to date`) %>%
  rename(deaths = `Number of COVID-19 confirmed deaths registered to date`) %>%
  mutate(new_deaths = deaths - replace_na(lag(deaths),0))%>%
  select(Date, deaths, new_deaths)


extend_case = extend_epi_records(base_data = tidy_case_data(case_df), extra_data = df_cases_tot, data_name='confirmed_cases')
extend_death = extend_epi_records(base_data = tidy_death_data(death_df), extra_data = df_deaths_tot, data_name='deaths')


scot_data = create_history(extend_case,tidy_pop_data(pop_df), upto = upto)$data
scot_death = create_history(extend_death,tidy_pop_data(pop_df), upto = upto, 
                            join = T, data_join = extend_case)$data


if(do.write) write.table(scot_data,paste0("./data/scot_data_",upto,".csv"),col.names = T,row.names = F,sep=",")
if(do.write) write.table(scot_death,paste0("./data/scot_deaths_",upto,".csv"),col.names = T,row.names = F,sep=",")

# create the age mixing matrices



polymoduk_norm=read_xlsx("./data/contact_matrices_152_countries/MUestimates_all_locations_2.xlsx",
                         sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                               "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"))

polymoduk_home=read_xlsx("./data/contact_matrices_152_countries/MUestimates_home_2.xlsx",
                         sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                               "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"))

polymoduk_other=read_xlsx("./data/contact_matrices_152_countries/MUestimates_other_locations_2.xlsx",
                          sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                                "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"))

polymoduk_work=read_xlsx("./data/contact_matrices_152_countries/MUestimates_work_2.xlsx",
                         sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                               "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"))

polymoduk_school=read_xlsx("./data/contact_matrices_152_countries/MUestimates_school_2.xlsx",
                           sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                                 "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"))


waifw_norm = create_mixmat(polymoduk_norm)
waifw_home = create_mixmat(polymoduk_home)
waifw_other = create_mixmat(polymoduk_other)
waifw_school = create_mixmat(polymoduk_school)
waifw_work = create_mixmat(polymoduk_work)
waifw_sdist = waifw_home + waifw_other

if(do.write) write.table(waifw_norm,"./data/waifw_norm.csv",col.names = F,row.names = F,sep=",")
if(do.write) write.table(waifw_home,"./data/waifw_home.csv",col.names = F,row.names = F,sep=",")
if(do.write) write.table(waifw_school,"./data/waifw_school.csv",col.names = F,row.names = F,sep=",")
if(do.write) write.table(waifw_work,"./data/waifw_work.csv",col.names = F,row.names = F,sep=",")
if(do.write) write.table(waifw_sdist,"./data/waifw_sdist.csv",col.names = F,row.names = F,sep=",")


# create the data with demographic data (age structure)

#source: https://www.scotlandscensus.gov.uk/ods-web/download/getDownloadFile.html?downloadFileIds=Scotland%20stds

census_age = read_csv(file ="./data/Scotland-HealthBoard-census-2011/DC1117SC.csv",
                      skip=5 , col_names = c("Health_Board","age","All people","Males","Females"), col_types = cols())

order_name = create_history(tidy_case_data(case_df),tidy_pop_data(pop_df), upto = upto)$name

scot_age <- create_age_pop(census_df = census_age, order_name)

if(do.write) write.table(scot_age, "./data/scot_age.csv",col.names = F,row.names = F,sep=",")

# create the data with age-structured parameters

p_h=c(0.143,	0.1141,	0.117,	0.102,	0.125,	0.2,	0.303,	0.114525) #https://www.ecdc.europa.eu/sites/default/files/documents/RRA-seventh-update-Outbreak-of-coronavirus-disease-COVID-19.pdf
cfr=c(0.001,	0.002,	0.002	,0.004,	0.013,	0.036,	0.114,	0.00525) #https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.12.2000256. NOTE: case defined as infected individuals (either symptomatic or asymptomatic)

if(do.write) write.table(create_age_epiparaam(p_h, cfr) , "./data/cfr_byage.csv",col.names = F,row.names = F,sep=",")
