#======================
# libraries
#======================
library(dplyr)
library(readr)
library(tidyr)
library(readxl)
library(lubridate)


#======================
# Settings
#======================
do.write = F


#======================
# Disease data
#======================
# Data published by Public Health England 
# https://www.gov.uk/government/publications/covid-19-track-coronavirus-cases

# Scottish Data from HPS
WATTY62DATA <- "https://raw.githubusercontent.com/watty62/Scot_covid19/master/data/processed/regional_cases.csv"
WATTY62POP <- "https://raw.githubusercontent.com/watty62/Scot_covid19/master/data/processed/HB_Populations.csv"
WATTY62DEATHS <- "https://raw.githubusercontent.com/watty62/Scot_covid19/master/data/processed/regional_deaths.csv"
scot_pop <- read_csv(file =WATTY62POP , col_types = cols()) %>% 
  mutate(hcw_per_rata = Population/sum(Population))
scot_tot_pop = as.data.frame(list(`Name` = "Scotland",`Population` = sum(scot_pop$Population)))%>% 
  mutate(hcw_per_rata = Population/sum(Population)) %>% 
  mutate_if(is.factor,as.character)
scot_pop <- bind_rows(scot_pop,scot_tot_pop)
scot_deaths <- read_csv(file = WATTY62DEATHS, col_types = cols()) 
scot_deaths[scot_deaths == "x"] <- -999
scot_deaths <- scot_deaths %>% 
  rename(`Scotland`= `Grand Total`) %>%
  mutate(Date = lubridate::dmy(Date))
scot_data <- read_csv(file = WATTY62DATA, col_types = cols()) %>%
  rename(`Scotland`= `Grand Total`) %>%
  mutate(Date = lubridate::dmy(Date),
         Area = "Scotland")

#formating  the case data (HPS) for model input
temp <- right_join(scot_pop %>% select(Name,Population),as.data.frame(t(
  as.matrix(
    scot_data %>% filter(Date <= max(scot_data$Date)) %>% 
      select(-Date,-Area)))) %>% tibble::rownames_to_column(), by = c("Name"="rowname")) %>% 
  select(-Name)
colnames(temp) <- (1:ncol(temp))-2

order_row= as.data.frame(t(
  as.matrix(
    scot_data %>% filter(Date <= max(scot_data$Date)) %>% 
      select(-Date,-Area)))) %>% tibble::rownames_to_column() %>% select(rowname)

if(do.write) write.table(temp,paste0("./data/scot_data_",max(scot_data$Date),".csv"),col.names = T,row.names = F,sep=",")

#formating  the death data (HPS) for model input
dtemp <- left_join(scot_data %>% select(Date),scot_deaths,by="Date") %>% 
  replace(is.na(.),0)

dtemp <- right_join(scot_pop %>% select(Name,Population),as.data.frame(t(
  as.matrix(
    dtemp %>% filter(Date <= max(scot_data$Date)) %>% 
      select(-Date)))) %>% tibble::rownames_to_column(), by = c("Name"="rowname")) %>% 
  select(-Name) %>% 
  mutate_if(is.factor,as.character) %>% 
  mutate_if(is.character,as.numeric) 
colnames(dtemp) <- (1:ncol(dtemp))-2

if(do.write) write.table(dtemp,paste0("./data/scot_deaths_",max(scot_data$Date),".csv"),col.names = T,row.names = F,sep=",")



#======================
# Age mixing matrices
#======================
# from https://doi.org/10.1371/journal.pcbi.1005697

#all location
polymoduk=read_xlsx("./data/contact_matrices_152_countries/MUestimates_all_locations_2.xlsx",
                    sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                          "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+")) %>%
  mutate(from = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"),
         group_age = ifelse(from %in% c("00-04","05-09","10-14","15-19"),"Under20",
                            ifelse(from %in% c("20-24","25-29"),"20-29",
                                   ifelse(from %in% c("30-34","35-39"),"30-39",
                                          ifelse(from %in% c("40-44","45-49"),"40-49",
                                                 ifelse(from %in% c("50-54","55-59"),"50-59",
                                                        ifelse(from %in% c("60-64","65-69"),"60-69",
                                                               ifelse(from %in% c("70-74","75+"),"Over70",NA)))))))) %>% 
  mutate(`Under20` = (`00-04`+`05-09`+`10-14`+`15-19`)/4,
         `20-29` = (`20-24`+`25-29`)/2,
         `30-39` =(`30-34`+`35-39`)/2,
         `40-49` = (`40-44`+`45-49`)/2,
         `50-59` = (`50-54`+`55-59`)/2,
         `60-69` = (`60-64`+`65-69`)/2,
         `Over70` = (`70-74`+`75+`)/2,
         `HCW`= (`25-29`+`30-34`+`35-39`+`40-44`+`45-49`+`50-54`)/6) %>% 
  ungroup() %>% 
  select(group_age,`Under20`,`20-29`,`30-39`,`40-49`,`50-59`,`60-69`,`Over70`,`HCW`) %>% 
  group_by(group_age) %>% 
  summarise_at(vars(-group_cols()),mean) %>% 
  ungroup() %>% 
  mutate(group_age = factor(group_age,levels = c("Under20","20-29","30-39","40-49","50-59","60-69","Over70","HCW"))) %>% 
  arrange(group_age) %>% 
  select(-group_age) 


polymoduk_HCW=read_xlsx("./data/contact_matrices_152_countries/MUestimates_all_locations_2.xlsx",
                        sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                              "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+")) %>%
  mutate(from = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"),
         group_age = ifelse(from %in% c("25-29","30-34","35-39","40-44","45-49","50-54"),"HCW",NA)) %>% 
  drop_na() %>% 
  mutate(`Under20` = (`00-04`+`05-09`+`10-14`+`15-19`)/4,
         `20-29` = (`20-24`+`25-29`)/2,
         `30-39` =(`30-34`+`35-39`)/2,
         `40-49` = (`40-44`+`45-49`)/2,
         `50-59` = (`50-54`+`55-59`)/2,
         `60-69` = (`60-64`+`65-69`)/2,
         `Over70` = (`70-74`+`75+`)/2,
         `HCW`= (`25-29`+`30-34`+`35-39`+`40-44`+`45-49`+`50-54`)/6) %>% 
  ungroup() %>% 
  select(group_age,`Under20`,`20-29`,`30-39`,`40-49`,`50-59`,`60-69`,`Over70`,`HCW`) %>% 
  group_by(group_age) %>% 
  summarise_at(vars(-group_cols()),mean) %>% 
  ungroup() %>% 
  mutate(group_age = factor(group_age,levels = c("Under20","20-29","30-39","40-49","50-59","60-69","Over70","HCW"))) %>% 
  arrange(group_age) %>% 
  select(-group_age) 

waifw_norm = bind_rows(polymoduk, polymoduk_HCW)

if(do.write) write.table(waifw_norm,"./data/waifw_norm.csv",col.names = F,row.names = F,sep=",")

#home only (i.e non-work, non-school)
polymoduk=read_xlsx("./data/contact_matrices_152_countries/MUestimates_home_2.xlsx",
                    sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                          "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+")) %>%
  mutate(from = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"),
         group_age = ifelse(from %in% c("00-04","05-09","10-14","15-19"),"Under20",
                            ifelse(from %in% c("20-24","25-29"),"20-29",
                                   ifelse(from %in% c("30-34","35-39"),"30-39",
                                          ifelse(from %in% c("40-44","45-49"),"40-49",
                                                 ifelse(from %in% c("50-54","55-59"),"50-59",
                                                        ifelse(from %in% c("60-64","65-69"),"60-69",
                                                               ifelse(from %in% c("70-74","75+"),"Over70",NA)))))))) %>% 
  mutate(`Under20` = (`00-04`+`05-09`+`10-14`+`15-19`)/4,
         `20-29` = (`20-24`+`25-29`)/2,
         `30-39` =(`30-34`+`35-39`)/2,
         `40-49` = (`40-44`+`45-49`)/2,
         `50-59` = (`50-54`+`55-59`)/2,
         `60-69` = (`60-64`+`65-69`)/2,
         `Over70` = (`70-74`+`75+`)/2,
         `HCW`= (`25-29`+`30-34`+`35-39`+`40-44`+`45-49`+`50-54`)/6) %>% 
  ungroup() %>% 
  select(group_age,`Under20`,`20-29`,`30-39`,`40-49`,`50-59`,`60-69`,`Over70`,`HCW`) %>% 
  group_by(group_age) %>% 
  summarise_at(vars(-group_cols()),mean) %>% 
  ungroup() %>% 
  mutate(group_age = factor(group_age,levels = c("Under20","20-29","30-39","40-49","50-59","60-69","Over70","HCW"))) %>% 
  arrange(group_age) %>% 
  select(-group_age) 

polymoduk_HCW=read_xlsx("./data/contact_matrices_152_countries/MUestimates_home_2.xlsx",
                        sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                              "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+")) %>%
  mutate(from = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"),
         group_age = ifelse(from %in% c("25-29","30-34","35-39","40-44","45-49","50-54"),"HCW",NA)) %>% 
  drop_na() %>% 
  mutate(`Under20` = (`00-04`+`05-09`+`10-14`+`15-19`)/4,
         `20-29` = (`20-24`+`25-29`)/2,
         `30-39` =(`30-34`+`35-39`)/2,
         `40-49` = (`40-44`+`45-49`)/2,
         `50-59` = (`50-54`+`55-59`)/2,
         `60-69` = (`60-64`+`65-69`)/2,
         `Over70` = (`70-74`+`75+`)/2,
         `HCW`= (`25-29`+`30-34`+`35-39`+`40-44`+`45-49`+`50-54`)/6) %>% 
  ungroup() %>% 
  select(group_age,`Under20`,`20-29`,`30-39`,`40-49`,`50-59`,`60-69`,`Over70`,`HCW`) %>% 
  group_by(group_age) %>% 
  summarise_at(vars(-group_cols()),mean) %>% 
  ungroup() %>% 
  mutate(group_age = factor(group_age,levels = c("Under20","20-29","30-39","40-49","50-59","60-69","Over70","HCW"))) %>% 
  arrange(group_age) %>% 
  select(-group_age) 


waifw_home = bind_rows(polymoduk, polymoduk_HCW)
if(do.write) write.table(waifw_home,"./data/waifw_home.csv",col.names = F,row.names = F,sep=",")


#home + other location (i.e non-work, non-school)
polymoduk=read_xlsx("./data/contact_matrices_152_countries/MUestimates_home_2.xlsx",
                    sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                          "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+")) + 
  read_xlsx("./data/contact_matrices_152_countries/MUestimates_other_locations_2.xlsx",
            sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                  "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+")) 
polymoduk=polymoduk %>%
  mutate(from = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"),
         group_age = ifelse(from %in% c("00-04","05-09","10-14","15-19"),"Under20",
                            ifelse(from %in% c("20-24","25-29"),"20-29",
                                   ifelse(from %in% c("30-34","35-39"),"30-39",
                                          ifelse(from %in% c("40-44","45-49"),"40-49",
                                                 ifelse(from %in% c("50-54","55-59"),"50-59",
                                                        ifelse(from %in% c("60-64","65-69"),"60-69",
                                                               ifelse(from %in% c("70-74","75+"),"Over70",NA)))))))) %>% 
  mutate(`Under20` = (`00-04`+`05-09`+`10-14`+`15-19`)/4,
         `20-29` = (`20-24`+`25-29`)/2,
         `30-39` =(`30-34`+`35-39`)/2,
         `40-49` = (`40-44`+`45-49`)/2,
         `50-59` = (`50-54`+`55-59`)/2,
         `60-69` = (`60-64`+`65-69`)/2,
         `Over70` = (`70-74`+`75+`)/2,
         `HCW`= (`25-29`+`30-34`+`35-39`+`40-44`+`45-49`+`50-54`)/6) %>% 
  ungroup() %>% 
  select(group_age,`Under20`,`20-29`,`30-39`,`40-49`,`50-59`,`60-69`,`Over70`,`HCW`) %>% 
  group_by(group_age) %>% 
  summarise_at(vars(-group_cols()),mean) %>% 
  ungroup() %>% 
  mutate(group_age = factor(group_age,levels = c("Under20","20-29","30-39","40-49","50-59","60-69","Over70","HCW"))) %>% 
  arrange(group_age) %>% 
  select(-group_age) 

polymoduk_HCW=read_xlsx("./data/contact_matrices_152_countries/MUestimates_home_2.xlsx",
                        sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                              "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"))  + 
  read_xlsx("./data/contact_matrices_152_countries/MUestimates_other_locations_2.xlsx",
            sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                  "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+")) 
polymoduk_HCW=polymoduk_HCW %>%
  mutate(from = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"),
         group_age = ifelse(from %in% c("25-29","30-34","35-39","40-44","45-49","50-54"),"HCW",NA)) %>% 
  drop_na() %>% 
  mutate(`Under20` = (`00-04`+`05-09`+`10-14`+`15-19`)/4,
         `20-29` = (`20-24`+`25-29`)/2,
         `30-39` =(`30-34`+`35-39`)/2,
         `40-49` = (`40-44`+`45-49`)/2,
         `50-59` = (`50-54`+`55-59`)/2,
         `60-69` = (`60-64`+`65-69`)/2,
         `Over70` = (`70-74`+`75+`)/2,
         `HCW`= (`25-29`+`30-34`+`35-39`+`40-44`+`45-49`+`50-54`)/6) %>% 
  ungroup() %>% 
  select(group_age,`Under20`,`20-29`,`30-39`,`40-49`,`50-59`,`60-69`,`Over70`,`HCW`) %>% 
  group_by(group_age) %>% 
  summarise_at(vars(-group_cols()),mean) %>% 
  ungroup() %>% 
  mutate(group_age = factor(group_age,levels = c("Under20","20-29","30-39","40-49","50-59","60-69","Over70","HCW"))) %>% 
  arrange(group_age) %>% 
  select(-group_age) 


waifw_sdist = bind_rows(polymoduk, polymoduk_HCW)
if(do.write) write.table(waifw_sdist,"./data/waifw_sdist.csv",col.names = F,row.names = F,sep=",")


#======================
# Populations structure
#======================

#population structure per health board (2011), from https://www.scotlandscensus.gov.uk/ods-web/data-warehouse.html
scot_age <- read_csv(file ="./data/HealthBoard_census/DC1117SC.csv",
                     skip=5 , col_names = c("Health_Board","age","All people","Males","Females"), col_types = cols()) %>% 
  mutate(group_age = ifelse(age %in% c("Under 1",as.character(1:20)),"Under20",
                            ifelse(age %in% as.character(21:29),"20-29",
                                   ifelse(age %in% as.character(30:39),"30-39",
                                          ifelse(age %in% as.character(40:49),"40-49",
                                                 ifelse(age %in% as.character(50:59),"50-59",
                                                        ifelse(age %in% as.character(60:69),"60-69",
                                                               ifelse(age %in% c(as.character(70:84),"85 to 89","90 to 94","95 and over"),"Over70",NA))))))),
         group_age = factor(group_age,levels= names(waifw_norm)[-8]))

scot_age <- scot_age %>% 
  drop_na(group_age) %>% 
  ungroup() %>% 
  mutate(Health_Board = gsub("&","and",Health_Board)) %>% 
  group_by(Health_Board,group_age) %>% 
  summarise(tot=sum(`All people`)) %>% 
  mutate(freq=tot/sum(tot))

if(do.write) {
  write.table(scot_age %>% 
                select(-tot) %>% 
                spread(group_age,freq) %>% 
                # filter(Health_Board != "Scotland") %>% 
                ungroup() %>% 
                right_join(order_row %>% rename(Health_Board = `rowname`), by = "Health_Board") %>% 
                select(-Health_Board) ,
              "./data/scot_age.csv",col.names = F,row.names = F,sep=",")
}


#======================
# age-structured parameters
#======================
#fixed age-structured probabilities
cfr_byage=data.frame(
  #assume HCW similar to 20-59
  agegroup=colnames(waifw_norm),
  p_h=c(0.143,	0.1141,	0.117,	0.102,	0.125,	0.2,	0.303,	0.114525), #https://www.ecdc.europa.eu/sites/default/files/documents/RRA-seventh-update-Outbreak-of-coronavirus-disease-COVID-19.pdf
  cfr=c(0.001,	0.002,	0.002	,0.004,	0.013,	0.036,	0.114,	0.00525) #https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.12.2000256. NOTE: case defined as infected individuals (either symptomatic or asymptomatic)
) %>% 
  mutate(p_d=cfr/p_h) #convert p(d) to p(d|h)

if(do.write){
  write.table(cfr_byage %>% 
                select(-agegroup) %>% 
                ungroup() ,
              "./data/cfr_byage.csv",col.names = F,row.names = F,sep=",")
}



       