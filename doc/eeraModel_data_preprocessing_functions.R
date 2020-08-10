#======================
# Settings
#======================
do.write = F

# Check that the current directory is repository

if(!file.exists(file.path(getwd(), "data", "contact_matrices_152_countries")))
{
  stop("This script must be run from the root of the EERA Model repository.")
}

# ----------------------------------
# FUNCTIONS
# ----------------------------------

tidy_case_data <- function(case_df){
  require(dplyr)
  scot_data <- case_df %>%
    rename(`Scotland`= `Grand Total`) %>%
    mutate(Date = lubridate::dmy(Date))
  
  return(scot_data)
  
}

tidy_death_data <- function(death_df){
  require(dplyr)
  
  death_df[death_df == "x"] <- -999   
  scot_deaths <- death_df %>% 
    rename(`Scotland`= `Grand Total`) %>%
    mutate(Date = lubridate::dmy(Date))
  
  return(scot_deaths)
  
}



tidy_pop_data <- function(pop_df){
  require(dplyr)
  
  pop_df <- pop_df %>% 
    mutate(prop = Population/sum(Population))
  
  pop_df_tot = as.data.frame(list(`Name` = "Scotland",`Population` = sum(pop_df$Population)))%>% 
    mutate(prop = Population/sum(Population)) %>% 
    mutate_if(is.factor,as.character)
  scot_pop <- bind_rows(pop_df,pop_df_tot)
  
  return(scot_pop)
  
}



# function merging epidemiological data reported before and after changes in recording practices
# base_data: pre change data
# extra_data: prost change data
# extra_data: name of the column of interest

extend_epi_records <- function(base_data, extra_data, data_name){

  # Scottish overall cases
  scot_cases_tot <- extra_data %>%
    mutate(Date = lubridate::ymd(Date))
  
  scot_data <- 
    base_data %>% 
    select(-Scotland,) %>% 
    full_join(scot_cases_tot %>% 
                select(Date,!!as.symbol(data_name)) %>% 
                rename(Scotland = !!as.symbol(data_name)), 
              by = "Date") %>% 
    mutate(Scotland = replace_na(Scotland,0)) %>% 
    replace(is.na(.), -999)
  
  return(scot_data)
  
}

# function creating the input data giving information of number of cases/deaths per day (in column) and for each 
# health boards (rows). the last row (15) is for the whole Scotland. 
# The first column is the population at risk.
# disease_df: data.frame / tibble providing the epidemiological data (cases / deaths) use in the model or all areas of interest (scottish health boards and Scotland)
# pop_df: data.frame / tibble providing the population size for all areas of interest (scottish health boards and Scotland)
# upto: date limit of interest (a date in the format "yyyy-mm-dd"). if NULL, upto = max(disease_df$Date)
# join: TRUE/FALSE. boolean to allow (or not) the joining of data.frame / tibble to ensure consistencies between data frame 
# data_join: disease data to match
create_history <- function(disease_df,pop_df, upto=NULL, join=F,data_join=NULL){
  require(dplyr)
  require(tidyr)
  require(lubridate)
  
  if(join & !is.null(data_join)){
    temp<- left_join(data_join %>% select(Date),disease_df,by="Date") %>% 
      replace(is.na(.),0)
    
  } else{
    temp = disease_df
    if(join & is.null(data_join)) message("Warning!! Base data not provided. Joining not implemented.")
  }
  
  #define record date limit
  if(is.null(upto)) upto = max(temp$Date)
  
  
  #formating  the case data (HPS) for model input
  rval_temp <- right_join(pop_df %>% 
                       select(Name,Population),
                     as.data.frame(
                       t(
                         as.matrix(
                           temp %>% 
                            filter(Date <= upto) %>% 
                            select(-Date)
                          )
                         )
                       ) %>% 
                        tibble::rownames_to_column(), 
                      by = c("Name"="rowname")) %>% 
                      select(-Name)
  colnames(rval_temp) <- (1:ncol(rval_temp))-2

  order_row= as.data.frame(
    t(
      as.matrix(
        temp %>% 
          filter(Date <= upto) %>% 
          select(-Date)
      )
    )
  ) %>% 
    tibble::rownames_to_column() %>% 
    select(rowname)
  
  df_return = list(
    data = rval_temp,
    name = order_row
  )
 return(df_return) 
  
}

#
#
#
#

# function formating the mixing matrices
# maxmatrix: age mixing matrix in its original format
create_mixmat <- function(mixmatrix){
  
  require(dplyr)
  
  polymoduk= mixmatrix %>%
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
  
  
  polymoduk_HCW= mixmatrix %>%
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
  
  waifw = bind_rows(polymoduk, polymoduk_HCW)
  
  return(waifw)
  
}

#
#
#
#
#

#population structure per health board (2011), from https://www.scotlandscensus.gov.uk/ods-web/data-warehouse.html
# census_df: census data frame giving number of individauls per age
# order_name: list of all administrative areas in order of scot_data
#
create_age_pop <- function(census_df, order_name){
  
  require(dplyr)


  scot_age <- census_df %>% 
    mutate(group_age = ifelse(age %in% c("Under 1",as.character(1:20)),"Under20",
                              ifelse(age %in% as.character(21:29),"20-29",
                                     ifelse(age %in% as.character(30:39),"30-39",
                                            ifelse(age %in% as.character(40:49),"40-49",
                                                   ifelse(age %in% as.character(50:59),"50-59",
                                                          ifelse(age %in% as.character(60:69),"60-69",
                                                                 ifelse(age %in% c(as.character(70:84),"85 to 89","90 to 94","95 and over"),"Over70",NA))))))),
           group_age = factor(group_age,levels= c("Under20", "20-29","30-39","40-49","50-59","60-69","Over70")))
  
  scot_age <- scot_age %>% 
    drop_na(group_age) %>% 
    ungroup() %>% 
    mutate(Health_Board = gsub("&","and",Health_Board)) %>% 
    group_by(Health_Board,group_age) %>% 
    summarise(tot=sum(`All people`)) %>% 
    mutate(freq=tot/sum(tot))


  scot_age = scot_age %>% 
    select(-tot) %>% 
    spread(group_age,freq) %>% 
    # filter(Health_Board != "Scotland") %>% 
    ungroup() %>% 
    right_join(order_name %>% rename(Health_Board = `rowname`), by = "Health_Board") %>% 
    select(-Health_Board) 
  
  return(scot_age)
}



#
# function formating the data with age-structured parameters
# p_h: probability of severe cases requiring hospitalisation
# cfr: case fatality ratio
#
create_age_epiparaam <- function(p_h, cfr){
  
  require(dplyr)
  
  #fixed age-structured probabilities
  cfr_byage=data.frame(
    #assume HCW similar to 20-59
    agegroup=colnames(waifw_norm),
    p_h=p_h,
    cfr=cfr
  ) %>% 
    mutate(p_d=cfr/p_h) #convert p(d) to p(d|h)
  
  cfr_byage = cfr_byage %>% select(-agegroup) %>% ungroup()
  
  return(cfr_byage)
  
}





# ----------------------------------
# create data use for model
# ----------------------------------

#create the case and death data

# cases pre change in recording practices
WATTY62DATA <- "https://raw.githubusercontent.com/watty62/Scot_covid19/master/data/processed/regional_cases.csv"
WATTY62POP <- "https://raw.githubusercontent.com/watty62/Scot_covid19/master/data/processed/HB_Populations.csv"
WATTY62DEATHS <- "https://raw.githubusercontent.com/watty62/Scot_covid19/master/data/processed/regional_deaths.csv"

pop_df <- read_csv(file =WATTY62POP , col_types = cols())
case_df = read_csv(file = WATTY62DATA, col_types = cols())
death_df <- read_csv(file = WATTY62DEATHS, col_types = cols()) 

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


upto = "2020-06-09"


extend_case = extend_epi_records(base_data = tidy_case_data(case_df), extra_data = df_cases_tot, data_name='confirmed_cases')
extend_death = extend_epi_records(base_data = tidy_case_data(death_df), extra_data = df_deaths_tot, data_name='deaths')


scot_data = create_history(extend_case,tidy_pop_data(pop_df), upto = upto)$data
scot_death = create_history(extend_death,tidy_pop_data(pop_df), upto = upto, 
                            join = T, data_join = extend_case)$data


if(do.write) write.table(scot_data,paste0("./data/scot_data_",upto,".csv"),col.names = T,row.names = F,sep=",")
if(do.write) write.table(scot_death,paste0("./data/scot_deaths_",upto,".csv"),col.names = T,row.names = F,sep=",")

# create the age mixing matrices

matrix_data_dir <- file.path(getwd(), "data", "contact_matrices_152_countries")

polymoduk_norm=read_xlsx(file.path(matrix_data_dir, "MUestimates_all_locations_2.xlsx"),
                    sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                          "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"))
  
polymoduk_home=read_xlsx(file.path(matrix_data_dir, "MUestimates_home_2.xlsx"),
                    sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                          "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"))

polymoduk_other=read_xlsx(file.path(matrix_data_dir, "MUestimates_other_locations_2.xlsx"),
                    sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                          "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"))

polymoduk_work=read_xlsx(file.path(matrix_data_dir, "MUestimates_work_2.xlsx"),
                          sheet="United Kingdom of Great Britain",col_names = c("00-04","05-09","10-14","15-19","20-24","25-29","30-34","35-39",
                                                                                "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"))

polymoduk_school=read_xlsx(file.path(matrix_data_dir, "MUestimates_school_2.xlsx"),
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

census_age = read_csv(file =file.path(getwd(), "data", "HealthBoard_census", "DC1117SC.csv"),
         skip=5 , col_names = c("Health_Board","age","All people","Males","Females"), col_types = cols())

order_name = create_history(tidy_case_data(case_df),tidy_pop_data(pop_df), upto = upto)$name

scot_age <- create_age_pop(census_df = census_age, order_name)

if(do.write) write.table(scot_age, file.path(getwd(), "data", "scot_age.csv"),col.names = F,row.names = F,sep=",")

# create the data with age-structured parameters

p_h=c(0.143,	0.1141,	0.117,	0.102,	0.125,	0.2,	0.303,	0.114525) #https://www.ecdc.europa.eu/sites/default/files/documents/RRA-seventh-update-Outbreak-of-coronavirus-disease-COVID-19.pdf
cfr=c(0.001,	0.002,	0.002	,0.004,	0.013,	0.036,	0.114,	0.00525) #https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.12.2000256. NOTE: case defined as infected individuals (either symptomatic or asymptomatic)

if(do.write) write.table(create_age_epiparaam(p_h, cfr) , file.path(getwd(), "data", "cfr_byage.csv"), col.names = F,row.names = F,sep=",")




