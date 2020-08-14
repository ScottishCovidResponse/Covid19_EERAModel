library(SCRCdataAPI)

data_root <- file.path("data", "example")
filename <- file.path(data_root, "data_for_scotland.h5")

files = list(
  age = file.path(data_root, "scot_age.csv"),
  data = file.path(data_root, "scot_data.csv"),
  deaths = file.path(data_root, "scot_deaths.csv")
)

age = read.csv(files$age, header=FALSE)
data = read.csv(files$data)
deaths = read.csv(files$deaths)

create_array(filename = filename,
              component = "age",
              array = as.matrix(age),
              dimension_names = list(rowvalue = 1:nrow(age), colvalue = 1:ncol(age))
              )
create_array(filename = filename,
             component = "data",
             array = as.matrix(data),
             dimension_names = list(rowvalue = 1:nrow(data), colvalue = colnames(data))
             )
create_array(filename = filename,
             component = "deaths",
             array = as.matrix(deaths),
             dimension_names = list(rowvalue = 1:nrow(deaths), colvalue = colnames(deaths))
             )