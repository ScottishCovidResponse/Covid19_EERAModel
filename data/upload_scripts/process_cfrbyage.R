##############################################################################
#                                                                            #
#                           POSTERIOR PARAMETERS                             #
#                                                                            #
#   The following script is used to generate a table to represent the        #
#   posterior parameters as read from the CSV file.                          #
#                                                                            #
#   @author : K. Zarebski                                                    #
#   @date   : last modified 2020-08-13                                       #
#                                                                            #
##############################################################################

library(SCRCdataAPI)
library(data.table)
library(magrittr)

date_accessed <- Sys.Date()
struct_version <- 0
dataset_version <- 0

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
    agegroup=c("Under20","20-29","30-39","40-49","50-59","60-69","Over70","HCW"),
    p_h=p_h,
    cfr=cfr
  ) %>% 
    dplyr::mutate(p_d=cfr/p_h) #convert p(d) to p(d|h)
  
  cfr_byage = cfr_byage %>% dplyr::select(-agegroup) %>% ungroup()
  
  return(cfr_byage)
  
}

p_h=c(0.143,	0.1141,	0.117,	0.102,	0.125,	0.2,	0.303,	0.114525) #https://www.ecdc.europa.eu/sites/default/files/documents/RRA-seventh-update-Outbreak-of-coronavirus-disease-COVID-19.pdf
cfr=c(0.001,	0.002,	0.002	,0.004,	0.013,	0.036,	0.114,	0.00525) #https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.12.2000256. NOTE: case defined as infected individuals (either symptomatic or asymptomatic)

input_table = data.table(create_age_epiparaam(p_h, cfr))

token_file <- file.path("data", "upload_scripts", "token.txt")
if (!file.exists(token_file)) {
  stop(paste("Failed to find file API token file at", token_file))
}
key <- read.table(token_file)

tmp <- as.Date(date_accessed, format = "%Y-%m-%d")

version_number <- paste(struct_version, gsub("-", "", tmp), dataset_version, sep = ".")
namespace <- "EERA"
product_name <- file.path("prob_hosp_and_cfr", "data_for_scotland")

# where is the data product saved? (locally, before being stored)
product_filename <- paste(version_number, "h5", sep = ".")

create_table(filename = product_filename,
             path = product_name,
             component = "cfr_byage",
             df = input_table)

# where is the data product stored?
product_storageRoot <- "boydorr"
product_path <- product_name

# data product storage root
product_storageRootId <- new_storage_root(name = product_storageRoot,
                                          root = "ftp://boydorr.gla.ac.uk/scrc/",
                                          key = key)
# namespace
namespaceId <- new_namespace(name = namespace,
                             key = key)


dataProductURIs <- upload_data_product(
  storage_root_id = product_storageRootId,
  name = product_name,
  processed_path = file.path(product_name, product_filename),
  product_path = paste(namespace, product_name, product_filename, sep = "/"),
  version = version_number,
  namespace_id = namespaceId,
  key = key)