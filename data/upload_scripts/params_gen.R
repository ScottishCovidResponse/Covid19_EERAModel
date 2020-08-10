##############################################################################
#                                                                            #
#                        EERA Model Parameter Upload                         #
#                                                                            #
#   This script adds new entries into the data registry for the EERA model   #
#   parameters and creates TOML files which should then be uploaded to the   #
#   relevant storage platform (ScottishCovidResponse/DataRepository on       #
#   GitHub at time of writing).                                              #
#                                                                            #
#   You will require an API token from the registry stored in:               #
    token_file <- file.path("data", "upload_scripts", "token.txt")           #
#                                                                            #
#   This script is based on the upload scripts by Sonia Mitchell.            #
#                                                                            #
#   @author K. Zarebski                                                      #
#   @date   last modified 2020-08-10                                         #
    global_version <- "0.1.0"                                                #
#                                                                            #
##############################################################################

library(SCRCdataAPI)    # Requires the latest SCRCdataAPI library

# Start constructing address locations
# two sets of data: 
#   - EERA/fixed-parameters
#   - EERA/prior-distributions

namespace <- "EERA"
prefix <- list(fixed="fixed-parameters", prior_dis="prior-distributions")

if(!file.exists(token_file))
{
    stop(paste("Failed to retrieve a data registry token file from", token_file))
}

key <- read.table(token_file)

# Assemble remote storage location objects
data_github_url <- "git@github.com:ScottishCovidResponse/DataRepository.git"
productStorageRoot <- "DataRepository"
storage_rootId <- new_storage_root(
  name = productStorageRoot,
  root = "https://raw.githubusercontent.com/ScottishCovidResponse/DataRepository/",
  key = key)

namespaceId <- new_namespace(name = namespace,
                             key = key)

# List of fixed parameters and their values, if a single parameter is updated
# without changing the others the version number should also be set
fixed <- list(
    list(name="K",         value=2000,   version=global_version),
    list(name="T_hos",     value=5,      version=global_version),
    list(name="T_inf",     value=1.5,    version=global_version),
    list(name="T_lat",     value=4,      version=global_version),
    list(name="T_rec",     value=11,     version=global_version),
    list(name="T_sym",     value=7,      version=global_version),
    list(name="day_shut",  value=19,     version=global_version),
    list(name="inf_asym",  value=1.0,    version=global_version),
    list(name="juvp_s",    value=0.1,    version=global_version),
    list(name="total_hcw", value=112974, version=global_version)
)

# List of prior distributions and their parameters, if a single distribution is
# updated without changing the others the version number should also be set
prior_dis <- list(
    list(name="rrd", type="gamma", params=list(k=1.0, theta=1.0), version=global_version),
    list(name="pinf", type="beta", params=list(alpha=3.0, beta=9.0), version=global_version),
    list(name="phcw", type="beta", params=list(alpha=3.0, beta=3.0), version=global_version),
    list(name="ps", type="beta", params=list(alpha=9.0, beta=3.0), version=global_version),
    list(name="q", type="beta", params=list(alpha=3.0, beta=3.0), version=global_version),
    list(name="lambda", type="uniform", params=list(a=1e-9, b=1e-6), version=global_version),
    list(name="d", type="beta", params=list(alpha=3.0, beta=3.0), version=global_version),
    list(name="chcw", type="poisson", params=list(lambda=42), version=global_version)
)

# Iterate through fixed parameters creating TOML objects and adding them
# as statements within the DataRegistry
for(param in fixed)
{
    name <- file.path(prefix$fixed, param$name)
    path <- paste("master", namespace, name, sep = "/")
    filename <- paste0(param$version, ".toml")
    component_name <- gsub("^.*/([^/]*)$", "\\1", name)

    create_estimate(filename = filename,
        path = file.path("data-raw", path),
        parameters = as.list(setNames(param$value, component_name)))
    upload_data_product(storage_root_id = storage_rootId,
                    name = name,
                    component_name = component_name,
                    processed_path = file.path("data-raw", path, filename),
                    product_path = file.path(path, filename),
                    version = param$version,
                    namespace_id = namespaceId,
                    key = key)
}

# Iterate through prior distibutions creating TOML objects and adding them
# as statements within the DataRegistry
for(dis in prior_dis)
{
    name <- file.path(prefix$prior_dis, dis$name)
    path <- paste("master", namespace, name, sep = "/")
    filename <- paste0(dis$version, ".toml")
    component_name <- gsub("^.*/([^/]*)$", "\\1", name)

    create_distribution(
        filename = filename,
        file.path("data-raw", path),
        name = dis$name,
        distribution = dis$type,
        parameters = dis$params
    )
    
    upload_data_product(storage_root_id = storage_rootId,
                    name = name,
                    component_name = component_name,
                    processed_path = file.path("data-raw", path, filename),
                    product_path = file.path(path, filename),
                    version = dis$version,
                    namespace_id = namespaceId,
                    key = key)

}