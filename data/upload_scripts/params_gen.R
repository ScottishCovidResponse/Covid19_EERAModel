library(SCRCdataAPI)

namespace <- "EERA"
prefix <- list(fixed="fixed-parameters", prior_dis="prior-distributions")
global_version <- "0.1.0"

key <- read.table("token.txt")
data_github_url <- "git@github.com:ScottishCovidResponse/DataRepository.git"
productStorageRoot <- "DataRepository"
storage_rootId <- new_storage_root(
  name = productStorageRoot,
  root = "https://raw.githubusercontent.com/ScottishCovidResponse/DataRepository/",
  key = key)
namespaceId <- new_namespace(name = namespace,
                             key = key)

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

prior_dis <- list(
    list(name="rrd", type="gamma", params=list(shape=1.0, scale=1.0), version=global_version)
)

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