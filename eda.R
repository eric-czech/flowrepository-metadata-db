
library(FlowRepositoryR)
library(tidyverse)
library(dbplyr)


FCS_BYTE_OFFSET <<- 65535

get_metadata <- function(dataset_id){
  data_set <- suppressWarnings(flowRep.get(dataset_id))
  
  res <- list()
  res$exp_id <- data_set@id[[1]]
  res$exp_name <- data_set@name[[1]]
  
  res$investigators <- data_set@primary.investigator %>% str_trim %>% str_c(collapse='|')
  res$investigators <- if (length(res$investigators) > 0) res$investigators else NA_character_
  
  res$researchers <- data_set@primary.researcher %>% str_trim %>% str_c(collapse='|')
  res$researchers <- if (length(res$researchers) > 0) res$researchers else NA_character_
  
  res$attachments <- data_set@attachments %>% map(~.@name) %>% unlist %>% str_trim %>% str_c(collapse='|')
  res$attachments <- if (length(res$attachments) > 0) res$attachments else NA_character_
  
  res$keywords <- unlist(data_set@keywords)
  res$fcs <- map(data_set@fcs.files[1:3], function(fcs_file){
    print(fcs_file@url)
    url_con <- url(fcs_file@url)
    open(url_con, "rb")
    con <- rawConnection(readBin(url_con, 'raw', FCS_BYTE_OFFSET))
    
    header <- flowCore:::readFCSheader(con)
    txt <- tryCatch(
      flowCore:::readFCStext(con, header, emptyValue=T),
      error = function(e) flowCore:::readFCStext(con, header, emptyValue=F)
    )
    
    r <- list()
    r$exp_id <- data_set@id[[1]]
    r$version <- fcs_file@fcs.version
    r$filename <- fcs_file@name 
    r$size <- fcs_file@size
    r$creator <- flowCore:::readFCSgetPar(txt, 'CREATOR', strict=F) %>% unname
    r$n_params <- as.integer(flowCore:::readFCSgetPar(txt, "$PAR"))
    r$param_channels <- flowCore:::readFCSgetPar(txt, paste("$P", 1:r$n_params, "N", sep="")) %>% str_trim
    r$param_names <- flowCore:::readFCSgetPar(txt, paste("$P", 1:r$n_params, "S", sep=""), strict=F) %>% str_trim
    
    close(url_con)
    close(con)
    r
  })
  res
}

##############
# Processing #
##############

process_datasets <- function(dataset_ids, max_failures=2){
  data <- list()
  failures <- 0
  for (dataset_id in dataset_ids){
    print(str_glue('Processing dataset {dataset_id}'))
    data[[dataset_id]] <- tryCatch(
      get_metadata(dataset_id),
      error = function(e){
        failures <<- failures + 1
        if (failures > max_failures)
          stop(str_glue('Number of dataset processing failures ({failures}) exceeds maximum ({max_failures})'))
        warning(str_glue('A failure occurred processing dataset {dataset_id}. Error message: {msg}"', msg=e$message))
        list()
      }
    )
  }
  data
}

#data_sets <- flowRep.ls()
# data_set <- flowRep.get("FR-FCM-ZZJ7")
data <- process_datasets(data_sets[1:3])


##############
# Extraction #
##############

## Extract experiment metadata
cs_exp <- # Build vector of scalar fields specific to experiments
  data %>% map(names) %>% unlist %>% unname %>% unique %>% 
  discard(~. %in% c('fcs', 'keywords'))
df_exp <- # Extract the above fields for each experiment
  data %>% map(`[`, cs_exp) %>% bind_rows 

## Extract experiment keywords
df_kwd <- 
  data %>% map('keywords') %>% enframe %>% unnest %>% 
  set_names('exp_id', 'keyword')

## Extract FCS file metadata
cs_fcs <- # Build vector of scalar fields specific to fcs files (but not channels)
  data %>% map('fcs') %>% map(~map(., names)) %>% 
  unlist %>% unname %>% unique %>% 
  discard(~. %in% c('param_channels', 'param_names'))
df_fcs <- # Extract the fields above for each fcs file
  data %>% map('fcs') %>% unlist(recursive=F) %>% map(`[`, cs_fcs) %>% bind_rows

## Extract FCS file channels
df_chl <- # Build single data frame mapping experiments to fcs files and their corresponding channels
  data %>% map('fcs') %>% unlist(recursive=F) %>% 
  map(`[`, c('param_channels', 'param_names', 'filename', 'exp_id')) %>% 
  map(as_tibble) %>% bind_rows



##########
# Export #
##########

db <- DBI::dbConnect(RSQLite::SQLite(), path = "data/flowrepository-metadata.sqlite")
# src_dbi(db)
# DBI::dbDisconnect(db_con)

copy_to(db, df_exp, 'experiments', indexes = list('exp_id'), overwrite=TRUE, temporary=FALSE)
copy_to(db, df_kwd, 'keywords', indexes = list('exp_id'), overwrite=TRUE, temporary=FALSE)
copy_to(db, df_fcs, 'fcs_files', indexes = list('exp_id'), overwrite=TRUE, temporary=FALSE)
copy_to(
  db, df_chl, 'fcs_channels', 
  indexes = list('exp_id', 'filename', 'param_names'),
  overwrite=TRUE, temporary=FALSE
)


tbl_exp <- tbl(db, "experiments")

# con <- url('https://flowrepository.org/experiments/1602/fcs_files/170708/download')
# #con <- file('/tmp/header.txt')
# open(con, "rb")
# a <- readBin(con, 'raw', 65535)
# c2 <- rawConnection(a)
# header <- flowCore:::readFCSheader(c2)
# txt <- flowCore:::readFCStext(c2, header, emptyValue = T)
# 
# file_size <- as.integer(flowCore:::readFCSgetPar(txt, "$TOT"))
# n_params <- as.integer(flowCore:::readFCSgetPar(txt, "$PAR"))
# param_name <- flowCore:::readFCSgetPar(txt, paste("$P", 1:n_params, "N", sep=""))
# param_desc <- flowCore:::readFCSgetPar(txt, paste("$P", 1:n_params, "S", sep=""))
# params <- set_names(unname(param_desc), param_name)
# flowCore:::readFCSgetPar(txt, 'CREATOR')
# adata <- flowCore:::makeFCSparameters(unname(param_name), txt, F, F, 0, 0)
# adata@data


