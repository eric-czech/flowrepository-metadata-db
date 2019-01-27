library(FlowRepositoryR)
library(tidyverse)
library(dbplyr)
library(logging)
basicConfig()


FCS_BYTE_OFFSET <<- 32768
FCS_MAX_FILES <<- 30

get_metadata <- function(dataset_id){
  data_set <- suppressWarnings(flowRep.get(dataset_id))
  
  res <- list()
  res$exp_id <- data_set@id[[1]]
  res$exp_name <- data_set@name[[1]]
  
  # Function to pull and concatenate string attributes and ensure
  # they are not empty vectors
  get_str <- function(prop) {
    v <- attr(data_set, prop) %>% str_trim %>% str_c(collapse='|')
    if (length(v) > 0) v else NA_character_
  }
  
  # Extract information about researchers and related attachments, all
  # of which could be empty, one, or multiple values
  res$investigators <- get_str('primary.investigator')
  res$researchers <- get_str('primary.researcher')
  res$dates <- get_str('experiment.dates')
  res$keywords <- unlist(data_set@keywords)
  # Avoid failing whole experiment parsing due to any potential attachment processing errors
  tryCatch({
      res$attachments <- data_set@attachments %>% map_chr(~.@name) %>% str_trim %>% str_c(collapse='|')
      res$attachments <- if (length(res$attachments) > 0) res$attachments else NA_character_
    }, error = function(e) {
      res$attachments <- '[Parser Error]'
  })
  
  
  
  # Limit the number of fcs files to process randomly, if necessary, 
  # and record the count before and after applying this restriction
  fcs_files <- data_set@fcs.files
  if (length(fcs_files) > FCS_MAX_FILES) 
    fcs_files <- fcs_files %>% base::sample(FCS_MAX_FILES)
  res$n_fcs_files_total <- length(data_set@fcs.files)
  res$n_fcs_files_parse <- length(fcs_files)
  
  # Process each fcs file entry separately where the focus is on extracting
  # channel names (but other file-global data is also collected)
  res$fcs <- map(fcs_files, function(fcs_file){
    
    # Read the first FCS_BYTE_OFFSET bytes into an array and then pass this array
    # to flowCore internal functions as if it was a standard file connection
    url_con <- url(fcs_file@url)
    open(url_con, "rb")
    con <- rawConnection(readBin(url_con, 'raw', FCS_BYTE_OFFSET))
    
    # Parse header and text segments (note that emptyValue corresponds to how
    # ambiguous double delimiters are treated and per the flowCore recommendations,
    # both are tried in the event of a failure)
    header <- flowCore:::readFCSheader(con)
    txt <- tryCatch(
      flowCore:::readFCStext(con, header, emptyValue=T),
      error = function(e) flowCore:::readFCStext(con, header, emptyValue=F)
    )
    
    # Extract file-level metadata
    r <- list()
    r$exp_id <- data_set@id[[1]]
    r$version <- fcs_file@fcs.version
    r$filename <- fcs_file@name 
    r$size <- fcs_file@size
    r$creator <- flowCore:::readFCSgetPar(txt, 'CREATOR', strict=F) %>% unname
    r$n_params <- as.integer(flowCore:::readFCSgetPar(txt, "$PAR"))
    
    # Extract non-scalar per-channel data
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

process_datasets <- function(dataset_ids, max_failures=50){
  data <- list()
  failures <- 0
  n <- length(dataset_ids)
  for (i in 1:n) {
    dataset_id <- dataset_ids[i]
    loginfo(str_glue('Processing dataset {dataset_id} ({i} of {n} | #failures = {failures})'))
    data[[dataset_id]] <- tryCatch(
      get_metadata(dataset_id),
      error = function(e){
        failures <<- failures + 1
        if (failures > max_failures)
          stop(str_glue('Number of dataset processing failures ({failures}) exceeds maximum ({max_failures})'))
        tryCatch({
          # logwarn (or loginfo) fail if the message string is too long
          msg <- str_sub(e$message, end=1000)
          logwarn(str_glue('A failure occurred processing dataset {dataset_id}. Error message: {msg}"'))
        }, error=function(e) {
          logwarn(str_glue('A failure occurred processing dataset {dataset_id}"'))
        })
        list()
      }
    )
  }
  loginfo(str_glue('Processing complete.  Number of failed datasets = {failures} (of {n})', n=length(dataset_ids)))
  data
}

# data_set <- flowRep.get("FR-FCM-ZZJ7")
data_sets <- flowRep.ls()
data_raw <- process_datasets(data_sets)
data_raw <- data

# Save parsed results to file immediately in case of downstream crash
saveRDS(data_raw, file = "data/build-data.rds")

# Remove list entries for failed datasets (returned as empty lists)
failed_data_set_ids <- data_raw %>% map(length) %>% unlist %>% keep(~. == 0) %>% names
loginfo(str_glue(
  'Removing the following {n} failed data set ids from results: {ids}', 
  n=length(failed_data_set_ids), ids=str_c(failed_data_set_ids, collapse=', ')
))
data <- names(data) %in% failed_data_set_ids %>% `!` %>% data_raw[.]
stopifnot(length(data) == length(data_raw) - length(failed_data_set_ids))

# d <- process_datasets('FR-FCM-ZYRT')

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
  data %>% map('keywords') %>% enframe %>% filter(!map_lgl(value, is.null)) %>% unnest %>% 
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
  map(as_tibble) %>% bind_rows %>%
  rename(param_name=param_names, param_channel=param_channels)

loginfo('Extraction complete')

##########
# Export #
##########

write_csv(df_exp, 'data/experiments.csv')
write_csv(df_kwd, 'data/keywords.csv')
write_csv(df_fcs, 'data/fcs_files.csv')
write_csv(df_chl, 'data/fcs_channels.csv')

# db <- DBI::dbConnect(RSQLite::SQLite(), path = "data/flowrepository-metadata.sqlite")
# copy_to(db, df_exp, 'experiments', indexes = list('exp_id'), overwrite=TRUE, temporary=FALSE)
# copy_to(db, df_kwd, 'keywords', indexes = list('exp_id'), overwrite=TRUE, temporary=FALSE)
# copy_to(db, df_fcs, 'fcs_files', indexes = list('exp_id'), overwrite=TRUE, temporary=FALSE)
# copy_to(
#   db, df_chl, 'fcs_channels', 
#   indexes = list('exp_id', 'filename', 'param_names'),
#   overwrite=TRUE, temporary=FALSE
# )

loginfo('Export complete')

