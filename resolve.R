library(tidyverse)

# Parameter/channel name normalization function with anything parenthetical 
# stripped out, trimmed and converted to upper case
prep_str <- function(x) x %>% str_replace_all('[\\(\\[\\{].*[\\)\\]\\}]', '') %>% str_trim %>% str_to_upper

# Function used to resolve HGNC symbol groups to a single preferred name where
# the shortest CD* name is used if present
get_preferred_symbol <- function(g){
  # Get shortest name matching CD prefix
  r <- g %>% map_int(str_length) %>% base::order() %>% g[.] %>% .[str_detect(., '^CD')]
  if (length(r) > 0) r[[1]] else g[[1]]
}

# Load HGNC symbol table with appended manual additions that are specific to flow/cytof
df_hgnc <- read_csv('data/param_map_hgnc.csv', col_types=list()) %>%
  
  # Remove any symbols shorter than 3 characters as these are likely to map incorrectly
  # to parts of parameter names.  Some of these also seem to be in error such as 'CD':
  #   symbol approved_symbol preferred_symbol symbol_lookup approved_lookup    id
  # <chr>  <chr>           <chr>            <chr>         <chr>           <int>
  # 1 CD     CELIAC2         CD               CD            CELIAC2          7097
  # 2 CD     CTLA4           CD               CD            CTLA4            7098
  # 3 CD     NOD2            CD               CD            NOD2             7099
  filter(str_length(symbol) >= 3 & str_length(approved_symbol) >= 3) %>%
  
  # Merge with manually specified mappings (may be any length symbols)
  bind_rows(read_csv('data/param_map_edits.csv', col_types=list())) %>%

  # Group by approved symbol and resolved to a "preferred" symbol which currently means finding 
  # the CD molecule name for the group, if one exists (otherwise the original approved name is used)
  group_by(approved_symbol) %>% 
  mutate(preferred_symbol=get_preferred_symbol(c(approved_symbol[1], symbol))) %>%
  ungroup %>%
    
  # Remove any NA values (there is only one at TOW)
  filter(!is.na(approved_symbol) & !is.na(preferred_symbol) & !is.na(symbol)) %>%
  
  # Normalize symbols except for "preferred" since this should always be used as-is (for display)
  mutate(symbol_lookup=prep_str(symbol), approved_lookup=prep_str(approved_symbol), id=1:nrow(.))

stopifnot(all(!is.na(df_hgnc)))

# Load FlowRepository meta tables
df_chl <- read_csv('data/fcs_channels.csv', col_types=list())

# Build data frame to use for hgnc name joins
idx_hgnc <- bind_rows(
  # Add approved symbol -> id mapping first, followed by alias symbol -> id and then
  # remove any duplicate symbols/keys (which will favor approved names)
  df_hgnc %>% select(key=approved_lookup, value=preferred_symbol),
  df_hgnc %>% select(key=symbol_lookup, value=preferred_symbol)
) %>% filter(!duplicated(key)) 

df_chl_rs <- df_chl %>% mutate(
  # Attempt exact match to channel description field with highest precedence followed by
  # splits on common separating characters to match first and last strings
  lookup1=prep_str(param_name),
  lookup2=lookup1 %>% str_split("\\s*[ _/.]\\s*") %>% map_chr(~head(., 1)),
  lookup3=lookup1 %>% str_split("\\s*[ _/.]\\s*") %>% map_chr(~tail(., 1)),
  
  # Do the same for the channel name
  lookup4=prep_str(param_channel),
  lookup5=lookup4 %>% str_split("\\s*[ _/.]\\s*") %>% map_chr(~head(., 1)),
  lookup6=lookup4 %>% str_split("\\s*[ _/.]\\s*") %>% map_chr(~tail(., 1)),
  
  # Lastly, attempt to match on channel or description using a split on hyphens with
  # lowest precedence as hyphens are more likely to be valid parts of symbols
  lookup7=lookup1 %>% str_split("\\s*[-]\\s*") %>% map_chr(~head(., 1)),
  lookup8=lookup1 %>% str_split("\\s*[-]\\s*") %>% map_chr(~tail(., 1)),
  lookup9=lookup4 %>% str_split("\\s*[-]\\s*") %>% map_chr(~head(., 1)),
  lookup10=lookup4 %>% str_split("\\s*[-]\\s*") %>% map_chr(~tail(., 1))
) %>% 
  left_join(idx_hgnc %>% rename(value1=value, lookup1=key), by='lookup1') %>%
  left_join(idx_hgnc %>% rename(value2=value, lookup2=key), by='lookup2') %>%
  left_join(idx_hgnc %>% rename(value3=value, lookup3=key), by='lookup3') %>%
  left_join(idx_hgnc %>% rename(value4=value, lookup4=key), by='lookup4') %>%
  left_join(idx_hgnc %>% rename(value5=value, lookup5=key), by='lookup5') %>%
  left_join(idx_hgnc %>% rename(value6=value, lookup6=key), by='lookup6') %>%
  left_join(idx_hgnc %>% rename(value7=value, lookup7=key), by='lookup7') %>%
  left_join(idx_hgnc %>% rename(value8=value, lookup8=key), by='lookup8') %>%
  left_join(idx_hgnc %>% rename(value9=value, lookup9=key), by='lookup9') %>%
  left_join(idx_hgnc %>% rename(value10=value, lookup10=key), by='lookup10') %>%
  mutate(
    param=coalesce(value1, value2, value3, value4, value5, value6, value7, value8, value9, value10, param_name, param_channel),
    resolution=case_when(
      !is.na(value1) ~ 'lookup',
      !is.na(value2) ~ 'lookup',
      !is.na(value3) ~ 'lookup',
      !is.na(value4) ~ 'lookup',
      !is.na(value5) ~ 'lookup',
      !is.na(value6) ~ 'lookup',
      !is.na(value7) ~ 'lookup',
      !is.na(value8) ~ 'lookup',
      !is.na(value9) ~ 'lookup',
      !is.na(value10) ~ 'lookup',
      !is.na(param_name) ~ 'original',
      TRUE ~ 'channel'
    )
  )
stopifnot(nrow(df_chl_rs) == nrow(df_chl))

# Export channels data with resolved parameter names
df_chl_rs %>% select(-starts_with('value'), -starts_with('lookup')) %>%
  write_csv('data/fcs_channels_resolved.csv')


