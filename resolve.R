#' Script intended for normalization of FCS parameter names.  
#' 
#' This process occurs in several steps:
#' 1. A table of HGNC symbols to names is loaded (expored from the HGNChelper package)
#' 2. Common symbols related to flow/cytof like DNA,FSC,SSC,etc are removed as well as any symbols less than 3 chars
#' 3. This table is appended to a list of manual mappings (see data/param_map_edits.csv)
#' 4. Parameter and channel names in flow metadata are loaded and split based on common naming schemes
#' 5. All split symbols (and originals) are matched against the lookup table
#' 6. Results conatining the raw channel/parameter names are written again but this time they also
#'    include the most likely target for the channel as well as what part of either the channel/parameter
#'    name was used to make the match and a "precedence" that indicates which parsing function resulted in the hit (if any)
library(tidyverse)

# Parameter/channel name normalization function with anything parenthetical 
# stripped out, trimmed and converted to upper case
#prep_str <- function(x) x %>% str_replace_all('[\\(\\[\\{].*[\\)\\]\\}]', '') %>% str_trim %>% str_to_upper
prep_str <- function(x) x %>% str_trim %>% str_to_upper

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
  
  # Remove gene symbols likely to counfound with flow specific param names
  filter(symbol %>%          str_detect('^SSC|^FSC|^CELL$|^RED|^APC|^PE$|^CYC$|^COMP$') %>% `!`) %>%
  filter(approved_symbol %>% str_detect('^SSC|^FSC|^CELL$|^RED|^APC|^PE$|^CYC$|^COMP$') %>% `!`) %>%
  
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


# Candidate symbol parsing functions, specified in order of precedence
parser_fns <- list()
# Return symbol for matching as-is (with trim)
parser_fns$f1=function(x) x %>% str_trim
# Split on any common top level delimiters
parser_fns$f2=function(x) parser_fns$f1(x) %>% str_split("\\s*[ _/(){}\\[\\]#]\\s*")
# Split on secondary delimiters more likely to be valid parts of names
parser_fns$f3=function(x) parser_fns$f2(x) %>% map(~str_split(., "\\s*[.-]\\s*") %>% unlist)

# Break symbols into individual candidates with associated precedence
get_candidates <- function(symbols){
  df <- list()
  ct <- 1
  nc <- length(symbols) * length(parser_fns)
  for (i in 1:length(symbols)){
    for (j in 1:length(parser_fns)){
      cands <- parser_fns[[j]](symbols[[i]]) %>% unlist %>% str_trim %>% str_to_upper
      # Assign precedence such that higher values are preferred
      df[[ct]] <- tibble(key=cands, precedence=nc-ct+1)
      ct <- ct + 1
    }
  }
  bind_rows(df) %>% 
    filter(!duplicated(key) & !is.na(key)) 
}

df_lkp <- df_chl %>% group_by(param_name, param_channel) %>% do({
    d <- .
    symbols <- c(d$param_name[[1]], d$param_channel[[1]])
    get_candidates(symbols) %>% mutate(param_name=symbols[[1]], param_channel=symbols[[2]])
  }) %>% ungroup %>%
  # Left join to lookup table so that all param name+channel combinations are preserved
  left_join(idx_hgnc, by='key') %>% 
  # Take first row with a match and highest precedence (or arbitrary row if there are no matches)
  group_by(param_name, param_channel) %>% slice(which.max(precedence * (!is.na(value)))) %>% ungroup %>%
  # Assign fields such that the term/symbol that led to the match is known as well as whether
  # or not a match was found or if the ultimate "param" field assigned is the original or channel name
  rename(term=key, param=value) %>%
  mutate(
    resolution=case_when(!is.na(param) ~ 'lookup', !is.na(param_name) ~ 'original', TRUE ~ 'channel'),
    term=case_when(!is.na(param) ~ term, TRUE ~ NA_character_)
  ) %>%
  mutate(param=coalesce(param, param_name, param_channel))

# df_lkp %>% group_by(resolution) %>% tally # At TOW:
# A tibble: 3 x 2
# resolution     n
# <chr>      <int>
# 1 channel     1007
# 2 lookup      4619
# 3 original    2397

# Ensure there is now one record per param name + channel combo
stopifnot(df_lkp %>% group_by(param_name, param_channel) %>% tally %>% pull(n) %>% max == 1)

# Join the resolved names to the original channels data frame (and ensure all records match)
df_chl_rs <- df_chl %>% inner_join(df_lkp, by=c('param_name', 'param_channel'), na_matches='na')
stopifnot(nrow(df_chl_rs) == nrow(df_chl))

# Export channels data with resolved parameter names
# > df_chl_rs
# A tibble: 100 x 8
# param_channel param_name   filename                                 exp_id      term  precedence param resolution
# <chr>         <chr>        <chr>                                    <chr>       <chr>      <dbl> <chr> <chr>     
# 1 SSC-A         NA           TS18 447 P3…Live.fcs                     FR-FCM-ZZFW NA             3 SSC-A channel   
# 2 SSC-A         NA           c5a.fcs                                  FR-FCM-ZZML NA             3 SSC-A channel   
# 3 (Er170)Di     CD3_Er170    IM_Exp6_14072014_cells_found.fcs         FR-FCM-ZZWC CD3            5 CD3   lookup    
# 4 FL1-H         CD63 FITC    963F12.008                               FR-FCM-ZZ93 CD63           5 CD63  lookup    
# 5 (Sm148)Di     CD16-148 (v) MP-062_1263-13_Week 4_viSNE_melanoma.fcs FR-FCM-ZYW8 CD16           4 CD16  lookup    
# 6 FITC-A        NKp46        9607_7_3_NKR.fcs                         FR-FCM-ZZ5N NKP46          6 CD335 lookup    
# 7 SSC-A         SSC-A        F08 85 t3 5000x.fcs                      FR-FCM-ZYBR NA             6 SSC-A original  
# 8 FSC-A         NA           Specimen_001_C6_C06_YER039C-A_0-2.fcs    FR-FCM-ZZEZ NA             3 FSC-A channel
df_chl_rs %>% select(-starts_with('value'), -starts_with('lookup')) %>%
  write_csv('data/fcs_channels_resolved.csv')


