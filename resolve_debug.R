#' Commands useful for debugging individual parameters following execution of the resolve.R script

# Estimate rate of matching
df_chl_rs %>% 
  group_by(exp_id, param, resolution) %>% tally %>% select(-n) %>% ungroup %>%
  group_by(resolution) %>% tally 
# # A tibble: 3 x 2
# resolution     n
# <chr>      <int>
# 1 channel     4178
# 2 lookup      8901
# 3 original    3550

# Find most common unmatched params
df_chl_rs %>% 
  filter(resolution == 'original') %>% 
  #filter(resolution == 'channel') %>% 
  group_by(exp_id, param) %>% tally %>% select(-n) %>% ungroup %>%
  group_by(param) %>% tally %>% 
  arrange(desc(n)) %>% print(n=100)

###### Search for individual param
#par_query <- 'CTLA[-]*4'
#par_query <- 'IFN'
#par_query <- 'GRANZ'
#par_query <- 'TCR'
#par_query <- 'TNF'
#par_query <- 'IL[-]*[\\d]+'
par_query <- 'FOXP'

# Search for param (w/ experiment count) and kind of resolution in matched records
df_chl_rs %>% filter(param %>% str_to_upper %>% str_detect(par_query)) %>% 
  group_by(exp_id, param, resolution) %>% tally %>% ungroup %>% select(-n) %>%
  group_by(param) %>% summarize(resolution=resolution[1], n=n()) %>% ungroup %>%
  arrange(desc(n)) %>% print(n=50)

# Look for any match in lookup table
df_hgnc_raw <- read_csv('data/param_map_hgnc.csv', col_types=list())
# df_hgnc_raw <- df_hgnc
df_hgnc_raw %>% gather %>% filter(value %>% str_detect(paste0('^', par_query, '$'))) # Exact
df_hgnc_raw %>% gather %>% filter(value %>% str_detect(par_query)) %>% # Approximate
  arrange(str_length(value)) %>% print(n=50)
df_hgnc_raw %>% filter(symbol %>% str_detect(paste0('^', par_query, '$'))) %>% print(n=50)
df_hgnc_raw %>% filter(approved_symbol %>% str_detect(paste0('^', par_query, '$'))) %>% print(n=50)

# Query raw channel/parameter names
df_chl_raw <- read_csv('data/fcs_channels.csv', col_types=list())
df_chl_raw %>% select(param_channel, param_name) %>% gather %>% mutate(value=str_to_upper(value)) %>%
  filter(value %>% str_detect(par_query)) %>% group_by(value) %>% tally %>% arrange(desc(n)) %>% print(n=100)

