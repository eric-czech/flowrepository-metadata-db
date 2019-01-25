
#db <- DBI::dbConnect(RSQLite::SQLite(), path = "data/flowrepository-metadata.sqlite")

# Load human symbol table exported from HGNChelper package
prep_str <- function(x) x %>% str_trim %>% str_to_upper
df_hgnc <- read_csv('data/hgnc_table.csv', col_types=list()) %>%
  filter(!is.na(approved_symbol)) %>%
  filter(!is.na(symbol)) %>%
  mutate(symbol_lookup=prep_str(symbol), approved_lookup=prep_str(approved_symbol), id=1:nrow(.))

# Load FlowRepository meta tables
df_exp <- read_csv('data/experiments.csv', col_types=list())
df_kwd <- read_csv('data/keywords.csv', col_types=list())
df_fcs <- read_csv('data/fcs_files.csv', col_types=list())
df_chl <- read_csv('data/fcs_channels.csv', col_types=list())

# Build data frame to use for hgnc name joins
idx_hgnc <- bind_rows(
    # Add approved symbol -> id mapping first, followed by alias symbol -> id and then
    # remove any duplicate symbols/keys (which will favor approved names)
    df_hgnc %>% select(key=approved_lookup, value=approved_symbol),
    df_hgnc %>% select(key=symbol_lookup, value=approved_symbol)
  ) %>% filter(!duplicated(key)) 

df_chl_rs <- df_chl %>% mutate(lookup=prep_str(param_name)) %>% 
  left_join(idx_hgnc %>% rename(param_name_hgnc=value, lookup=key), by='lookup') %>%
  mutate(param=coalesce(param_name_hgnc, param_name, param_channel))

top_params <- df_chl_rs %>% 
  group_by(exp_id, param) %>% tally %>% ungroup %>% select(-n) %>%
  group_by(param) %>% tally %>% ungroup %>% arrange(desc(n)) %>% head(500)

df_chl_rs_pivot <- df_chl_rs %>% 
  group_by(exp_id, param) %>% tally %>% ungroup %>% filter(param %in% pull(top_params, param)) %>%
  spread(param, n)

m <- as.matrix(df_chl_rs_pivot %>% select(-exp_id))
m <- ifelse(!is.na(m) & (m > 0), 1, 0)
rownames(m) <- pull(df_chl_rs_pivot, exp_id)
#heatmap.2(m[1:100, 1:100], trace='none', cexCol=.2, cexRow=.2)
heatmap.2(m, trace='none', cexCol=.3, cexRow=.3)




# https://stackoverflow.com/questions/33777165/correlation-using-funs-in-dplyr
dft <- data.frame(g=c('g1', 'g1', 'g1', 'g2', 'g2'), c1=c(1, 1, 0, 0, 0), c2=c(1, 0, 1, 0, 1), c3=c(0, 0, 1, 0, 1))
dft %>% select(-g) %>% rownames_to_column %>% gather(channel, value, starts_with('c')) %>%
  full_join(., ., by=c('rowname')) %>% arrange(channel.x, channel.y)


