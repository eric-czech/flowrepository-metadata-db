
library(tidyverse)

# Load FlowRepository meta tables
df_exp <- read_csv('data/experiments.csv', col_types=list())
df_kwd <- read_csv('data/keywords.csv', col_types=list())
df_fcs <- read_csv('data/fcs_files.csv', col_types=list())

valid_exp_ids <- df_exp %>% 
  # Remove duplicated "Twins TS" experiments
  filter((exp_name %>% str_detect('Twins TS') %>% `!`) | (exp_name == 'Twins TS01')) %>% 
  pull(exp_id)
df_chl <- read_csv('data/fcs_channels_resolved.csv', col_types=list()) %>%
  filter(exp_id %in% valid_exp_ids)


top_params <- df_chl %>% 
  filter(resolution == 'lookup') %>% 
  group_by(exp_id, param) %>% tally %>% ungroup %>% select(-n) %>%
  group_by(param) %>% tally %>% ungroup %>% arrange(desc(n)) %>% head(400)

# top_params %>% print(n=100)

df_chl_pivot <- df_chl %>% 
  group_by(exp_id, param) %>% tally %>% ungroup %>% select(-n) %>%
  filter(param %in% pull(top_params, param)) %>%
  group_by(exp_id) %>% mutate(n_params=length(param)) %>% ungroup %>%
  filter(n_params > 10) %>%
  spread(param, n_params)
# dim(df_chl_pivot)

m <- as.matrix(df_chl_pivot %>% select(-exp_id))
m <- ifelse(!is.na(m) & (m > 0), 1, 0)
rownames(m) <- pull(df_chl_pivot, exp_id)
#heatmap.2(m[1:100, 1:100], trace='none', cexCol=.2, cexRow=.2)
#gplots::heatmap.2(m, trace='none', cexCol=.3, cexRow=.3)


library(plotly)
y <- hclust(dist(m, method='binary'))
x <- hclust(dist(t(m), method='binary'))
mo <- m[order.dendrogram(as.dendrogram(y)), order.dendrogram(as.dendrogram(x))]
plot_ly(z=mo, y=rownames(mo), x=colnames(mo), type='heatmap') %>% plotly::layout(margin=list(l=100,b=75))
#plot_ly(z=m[y$order, x$order], type='heatmap')


#colnames(mo) %>% (function(v) which(v == 'CD335'))
#colnames(mo) %>% (function(v) v[which(v == 'CD107A'):which(v == 'CD62L')]) %>% str_c(collapse=' ')



# Counts of params
top_params %>% filter(param %in% colnames(mo)) %>% head(100) %>%
  mutate(param=factor(param, levels=param)) %>%
  plot_ly(x=~param, y=~n, type='bar') %>%
  plotly::layout(margin=list(b=200))

# Biggest experiments per panel size
# df_chl %>% 
#   filter(resolution == 'lookup') %>% 
#   group_by(exp_id, filename, param) %>% tally %>% select(-n) %>% ungroup %>%
#   group_by(exp_id, filename) %>% tally %>% ungroup %>%
#   group_by(exp_id) %>% summarize(n_max=max(n)) %>%
#   arrange(desc(n_max)) %>% data.frame %>% head(30)

df_chl %>% 
  filter(resolution == 'lookup') %>% 
  group_by(exp_id, filename, param) %>% tally %>% select(-n) %>% ungroup %>%
  group_by(exp_id, filename) %>% do({
    d <- .
    is_lympho <- d$param %>% str_detect('^CD3$') %>% any
    data.frame(n=nrow(d), is_lympho=is_lympho)
  }) %>% ungroup %>%
  filter(is_lympho) %>% 
  group_by(exp_id) %>% summarize(n_max=max(n)) %>%
  arrange(desc(n_max)) %>% data.frame %>% head(30)

df_chl %>% filter(exp_id == 'FR-FCM-ZYWY' & resolution == 'lookup') %>%
  group_by(param) %>% tally %>% print(n=150)


# https://stackoverflow.com/questions/33777165/correlation-using-funs-in-dplyr
# dft <- data.frame(g=c('g1', 'g1', 'g1', 'g2', 'g2'), c1=c(1, 1, 0, 0, 0), c2=c(1, 0, 1, 0, 1), c3=c(0, 0, 1, 0, 1))
# dft %>% select(-g) %>% rownames_to_column %>% gather(channel, value, starts_with('c')) %>%
#   full_join(., ., by=c('rowname')) %>% arrange(channel.x, channel.y)


