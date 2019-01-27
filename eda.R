library(FlowRepositoryR)
library(tidyverse)

con <- url('https://flowrepository.org/experiments/1602/fcs_files/170708/download')
#con <- file('/tmp/header.txt')
open(con, "rb")
a <- readBin(con, 'raw', 65535)
c2 <- rawConnection(a)
header <- flowCore:::readFCSheader(c2)
txt <- flowCore:::readFCStext(c2, header, emptyValue = T)

file_size <- as.integer(flowCore:::readFCSgetPar(txt, "$TOT"))
n_params <- as.integer(flowCore:::readFCSgetPar(txt, "$PAR"))
param_name <- flowCore:::readFCSgetPar(txt, paste("$P", 1:n_params, "N", sep=""))
param_desc <- flowCore:::readFCSgetPar(txt, paste("$P", 1:n_params, "S", sep=""))
params <- set_names(unname(param_desc), param_name)
flowCore:::readFCSgetPar(txt, 'CREATOR')
adata <- flowCore:::makeFCSparameters(unname(param_name), txt, F, F, 0, 0)
adata@data

# FCS EDA

library(flowCore)

fr <- flowCore::read.FCS('~/Downloads/AUX6 fresh_cct.fcs') 
df <- fr %>% exprs %>% as_tibble %>% set_names(fr@parameters@data$desc %>% make.names(unique=T)) 
var <- 'gdTCR'
trans <- function(v)asinh(v)
#df %>% select(!!!var) %>% mutate_all(funs(trans)) %>% ggplot(aes_string(x=var)) + geom_histogram(bins=64)

devtools::source_gist('bf63a4833cf991a1595e3fc503856c4f')
df %>% select(CD45RA, CCR7) %>% mutate_all(funs(trans)) %>% mutate(density=get_density(CD45RA, CCR7)) %>% 
  ggplot(aes(x=CD45RA, y=CCR7, color=density)) + geom_point(size=.3, alpha=.5) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))

