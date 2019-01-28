# flowrepository-metadata-db

This project contains a script and the resulting data for collecting and analyzing metadata from [FlowRepository](https://flowrepository.org/).  The primary purpose for this is to collect FCS metadata that isn't otherwise readily available.  For example, the data and code here is good for answering questions like:

1. What datasets have at least X of N possible intracellular or surface markers?
2. What FCS files available across the whole of FlowRepository have the largest panels?
3. What published datasets are most similar to another dataset of interest?

The logic for data collection exists in [build.R](build.R), which downloads the first 32k of each FCS file, parses the header, and stores channel names along with other dataset and file level metadata.  The biggest challenge in using this information is normalizing all the
different parameter naming conventions into a common nomenclature.  An effort is made in this project to do this by using an [HGNC](https://www.genenames.org/) lookup table (from [HGNChelper](https://cran.r-project.org/web/packages/HGNChelper/index.html)) as well as parsing and matching logic in [resolve.R](resolve.R).  At TOW, this resolves cases like this:


| Channel or Parameter Name  | Result |
| ------------- | ------------- |
| 158Gd_CD33  | CD33  |
| PE CD203c  | CD203c  |
| CXCR5  | CD185  |
| IL-15Ra | CD215 |
| CD235A BIOTIN/QD 585  | CD235A  |
| CD62L-APC-Cy7 | CD62L |
| 162Dy_TIM3_Dy162Di|CD366|
| 152Sm_CD95_-_FAS|CD95|
| CD3(Er170)Di|CD3|
| 173Yb_Tbet_Yb173Di|TBX21|
| CD14.PerCP/ PI-A|CD14|

You get the picture.  The matching isn't perfect but using typical name delimiters like ```[./ -#]``` and joining the resulting terms to a master lookup table resolves most cases, and gives preference to CD designations for the sake of cross-dataset comparison.  There is also a manual list of mappings in [data/param_map_edits.csv]([data/param_map_edits.csv]) that can be altered to remap with new manual additions (or with changes made to resolve.R).  This list was built with a focus on T cell markers so more entries might be necessary to reach a good normalization for other kinds of studies.

Analyzing the data requires nothing more than tidyverse since all data collected is stored as csvs.  See the [data](data) folder for results as of Jan 2019.  See [overview.Rmd](overview.md) for table content.


