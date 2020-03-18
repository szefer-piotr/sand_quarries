# Species accumulation curves

# Upload data and necessary packages

#https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html


source("data/data_processing.R")
library(iNEXT)

out <- iNEXT(stage1, q = c(0,1,2), datatype = "abundance")
ggiNEXT(out, type=3, facet.var = "site", color.var = "site")
