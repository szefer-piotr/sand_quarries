# Ind val

source("data/data_processing.R")
library(labdsv)

# Indicator values
ivs <- indval(asc, clustering = stages$succession)

dfr <-data.frame(Cluster = ivs$maxcls,
                 IndicatorValue = ivs$indcls,
      Probability = ivs$pval)

dfr <- dfr[dfr$Probability <= 0.05,]
dfr <- dfr[order(dfr$Cluster), ]
dfr <- with(dfr, dfr[order(-Cluster, IndicatorValue, decreasing = TRUE),])

                     