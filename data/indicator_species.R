# Ind val

source("data/data_processing.R")
library(labdsv)

set.seed(123)
# Indicator values
ivs <- indval(asc, clustering = stages$succession, numitr = 9999)

dfr <-data.frame(Cluster = ivs$maxcls,
                 IndicatorValue = ivs$indcls,
      Probability = ivs$pval)

dfr <- dfr[dfr$Probability <= 0.05,]
dfr <- dfr[order(dfr$Cluster), ]
dfr <- with(dfr, dfr[order(-Cluster, IndicatorValue, decreasing = TRUE),])

# Contingency table
gps <- cbind(colnames(asc),groups)
rownames(gps) <- gps[,1]

dfr$group <- gps[rownames(dfr),2]
ct <- table(dfr$Cluster, dfr$group)
rownames(ct) <- c("SI","SII","SIII")

# Decomposition of 3x3 contingency table ----

# Indicator species for trophic groups
ivtable.entries <- c(ct[1,],  ct[2,],  ct[3,])
ivtable <- t(as.table(matrix(ivtable.entries, nrow = 3, byrow = TRUE, dimnames = list(Stage = c('Early', 'Mid', 'Late'), Group = c('Herb', 'Klep', 'Pred')))))
addmargins(ivtable)

# Independence between Stage and Group
library(MASS)
ivtable.loglm <- loglm( ~ Stage + Group, data = ivtable) 
ivtable.loglm

# Reject the hypothesis that Stage and Group are independent. There is a strong relationship between these two.

# Subtables analysis

# Table 3.3 page 83 Agresti book

# Two rows Herb and Klep vs mid late that have similar proportions
ivtable.prop.mar.1.rows.12.table <- prop.table(ivtable[1:2,], margin = 1)
# transform to percentages
ivtable.percent.mar.1.rows.12.table <- 100*ivtable.prop.mar.1.rows.12.table
ivtable.A <- ivtable[1:2, 2:3]
ivtable.A.loglm <- loglm( ~ Group + Stage, data = ivtable.A) 
ivtable.A.loglm
# Homogenous

margin.table(ivtable.A, margin = 1)
ivtable.B <- as.table(matrix(c(ivtable[1:2, 1],margin.table(ivtable.A, margin = 1)),
                            nrow = 2, 
                            byrow = FALSE, 
                            dimnames = list(Group = c('Herb','Klep'), 
                                            Stage = c('Early', 'Mid+Late'))))  
ivtable.B
ivtable.B.loglm <- loglm( ~  Stage + Group, data = ivtable.B) 
ivtable.B.loglm
# Homogenous (but with zero in one cell)

#We could start with 2x3 table 
ivtable.Herb.Klep.loglm <- loglm( ~ Group + Stage, data = ivtable[1:2, ])
ivtable.Herb.Klep.loglm
# Seems like first two rows are homogenous. with G^2 = 1.2431761
# 0.02009123 + 1.2230849

# Since Ag.3.6.A is homogenous we could compute column marginal totals
ivtable.C <- as.table(matrix(c(margin.table(ivtable.A, margin = 2), 
                               ivtable[3, 2:3]),
                            nrow = 2, 
                            byrow = TRUE, 
                            dimnames = list(Group = c('Herb+Klep','Pred'), 
                                            Stage = c('Mid', 'Late')))) 
ivtable.C.loglm <- loglm( ~  Stage + Group, data = ivtable.C) 
ivtable.C.loglm
# Nonhomogenous - Group and Stage are related

# slopy code
ivtable.D <- as.table(matrix(c(margin.table(ivtable.B, margin = 2)[1],
                               ivtable[3,1],
                              margin.table(ivtable.C, margin = 1)),
                            nrow = 2, 
                            byrow = FALSE, 
                            dimnames = list(Group = c('Herb+Klep','Pred'), 
                                            Stage = c('Early', 'Mid+Late'))))

ivtable.D
ivtable.D.loglm <- loglm( ~  Stage + Group, data = ivtable.D) 
ivtable.D.loglm

# Test for 12 decimal places
round(ivtable.loglm$lrt,12) == round((ivtable.A.loglm$lrt+ivtable.B.loglm$lrt+ivtable.C.loglm$lrt+ivtable.D.loglm$lrt),12)

ivtable.A
ivtable.B
ivtable.C
ivtable.D
