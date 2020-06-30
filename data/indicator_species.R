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

# Decomposition of 3x3 contingency table
# https://rstudio-pubs-static.s3.amazonaws.com/65435_a8b26773b5d64138b6411a5aa306a5b9.html
Ag.3.5.table.entries <- c(ct[1,],  ct[2,],  ct[3,])
Ag.3.5 <- as.table(matrix(Ag.3.5.table.entries, nrow = 3, byrow = TRUE, dimnames = list(Stage = c('Early', 'Mid', 'Late'), Group = c('Herb', 'Klep', 'Pred'))))
Ag.3.5 <- t(Ag.3.5)
addmargins(Ag.3.5)

# Independence between Stage and Group
library(MASS)
Ag.3.5.loglm <- loglm( ~ Stage + Group, data = Ag.3.5) 
Ag.3.5.loglm

# Reject the hypothesis that Stage and Group are independent. There is a strong relationship between these two.

# Subtables analysis

# Table 3.3 page 83 Agresti book

# Two rows Herb and Klep vs mid late that have similar proportions
Ag.3.5.prop.mar.2.rows.12.table <- prop.table(Ag.3.5[1:2,], margin = 2)
# transform to percentages
Ag.3.5.percent.mar.2.rows.12.table <- 100*Ag.3.5.prop.mar.2.rows.12.table
Ag.3.6.A <- Ag.3.5[1:2, 2:3]

Ag.3.5.prop.mar.1.cols.23.table <- prop.table(Ag.3.5[,2:3], margin = 1)
# transform to percentages
Ag.3.5.percent.mar.1.cols.23.table <- 100*Ag.3.5.prop.mar.1.cols.23.table
Ag.3.6.A <- Ag.3.5[1:2, 2:3]
Ag.3.6.A.loglm <- loglm( ~ Group + Stage, data = Ag.3.6.A) 
Ag.3.6.A.loglm
# Homogenous

margin.table(Ag.3.6.A, margin = 1)
Ag.3.6.B <- as.table(matrix(c(Ag.3.5[1:2, 1],margin.table(Ag.3.6.A, margin = 1)),
                            nrow = 2, 
                            byrow = FALSE, 
                            dimnames = list(Group = c('Herb','Klep'), 
                                            Stage = c('Early', 'Mid+Late'))))  
Ag.3.6.B
Ag.3.6.B.loglm <- loglm( ~  Stage + Group, data = Ag.3.6.B) 
Ag.3.6.B.loglm
# Homogenous

#We could start with 2x3 table 
Ag.3.5.Herb.Klep.loglm <- loglm( ~ Group + Stage, data = Ag.3.5[1:2, ])
Ag.3.5.Herb.Klep.loglm
# Seems like first two rows are homogenous. with G^2 = 1.2431761
# 0.02009123 + 1.2230849

# Since Ag.3.6.A is homogenous we could compute column marginal totals
Ag.3.6.C <- as.table(matrix(c(margin.table(Ag.3.6.A, margin = 2), 
                              Ag.3.5[3, 2:3]),
                            nrow = 2, 
                            byrow = TRUE, 
                            dimnames = list(Group = c('Herb+Klep','Pred'), 
                                            Stage = c('Mid', 'Late')))) 
Ag.3.6.C.loglm <- loglm( ~  Stage + Group, data = Ag.3.6.C) 
Ag.3.6.C.loglm
# Nonhomogenous - Group and Stage are related

# slopy code
Ag.3.6.D <- as.table(matrix(c(margin.table(Ag.3.6.B, margin = 2)[1],
                              Ag.3.5[3,1],
                              margin.table(Ag.3.6.C, margin = 1)),
                            nrow = 2, 
                            byrow = FALSE, 
                            dimnames = list(Group = c('Herb+Klep','Pred'), 
                                            Stage = c('Early', 'Mid+Late'))))

Ag.3.6.D
Ag.3.6.D.loglm <- loglm( ~  Stage + Group, data = Ag.3.6.D) 
Ag.3.6.D.loglm

# Test for 12 decimal places
round(Ag.3.5.loglm$lrt,12) == round((Ag.3.6.A.loglm$lrt+Ag.3.6.B.loglm$lrt+Ag.3.6.C.loglm$lrt+Ag.3.6.D.loglm$lrt),12)

Ag.3.6.A
Ag.3.6.B
Ag.3.6.C
Ag.3.6.D


###### PREVIOUS DECOMPOSITION
# # Two columns Early and Late that have similar proportions
# Ag.3.5.prop.mar.1.cols.13.table <- prop.table(Ag.3.5[,c(1,3)], margin = 1)
# # transform to percentages
# Ag.3.5.percent.mar.1.cols.13.table <- 100*Ag.3.5.prop.mar.1.cols.13.table
# Ag.3.6.A <- Ag.3.5[c(1,3), c(1,3)]
# Ag.3.6.A
# Ag.3.6.A.loglm <- loglm( ~ Stage + Group, data = Ag.3.6.A) 
# Ag.3.6.A.loglm
# 
# # Ho is true, this is homogenous: it is true that its row and column marginal totals contain the information that lie within the individual cells
# margin.table(Ag.3.6.A, margin = 2)
# Ag.3.6.B <- as.table(matrix(c(margin.table(Ag.3.6.A, margin = 2), 
#                               Ag.3.5[2,c(1,3)]), 
#                             nrow = 2, 
#                             byrow = TRUE, 
#                             dimnames = list(Group = c('Herb+Pred', 'Klep'), 
#                                             Stage = c('Early', 'Late'))))  
# Ag.3.6.B
# 
# Ag.3.6.B.loglm <- loglm( ~  Stage + Group, data = Ag.3.6.B) 
# Ag.3.6.B.loglm
# # H0 rejected - nonhomogenous
# 
# 
# # present ‘percentages’ table
# round(Ag.3.5.percent.mar.1.cols.12.table, 1)
# Ag.3.6.A <- Ag.3.5[1:2, 1:2]
# Ag.3.6.A
# Ag.3.6.A.loglm <- loglm( ~ Stage + Group, data = Ag.3.6.A) 
# Ag.3.6.A.loglm
# # Herb and Klep are homogenous in Early and Mid stage
# margin.table(Ag.3.6.A, margin = 1)
